/*
*
*      Copyright (c)  2011-2026  Ivan Balaban
*      ivvaan@gmail.com
*
This file is part of Seg_int library.

Seg_int is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Seg_int is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Seg_int.  If not, see <http://www.gnu.org/licenses/>.
*/
// seg_int.cpp : Defines the entry point for the console application.
//
// seg_int.cpp : Defines the entry point for the console application.
//
#include "test_coll.h"
#include "options_base.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

#if defined(_WIN32)
#include <limits>
#endif
// Page template: svg_player.html is a C++ raw string literal, so it is #included as is and
// also opens in a browser by itself (with demo data), which is handy for debugging the player.
// MSVC limits a single string literal to 16380 bytes (C2026): keep svg_player.html below that.
const char* STYLE_temlate =
#include "svg_player.html"
;

const char* to_insert_SVG = "<!--inserthere-->";

#if defined(_WIN32) && !defined(_DEBUG)
#include <windows.h>

static void apply_bench_settings(long long affinity_cpu, int pri_preset,
  DWORD& process_pri_save, int& thread_pri_save, DWORD_PTR& system_affinity_mask)
{
  process_pri_save = GetPriorityClass(GetCurrentProcess());
  thread_pri_save = GetThreadPriority(GetCurrentThread());

  DWORD process_pri = process_pri_save;
  int thread_pri = thread_pri_save;

  switch (pri_preset) {
  case 1: process_pri = ABOVE_NORMAL_PRIORITY_CLASS; thread_pri = THREAD_PRIORITY_ABOVE_NORMAL; break;
  case 2: process_pri = HIGH_PRIORITY_CLASS; thread_pri = THREAD_PRIORITY_HIGHEST; break;
  case 3: process_pri = REALTIME_PRIORITY_CLASS; thread_pri = THREAD_PRIORITY_TIME_CRITICAL; break;
  case 0:
  default: break;
  }

  if (pri_preset)
    SetPriorityClass(GetCurrentProcess(), process_pri);

  if (pri_preset)
    SetThreadPriority(GetCurrentThread(), thread_pri);

  system_affinity_mask = 0;

  if (affinity_cpu >= 0) {
    DWORD_PTR proc_mask = 0, sys_mask = 0;
    if (GetProcessAffinityMask(GetCurrentProcess(), &proc_mask, &sys_mask) != 0) {
      system_affinity_mask = sys_mask;
      if (affinity_cpu < (long long)(sizeof(DWORD_PTR) * 8)) {
        DWORD_PTR mask = (DWORD_PTR)1 << affinity_cpu;
        if (mask & proc_mask)
          SetThreadAffinityMask(GetCurrentThread(), mask);
      }
    }
  }
}

static void restore_bench_settings(DWORD process_pri_save, int thread_pri_save, DWORD_PTR system_affinity_mask)
{
  SetThreadPriority(GetCurrentThread(), thread_pri_save);
  SetPriorityClass(GetCurrentProcess(), process_pri_save);
  (void)system_affinity_mask;
}

#else
static void apply_bench_settings(long long, int, unsigned long&, int&, unsigned long&) {}
static void restore_bench_settings(unsigned long, int, unsigned long) {}
#endif
#undef small

double benchmark(const Options& opt, PSeg seg_coll, int4 alg, double& res)
{
  double timeit, tottime = 0;
  uint4 n_call = 0;
  using namespace std::chrono;
  double mint = 1.0e10;

#if defined(_WIN32) && !defined(_DEBUG)
  DWORD process_pri_save = 0;
  int thread_pri_save = THREAD_PRIORITY_NORMAL;
  DWORD_PTR system_affinity_mask = 0;
  apply_bench_settings(opt.affinity_cpu, opt.pri_preset, process_pri_save, thread_pri_save, system_affinity_mask);
#endif

  do
  {
    auto find_intersections = get_find_intersections_func(opt.reg_stat);
    auto start = high_resolution_clock::now();
    res = find_intersections(opt, seg_coll, alg);
    tottime += timeit = static_cast<duration<double>>(high_resolution_clock::now() - start).count();

    if (timeit < mint)
      mint = timeit;
    n_call++;
  }
#ifdef NDEBUG
  while ((tottime < 3) || (n_call < opt.n_repeat));
#else
  while (false);
#endif

#if defined(_WIN32) && !defined(_DEBUG)
  restore_bench_settings(process_pri_save, thread_pri_save, system_affinity_mask);
#endif

  return mint;
}

void perform_tests(const Options& opt, PSeg seg_coll)
{
  double exec_time[33], nInt[33];
  const char* ss = "Llagi";
  const char* sd = "rlmspc";
  const char* stat_names[] = { "inters. numb","max inters. per segm","inters. numb","inters. numb" ,"inters. numb" };

  if (!opt.print_less)printf("\ntesting... ************************************************\n");

  for (int4 a = sizeof(alg_list) / sizeof(alg_list[0]) - 1; a > -1; a--) {
    if (opt.alg & alg_list[a]) {
      exec_time[a] = -1;
      exec_time[a] = benchmark(opt, seg_coll, alg_list[a], nInt[a]);

      if (opt.print_less)
        printf("a%i;s%c;d%c;S%i;n%i;i%13.0f;t%6.5f;p%f\n", alg_list[a], ss[opt.seg_type], sd[opt.distr_type], opt.random_seed, opt.n, nInt[a], exec_time[a], opt.distr_param);
      else
        printf("alg=%s; %s=%13.0f; exec time=%6.5f\n", alg_names[a], stat_names[opt.reg_stat], nInt[a], exec_time[a]);
    }
  };

  if (opt.rtime_printout && (!opt.print_less) && (opt.reg_stat != _Registrator::per_segm_reg_max_per_segm_stat)) {
    double check_time = -1;
    if ((opt.alg & triv)) {
      double n_checks = ((double)opt.n) * 0.5 * (opt.n - 1);
      check_time = exec_time[0] / n_checks;
      printf("\ntrivial alg made %13.0f intersection checks\n", n_checks);
      printf("one intersection check takes %6.3f ns, let's use intersection check time (ICT) as a time unit \n\n", 1E+09 * check_time);
    }
    else if (opt.ICT > 0) {
      printf("\none intersection check presumably takes %6.3f ns, let's use this value as a time unit (ICT) \n\n", opt.ICT);
      check_time = 1E-09 * opt.ICT;
    };
    if (check_time > 0)
      for (int4 a = 0; a < sizeof(alg_list) / sizeof(alg_list[0]); a++) {
        if (opt.alg & alg_list[a]) {
          if (exec_time[a] > 0)
            printf("%s finds one intersection in %6.3f ICT \n", alg_names[a], exec_time[a] / check_time / nInt[a]);
        }
      };
  }
};

void WriteSvgHtml(const Options& opt, PSeg seg_coll)
{
  if(!opt.fname)
    return;
  const char* insert_pos = strstr(STYLE_temlate, to_insert_SVG);
  assert(insert_pos != nullptr);
  if (!insert_pos)
    return;

  std::ostringstream svg;
  write_SVG(&svg, opt, seg_coll, opt.alg);

  std::ofstream svgf(opt.fname);
  if (!svgf.is_open())
    return;
  svgf.write(STYLE_temlate, insert_pos - STYLE_temlate);
  svgf << svg.str();
  svgf << insert_pos;
};

int main(int argc, char* argv[])
{
  Options opt;

#ifdef _DEBUG
  _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF);
#endif

  if(!opt.ParseArgs(argc, argv))
  return 0;

  opt.CheckOptions();
  opt.PrintActual();

  CRandomValueGen rv(opt.random_seed);
  auto seg_coll = create_test_collection(opt, rv);

  perform_tests(opt, seg_coll);

  WriteSvgHtml(opt, seg_coll);

  if (opt.wait) { printf("\npress 'Enter' to continue"); getchar(); }
  return 0;
}









