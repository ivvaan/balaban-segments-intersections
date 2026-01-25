/*
*
*      Copyright (c)  2011-2020  Ivan Balaban 
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

#include "test_coll.h"
#include "utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

const char* STYLE_temlate = R"WYX(
<!DOCTYPE html>
<html>
<style>
line{stroke:blue;
stroke-width:1px;
vector-effect: non-scaling-stroke;}
path{stroke:blue;
stroke-width:1px;
fill:none;
vector-effect: non-scaling-stroke;}
.int_line{stroke:green;}
.line1{stroke:cadetblue;}
.line2{stroke:darkblue;}
.edge{stroke:black;}
.arc{stroke:teal;}
circle{r:0.001;fill: gold;}
.triv {r:0.0010;fill: purple;}
.ssw {r:0.0016;fill: sandybrown;}
.fast {r:0.0020;fill: orange;}
.optimal {r:0.0024;fill: khaki;}
.parallel {r:0.0028;fill: darkorange;}
.mem_save {r:0.0032;fill: salmon;}
</style>
<script>
var state=1;
var dir=1;
var zoom=function(){
svg=document.getElementsByTagName('svg')[0];
mul=1.25;
if(dir<0) {mul=0.8;}
var h=parseFloat(svg.getAttribute('height'));
svg.setAttribute('height',h*mul+'%');
var w=parseFloat(svg.getAttribute('width'));
svg.setAttribute('width',w*mul+'%');
state+=dir;
if(state>=8){dir=-1;}
if(state<=1){dir=1;}
}
</script>
<body ondblclick='zoom()'>
<!--inserthere-->
</body>
</html>
)WYX";

const char* to_insert_SVG = "<!--inserthere-->";


#if defined(_WIN32) && !defined(_DEBUG)
#include <windows.h>

#define ON_BENCH_BEG \
auto hThread=GetCurrentThread();\
auto hProcess=GetCurrentProcess();\
auto dwProcessPriSave = GetPriorityClass(hProcess);\
SetPriorityClass(hProcess, REALTIME_PRIORITY_CLASS);\
auto dwThreadPriSave = GetThreadPriority(GetCurrentThread()); \
if(thread_affinity_mask){SetThreadPriority(hThread, THREAD_PRIORITY_TIME_CRITICAL);\
SetThreadAffinityMask(hThread, thread_affinity_mask);}

#define ON_BENCH_END \
SetThreadPriority(hThread, dwThreadPriSave);\
SetPriorityClass(hProcess, dwProcessPriSave);
#undef small
#else
#define ON_BENCH_BEG
#define ON_BENCH_END
#endif

//#include <tchar.h>

#ifdef _DEBUG
#include <crtdbg.h>
#endif

struct Options
{
  int4 alg = 31;
  int4 seg_type = 2;
  int4 n = 20000;
  int4 distr_type = 0;
  REAL distr_param = 1.0;

  BOOL print_less = FALSE;
  BOOL wait = FALSE;
  BOOL rtime_printout = FALSE;
  bool use_counters = false;

  uint4 reg_stat = 2;
  const char* fname = nullptr;

  double ICT = -1;
  unsigned random_seed = 317;

  int4 range_for_int_seg = full_int_range;
};

double _benchmark_new(const Options& opt, PSeg seg_coll, int4 alg, double& res)
{
  double timeit, tottime = 0;
  int4 n_call = 0;
  using namespace std::chrono;
  double mint = 1.0e10;
  unsigned long thread_affinity_mask = alg & fast_parallel ? 0 : 8;

  ON_BENCH_BEG
  do
  {
    auto find_intersections = get_find_intersections_func(opt.reg_stat);
    auto start = high_resolution_clock::now();
    res = find_intersections(opt.seg_type, opt.n, seg_coll, alg, opt.reg_stat, opt.range_for_int_seg);
    tottime += timeit = static_cast<duration<double>>(high_resolution_clock::now() - start).count();

    if (timeit < mint)
      mint = timeit;
    n_call++;
  }
#ifdef NDEBUG
  while ((tottime < 3) || (n_call < 3));
#else
  while (false);
#endif
  ON_BENCH_END

  return mint; //tottime / n_call;
}

void perform_tests(const Options& opt, PSeg seg_coll)
{
  double exec_time[33], nInt[33];
  const char* ss = "Llag";
  const char* sd = "rlmspc";
  char counters_string[256], * counters_mute;
  const char* stat_names_new[] = { "inters. numb","max inters. per segm","inters. numb","inters. numb" };
  const char** stat_names = stat_names_new;

  counters_string[0] = 0;

#ifdef COUNTERS_ON
  if (opt.use_counters)strcpy(counters_string, "\ncounters:");
#endif
  counters_mute = counters_string + strlen(counters_string);
  if (!opt.print_less)printf("\nnew implementation testing... ************************************************\n");

  for (int4 a = sizeof(alg_list) / sizeof(alg_list[0]) - 1; a > -1; a--)
  {
    if (opt.alg & alg_list[a])
    {
      exec_time[a] = -1;
      exec_time[a] = _benchmark_new(opt, seg_coll, alg_list[a], nInt[a]);

      if (opt.print_less)
        printf("a%i;s%c;d%c;S%i;n%i;i%13.0f;t%6.5f;p%f;%s\n", alg_list[a], ss[opt.seg_type], sd[opt.distr_type], opt.random_seed, opt.n, nInt[a], exec_time[a], opt.distr_param, counters_mute);
      else
        printf("alg=%s; %s=%13.0f; exec time=%6.5f;%s\n", alg_names[a], stat_names[opt.reg_stat], nInt[a], exec_time[a], counters_string);
    }
  };

  if (opt.rtime_printout && (!opt.print_less) && (opt.reg_stat != 1)) {
    if ((opt.alg & triv)) {
      double n_checks = ((double)opt.n) * 0.5 * (opt.n - 1);
      double check_time = exec_time[0] / n_checks;
      printf("\ntrivial alg made %13.0f intersection checks\n", n_checks);
      printf("one intersection check takes %6.3f ns, let's use intersection check time (ICT) as a time unit \n\n", 1E+09 * check_time);
      for (int4 a = 0; a < sizeof(alg_list) / sizeof(alg_list[0]); a++)
      {
        if (opt.alg & alg_list[a])
        {
          if (exec_time[a] > 0)
            printf("%s finds one intersection in %6.3f ICT \n", alg_names[a], exec_time[a] / check_time / nInt[a]);
        }
      };
    }
    else if (opt.ICT > 0) {
      printf("\none intersection check presumably takes %6.3f ns, let's use this value as a time unit (ICT) \n\n", opt.ICT);
      double check_time = 1E-09 * opt.ICT;
      for (int4 a = 0; a < sizeof(alg_list) / sizeof(alg_list[0]); a++)
      {
        if (opt.alg & alg_list[a])
        {
          if (exec_time[a] > 0)
            printf("%s finds one intersection in %6.3f ICT \n", alg_names[a], exec_time[a] / check_time / nInt[a]);
        }
      };

    }
  }
};


static void TryWriteSvgHtml(const Options& opt, PSeg seg_coll)
{
  if (!opt.fname || (opt.n > max_SVG_items))
    return;


  std::ostringstream svg;
  write_SVG(&svg, opt.seg_type, opt.n, seg_coll, opt.alg, opt.reg_stat,opt.range_for_int_seg);

  const char* insert_pos = strstr(STYLE_temlate, to_insert_SVG);
  if (!insert_pos)
    return;
  std::ofstream svgf(opt.fname);
  if (!svgf.is_open())
    return;
  svgf.write(STYLE_temlate, insert_pos - STYLE_temlate);
  svgf << svg.str();
  svgf << "</svg>\n";
  svgf << insert_pos;
}
int main(int argc, char* argv[])
{
  Options opt;
  const char* ss = "Llagi", * sd = "rlmspc", * sr = "pPcrC";
  int4 range_for_int_seg;

#ifdef _DEBUG
  _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF);
#endif


  if (argc == 1)
  {
#ifdef COUNTERS_ON
    printf("usage: seg_int -aA -sS -SR -dD -nN -pP -rR -m -eE -w -SN -fhtmfile -c\n");
#else    
    printf("usage: seg_int -aA -sS -SR -dD -nN -pP -rR -m -eE -w -SN -fhtmfile\n");
#endif
    printf(
      R"WYX(example: seg_int -a14 -sa -dp -n20000 -p5.5
-aA: type of algorithms tested
 A=1: trivial
 A=2: simple sweep
 A=4: balaban fast
 A=8: balaban optimal
 A=16: balaban fast parallel for 6 treads
 if you want several algorithms to test just sum up their values
 i.e. A=31 (=1+2+4+8+16) all algorithms
-sS: type of segments
 S=l: line segments representation y=a*x+b,x1<=x<=x2; a,b,x1,x2 - reals
 S=L: line segments representation r=b+a*t,0<=t<=1; b,a - vectors
 S=a: arcs
 S=g: graph: test example is simple circle of N vertices and N edges
 S=i: same as L but vectors a and b have integer coords; if a number is appended, e.g. -si100, coordinates are in range -100 to 100; if not, range is 1/2 of int range
-nN: number of segments
-pP: parameter value, must be positive
-dD: type of distribution
 D=r: random segments
 D=l: random length almost x parallel segments, the bigger
  distr_param/N the less parallel segments
 D=m: mixed random length almost x parallel 'long' segments
  (33%%) and 'small' segments(67%%), the bigger distr_param/N
  the less parallel 'long' and longer 'short' segments
 D=s: short segments: random segment with  length multiplied by distr_param/N
 D=distr_param: random segment with  length multiplied by distr_param
 D=c: segments ends are on the opposite sides of unit circle, each
  segment intesect each
-rR: type of registrar used in new implementation and type of result statistic
 R=c: total intersection counting registrator; total count statistic
 R=C: total intersection counting registrator; total count statistic, but it pretend it register intersection points
 R=p: total intersection counting and per segment intersection counting registrator; 
  total count statistic
 R=P: total intersection counting and per segment intersection counting registrator; 
  max intersections pre segment statistic
 R=r: really storing pairs and intersections registrator (be carefull with memory!!!); 
  total count statistic. As we can have O(N^2) int. the option is limited to N=20000 max
-SR: capital S for random seed; R - random seed value; if R=0 - non pseudo random generator used
-fhtmfile: if presented, program writes SVG picture to htmfile. For examle -fC:/tmp/res.htm
  To limit resulting file size option works only for 5000 segments and less, also only first 
  150000 intersections are drawn.
-m: print truncated information in one row (to make a table)
-w: stop and wait for an input befor exit
)WYX"
    );

#ifdef COUNTERS_ON
    printf("-c: counters are printed, if presented\n");
#endif    
    return 0;
  }
  else
  {
    for (int i = 1; i < argc; i++)
      if (argv[i][0] == '-')
      {
        switch (argv[i][1])
        {
        case 'a':
        {
          opt.alg = atoi(argv[i] + 2);
          if ((opt.alg < 1) || (opt.alg > 31))
          {
            opt.alg = 4; printf("some error in -a param. 4 used instead.\n");
          }
        }
        break;
        case 's':
        {
          switch (argv[i][2])
          {
          case 'L':opt.seg_type = _Segment::line1; break;
          case 'l':opt.seg_type = _Segment::line2; break;
          case 'a':opt.seg_type = _Segment::arc; break;
          case 'g':opt.seg_type = _Segment::graph; break;
          case 'i': {
            opt.seg_type = _Segment::intline;
            opt.range_for_int_seg = atoi(argv[i] + 3);
            if (opt.range_for_int_seg == 0) opt.range_for_int_seg = full_int_range;
            if (opt.range_for_int_seg < 5 || opt.range_for_int_seg > full_int_range) {
              printf("some error in -si param. full_int_range used as a range.\n");
              opt.range_for_int_seg = full_int_range;
            };
          };
                  break;
          default:
          {
            printf("some error in -s param. l used instead.\n");
            opt.seg_type = _Segment::line2;
          }
          };
        };
        break;
        case 'd':
        {
          switch (argv[i][2])
          {
          case 'r':opt.distr_type = _Distribution::random; break;
          case 'l':opt.distr_type = _Distribution::parallel; break;
          case 'm':opt.distr_type = _Distribution::mixed; break;
          case 's':opt.distr_type = _Distribution::small; break;
          case 'p':opt.distr_type = _Distribution::param_defined; break;
          case 'c':opt.distr_type = _Distribution::circle; break;
          default:
          {
            printf("some error in -distr_type param. r used instead.\n");
            opt.distr_type = _Distribution::random;
          }
          };
        };
        break;
        case 'r':
        {
          switch (argv[i][2])
          {
          case 'c':opt.reg_stat = _Registrator::just_count; break;
          case 'C':opt.reg_stat = _Registrator::just_count_require_intersections; break;
          case 'p':opt.reg_stat = _Registrator::per_segm_reg_just_count_stat; break;
          case 'P':opt.reg_stat = _Registrator::per_segm_reg_max_per_segm_stat; break;
          case 'r':opt.reg_stat = _Registrator::store_pairs_and_ints_just_count_stat; break;
          default:
            opt.reg_stat = 2;
          };
        };
        break;
        case 'n':
        {
          opt.n = atoi(argv[i] + 2);
          if ((opt.n < 2))
          {
            opt.n = 10000; printf("some error in -n param. 10000 used instead.\n");
          }
        }
        break;
        case 'S':
        {
          opt.random_seed = atoi(argv[i] + 2);
        }
        break;
        case 'p':
        {
          opt.distr_param = fabs(atof(argv[i] + 2));
          if (opt.distr_param <= 0)
          {
            opt.distr_param = 1.0; printf("some error in -distr_param param. 1.0 used instead.\n");
          }
        }
        break;
        case 'm':
        {
          opt.print_less = TRUE;
        }
        break;
        case 'e':
        {
          opt.ICT = fabs(atof(argv[i] + 2));
          opt.rtime_printout = TRUE;
        }
        break;
        case 'w':
        {
          opt.wait = TRUE;
        }
        break;
        case 'f':
        {
          opt.fname = argv[i] + 2;
        }
        break;
#ifdef COUNTERS_ON
        case 'c':
        {
          opt.use_counters = true;
        }
        break;
#endif
        }
      }
  }

  if ((opt.seg_type == _Segment::graph) && (opt.distr_type != _Distribution::random)) { printf("-sg  is compartible only with -dr!\n"); if (opt.wait) { printf("\npress 'Enter' to continue"); getchar(); } return 0; }
  if ((opt.reg_stat == _Registrator::store_pairs_and_ints_just_count_stat) && (opt.n > max_truereg_items)) {
    printf("Too many segments (and possible intersections!) for true registration! Just counting used instead.\n");
    opt.reg_stat = _Registrator::just_count;
  }

  if (!opt.print_less) {
    if (opt.seg_type != _Segment::intline)
      printf("params are: -a%i -s%c -d%c -r%c -n%i -S%i -p%f\n", opt.alg, ss[opt.seg_type], sd[opt.distr_type], sr[opt.reg_stat], opt.n, opt.random_seed, opt.distr_param);
    else
      printf("actual params is: -a%i -si%i -d%c -r%c -n%i -S%i -p%f -e%f\n", opt.alg, opt.range_for_int_seg, sd[opt.distr_type], sr[opt.reg_stat], opt.n, opt.random_seed, opt.distr_param, opt.ICT);
  }
  CRandomValueGen rv(opt.random_seed);
  PSeg seg_coll = create_test_collection(opt.seg_type, opt.n, opt.distr_type, opt.distr_param, rv, nullptr);

  TryWriteSvgHtml(opt, seg_coll);

  perform_tests(opt, seg_coll);
  delete_test_collection(opt.seg_type, seg_coll, nullptr);

  if (opt.wait) { printf("\npress 'Enter' to continue"); getchar(); }
  return 0;
}



