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
#include <strstream>

const char* STYLE_temlate = R"WYX(
<!DOCTYPE html>
<html>
<style>
line{stroke:blue;
stroke-width:0.001;}
path{stroke:blue;
stroke-width:0.001;
fill:none;}
.line1{stroke:cadetblue;}
.line2{stroke:darkblue;}
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
SetThreadPriority(hThread, THREAD_PRIORITY_TIME_CRITICAL);

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

double ICT = -1;

unsigned random_seed = 317;

void search_problem(uint4 n) {
  double nfast = 0;
  double ntriv = 0;
  uint4 r_seed = 1;
  CRandomValueGen rv(r_seed);
  auto find_intersections = get_find_intersections_func(_Registrator::just_count);
  while (nfast == ntriv) {
    rv.SetSeeed(++r_seed);
    auto seg_coll = create_test_collection(_Segment::arc, n, _Distribution::random, 1.0, rv, nullptr);
    nfast = find_intersections(_Segment::arc, n, seg_coll, _Algorithm::fast, _Registrator::just_count);
    ntriv = find_intersections(_Segment::arc, n, seg_coll, _Algorithm::simple_sweep, _Registrator::just_count);
    if (r_seed % 1000000 == 0)printf("%i ", r_seed/1000000);
    
    delete_test_collection(_Segment::arc, seg_coll, nullptr);

  }
  printf("\nfound %i \n", r_seed);
  {printf("\npress 'Enter' to continue"); getchar(); }
};


double _benchmark_old(char* counters_string, int4 n, PSeg* seg_ptr_coll, int4 s_type, int4 alg, double& res, bool dont_need_ip = false)
  {
  double timeit,tottime=0;
  double c[8] = {0,0,0,0,0,0,0,0};
  int4 n_call=0;
  using namespace std::chrono;
  double mint = 1.0e10;
  do
  {
      auto start = high_resolution_clock::now();
      res = find_intersections(s_type, n, seg_ptr_coll, alg, c, dont_need_ip);
	  tottime+=timeit = static_cast<duration<double>>(high_resolution_clock::now() - start).count();
      if (timeit < mint)
          mint = timeit;
	  n_call++;
  }
#ifndef _DEBUG  
  while ((tottime < 3)||(n_call<3));
#else
  while (false);
#endif
  if(counters_string)
   sprintf(counters_string,"%11.2f;%11.2f;%11.2f;%11.2f;%11.2f;%11.2f;%11.2f;%11.2f",c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);
   //printf("avg Size %11.2f\n",c[3]/(c[4]+1));

  //return (stop - start) / (n_call*CLOCKS_PER_SEC);
  return mint; //tottime / n_call;
  }

double _benchmark_new(int4 n, PSeg seg_coll, int4 s_type, int4 alg, double& res,uint4 reg_stat)
  {
  double timeit, tottime = 0;
  int4 n_call = 0;
  using namespace std::chrono;
  double mint = 1.0e10;
  do
  {
    auto find_intersections = get_find_intersections_func(reg_stat);
    auto start = high_resolution_clock::now();
    res = find_intersections(s_type, n, seg_coll,alg,reg_stat);
    tottime += timeit = static_cast<duration<double>>(high_resolution_clock::now() - start).count();
    set_SVG_stream(nullptr);
    if (timeit < mint)
      mint = timeit;
    n_call++;
  } 
#ifdef NDEBUG
  while ((tottime < 3) || (n_call < 3));
#else
  while (false);
#endif
  //return (stop - start) / (n_call*CLOCKS_PER_SEC);
  return mint; //tottime / n_call;
}

enum _Implementation
{
  impl_old = 1,
  impl_new = 2
};

void perform_tests(bool use_counters,int4 impl,int4 alg,int4 seg_type,int4 distr_type,REAL distr_param,bool print_less,bool rtime_printout,bool dont_need_ip,int4 n, PSeg seg_coll,PSeg *seg_ptr_coll,uint4 reg_stat) {
  double exec_time[33], nInt[33];
  const char *ss = "Llag", *sd = "rlmspc";
  char counters_string[256],*counters_mute;
  int4 alg_list[] = { triv, simple_sweep, fast, optimal, fast_parallel, bentley_ottmann,fast_no_ip,mem_save };
  const char *alg_names[] = { "trivial","simple_sweep","fast","optimal","fast_parallel","bentley_ottmann","fast no inters points","fast 'no R'" };
  const char *stat_names_new[] = { "inters. numb","max inters. per segm","inters. numb","inters. numb" };
  const char *stat_names_old[] = { "inters. numb","inters. numb","inters. numb","inters. numb" };
  const char **stat_names= stat_names_new;

  counters_string[0] = 0;

#ifdef COUNTERS_ON
  if(use_counters)strcpy(counters_string, "\ncounters:");
#endif
  counters_mute = counters_string + strlen(counters_string);
  if (impl == _Implementation::impl_old)
  {
    if (!print_less)printf("\nold implementation testing... ************************************************\n");
    stat_names = stat_names_old;
  }
  else
    if (!print_less)printf("\nnew implementation testing... ************************************************\n");

  ON_BENCH_BEG

    if ((impl == impl_old)&& (seg_type == graph)) { if (!print_less)printf("old implementation can't deal with graph segments\n"); return; }
    auto SVG_str = get_SVG_stream();
    for (int4 a = sizeof(alg_list) / sizeof(alg_list[0]) - 1; a>-1; a--)
      //for (int4 a = 0;a<sizeof(alg_list) / sizeof(alg_list[0]); a++)
    {
      if (alg&alg_list[a])
      {
        set_SVG_stream(SVG_str);
        exec_time[a]=-1;
        if ((alg_list[a] == fast_no_ip) && (seg_type == arc)) { if (!print_less)printf("fast no inters. points algorithm can handle only line segments\n"); continue; }
        if ((alg_list[a] == bentley_ottmann) && (impl == impl_new)) { if (!print_less)printf("new implementation does't support segments functions for Bentley & Ottman algorithm\n"); continue; }
        if ((alg_list[a] == fast_no_ip) && (impl == impl_new)) { if (!print_less)printf("new implementation does't support segments functions for no inters. points algorithm \n"); continue; }
        if (impl == impl_old)
          exec_time[a] = _benchmark_old(use_counters?counters_mute:NULL, n, seg_ptr_coll, seg_type, alg_list[a], nInt[a], dont_need_ip);
        else
          exec_time[a] = _benchmark_new(n, seg_coll, seg_type, alg_list[a], nInt[a],reg_stat);
        
        if (print_less)
          printf("I%i;a%i;s%c;d%c;S%i;n%i;i%13.0f;t%6.5f;p%f;%s\n", impl, alg_list[a], ss[seg_type], sd[distr_type], random_seed, n, nInt[a], exec_time[a], distr_param, counters_mute);
        else
          printf("alg=%s; %s=%13.0f; exec time=%6.5f;%s\n", alg_names[a],stat_names[reg_stat],  nInt[a], exec_time[a], counters_string);
      }
    };
  ON_BENCH_END

    if (rtime_printout && (!print_less) && (reg_stat != 1)) {
      if ((alg & triv) ){
        double n_checks = n * 0.5 * (double)(n - 1);
        double check_time = exec_time[0] / n_checks;
        printf("\ntrivial alg made %13.0f intersection checks\n", n_checks);
        printf("one intersection check takes %6.3f ns, let's use intersection check time (ICT) as a time unit \n\n", 1E+09 * check_time);
        for (int4 a = 0; a < sizeof(alg_list) / sizeof(alg_list[0]); a++)
        {
          if (alg & alg_list[a])
          {
            if (exec_time[a] > 0)
              printf("%s finds one intersection in %6.3f ICT \n", alg_names[a], exec_time[a] / check_time / nInt[a]);
          }
        };
      }
      else if(ICT>0){
        printf("\none intersection check presumably takes %6.3f ns, let's use this value as a time unit (ICT) \n\n", ICT);
        double check_time = 1E-09 * ICT ;
        for (int4 a = 0; a < sizeof(alg_list) / sizeof(alg_list[0]); a++)
        {
          if (alg & alg_list[a])
          {
            if (exec_time[a] > 0)
              printf("%s finds one intersection in %6.3f ICT \n", alg_names[a], exec_time[a] / check_time / nInt[a]);
          }
        };

      }
    }

};


int main(int argc, char* argv[])
  {
  int4 alg=255, seg_type=2, n = 20000, distr_type=0;
  REAL distr_param=1.0;
  BOOL print_less=FALSE,print_counters=FALSE,wait=FALSE, rtime_printout=FALSE;
  bool use_counters = false;
  char rpar[]="-r";
  int4 impl = _Implementation::impl_old + _Implementation::impl_new;
  const char *ss = "Llag", *sd = "rlmspc",*sr="pPcr";
  const char *alg_names[] = { "trivial","simple_sweep","fast","optimal","fast_parallel","bentley_ottmann","fast no inters points","fast 'no R'" };
  uint4 reg_stat = 2;
  const char* fname = nullptr;

#ifdef _DEBUG
  _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF);
#endif
  

  if (argc == 1)
  {
#ifdef COUNTERS_ON
    printf("usage: seg_int -aA -sS -dD -nN -pP -iI -rR -m -w -c\n");
#else    
    printf("usage: seg_int -aA -sS -dD -nN -pP -rR -m -e -w -SN -fhtmfile\n");
#endif    
    printf("example: seg_int -a14 -sa -dp -n20000 -p5.5\n");
    printf("-aA: type of algorithms tested\n");
    printf(" A=1: trivial\n");
    printf(" A=2: simple sweep\n");
    printf(" A=4: balaban fast\n");
    printf(" A=8: balaban optimal\n");
    printf(" A=16: balaban fast parallel for 6 treads\n");
    printf(" A=32: bentley & ottmann\n");
    printf(" A=64: balaban fast;  intersection points aren't reported\n");
    printf("  (only intersecting pairs)\n");
    printf(" A=128: balaban fast; memory save algorithm 4 bytes less per segment\n");
    printf(" if you want several algorithms to test just sum up their values\n");
    printf(" i.e. A=255 (=1+2+4+8+16+32+64+128) all algorithms\n");
    printf("-sS: type of segments\n");
    printf(" S=l: line segments representation y=a*x+b,x1<=x<=x2; a,b,x1,x2 - reals\n");
    printf(" S=L: line segments representation r=b+a*t,0<=t<=1; b,a - vectors\n");
    printf(" S=a: arcs\n");
    printf(" S=g: graph: test example is simple circle of N vertices and N edges\n");
    printf("-iI: implementation\n");
    printf(" I=1: 'old' function pointers segment collection interface\n");
    printf(" I=2: 'new' template based segment collection interface\n");
    printf(" I=3: both\n");
    printf("-nN: number of segments\n");
    printf("-pP: parameter value, must be positive\n");
    printf("-dD: type of distribution\n");
    printf(" D=r: random segments\n");
    printf(" D=l: random length almost x parallel segments, the bigger\n");
    printf("  distr_param/N the less parallel segments\n");
    printf(" D=m: mixed random length almost x parallel 'long' segments\n");
    printf("  (33%%) and 'small' segments(67%%), the bigger distr_param/N\n");
    printf("  the less parallel 'long' and longer 'short' segments\n");
    printf(" D=s: short segments: random segment with  length multiplied by distr_param/N\n");
    printf(" D=distr_param: random segment with  length multiplied by distr_param\n");
    printf(" D=c: segments ends are on the opposite sides of unit circul, each\n");
    printf("  segment intesect each\n");

    printf("-rR: type of registrator used in new implementation and type of result statistic\n");
    printf(" R=c: total intersection counting registrator; total count statistic\n");
    printf(" R=p: total intersection counting and per segment intersection counting registrator; total count statistic\n");
    printf(" R=P: total intersection counting and per segment intersection counting registrator; max intersections pre segment statistic\n");
    printf(" R=r: really storing pairs and intersections registrator (be carefull with memory!!!); total count statistic\n");
    printf("  as we can have O(N^2) intersections this option is limited to N=20000 max\n");

    printf("-SR: capital S for random seed\n");
    printf(" R=0: no seed, randomized generation\n");
    printf(" R!=0: R is the seed\n");

    printf("-m: if presented, means 'print less' mode - few information printed\n");

    printf("-eE: if presented, for each alg prints out relative (compared to checking two segments' intersection ICT) \n");
    printf(" time to find one intersection (works only  print_less off)\n");
    printf(" E: assumed ICT (nanoseconds) - intersection of two segments check (and register if found) time\n");

    printf("-w: if presented, program wait until 'Enter' pressed before closing\n");
    
    printf("-fhtmfile: if presented, program writes SVG picture to htmfile. For examle -fC:/tmp/res.htm\n");
    printf(" intersection points are drawn only with -rr option (otherwise they are not stored)\n");
    printf(" to limit resulting file size only first %i intersections are drawn\n", max_SVG_points);
    
#ifdef COUNTERS_ON
    printf("-c: counters are printed, if presented\n");
#endif    
    printf("Important! -sa is not compartible with -dm, -dl and -r options!\n");
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
            alg = atoi(argv[i] + 2);
            if ((alg < 1) || (alg > 255))
            {
              alg = 4; printf("some error in -a param. 4 used instead.\n");
            }
          }
          break;
          case 'i':
          {
            impl = atoi(argv[i] + 2);
            if ((impl < 1) || (impl > 3))
            {
              impl = 3; printf("some error in -i param. 3 used instead.\n");
            }
          }
          break;
          case 's':
          {
            switch (argv[i][2])
            {
            case 'L':seg_type = line1; break;
            case 'l':seg_type = line2; break;
            case 'a':seg_type = arc; break;
            case 'g':seg_type = graph; break;
            default:
            {
              printf("some error in -s param. l used instead.\n");
              seg_type = line2;
            }
            };
          };
          break;
          case 'd':
          {
            switch (argv[i][2])
            {
            case 'r':distr_type = random; break;
            case 'l':distr_type = parallel; break;
            case 'm':distr_type = mixed; break;
            case 's':distr_type = small; break;
            case 'p':distr_type = param_defined; break;
            case 'c':distr_type = circul; break;
            default:
            {
              printf("some error in -distr_type param. r used instead.\n");
              distr_type = random;
            }
            };
          };
          break;
          case 'r':
          {
            switch (argv[i][2])
            {
            case 'c':reg_stat = _Registrator::just_count; break;
            case 'p':reg_stat = _Registrator::per_segm_reg_just_count_stat; break;
            case 'P':reg_stat = _Registrator::per_segm_reg_max_per_segm_stat; break;
            case 'r':reg_stat = _Registrator::store_pairs_and_ints_just_count_stat; break;
            default:
              reg_stat = 2;
            };
          };
          break;
          case 'n':
          {
            n = atoi(argv[i] + 2);
            if ((n < 2))
            {
              n = 10000; printf("some error in -n param. 10000 used instead.\n");
            }
          }
          break;
          case 'S':
          {
            random_seed = atoi(argv[i] + 2);
          }
          break;
          case 'p':
          {
            distr_param = fabs(atof(argv[i] + 2));
            if (distr_param <= 0)
            {
              distr_param = 1.0; printf("some error in -distr_param param. 1.0 used instead.\n");
            }
          }
          break;
          case 'm':
          {
            print_less = TRUE;
          }
          break;
          case 'e':
          {
            ICT = fabs(atof(argv[i] + 2));
            rtime_printout = TRUE;
          }
          break;
          case 'w':
          {
            wait = TRUE;
          }
          break;
          case 'f':
          {
            fname = argv[i] + 2;;
          }
          break;
#ifdef COUNTERS_ON
          case 'c':
          {
            use_counters = true;
          }
          break;
  #endif
          }
      }
  }
  if((seg_type==graph)&&(distr_type!=random)) {printf("-sg  is compartible only with -dr!\n"); if (wait) { printf("\npress 'Enter' to continue"); getchar(); } return 0;}
  if((seg_type==arc)&&(distr_type==mixed)) {printf("-sa -dm is not compartible!\n"); if (wait) { printf("\npress 'Enter' to continue"); getchar(); } return 0;}                  
  if((seg_type==arc)&&(distr_type==parallel)) {printf("-sa -dl is not compartible!\n"); if (wait) { printf("\npress 'Enter' to continue"); getchar(); } return 0;}
  if ((reg_stat == _Registrator::store_pairs_and_ints_just_count_stat) && (n > 20000)) {
    printf("Too many segments (and possible intersections!) for true registration! Just counting used instead.\n"); 
    reg_stat == _Registrator::just_count;
    }

  if (!print_less)printf("actual params is: -a%i -s%c -d%c -r%c -i%i -n%i -S%i -p%f -e%f\n", alg, ss[seg_type], sd[distr_type],sr[reg_stat], impl, n, random_seed, distr_param,ICT);
  //search_problem(n);
  //return 0;
  PSeg *seg_ptr_coll = nullptr;
  PSeg seg_coll;
  CRandomValueGen rv(random_seed);

  if(impl&_Implementation::impl_old)
    seg_coll = create_test_collection(seg_type, n, distr_type, distr_param,rv, &seg_ptr_coll);
  else
    seg_coll = create_test_collection(seg_type, n, distr_type, distr_param,rv, nullptr);

  bool dont_need_ip = (alg & fast_no_ip) != 0;
  std::ostrstream SVG_str;
  std::ofstream svgf;
  bool stream_ok=fname && (svgf.open(fname), svgf.is_open());
  if (stream_ok) {
    coll_to_SVG(seg_coll, seg_type, n, &SVG_str);
    set_SVG_stream(&SVG_str);
  }
  if (impl&_Implementation::impl_old)
    perform_tests(use_counters, _Implementation::impl_old, alg, seg_type, distr_type, distr_param, print_less, rtime_printout, dont_need_ip, n, seg_coll, seg_ptr_coll,reg_stat);
  if (impl&_Implementation::impl_new)
    perform_tests(use_counters, _Implementation::impl_new, alg, seg_type, distr_type, distr_param, print_less, rtime_printout, dont_need_ip, n, seg_coll, seg_ptr_coll,reg_stat);

  
  if (stream_ok) {
    auto insert_pos = strstr(STYLE_temlate, to_insert_SVG);
    {
      std::ostream_iterator<char> it(svgf);
      std::copy(STYLE_temlate, insert_pos, it);
    }
    SVG_str << std::ends;
    svgf << SVG_str.str();
    svgf << "</svg>\n" << insert_pos;
    svgf.close();
  }
  delete_test_collection(seg_type, seg_coll, seg_ptr_coll);

  if(wait){printf("\npress 'Enter' to continue"); getchar();}
  return 0;
  }


