#pragma once
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

#include "test_coll.h"

inline constexpr const char* help_msg = R"WYX(
example: seg_int -a14 -sa -dp -n20000 -p5.5
-aA: type of algorithms tested
 A=1: trivial
 A=2: simple sweep
 A=4: balaban fast
 A=8: balaban optimal
 A=16: balaban fast parallel (threads via -T, default 8)
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
-eE: prints intersection find time compare to trivial intersection check time (ICT) for all tested algorithms 
  E - ICT estimation in nanoseconds (you can take it from small segment sets)
  if E not presented (just -e), program calculates ICT from trivial algorithm run. In the later 
  case trivial algorithm must be selected for testing and for the large N it can take quite a time.
-SR: capital S for random seed; R - random seed value; if R=0 - non pseudo random generator used
-Tn: number of threads for the parallel algorithm (A=16). Default 8.
-fhtmfile: if presented, program writes SVG picture to htmfile. For examle -fC:/tmp/res.htm
  To limit resulting file size option works only for 5000 segments and less, also only first 
  150000 intersections are drawn.
-m: print truncated information in one row (to make a table)
-w: stop and wait for an input befor exit
-An: bind current thread to logical CPU n (0-based). -aff-1 (default) disables affinity.
-Pp: thread priority preset:
 p=0: don't change priorities (default)
 p=1: ABOVE_NORMAL process + THREAD_PRIORITY_ABOVE_NORMAL
 p=2: HIGH process + THREAD_PRIORITY_HIGHEST
 p=3: REALTIME process + THREAD_PRIORITY_TIME_CRITICAL
)WYX";

struct Options : public SegmentsAndRegOptions
{
  uint4 alg = 31;

  BOOL print_less = FALSE;
  BOOL wait = FALSE;
  BOOL rtime_printout = FALSE;
  bool use_counters = false;

  const char* fname = nullptr;

  double ICT = -1;
  unsigned random_seed = 317;

  long long affinity_cpu = -1;
  int pri_preset = 0;

  bool ParseArgs(int argc, char* argv[])
  {
    if (argc == 1) {
      printf("usage: seg_int -aA -sS -SR -dD -nN -pP -rR -m -eE -w -SN -Tn -fhtmfile -An -Pp\n");
      PrintHelp();
      return false;
    }
    for (int i = 1; i < argc; i++)
      if (argv[i][0] == '-') {
        switch (argv[i][1])
        {
        case 'a':
        {
          alg = atoi(argv[i] + 2);
          if ((alg < 1) || (alg > 31)){
            alg = 4; printf("some error in -a param. 4 used instead.\n");
          }
        }
        break;
        case 's':
        {
          switch (argv[i][2])
          {
            case 'L':seg_type = _Segment::line1; break;
            case 'l':seg_type = _Segment::line2; break;
            case 'a':seg_type = _Segment::arc; break;
            case 'g':seg_type = _Segment::graph; break;
            case 'i': {
              seg_type = _Segment::intline;
              range_for_int_seg = atoi(argv[i] + 3);
              if (range_for_int_seg == 0) range_for_int_seg = full_int_range;
              if (range_for_int_seg < 5 || range_for_int_seg > full_int_range) {
                printf("some error in -si param. full_int_range used as a range.\n");
                range_for_int_seg = full_int_range;
              };
            };
            break;
            default:
            {
              printf("some error in -s param. l used instead.\n");
              seg_type = _Segment::line2;
            }
          };
        };
        break;
        case 'd':
        {
            switch (argv[i][2])
            {
            case 'r':distr_type = _Distribution::random; break;
            case 'l':distr_type = _Distribution::parallel; break;
            case 'm':distr_type = _Distribution::mixed; break;
            case 's':distr_type = _Distribution::small; break;
            case 'p':distr_type = _Distribution::param_defined; break;
            case 'c':distr_type = _Distribution::circle; break;
            default:
            {
              printf("some error in -distr_type param. r used instead.\n");
              distr_type = _Distribution::random;
            }
          };
        };
        break;
        case 'r':
        {
          switch (argv[i][2])
          {
          case 'c':reg_stat = _Registrator::just_count; break;
          case 'C':reg_stat = _Registrator::just_count_require_intersections; break;
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
          if ((n < 2))      {
            n = 10000; printf("some error in -n param. 10000 used instead.\n");
          }
        }
        break;
        case 'S':
        {
          random_seed = atoi(argv[i] + 2);
        }
        break;
        case 'T':
        {
          auto n_threads = atoi(argv[i] + 2);
          if (n_threads < 1) {
            n_threads = n_threads_default;
            printf("some error in -T param. %u used instead.\n", n_threads_default);
          }
          else
            n_threads = static_cast<uint4>(n_threads);
        }
        break;
        case 'p':
        {
          distr_param = fabs(atof(argv[i] + 2));
          if (distr_param <= 0) {
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
          fname = argv[i] + 2;
        }
        break;
        case 'A':
        {
          affinity_cpu = _strtoi64(argv[i] + 2, nullptr, 10);
        }
        break;
        case 'P':
        {
          pri_preset = atoi(argv[i] + 2);
          if (pri_preset < 0) pri_preset = 0;
          if (pri_preset > 3) pri_preset = 3;
        }
        break;
        }
      }
      return true;
  }

  void CheckOptions()
  {
    if ((seg_type == _Segment::graph) && (distr_type != _Distribution::random)) {
      printf("-sg  is compartible only with -dr!\n");
      distr_type = _Distribution::random;
    }
    if ((reg_stat == _Registrator::store_pairs_and_ints_just_count_stat) && (n > max_truereg_items)) {
      printf("Too many segments (and possible intersections!) for true registration! Just counting used instead.\n");
      reg_stat = _Registrator::just_count;
    }
    if ((fname) && (n > max_SVG_items)) {
      printf("Too many segments for SVG output, file will not be created.\n");
      fname = nullptr;
    }
  }

  void PrintHelp() const
  {
    printf("%s", help_msg);
  }

  void PrintActual() const
  {
    const char* ss = "Lligi", * sd = "rlmspc", * sr = "pPcrC";
    if (!print_less) {
      if (seg_type != _Segment::intline)
        printf("actual params are: -a%i -s%c -d%c -r%c -n%i -S%i -T%u -p%f -A%lld -P%i\n", alg, ss[seg_type], sd[distr_type], sr[reg_stat], n, random_seed, n_threads, distr_param, affinity_cpu, pri_preset);
      else
        printf("actual params are: -a%i -si%i -d%c -r%c -n%i -S%i -T%u -p%f -e%f -A%lld -P%i\n", alg, range_for_int_seg, sd[distr_type], sr[reg_stat], n, random_seed, n_threads, distr_param, ICT, affinity_cpu, pri_preset);
    }
  }
};
