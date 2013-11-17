/*
*
*      Copyright (c)  2011-2013  Ivan Balaban 
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>


//#include <tchar.h>

#ifdef _DEBUG
#include <crtdbg.h>
#endif



double _benchmark(char *counters_string,int4 n,PSeg *coll,int4 s_type,int4 alg,BOOL dont_need_ip,double &res)
  {
  double start,stop;
  double c[8];
  int4 n_call=0;
  start=clock();
  do
  {
  res=find_intersections(s_type,n,coll,alg,dont_need_ip,c);
  stop=clock();
  n_call++;
  }
  while ((stop-start)<CLOCKS_PER_SEC);
  if(counters_string)
   sprintf(counters_string,"%11.2f;%11.2f;%11.2f;%11.2f;%11.2f;%11.2f;%11.2f;%11.2f",c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);
   //printf("avg Size %11.2f\n",c[3]/(c[4]+1));
  return (stop-start)/(n_call*CLOCKS_PER_SEC);
  }



int main(int argc, char* argv[])
  {
  int4 alg=63, seg_type=2, n = 20000, d=0;
  REAL p=1.0;
  double exec_time[33],nInt[33];
  char *ss="Lla",*sd="rlmsp";
  BOOL mute=FALSE,print_counters=FALSE,wait=FALSE,dont_need_ip=FALSE;
  char counters_string[256],*pcs=NULL,*counters_mute,rpar[]="-r";
  int4 alg_list[]={ triv, simple_sweep, fast, optimal, fast_parallel, bentley_ottmann};
  char *alg_names[]={"trivial","simple_sweep","fast","optimal","fast_parallel","bentley_ottmann"};

#ifdef _DEBUG
  _CrtSetDbgFlag(_CRTDBG_CHECK_ALWAYS_DF);
#endif

#ifdef COUNTERS_ON
strcpy(counters_string,"\ncounters:");
counters_mute=counters_string+strlen(counters_string);
#else
counters_string[0]=0;
counters_mute=counters_string;
#endif

  

  if (argc==1)
    {
#ifdef COUNTERS_ON
    printf("usage: seg_int -aA -sS -dD -nN -pP -r -m -w -c\n");
#else    
    printf("usage: seg_int -aA -sS -dD -nN -pP -r -m -w\n");
#endif    
    printf("example: seg_int -a14 -sa -dp -n20000 -p5.5\n");
    printf("-aA: type of algorithms tested\n");
    printf(" A=1: trivial\n");
    printf(" A=2: simple sweep\n");
    printf(" A=4: balaban fast\n");
    printf(" A=8: balaban optimal\n");
    printf(" A=16: balaban fast parallel for 2 cores\n");
    printf(" A=32: bentley & ottmann\n");
    printf(" if you want several algorithms to test just sum up their values\n");
    printf(" i.e. A=63 (=1+2+4+8+16+32) all algorithms\n");
    printf("-sS: type of segments\n");
    printf(" S=l: line segments representation y=a*x+b,x1<=x<=x2; a,b,x1,x2 - reals\n");
    printf(" S=L: line segments representation r=b+a*t,0<=t<=1; b,a - vectors\n");
    printf(" S=a: arcs\n");
    printf("-nN: number of segments\n");
    printf("-pP: parameter value, must be positive\n");
    printf("-dD: type of distribution\n");
    printf(" D=r: random segments\n");
    printf(" D=l: random length almost x parallel segments, the bigger p/N the less parallel segments\n");
    printf(" D=m: mixed random length almost x parallel 'long' segments (33%) and 'small' segments (67%), the bigger p/N the less parallel 'long' and longer 'short' segments\n");
    printf(" D=s: short segments: random segment with  length multiplied by p/N\n");
    printf(" D=p: random segment with  length multiplied by p\n");
    printf("-r: if presented, fast algorithm doesn't find intersection points but only intersecting pairs\n");
    printf("-m: if presented, means 'mute' mode - few information printed\n");
    printf("-w: if presented, program wait until 'Enter' pressed before closing\n");
#ifdef COUNTERS_ON
    printf("-c: counters are printed, if presented\n");
#endif    
    printf("Important! -sa is not compartible with -dm, -dl and -r options!\n");
    return 0;
    }
  else
    {
    for(int i =1;i<argc;i++)
      if (argv[i][0]=='-') 
        {
        switch(argv[i][1])
          {
          case 'a':
            {
            alg=atoi(argv[i]+2);
            if((alg<1)||(alg>63))
              {alg=4; printf("some error in -a param. 4 used instead.\n");}
            }
            break;
          case 's':
            {
            switch(argv[i][2])
              {
              case 'L':seg_type=line1;break;
              case 'l':seg_type=line2;break;
              case 'a':seg_type=arc;break;
              default:
                {
                printf("some error in -s param. l used instead.\n");
                seg_type=line2; 
                }
              };
            };
            break;
          case 'd':
            {
            switch(argv[i][2])
              {
              case 'r':d=random;break;
              case 'l':d=parallel;break;
              case 'm':d=mixed;break;
              case 's':d=small;break;
              case 'p':d=param_defined;break;
              default:
                {
                printf("some error in -d param. r used instead.\n");
                d=random; 
                }
              };
            };
            break;
          case 'n':
            {
            n=atoi(argv[i]+2);
            if((n<2))
              {n=10000; printf("some error in -n param. 10000 used instead.\n");}
            }
            break;
          case 'p':
            {
            p=fabs(atof(argv[i]+2));
            if(p<=0)
              {p=1.0; printf("some error in -p param. 1.0 used instead.\n");}
            }
            break;
          case 'r':
            {
            dont_need_ip=TRUE;
            }
            break;
          case 'm':
            {
            mute=TRUE;
            }
            break;
          case 'w':
            {
            wait=TRUE;
            }
            break;
#ifdef COUNTERS_ON
          case 'c':
            {
              pcs=counters_mute;
            }
            break;
#endif
          }
        }
    }
  if(!dont_need_ip)rpar[0]=0;
  if(!mute)printf("actual params is: -a%i -s%c -d%c -n%i -p%f %s\n", alg, ss[seg_type], sd[d], n, p,rpar);
  if((seg_type==arc)&&(d==mixed)) {printf("-sa -dm is not compartible!/n"); return 0;}
  if((seg_type==arc)&&(d==parallel)) {printf("-sa -dl is not compartible!/n"); return 0;}
  PSeg *coll=create_test_collection(seg_type,n,d,p);
  for(int4 a=0;a<sizeof(alg_list)/sizeof(alg_list[0]);a++)
    {
      if(alg&alg_list[a])
        {
        exec_time[a]=_benchmark(pcs,n,coll,seg_type,alg_list[a],dont_need_ip,nInt[a]);
        if(mute)
          printf("a%i;s%c;d%c;r%i;n%i;i%13.0f;t%6.5f;p%f;%s\n",alg_list[a],ss[seg_type], sd[d],dont_need_ip,n,nInt[a],exec_time[a],p,counters_mute);
        else
          printf("%s intersections=%13.0f time=%6.5f%s\n",alg_names[a],nInt[a],exec_time[a],counters_string);
        }
     };
  //printf("ratio fast  =%6.3f\n",0.5*(n*exec_time[2]*(n-1))/(exec_time[0]*nInt[2]));

  delete_test_collection(seg_type,coll);
  if(wait){printf("press 'Enter' to continue"); getchar();}
  return 0;
  }


