/*
*
*      Copyright (c)  2011  Ivan Balaban 
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



double _benchmark(char *counters_string,int4 n,PSeg *coll,int4 s_type,int4 alg,double &res)
  {
  double start,stop;
  double c[8];
  start=clock();
  res=find_intersections(s_type,n,coll,alg,c);
  stop=clock();
  if(counters_string)
   sprintf(counters_string,"%13.0f;%13.0f;%13.0f;%13.0f;%13.0f;%13.0f;%13.0f;%13.0f",c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);
  return (stop-start)/CLOCKS_PER_SEC;
  }



int main(int argc, char* argv[])
  {
  int4 alg=63, seg_type=2, n = 20000, d=0;
  REAL p=1.0;
  double exec_time[33],nInt[33];
  char *ss="Lla",*sd="rmsp";
  BOOL mute=FALSE,print_counters=FALSE;
  char counters_string[256],*pcs=NULL,*counters_mute;

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
    printf("usage: seg_int -aA -sS -dD -nN -pP -m -c\n");
    printf("example: seg_int -t14 -sa -dp -n20000 -p5.5\n");
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
    printf("-dD: type of distribution\n");
    printf(" D=r: random segments\n");
    printf(" D=m: mixed long and short parallel segments\n");
    printf(" D=s: short segments\n");
    printf(" D=p: random segment with parameter defined length\n");
    printf("Important! -sa -dm is not compartible!\n");
    printf("-nN: number of segments\n");
    printf("-pP: parameter value, must be positive\n");
    printf("-m: if presented means 'mute' mode - few information printed\n");
    printf("-c: counters are printed if presented\n");

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
                printf("some error in -s param. L used instead.\n");
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
              case 'm':d=mixed;break;
              case 's':d=small;break;
              case 'p':d=param_defined;break;
              default:
                {
                printf("some error in -d param. r used instead.\n");
                seg_type=random; 
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
            if(p==0)
              {p=1.0; printf("some error in -p param. 1.0 used instead.\n");}
            }
            break;
          case 'm':
            {
            mute=TRUE;
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
  if(!mute)printf("actual param string is: -a%i -s%c -d%c -n%i -p%f\n", alg, ss[seg_type], sd[d], n, p);
  if((seg_type==arc)&&(d==mixed)) {printf("-sa -dm is not compartible!/n"); return 0;}
  
  PSeg *coll=create_test_collection(seg_type,n,d,p);
   
  if(alg&fast)
    {
    exec_time[fast]=_benchmark(pcs,n,coll,seg_type,fast,nInt[fast]);
    if(mute)
      printf("a4;s%c;d%c;n%i;i%13.0f;t%8.2f;p%f;%s\n",ss[seg_type], sd[d],n,nInt[fast],exec_time[fast],p,counters_mute);
    else
      printf("fast   intersections=%13.0f time=%8.2f%s\n",nInt[fast],exec_time[fast],counters_string);
    }

  if(alg&optimal)
    {
    exec_time[optimal]=_benchmark(pcs,n,coll,seg_type,optimal,nInt[optimal]);
    if(mute)
      printf("a8;s%c;d%c;n%i;i%13.0f;t%8.2f;p%f;%s\n",ss[seg_type], sd[d],n,nInt[optimal],exec_time[optimal],p,counters_mute);
    else
      printf("optimal  intersections=%13.0f time=%8.2f%s\n",nInt[optimal],exec_time[optimal],counters_string);
    }

  if(alg&fast_parallel)
    {
    exec_time[fast_parallel]=_benchmark(pcs,n,coll,seg_type,fast_parallel,nInt[fast_parallel]);
    if(mute)
      printf("a16;s%c;d%c;n%i;i%13.0f;t%8.2f;p%f;%s\n",ss[seg_type], sd[d],n,nInt[fast_parallel],exec_time[fast_parallel],p,counters_mute);
    else
      printf("fast_parallel intersections=%13.0f time=%8.2f%s\n",nInt[fast_parallel],exec_time[fast_parallel],counters_string);
    }

  if(alg&bentley_ottmann)
    {  
    exec_time[bentley_ottmann]=_benchmark(pcs,n,coll,seg_type,bentley_ottmann,nInt[bentley_ottmann]);
    if(mute)
      printf("a32;s%c;d%c;n%i;i%13.0f;t%8.2f;p%f;%s\n",ss[seg_type], sd[d],n,nInt[bentley_ottmann],exec_time[bentley_ottmann],p,counters_mute);
    else
      printf("bentley_ottmann intersections=%13.0f time=%8.2f%s\n",nInt[bentley_ottmann],exec_time[bentley_ottmann],counters_string);
    }


  if(alg&simple_sweep)
    {  
    exec_time[simple_sweep]=_benchmark(pcs,n,coll,seg_type,simple_sweep,nInt[simple_sweep]);
    if(mute)
      printf("a2;s%c;d%c;n%i;i%13.0f;t%8.2f;p%f;%s\n",ss[seg_type], sd[d],n,nInt[simple_sweep],exec_time[simple_sweep],p,counters_mute);
    else
      printf("simple_sweep intersections=%13.0f time=%8.2f%s\n",nInt[simple_sweep],exec_time[simple_sweep],counters_string);
    }

  if(alg&triv)
    {  
    exec_time[triv]=_benchmark(pcs,n,coll,seg_type,triv,nInt[triv]);
    if(mute)
      printf("a1;s%c;d%c;n%i;i%13.0f;t%8.2f;p%f;%s\n",ss[seg_type], sd[d],n,nInt[triv],exec_time[triv],p,counters_mute);
    else
      printf("triv       intersections=%13.0f time=%8.2f%s\n",nInt[triv],exec_time[triv],counters_string);
    }

  /*printf("ratio fast  =%6.3f\n",0.5*(n*exec_time[fast]*(n-1))/(exec_time[triv]*nInt[fast]));*/

  delete_test_collection(seg_type,coll);
  //getchar();
  return 0;
  }


