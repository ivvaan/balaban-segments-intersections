// robust_seg_int.cpp : Defines the entry point for the console application.
//


#include "seg_coll.h"
#include "robust_int_finder.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

double _benchmark(CRobustSegmentCollection *scoll,int4 alg,double &res)
  {
  double start,stop;
  int4 n_call=0;
  CIntersectionFinder int_finder;
  start=clock();
  switch(alg)
  {
  case 0:int_finder.trivial(scoll);break;
  case 1:int_finder.simple_sweep(scoll);break;
  case 2:int_finder.balaban_fast(scoll);break;
  };
  res=get_int_numb_and_reset(scoll);
  stop=clock();
  return (stop-start)/(CLOCKS_PER_SEC);
  }


int main(int argc, char* argv[])
{
 int4 nv=2000;
 int4 ne=2*nv;
 int4 alg=7;
 double exec_time[33],nInt[33];
 bool wait=FALSE;
 int4 alg_list[]={ 1, 2, 4};
 char *alg_names[]={"trivial","simple_sweep","fast"};

   if (argc==1)
    {

    printf("usage: robust_seg_int -aA -vV -eE -w\n");
    printf("example: robust_seg_int -a7 -v2000 -e4000\n");
    printf("-aA: type of algorithms tested\n");
    printf(" A=1: trivial\n");
    printf(" A=2: simple sweep\n");
    printf(" A=4: balaban fast\n");
    printf(" if you want several algorithms to test just sum up their values\n");
    printf(" i.e. A=7 (=1+2+4) all algorithms\n");
    printf("-vV: number of randomly distributed vertices\n");
    printf("-eE: number of edges randomly connects the vertices\n");
    printf("-w: if presented, program wait until 'Enter' pressed before closing\n");
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
            if((alg<1)||(alg>7))
              {alg=4; printf("some error in -a param. 4 used instead.\n");}
            }
            break;
          case 'v':
            {
            nv=atoi(argv[i]+2);
            if((nv<2))
              {nv=2000; printf("some error in -v param. 2000 used instead.\n");}
            }
            break;
          case 'e':
            {
            ne=atoi(argv[i]+2);
            if(ne<2)
              {nv=4000; printf("some error in -e param. 4000 used instead.\n");}
            }
            break;
          case 'w':
            {
            wait=true;
            }
            break;
          }
        }
    }

 CRobustSegmentCollection *scoll=NULL;
 
 scoll=create_segment_collection (nv, ne);

   for(int4 a=0;a<sizeof(alg_list)/sizeof(alg_list[0]);a++)
    {
      if(alg&alg_list[a])
        {
        exec_time[a]=_benchmark(scoll,a,nInt[a]);
          printf("%i vertices %i edges; %s int=%12.0f time=%6.4f%s\n",nv,ne,alg_names[a],nInt[a],exec_time[a]);
        }
     };
 delete_segment_collection (scoll);
 if(wait){printf("press 'Enter' to continue"); getchar();}
	return 0;
}

