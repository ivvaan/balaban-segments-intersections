/*																 fou
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
#include "test_coll.h"
#include <stdio.h>
#include "registrator.h"

#include "segments.h"
#include "lines1_collection.h"
#include "lines2_collection.h"
#include "arcs_collection.h"
#include "graph_collection.h"

#include "new_int_finder.h"
#include "fast_finder.h"
#include "mem_save_finder.h"
#include "optimal_finder.h"
#include "intersection_finder.h"


constexpr int4 n_threads = 6;


// creating and deleting collections
PSeg first_segment_ptr = NULL;
template <class T>
minmaxrect get_minmax(int4 n, T c[])
{
  TPlaneVect bp = c[0].BegPoint();
  TPlaneVect ep = c[0].EndPoint();
  REAL xmin = bp.x;
  REAL xmax = ep.x;
  REAL ymin = MIN(bp.y, ep.y);
  REAL ymax = MAX(bp.y, ep.y);
  for (int4 i = 1; i < n; ++i) {
    bp = c[i].BegPoint();
    ep = c[i].EndPoint();
    xmin = MIN(xmin, bp.x);
    xmax = MAX(xmax, ep.x);
    ymin = std::min({ ymin,bp.y,ep.y });
    ymax = std::max({ ymax,bp.y,ep.y });
  }

  return { {xmin,ymin},{xmax,ymax} };
};



PSeg create_test_collection(int4 seg_type, int4 n, int4 distr, REAL par, CRandomValueGen& random_gen, PSeg** seg_ptr_coll_ptr)
{
  PSeg* colls = nullptr;
  int4 i;
  constexpr int4 sentinels = 2;
  if (seg_ptr_coll_ptr != nullptr)colls = new PSeg[n];
  if (distr != param_defined)par /= n / 33.0;
  if (distr == parallel)par /= 2.4;
  PSeg result;
  switch (seg_type)
  {
  case line1:
  {
    TLineSegment1* Colls = new TLineSegment1[n + sentinels];
    if (seg_ptr_coll_ptr != nullptr)for (i = 0; i < n; i++)colls[i] = Colls + i;
    for (i = 0; i < n; i++)
      Colls[i].InitRandom(random_gen, i, distr, par);
    auto mmr = get_minmax(n, Colls).get_scaled(2.0);
    Colls[n] = TLineSegment1(mmr.ld, { mmr.rt.x,mmr.ld.y });
    Colls[n + 1] = TLineSegment1({ mmr.ld.x,mmr.rt.y }, mmr.rt);
    result = Colls;
    first_segment_ptr = Colls;
  }; break;
  case line2:
  {
    TLineSegment2* Colls = new TLineSegment2[n + sentinels];
    if (seg_ptr_coll_ptr != nullptr)for (i = 0; i < n; i++)colls[i] = Colls + i;
    for (i = 0; i < n; i++)
      Colls[i].InitRandom(random_gen, i, distr, par);
    auto mmr = get_minmax(n, Colls).get_scaled(2.0);
    Colls[n].Init({ mmr.ld, { mmr.rt.x,mmr.ld.y } });
    Colls[n + 1].Init({ { mmr.ld.x,mmr.rt.y }, mmr.rt });
    result = Colls;
    first_segment_ptr = Colls;
  }; break;
  case arc:
  {
    TArcSegment* Colls = new TArcSegment[n];
    if (seg_ptr_coll_ptr != nullptr)for (i = 0; i < n; i++)colls[i] = Colls + i;
    for (i = 0; i < n; i++)
      Colls[i].InitRandom(random_gen, i, distr, par);
    result = Colls;
    first_segment_ptr = Colls;
  }; break;
  case graph:
  {
    TPlaneVect* Colls = new TPlaneVect[n];
    for (i = 0; i < n; i++)
      Colls[i].InitRandom(random_gen);
    result = Colls;
    first_segment_ptr = Colls;
  }; break;

  }

  if (seg_ptr_coll_ptr != nullptr)*seg_ptr_coll_ptr = colls;
  return result;
};


void  delete_test_collection(int4 seg_type, PSeg seg_coll, PSeg* seg_ptr_coll)
{
  switch (seg_type)
  {
  case line1:delete[](TLineSegment1*)seg_coll; break;
  case line2:delete[](TLineSegment2*)seg_coll; break;
  case arc:delete[](TArcSegment*)seg_coll; break;
  }
  if (seg_ptr_coll)delete[]  seg_ptr_coll;
};


//for old implementation

template <class SEGMENT>
class SegmentFunctions
  {
  public:
    static void register_intersection(PRegObj r,PSeg s1,PSeg s2,int4 nInt,TPlaneVect *points)
      {
#ifdef PRINT_SEG_AND_INT
      int4 n1=(SEGMENT*)s1-(SEGMENT*)first_segment_ptr,n2=(SEGMENT*)s2-(SEGMENT*)first_segment_ptr;
      if (n1<n2)
		  printf("found int %i, %i\n",n1,n2);
      else
      printf("found int %i, %i\n", n2, n1);
#endif
      (*(double*)r) += nInt;
      };

    static REAL __YAtX(PSeg s,REAL X)
      {
        return ((SEGMENT *)s)->YAtX(X);
      };

    static bool __Below(REAL x,PSeg s1,PSeg s2)
      {
      return below(x,(SEGMENT*)s1,(SEGMENT*)s2);
      };
    template<bool  _reg_ip>
    static  int4 __FindAndRegIPoints(PSeg s1,PSeg s2,PRegObj intersection_registrator)
      {
      TPlaneVect P[2];                     
      int4 n=IntPoint<_reg_ip>((SEGMENT*)s1,(SEGMENT*)s2,P);
      if(n)
		  register_intersection(intersection_registrator,s1,s2,n, _reg_ip?P:NULL);
      return n; 
      };

    static  int4 __FindAndRegIPointOnRightOfSWL(REAL swl_x,PSeg s1,PSeg s2,TPlaneVect *p1,TPlaneVect *p2)
      {
          TPlaneVect p[2];
          int4 n=IntPoint<true>((SEGMENT*)s1,(SEGMENT*)s2,p);
//              if(n&&(p->x>swl_x)) *p1=*p;
              
          if(n==0)return 0;
          TPlaneVect *pt=p1;
          int4 res=0;
          int4 f=0,l=1;
          if((n==2)&&(p[0].x>p[1].x)){f=1;l=0;};
          if(p[f].x>swl_x)
            {
                *pt=p[f];
                pt=p2;
                res++;
            }
          if((n==2)&&(p[l].x>swl_x))
            {
                *pt=p[l];
                res++;
            } 
            return res;   
      };
    
    template<bool  _reg_ip>
    static int4 __FindAndRegIPointsInStripe(REAL b,REAL e, PSeg s1,PSeg s2,PRegObj intersection_registrator)
      {
      TPlaneVect P[2];
      int4 n=StripePoint<_reg_ip>(b,e,(SEGMENT*)s1,(SEGMENT*)s2,P);
      if(n)register_intersection(intersection_registrator,s1,s2,n, _reg_ip ? P : NULL);
      return n; 
      };
    static void __BegPoint(PSeg s,REAL &x,REAL &y)
      {
      ((SEGMENT*)s)->BegPoint(x,y);
      };
    static  void __EndPoint(PSeg s,REAL &x,REAL &y)
      {
      ((SEGMENT*)s)->EndPoint(x,y);
      };
    static  int4 __Under(PSeg s1,PSeg s2)
      {
      return ((SEGMENT*)s1)->under(((SEGMENT*)s2)->BegPoint());
      };

    static int4 __IntInside(REAL b,REAL e,PSeg s1,PSeg s2)
      {
      return IntInside(b,e,(SEGMENT*)s1,(SEGMENT*)s2);
      };  

  };


const int4 reg_obj_margin = 64;// for reg objects to be in different CPU cash blocks
double reg_objects[reg_obj_margin*n_threads];

double find_intersections(int4 seg_type, int4 SN, PSeg* colls, int4 alg, double* counters, bool dont_need_ip )
  {
  if(colls==NULL) return 0;
  double int_numb=0;
  BOOL is_line;
  if ((seg_type == line1) || (seg_type == line2))
	  {
		  CLineFinder intersection_finder;
		  switch (seg_type)
			  {
			  case line1:
				  intersection_finder.set_segm_fuctions(
					  SegmentFunctions<TLineSegment1>::__Below,
					  dont_need_ip ? SegmentFunctions<TLineSegment1>::__FindAndRegIPoints<false> : SegmentFunctions<TLineSegment1>::__FindAndRegIPoints<true>,
					  dont_need_ip ? SegmentFunctions<TLineSegment1>::__FindAndRegIPointsInStripe<false> : SegmentFunctions<TLineSegment1>::__FindAndRegIPointsInStripe<true>,
					  SegmentFunctions<TLineSegment1>::__IntInside,
					  SegmentFunctions<TLineSegment1>::__BegPoint,
					  SegmentFunctions<TLineSegment1>::__EndPoint,
					  SegmentFunctions<TLineSegment1>::__Under,
					  SegmentFunctions<TLineSegment1>::__FindAndRegIPointOnRightOfSWL,
					  SegmentFunctions<TLineSegment1>::__YAtX,
					  SegmentFunctions<TLineSegment1>::register_intersection,
					  &int_numb,
					  is_line = TLineSegment1::is_line
				  );
				  break;
			  case line2:
				  intersection_finder.set_segm_fuctions(
					  SegmentFunctions<TLineSegment2>::__Below,
					  dont_need_ip ? SegmentFunctions<TLineSegment2>::__FindAndRegIPoints<false> : SegmentFunctions<TLineSegment2>::__FindAndRegIPoints<true>,
					  dont_need_ip ? SegmentFunctions<TLineSegment2>::__FindAndRegIPointsInStripe<false> : SegmentFunctions<TLineSegment2>::__FindAndRegIPointsInStripe<true>,
					  SegmentFunctions<TLineSegment2>::__IntInside,
					  SegmentFunctions<TLineSegment2>::__BegPoint,
					  SegmentFunctions<TLineSegment2>::__EndPoint,
					  SegmentFunctions<TLineSegment2>::__Under,
					  SegmentFunctions<TLineSegment2>::__FindAndRegIPointOnRightOfSWL,
					  SegmentFunctions<TLineSegment2>::__YAtX,
					  SegmentFunctions<TLineSegment2>::register_intersection,
					  &int_numb,
					  is_line = TLineSegment2::is_line
				  );
				  break;
			  }
		  switch (alg)
			  {
				  case triv:intersection_finder.trivial(SN, colls); break;
				  case simple_sweep:intersection_finder.simple_sweep(SN, colls); break;
				  case fast: intersection_finder.balaban_fast(SN, colls); break;
				  case mem_save: intersection_finder.balaban_memory_save(SN, colls); break;
				  case fast_no_ip: intersection_finder.balaban_no_ip(SN, colls); break;
				  case optimal:intersection_finder.balaban_optimal(SN, colls); break;
				  case fast_parallel:
					  {
						  PRegObj additional_reg_obj[n_threads];
						  for (int i = 0; i < n_threads - 1; i++)
						  {
							  reg_objects[i*reg_obj_margin] = 0;
							  additional_reg_obj[i] = (PRegObj *)reg_objects + i*reg_obj_margin;
						  }
						  intersection_finder.fast_parallel(SN, colls, n_threads, additional_reg_obj);
						  for (int i = 0; i < n_threads - 1; i++)int_numb += *(double*)additional_reg_obj[i];
					  }
				  break;
				  case bentley_ottmann:intersection_finder.bentley_ottmann(SN, colls); break;
			  };
		  memcpy(counters, intersection_finder.my_counter, sizeof(intersection_finder.my_counter));
      }

  if (seg_type == arc)
  {
	  CCommonFinder intersection_finder;
	  intersection_finder.set_segm_fuctions(
		  SegmentFunctions<TArcSegment>::__Below,
		  SegmentFunctions<TArcSegment>::__FindAndRegIPoints<true>,
		  SegmentFunctions<TArcSegment>::__FindAndRegIPointsInStripe<true>,
		  SegmentFunctions<TArcSegment>::__IntInside,
		  SegmentFunctions<TArcSegment>::__BegPoint,
		  SegmentFunctions<TArcSegment>::__EndPoint,
		  SegmentFunctions<TArcSegment>::__Under,
		  SegmentFunctions<TArcSegment>::__FindAndRegIPointOnRightOfSWL,
		  SegmentFunctions<TArcSegment>::__YAtX,
		  SegmentFunctions<TArcSegment>::register_intersection,
		  &int_numb,
		  is_line = TArcSegment::is_line
	  );
	  switch (alg)
	  {
		  case triv:intersection_finder.trivial(SN, colls); break;
		  case simple_sweep:intersection_finder.simple_sweep(SN, colls); break;
		  case fast: intersection_finder.balaban_fast(SN, colls); break;
      case mem_save:intersection_finder.balaban_memory_save(SN, colls);break;
		  case fast_no_ip:  break;//incompatible
		  case optimal:intersection_finder.balaban_optimal(SN, colls); break;
		  case fast_parallel:
		  {
			  PRegObj additional_reg_obj[n_threads];
			  for (int i = 0; i < n_threads - 1; i++)
			  {
				  reg_objects[i*reg_obj_margin] = 0;
				  additional_reg_obj[i] = (PRegObj *)reg_objects + i*reg_obj_margin;
			  }
			  intersection_finder.fast_parallel(SN, colls, n_threads, additional_reg_obj);
			  for (int i = 0; i < n_threads - 1; i++)int_numb += *(double*)additional_reg_obj[i];
		  }
		  break;
		  case bentley_ottmann:intersection_finder.bentley_ottmann(SN, colls); break;
	  }

	  memcpy(counters, intersection_finder.my_counter, sizeof(intersection_finder.my_counter));
  }

  return int_numb;
  };


  template<template<class>class SegColl,class Counter>
  double find_int( int4 n, SegColl<Counter> &coll,int4 alg,uint4 stat)
  {
    switch (alg) {
        case triv: {
            CTrivialIntFinder fi;
            fi.find_intersections(coll);
        } break;
        case simple_sweep:{
            CSimpleSweepIntFinder fi;
            fi.find_intersections(coll);
        };
            break;
        case fast: {
            CFastIntFinder fi;
            fi.prepare_ends(coll);
            fi.find_intersections(coll);
        } break;
        case mem_save:{
            CMemSaveIntFinder fi;
            fi.prepare_ends(coll);
            fi.find_intersections(coll);
        }
            break;
        case optimal: {
          COptimalIntFinder fi;
          fi.prepare_ends(coll);
          fi.find_intersections(coll);
        }
            break;
        case fast_parallel: {
          constexpr uint4 cache_line_size = 64;
          constexpr uint4 reg_margin = 2 + cache_line_size / sizeof(Counter);// for reg objects to be in different CPU cache lines
          Counter reg_objects[reg_margin * n_threads];
          Counter* additional_reg_obj[n_threads];
            for (int i = 0; i < n_threads - 1; i++) {
                additional_reg_obj[i] = reg_objects + i * reg_margin;
            }
            CFastIntFinder fi;
            fi.prepare_ends(coll);
            fi.find_intersections(coll, n_threads, additional_reg_obj);
            coll.GetRegistrator()->combine_reg_data(n_threads, additional_reg_obj);

        } break;

    }

    return coll.GetRegistrator()->get_stat(stat);
  };

  template<class Counter>
  double _find_int(int4 seg_type, int4 n, PSeg segs, int4 alg,uint4 stat)
  {
    Counter reg;
    switch (seg_type) {
      case line1: {
        CLine1SegmentCollection<Counter> coll(n, segs, &reg);
        return find_int(n, coll, alg, stat);
      };
      case line2: {
        CLine2SegmentCollection<Counter> coll(n, segs, &reg);
        return find_int(n, coll, alg, stat); 
      };
      case arc:{
        CArcSegmentCollection<Counter>  coll(n, segs, &reg);
        return find_int(n, coll, alg, stat);
      };
      case graph: {
        CGraphSegmentCollection<Counter> coll(n, segs, &reg);
        return find_int(n, coll, alg, stat);
      };
    }
    return 0;
  };

  void write_SVG(chostream* SVG_stream, int4 seg_type, int4 n, PSeg segs, int4 algs, uint4 stat)
  {
    auto apply_algs = [=](auto& coll) {
      coll.coll_to_SVG(SVG_stream);
      for (int4 a = sizeof(alg_list) / sizeof(alg_list[0]) - 1; a > -1; --a)
        if (algs & alg_list[a]) {
          auto alg = alg_list[a];
          TrueRegistrator reg;
          coll.SetRegistrator(&reg);
          find_int(n, coll, alg, stat);
          reg.write_SVG(alg, SVG_stream);
        }
    };

    switch (seg_type) {
      case line1: {
        CLine1SegmentCollection<TrueRegistrator> coll(n, segs, nullptr);
        apply_algs(coll);
      }; break;
      case line2: {
        CLine2SegmentCollection<TrueRegistrator> coll(n, segs, nullptr);
        apply_algs(coll);
      }; break;
      case arc: {
        CArcSegmentCollection<TrueRegistrator>  coll(n, segs, nullptr);
        apply_algs(coll);
      }; break;
      case graph: {
        CGraphSegmentCollection<TrueRegistrator> coll(n, segs, nullptr);
        apply_algs(coll);
      }; break;
    }
  };


  find_intersections_func get_find_intersections_func(uint4 reg_type)
  {
    if (reg_type == _Registrator::just_count)return _find_int<SimpleCounter>;
    if (reg_type == _Registrator::store_pairs_and_ints_just_count_stat)
      return _find_int<TrueRegistrator>;
    return _find_int<PerSegmCounter>;
  };


 