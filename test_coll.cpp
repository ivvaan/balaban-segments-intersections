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
#include "int_lines_coll.h"

#include "new_int_finder.h"
#include "fast_finder.h"
#include "mem_save_finder.h"
#include "optimal_finder.h"
#include "intersection_finder.h"


constexpr int4 n_threads = 6;


// creating and deleting collections
PSeg first_segment_ptr = NULL;



PSeg create_test_collection(int4 seg_type, int4 n, int4 distr, REAL par, CRandomValueGen& random_gen, PSeg** seg_ptr_coll_ptr)
{
  PSeg* colls = nullptr;
  int4 i;
  constexpr int4 sentinels = 2;
  if (seg_ptr_coll_ptr != nullptr)colls = new PSeg[n];
  if (distr != _Distribution::param_defined)par /= n / 33.0;
  if (distr == _Distribution::parallel)par /= 2.4;
  PSeg result=nullptr;
  switch (seg_type)
  {
  case _Segment::line1:
  case _Segment::intline:
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
  case _Segment::line2:
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
  case _Segment::arc:
  {
    TArcSegment* Colls = new TArcSegment[n];
    if (seg_ptr_coll_ptr != nullptr)for (i = 0; i < n; i++)colls[i] = Colls + i;
    for (i = 0; i < n; i++)
      Colls[i].InitRandom(random_gen, i, distr, par);
    result = Colls;
    first_segment_ptr = Colls;
  }; break;
  case _Segment::graph:
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
  case _Segment::line1:delete[](TLineSegment1*)seg_coll; break;
  case _Segment::line2:delete[](TLineSegment2*)seg_coll; break;
  case _Segment::arc:delete[](TArcSegment*)seg_coll; break;
  }
  if (seg_ptr_coll)delete[]  seg_ptr_coll;
};



template<template<class>class SegColl,class Counter>
  auto find_int( int4 n, SegColl<Counter> &coll,int4 alg,uint4 stat)
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
        case optimal: {
          COptimalIntFinder fi;
          fi.prepare_ends(coll);
          fi.find_intersections(coll);
        }
            break;
        case fast_parallel: {
          constexpr uint4 cache_line_size = std::hardware_constructive_interference_size;
          constexpr uint4 reg_margin = (cache_line_size + sizeof(Counter) - 1) / sizeof(Counter);// for reg objects to be in different CPU cache lines
          Counter reg_objects[reg_margin * n_threads];
          Counter* additional_reg_obj[n_threads];
            for (int i = 0; i < n_threads - 1; i++) {
                additional_reg_obj[i] = reg_objects + i * reg_margin;
                additional_reg_obj[i]->Alloc(n);
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
  auto _find_int(int4 seg_type, int4 n, PSeg segs, int4 alg,uint4 stat,int4 range)
  {
    Counter reg;
    reg.Alloc(n);
    switch (seg_type) {
      case _Segment::line1: {
        CLine1SegmentCollection<Counter> coll(n, segs, &reg);
        return find_int(n, coll, alg, stat);
      };
      case intline: {
        CIntegerSegmentCollection<Counter> coll(n, (TLineSegment1*) segs, &reg, range);
        return find_int(n, coll, alg, stat);
      };
      case line2: {
          CLine2SegmentCollection<Counter> coll(n, segs, &reg);
          return find_int(n, coll, alg, stat); 
      };
      case _Segment::arc:{
        CArcSegmentCollection<Counter>  coll(n, segs, &reg);
        return find_int(n, coll, alg, stat);
      };
      case _Segment::graph: {
        CGraphSegmentCollection<Counter> coll(n, segs, &reg);
        return find_int(n, coll, alg, stat);
      };
    }
    return 0ULL;
  };

  void write_SVG(std::ostream* svg_stream, int4 seg_type, int4 n, PSeg segs, int4 algs, uint4 stat,int4 range)
  {
    if (!svg_stream)
      return;

    auto apply_algs = [=](auto& coll) {
      coll.coll_to_SVG(svg_stream);
      for (int4 a = sizeof(alg_list) / sizeof(alg_list[0]) - 1; a > -1; --a)
        if (algs & alg_list[a]) {
          auto alg = alg_list[a];
          TrueRegistrator reg;
          coll.SetRegistrator(&reg);
          find_int(n, coll, alg, stat);
          reg.write_SVG(alg, svg_stream);
        }
    };

    switch (seg_type) {
    case _Segment::intline: {
      CIntegerSegmentCollection<SimpleCounter> coll(n,(TLineSegment1 *) segs, nullptr,range);
      coll.coll_to_SVG(svg_stream);
    }; break;
    case _Segment::line1: {
      CLine1SegmentCollection<TrueRegistrator> coll(n, segs, nullptr);
      apply_algs(coll);
    }; break;
    case _Segment::line2: {
      CLine2SegmentCollection<TrueRegistrator> coll(n, segs, nullptr);
      apply_algs(coll);
    }; break;
    case _Segment::arc: {
      CArcSegmentCollection<TrueRegistrator>  coll(n, segs, nullptr);
      apply_algs(coll);
    }; break;
    case _Segment::graph: {
      CGraphSegmentCollection<TrueRegistrator> coll(n, segs, nullptr);
      apply_algs(coll);
    }; break;
    }
  };


  find_intersections_func get_find_intersections_func(uint4 reg_type)
  {
    
    if (reg_type == _Registrator::just_count)return _find_int<SimpleCounter>;
    if (reg_type == _Registrator::just_count_require_intersections)
      return _find_int<SimpleCounter2>;
    if (reg_type == _Registrator::store_pairs_and_ints_just_count_stat)
      return _find_int<TrueRegistrator>;
    return _find_int<PerSegmCounter>;
  };