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
#include "optimal_finder.h"


constexpr int4 n_threads = 6;


// creating and deleting collections
PSeg first_segment_ptr = NULL;

SegCollection::SegCollection(int4 seg_type, PSeg segs) noexcept
  : seg_type_(seg_type),
  segs_(segs)
{
}

SegCollection::SegCollection(SegCollection&& other) noexcept
  : seg_type_(other.seg_type_),
  segs_(other.segs_)
{
  other.seg_type_ = -1;
  other.segs_ = nullptr;
}

SegCollection& SegCollection::operator=(SegCollection&& other) noexcept
{
  if (this != &other)
  {
    reset();
    seg_type_ = other.seg_type_;
    segs_ = other.segs_;
    other.seg_type_ = -1;
    other.segs_ = nullptr;
  }
  return *this;
}

SegCollection::~SegCollection()
{
  reset();
}

void SegCollection::reset() noexcept
{
  if (!segs_)
    return;

  switch (seg_type_)
  {
  case _Segment::line1:
  case _Segment::intline:
    delete[](TLineSegment1*)segs_;
    break;
  case _Segment::line2:
    delete[](TLineSegment2*)segs_;
    break;
  case _Segment::arc:
    delete[](TArcSegment*)segs_;
    break;
  case _Segment::graph:
    delete[](TPlaneVect*)segs_;
    break;
  default:
    break;
  }

  segs_ = nullptr;
  seg_type_ = -1;
}

SegCollection create_test_collection(const SegmentsAndRegOptions& opt, CRandomValueGen& random_gen)
{
  int4 i;
  constexpr int4 sentinels = 2;

  auto par = opt.distr_param;
  if (opt.distr_type != _Distribution::param_defined)par /= opt.n / 33.0;
  if (opt.distr_type == _Distribution::parallel)par /= 2.4;

  PSeg result = nullptr;
  switch (opt.seg_type) {
    case _Segment::line1:
    case _Segment::intline:
    {
      TLineSegment1* Colls = new TLineSegment1[opt.n + sentinels];
      for (i = 0; i < opt.n; i++)
        Colls[i].InitRandom(random_gen, i, opt.distr_type, par);
      auto mmr = get_minmax(opt.n, Colls).get_scaled(2.0);
      Colls[opt.n] = TLineSegment1(mmr.ld, { mmr.rt.x,mmr.ld.y });
      Colls[opt.n + 1] = TLineSegment1({ mmr.ld.x,mmr.rt.y }, mmr.rt);
      result = Colls;
      first_segment_ptr = Colls;
    }; break;
    case _Segment::line2:
    {
      TLineSegment2* Colls = new TLineSegment2[opt.n + sentinels];
      for (i = 0; i < opt.n; i++)
        Colls[i].InitRandom(random_gen, i, opt.distr_type, par);
      auto mmr = get_minmax(opt.n, Colls).get_scaled(2.0);
      Colls[opt.n].Init({ mmr.ld, { mmr.rt.x,mmr.ld.y } });
      Colls[opt.n + 1].Init({ { mmr.ld.x,mmr.rt.y }, mmr.rt });
      result = Colls;
      first_segment_ptr = Colls;
    }; break;
    case _Segment::arc:
    {
      TArcSegment* Colls = new TArcSegment[opt.n];
      for (i = 0; i < opt.n; i++)
        Colls[i].InitRandom(random_gen, i, opt.distr_type, par);
      result = Colls;
      first_segment_ptr = Colls;
    }; break;
    case _Segment::graph:
    {
      TPlaneVect* Colls = new TPlaneVect[opt.n];
      for (i = 0; i < opt.n; i++)
        Colls[i].InitRandom(random_gen);
      result = Colls;
      first_segment_ptr = Colls;
    }; break;
  };
  return SegCollection(opt.seg_type, result);
};





template<template<class>class SegColl, class Counter>
auto find_int(int4 n, SegColl<Counter>& coll, int4 alg, uint4 stat)
{
  switch (alg) {
    case triv: {
      CTrivialIntFinder fi;
      fi.find_intersections(coll);
    } break;
    case simple_sweep: {
      CSimpleSweepIntFinder fi;
      fi.find_intersections(coll);
    }; break;
    case fast: {
      CFastIntFinder fi;
      fi.prepare_ends(coll);
      fi.find_intersections(coll);
    } break;
    case optimal: {
      COptimalIntFinder fi;
      fi.prepare_ends(coll);
      fi.find_intersections(coll);
    } break;
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
  };
  return coll.GetRegistrator()->get_stat(stat);
};


template<class Counter>
auto _find_int(const SegmentsAndRegOptions& opt, PSeg segs, int4 alg)
{
  Counter reg;
  reg.Alloc(opt.n);
  switch (opt.seg_type) {
    case _Segment::line1: {
      CLine1SegmentCollection<Counter> coll(opt.n, segs, &reg);
      return find_int(opt.n, coll, alg, opt.reg_stat);
    };
    case _Segment::intline: {
      CIntegerSegmentCollection<Counter> coll(opt.n, (TLineSegment1*)segs, &reg, opt.range_for_int_seg);
      return find_int(opt.n, coll, alg, opt.reg_stat);
    };
    case _Segment::line2: {
      CLine2SegmentCollection<Counter> coll(opt.n, segs, &reg);
      return find_int(opt.n, coll, alg, opt.reg_stat);
    };
    case _Segment::arc: {
      CArcSegmentCollection<Counter>  coll(opt.n, segs, &reg);
      return find_int(opt.n, coll, alg, opt.reg_stat);
    };
    case _Segment::graph: {
      CGraphSegmentCollection<Counter> coll(opt.n, segs, &reg);
      return find_int(opt.n, coll, alg, opt.reg_stat);
    };
  };
  return 0ULL;
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

void write_SVG(std::ostream* svg_stream, const SegmentsAndRegOptions& opt, PSeg segs, int4 algs)
{
  if (!svg_stream)
    return;

  auto apply_algs = [&](auto& coll) {
    coll.coll_to_SVG(svg_stream);
    for (int4 a = sizeof(alg_list) / sizeof(alg_list[0]) - 1; a > -1; --a)
      if (algs & alg_list[a]) {
        auto alg = alg_list[a];
        TrueRegistrator reg;
        coll.SetRegistrator(&reg);
        find_int(opt.n, coll, alg, 0);
        reg.write_SVG(alg, svg_stream);
      }
    };

  switch (opt.seg_type) {
    case _Segment::intline: {
      CIntegerSegmentCollection<SimpleCounter> coll(opt.n, (TLineSegment1*)segs, nullptr, opt.range_for_int_seg);
      coll.coll_to_SVG(svg_stream);
    }; break;
    case _Segment::line1: {
      CLine1SegmentCollection<TrueRegistrator> coll(opt.n, segs, nullptr);
      apply_algs(coll);
    }; break;
    case _Segment::line2: {
      CLine2SegmentCollection<TrueRegistrator> coll(opt.n, segs, nullptr);
      apply_algs(coll);
    }; break;
    case _Segment::arc: {
      CArcSegmentCollection<TrueRegistrator>  coll(opt.n, segs, nullptr);
      apply_algs(coll);
    }; break;
    case _Segment::graph: {
      CGraphSegmentCollection<TrueRegistrator> coll(opt.n, segs, nullptr);
      apply_algs(coll);
    }; break;
  };
};

