#pragma once
/*
*
*      Copyright (c)  2011-2026  Ivan Balaban
*      ivvaan@gmail.com
*
This file is part of Seg_int library.

SPDX-License-Identifier: BSL-1.0
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
https://www.boost.org/LICENSE_1_0.txt)
*/
// NEW IMPLEMENTATION

#include "utils.h"
#include <math.h>
#include <memory.h>
#include <time.h>
#include <cassert>
//#include <array>
//#include <memory>
//#include "test_coll.h"



class CTrivialIntFinder
{
public:
  template<TrivialSegColl SegmentsColl>
  void find_intersections(SegmentsColl &segments) 
  {
    uint4 N = segments.GetSegmNumb();
    for (uint4 i = 0; i + 1 < N; ++i)
    {
      segments.SetCurSeg(i);
      for (uint4 j = i + 1; j < N; ++j)
        segments.TrivCurSegIntWith(j);

    }
    segments.IntersectionsFindingDone();
  };

};



class CSimpleSweepIntFinder {
public:
  template <SweepSegColl SegmentsColl>
    void find_intersections(SegmentsColl& segments)
    {
        uint4 n = segments.GetSegmNumb();
        uint4* sgm;
        DECL_RAII_ARR(sgm, n * 3);
        uint4 *Ends = sgm+n;
        auto nEvents=segments.PrepareEndpointsSortedList(Ends);

        auto pos = sgm;
        for (auto i = Ends,l = i + nEvents; i != l; ++i) {
            auto s = SegmentsColl::get_segm(*i);
            if (SegmentsColl::is_last(*i)) 
            {
              auto last = std::remove(sgm, pos, s);
              assert(last + 1 == pos);
              pos=last;
            } 
            else 
            {
                segments.SetCurSeg(s);
                for (auto j = sgm; j != pos; ++j)
                    segments.SSCurSegIntWith(*j);
                *pos = s;
                ++pos;
            }
        }
        segments.IntersectionsFindingDone();
    };
};
