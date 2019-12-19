#pragma once
/*
*
*      Copyright (c)  2011-2017  Ivan Balaban
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

#include "utils.h"
#include <math.h>
#include <memory.h>
#include <time.h>
//#include <array>
//#include <memory>
//#include "test_coll.h"

#define MY_FREE_ARR_MACRO(a) \
    if (a!=nullptr) {                 \
        delete[] a;          \
        a = nullptr;            \
    }

template<class SegmentsColl>
class CTrivialIntFinder
{
public:

  void find_intersections(SegmentsColl *segments) 
  {
    uint4 N = segments->GetSegmNumb();
    for (uint4 i = 0; i < N; ++i)
    {
      segments->SetCurSeg(i);
      for (uint4 j = i + 1; j < N; ++j)
        segments->TrivCurSegIntWith(j);

    }
  };

};


template <class SegmentsColl>
class CSimpleSweepIntFinder {
public:
    void find_intersections(SegmentsColl* segments)
    {
        uint4 n = segments->GetSegmNumb();
        uint4* sgm = new uint4[n * 3];
        uint4 *Ends = sgm+n;
        segments->PrepareEndpointsSortedList(Ends);

        auto pos = sgm;
        for (auto i = Ends,l=i+2*n; i < l; ++i) {
            auto s = SegmentsColl::get_segm(*i);
            if (SegmentsColl::is_last(*i)) 
            {
                auto j = --pos;
                auto cur = *j;
                while (cur != s) {
                    --j;
                    auto buf = cur;
                    cur = *j;
                    *j = buf;
                }
            } 
            else 
            {
                segments->SetCurSeg(s);
                for (auto j = sgm; j < pos; ++j)
                    segments->SSCurSegIntWith(*j);
                *pos = s;
                ++pos;
            }
        }
        delete[] sgm;
    };
};


  struct ProgramStackRec //structure to use program stack to store list of starcases above current call
  {
    ProgramStackRec *prev=nullptr; //link to program stack record storing parent staircase information 
    int4 Q_pos; //starting position and 
    uint4 right_bound;//right boind of current staircase
    inline ProgramStackRec(int4 qp) :Q_pos(qp) {};
    inline ProgramStackRec(int4 qp, uint4 rb) : Q_pos(qp), right_bound(rb) {};
    inline ProgramStackRec *Set(ProgramStackRec *p, uint4 rb) { prev = p, right_bound = rb; return this; };
  };

template <class SegmentsColl>
class CommonImpl 
{
protected:
    static const int4 max_call = 5; //max number of sequential recursive call (opt)FindR before dividing current strip

    uint4 *SegL = nullptr, *SegR = nullptr, *ENDS = nullptr;
    int4 *Q = nullptr;
    void AllocMem(uint4 N)
    {
        SegL = new uint4[N];
        SegR = new uint4[N];
        ENDS = new uint4[2 * N];
    };

    void FreeMem()
    {
        MY_FREE_ARR_MACRO(SegL);
        MY_FREE_ARR_MACRO(SegR);
        MY_FREE_ARR_MACRO(ENDS);
    };

    void prepare_ends(SegmentsColl* segments)
    {
        auto N = segments->GetSegmNumb();
        AllocMem(N);
        segments->PrepareEndpointsSortedList(ENDS);
        for (uint4 i = 0; i < 2 * N; i++)
            if (SegmentsColl::is_last(ENDS[i]))
                SegR[SegmentsColl::get_segm(ENDS[i])] = i;
            else
                SegL[SegmentsColl::get_segm(ENDS[i])] = i;
    };

    //functions for fast algorithm
    void FindInt(SegmentsColl* segments, int4 qb, int4 qe, int4 l)
    {
        int4 c = l;
        while ((c > qb) && (segments->FindCurSegIntWith(Q[c]))) //first get intersections below
            --c;
        if (SegmentsColl::is_line_segments && (c != l))
            return; //if found and segment is line it can't be any more
        c = l + 1;
        while ((c <= qe) && (segments->FindCurSegIntWith(Q[c]))) // get intersections above
            ++c;
    };

    void FindIntI(SegmentsColl* segments, uint4 r_index, ProgramStackRec* stack_pos)
    {
        while (stack_pos->right_bound <= r_index)
            stack_pos = stack_pos->prev; // go from bottom to top and find staircase to start
        int4 l, r, m, QB, QE = stack_pos->Q_pos;
        for (stack_pos = stack_pos->prev; stack_pos!=nullptr; stack_pos = stack_pos->prev) {
            l = QB = stack_pos->Q_pos;
            r = QE + 1;
            while ((r - l) > 1) // binary search
            {
                m = (r + l) >> 1; //        m=(r+l)/2;
                if (segments->UnderCurPoint(Q[m]))
                    l = m;
                else
                    r = m;
            }
            FindInt(segments, QB, QE, l);
            QE = QB; // move staircase bound to parent staircase
        };
    };
};




//#undef MY_FREE_ARR_MACRO