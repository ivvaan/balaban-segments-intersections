#ifndef NEW_INT_FINDER_FOR_SEGMENT_INTERSECTION
#define NEW_INT_FINDER_FOR_SEGMENT_INTERSECTION
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
  template<class SegmentsColl>
  void find_intersections(SegmentsColl &segments) 
  {
    uint4 N = segments.GetSegmNumb();
    for (uint4 i = 0; i < N; ++i)
    {
      segments.SetCurSeg(i);
      for (uint4 j = i + 1; j < N; ++j)
        segments.TrivCurSegIntWith(j);

    }
  };

};



class CSimpleSweepIntFinder {
public:
  template <class SegmentsColl>
    void find_intersections(SegmentsColl& segments)
    {
        uint4 n = segments.GetSegmNumb();
        uint4* sgm = new uint4[n * 3];
        uint4 *Ends = sgm+n;
        segments.PrepareEndpointsSortedList(Ends);

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
                segments.SetCurSeg(s);
                for (auto j = sgm; j < pos; ++j)
                    segments.SSCurSegIntWith(*j);
                *pos = s;
                ++pos;
            }
        }
        delete[] sgm;
    };
};


  struct ProgramStackRec //structure to use program stack to store list of starcases above current call
  {
    int4 Q_pos; //starting position and 
    uint4 right_bound;//right boind of current staircase
    ProgramStackRec* prev = nullptr; //link to program stack record storing parent staircase information 
    ProgramStackRec(int4 qp) :Q_pos(qp) {};
    ProgramStackRec(int4 qp, uint4 rb) : Q_pos(qp), right_bound(rb) {};
    ProgramStackRec(int4 qp, uint4 rb, ProgramStackRec* pr) : Q_pos(qp), right_bound(rb), prev(pr) {};
    ProgramStackRec* Set(ProgramStackRec* p, uint4 rb) { prev = p, right_bound = rb; return this; };

  };


class CommonImpl 
{
public:
    bool dont_split_stripe;

    template <class SegmentsColl>
    void prepare_ends(SegmentsColl& segments)
    {
      auto N= segments.GetSegmNumb();
      ENDS = new uint4[2 * N];
      segments.PrepareEndpointsSortedList(ENDS);
      SegL = new uint4[N];
      SegR = new uint4[N];
      uint4 max_segm_on_vline=0,nsegm_on_vline = 0;
      for (uint4 i = 0; i < 2 * N; i++)
        if (SegmentsColl::is_last(ENDS[i]))
        {
          SegR[SegmentsColl::get_segm(ENDS[i])] = i;
          --nsegm_on_vline;
        }
        else
        {
          SegL[SegmentsColl::get_segm(ENDS[i])] = i;
          ++nsegm_on_vline;
          if (nsegm_on_vline > max_segm_on_vline)
            max_segm_on_vline = nsegm_on_vline;
        }
      LR_len = max_segm_on_vline+1;
      nTotSegm = N;
          
    };

    ~CommonImpl() { FreeMem(); };
protected:

  uint4 LR_len = 0;
  uint4 nTotSegm = 0;
  uint4 len_of_Q = 0;
  uint4 L_size = 1;

  uint4* SegL = nullptr, * SegR = nullptr, * ENDS = nullptr;
  int4* Q = nullptr;
  int4* L = nullptr;

  void FreeMem()
  {
    MY_FREE_ARR_MACRO(SegL);
    MY_FREE_ARR_MACRO(SegR);
    MY_FREE_ARR_MACRO(ENDS);
  };

  template <class IntersectionFinder, class SegmentsColl>
  static  void _FindR(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos, int4 call_numb, int4 max_call = 30)
  {
    do {
      if (i_f.L_size == 0)break;
      auto Q_pos = i_f.Split(segments, i_f.Q + ladder_start_index, interval_right_index);
      if (Q_pos == 0)break;
      Q_pos += ladder_start_index;
      ProgramStackRec stack_rec(Q_pos, interval_right_index, stack_pos);
      if (i_f.dont_split_stripe && (call_numb < max_call)) { //if found a lot of intersections repeat FindR
        _FindR(i_f, segments, Q_pos, interval_left_index, interval_right_index, &stack_rec, call_numb + 1, max_call);
        i_f.Merge(segments, interval_left_index, ladder_start_index, Q_pos);
        return;
      }
      uint4 m = (interval_left_index + interval_right_index) >> 1;
      if ((interval_left_index + 4 > interval_right_index)|| (call_numb > 1)&& (i_f.L_size != 0)) {
        // if L contains a lot of segments then cut on two parts
        max_call -= 2;
        FindR(i_f, segments, Q_pos, interval_left_index, m, &stack_rec, 0, max_call);
        i_f.InsDel(segments, m, &stack_rec);
        FindR(i_f, segments, Q_pos, m, interval_right_index, &stack_rec, 0, max_call);
      } else {// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
        max_call -= 4;
        uint4 q = (interval_left_index + m) >> 1;
        FindR(i_f, segments, Q_pos, interval_left_index, q, &stack_rec, 0, max_call);
        i_f.InsDel(segments, q, &stack_rec);
        FindR(i_f, segments, Q_pos, q, m, &stack_rec, 0, max_call);
        i_f.InsDel(segments, m, &stack_rec);
        q = (interval_right_index + m) >> 1;
        FindR(i_f, segments, Q_pos, m, q, &stack_rec, 0, max_call);
        i_f.InsDel(segments, q, &stack_rec);
        FindR(i_f, segments, Q_pos, q, interval_right_index, &stack_rec, 0, max_call);
      }
        //actually works without this line, but it simplifies segment collection class
      //protocol
      segments.SetCurStripeLeft(i_f.ENDS[interval_left_index]);

      i_f.Merge(segments, interval_left_index, ladder_start_index, Q_pos);
      return;
    } while (false);
    //if L or Q empty cut into 8 parts
    constexpr const int4 divide_into = 8;
    max_call -= 6;
    double step = MAX(1.0, static_cast<double>(interval_right_index - interval_left_index) / divide_into);
    double rb = 0.25 + interval_left_index;
    uint4 right_bound = interval_left_index;
    while (true) {
      auto left_bound = right_bound;
      right_bound = rb += step;
      FindR(i_f, segments, ladder_start_index, left_bound, right_bound, stack_pos, 0, max_call);
      if (right_bound == interval_right_index)
        return;
      i_f.InsDel(segments, right_bound, stack_pos);
    };
  };

  template <class IntersectionFinder, class SegmentsColl>
  static  void FindR(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos, int4 call_numb, int4 max_call = 30)
  {
    segments.SetCurStripe(i_f.ENDS[interval_left_index], i_f.ENDS[interval_right_index]);
    if (interval_right_index - interval_left_index == 1) {
      if (i_f.L_size > 1)
        i_f.SearchInStrip(segments, ladder_start_index);
      return;
    }
    _FindR(i_f, segments, ladder_start_index, interval_left_index, interval_right_index, stack_pos, call_numb, max_call);

  }

  //functions for fast algorithm
  template <class SegmentsColl>
  static void FindInt(SegmentsColl& segments, int4* const qb, int4* const qe, int4* l)
  {
    auto c = segments.FindCurSegIntDownWith(l,qb);
    if ((SegmentsColl::is_line_segments && (c != l)))
      return; //if found and segment is line or no stair above it can't be any more
    segments.FindCurSegIntUpWith(++l, qe);
  };

  template <class SegmentsColl>
  void FindIntI(SegmentsColl& segments, uint4 r_index, ProgramStackRec* stack_pos) const
  {
    while (stack_pos->right_bound <= r_index)
      stack_pos = stack_pos->prev; // go from bottom to top and find staircase to start
    auto r = stack_pos->Q_pos + 1;
    int4 l, m, qb;
    for (stack_pos = stack_pos->prev; stack_pos != nullptr; stack_pos = stack_pos->prev) {
      l = qb = stack_pos->Q_pos;
      auto Qe = Q + r;
      while (r > l + 1) // binary search
      {
        m = (r + l) / 2; // 
        if (segments.UnderCurPoint(Q[m]))
          l = m;
        else
          r = m;
      }
      FindInt(segments, Q+qb, Qe, Q+l);
      r = qb+1; // move staircase bound to parent staircase
    };
  };

  template<class SegmentsColl>
  void SearchInStripNonLineSeg(SegmentsColl& segments, int4* _L, int4* _Q)//simplified version for SearchInStrip
  {
    auto Size = L_size;
    do
    {
      auto Q_tail = Q + len_of_Q;
      int4 new_L_size = 0;
      int4 step_index = 1;
      _Q[1] = _L[0];
      for (uint4 cur_L_pos = 1; cur_L_pos < Size; cur_L_pos++)
      {
        auto step = step_index;
        auto cur_seg = _L[cur_L_pos];
        segments.SetCurSegCutBE(cur_seg);
        while ((step) && (segments.FindCurSegIntDownWith(_Q[step])))
          --step;

        if (step_index == step)
          _Q[++step_index] = cur_seg;
        else
        {
          _L[new_L_size++] = cur_seg;
          *(--Q_tail) = step_index;
        }
      }
      Q_tail = Q + len_of_Q;
      for (int4 i = 0; i < new_L_size; ++i)
      {
        auto c = *(--Q_tail);
        if (c >= step_index) break;
        segments.SetCurSegCutBE(_L[i]);
        while ((c < step_index) && (segments.FindCurSegIntDownWith(_Q[++c])));
      }
      _Q += step_index;
      Size = new_L_size;
    } while (Size);

  };

  template<class SegmentsColl>
  void SearchInStripLineSeg(SegmentsColl& segments, int4* L_) {
    auto Size = L_size;
    auto _L = L_ - 1;
    for (uint4 i = 1; i < Size; i++)
    {
      auto sn = L_[i];
      segments.SetCurSegCutBE(sn);
      if (segments.FindCurSegIntDownWith(_L[i])) {
        L_[i] = _L[i];
        uint4 j = i - 1;
        for (; (j) && (segments.FindCurSegIntDownWith(_L[j])); --j)
          L_[j] = _L[j];
        L_[j] = sn;
      }
    }
  }

};


uint4 get_max_call(uint4 N) {
    uint4 max_call = 8;
    for (; N; N >>= 1) max_call += 2;
    return max_call;
};

//#undef MY_FREE_ARR_MACRO
#endif