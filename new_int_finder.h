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
    for (uint4 i = 0; i + 1 < N; ++i)
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
        uint4* sgm;
        DECL_RAII_ARR(sgm, n * 3);
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
    bool isnot_top() { return prev!=nullptr; };
  };


class CommonImpl 
{
public:
    bool dont_cut_stripe;

    template <class SegmentsColl>
    void prepare_ends(SegmentsColl& segments)
    {
      // moved to SegmentsColl::Prepare() - it computes ENDS, SegL, SegR and statistics.
      auto pr = segments.Prepare();
      // set internal common statistics derived from collection Prepare result
      nTotSegm = segments.GetSegmNumb();
      nTotX = pr.ends_len;
      LR_len = pr.max_segm_on_vline + 3; // +3 for sentinels/guards
      avr_segm_on_vline = pr.avr_segm_on_vline;
    };

    template <class SegmentsColl>
    void FindIntI(SegmentsColl& segments, int4 sn, ProgramStackRec* stack_pos) const
    {
      auto r_index = segments.GetSegR(sn);
      while (stack_pos->right_bound <= r_index)
        stack_pos = stack_pos->prev; // go from bottom to top and find staircase to start
      int4  m, len;
      int4* Qb, * Ql, * Qe = Q + (stack_pos->Q_pos + 1);
      for (stack_pos = stack_pos->prev; stack_pos != nullptr; stack_pos = stack_pos->prev) {
        Qb = Ql = Q + stack_pos->Q_pos;
        len = Qe - Qb;
        while (len > 1) // binary search
        {
          m = len / 2; //
          if (segments.UnderCurPoint(Ql[m]))
            Ql += m;
          len -= m;
        }
        FindInt(segments, Qb, Qe, Ql);
        Qe = Qb + 1; // move staircase bound to parent staircase
      };
    };



    ~CommonImpl() { FreeMem(); };
protected:
  constexpr static uint4 cut_margin = 24;
  constexpr static uint4 min_strip_width = 6;// at least 2

  int4* Q = nullptr;
  int4* L = nullptr;
  double avr_segm_on_vline = 0;
  uint4 LR_len = 0;
  uint4 nTotSegm = 0;
  uint4 len_of_Q = 0;
  uint4 L_size = 0;
  uint4 nTotX = 0;


  template<class SegmentsColl, bool HasSentinels = has_sentinels<SegmentsColl> >
  struct CSentinel {
    CSentinel(SegmentsColl& coll, int4& ini) {};
  };

  template<class SegmentsColl>
  struct CSentinel <SegmentsColl, true> {
    CSentinel(SegmentsColl& coll, int4& ini) : to_restore(ini) {
      prev_val = to_restore;
      to_restore = coll.get_sentinel(false);
    };
    ~CSentinel() {
      to_restore = prev_val;
    };
  private:
    int4& to_restore;
    int4 prev_val;
  };

  auto GetQTail() {
    return Q + len_of_Q;
  }

  void FreeMem()
  {
  };

  template <class IntersectionFinder, class SegmentsColl>
  static void SISFindR(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    for (auto i = interval_left_rank + 1; i != interval_right_rank; ++i) {
      if (i_f.L_size > 1) {
        segments.SetCurStripe(i - 1, i);
        i_f.SearchInStrip(segments, ladder_start_index);
      }
      i_f.InsDel(segments, i, stack_pos);
    }
    //the stripe right bound needs to be installed correctly even if i_f.L_size <= 1
    segments.SetCurStripe(interval_right_rank - 1, interval_right_rank);
    if (i_f.L_size > 1)i_f.SearchInStrip(segments, ladder_start_index);
  }


  template <class IntersectionFinder, class SegmentsColl>
  static void FindRNoChecks(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    int4 Q_pos;
    //if L or Q empty cut the strip into 2^divide_pow substrips
    constexpr const uint4 divide_pow = 4;
    if ((i_f.L_size == 0)||(0 == (Q_pos = i_f.Split(segments, i_f.Q + ladder_start_index, interval_right_rank))))
      return MultipleCutting(i_f, segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos, divide_pow);
    //normal processing otherwise
    Q_pos += ladder_start_index;
    ProgramStackRec stack_rec(Q_pos, interval_right_rank, stack_pos);// creating new stack list entry for "bubbling up" 
    if (i_f.dont_cut_stripe ) { //if found a lot of intersections repeat FindR
      FindRNoChecks(i_f, segments, Q_pos, interval_left_rank, interval_right_rank, &stack_rec);
    }else{// cut at the middle into two substrips 
      uint4 middle = (interval_left_rank + interval_right_rank) / 2;
      FindR(i_f, segments, Q_pos, interval_left_rank, middle, &stack_rec);
      i_f.InsDel(segments, middle, &stack_rec);
      FindR(i_f, segments, Q_pos, middle, interval_right_rank, &stack_rec);
      //actually works without SetCurStripeLeft, but it simplifies segment collection class protocol
      //segments.SetCurStripeLeft(i_f.interval_left_rank);
    }
    i_f.Merge(segments, interval_left_rank, ladder_start_index, Q_pos);
  };

  template <class IntersectionFinder, class SegmentsColl>
  static void FindR(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    if (interval_right_rank < interval_left_rank + min_strip_width) {//end of recursion for narrow stripes
      SISFindR(i_f, segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos);
    } else {// normal step
      segments.SetCurStripe(interval_left_rank, interval_right_rank);
      FindRNoChecks(i_f, segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos);
    }
  }

  uint4 GetDivPow(uint4 l) {
    if (l < 64) return 2;
    int4 res = 0;
    for (uint4 target = 2.0*l / MAX(avr_segm_on_vline, 8.0);target; target >>= 1)
      ++res;
    return MAX(res,4);
  };

  //cuts the strip into 2^divide_pow substrips
  template <class IntersectionFinder, class SegmentsColl>
  static void MultipleCutting(IntersectionFinder& i_f,SegmentsColl& segments, int4 ladder_start_index, uint4 left_bound, uint4 interval_right_rank, ProgramStackRec* stack_pos, uint4 divide_pow)
  {
    if (interval_right_rank < left_bound + (min_strip_width << divide_pow)) // end of recursion for narrow stripes
      return SISFindR(i_f, segments, ladder_start_index, left_bound, interval_right_rank, stack_pos);  //if strip narrow just apply SISFindR

    uint8 step = interval_right_rank - left_bound;
    uint8 rb = (uint8)left_bound << divide_pow;
    uint4 right_bound = (rb += step) >> divide_pow;
    segments.SetCurStripe(left_bound, right_bound);
    FindRNoChecks(i_f, segments, ladder_start_index, left_bound, right_bound, stack_pos);
    do {
      left_bound = right_bound;
      right_bound = (rb += step) >> divide_pow;
      i_f.InsDel(segments, left_bound, stack_pos);
      segments.SetCurStripe(left_bound, right_bound);
      FindRNoChecks(i_f, segments, ladder_start_index, left_bound, right_bound, stack_pos);
    } while (right_bound != interval_right_rank);
  }

  //functions for fast algorithm
  template <class SegmentsColl>
  static void FindInt(SegmentsColl& segments, int4* const qb, int4* const qe, int4* l)
  {
    auto c = l;
    while ((c != qb) && segments.FindCurSegIntDownWith(*c)) //first get intersections below
      --c;

    constexpr bool line_seg = (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) == _Coll_flag_state::state_true);

    if ((line_seg && (c != l)))
      return; //if found and segment is line or no stair above it can't be any more
    do {
      ++l;
    } while ((l != qe) && segments.FindCurSegIntUpWith(*l)); // get intersections above
  };

  template<class SegmentsColl>
  void SearchInStripNonLineSeg(SegmentsColl& segments, int4* _L, int4* _Q)//simplified version for SearchInStrip
  {
    auto last_L = _L + L_size;
    auto first_L = _L + 1;
    auto _Q_pos = _Q;
    auto T_pos = GetQTail();
    do
    {
      auto new_L_pos = _L;
      _Q_pos = _Q;
      segments.SetCurSegCutBE(*++_Q_pos = _L[0]);
      auto Q_tail = T_pos;
      auto cur_L = first_L;
      for (CSentinel sentinel(segments , *_Q); cur_L < last_L; ++cur_L)
      {
        segments.SetCurSegCutBE(*cur_L);
        
        //auto step = _Q_pos;
        //while ((step!=_Q) && (segments.FindCurSegIntDownWith(*step))) --step;
        //if (_Q_pos == step)

        if (_Q_pos == segments.FindCurSegIntDownWith(_Q_pos, _Q))
          *++_Q_pos = *cur_L;
        else
        {
          *new_L_pos++ = *cur_L;
          *--Q_tail = _Q_pos-_Q;
        }
      }
      if (new_L_pos == _L)// L is empty
        return;
      *--Q_tail = _Q_pos - _Q;//place a guard for the loop below
      ++_Q; ++_Q_pos;
      if constexpr (has_sentinels<SegmentsColl>)
        *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel
      
      Q_tail = T_pos;
      auto cur_Q = _Q + *--Q_tail;
      for (cur_L = _L ; cur_Q != _Q_pos; ++cur_L)
      {
        segments.SetCurSegCutBE(*cur_L);
        segments.FindCurSegIntUpWith(cur_Q, _Q_pos);
        cur_Q = _Q + *--Q_tail;
      }
      _Q = --_Q_pos;
      last_L = new_L_pos;
    } while (last_L > first_L);//L contains >1 segment
//    if (last_L == first_L)
      *++_Q = _L[0];
  };

  template<class SegmentsColl>
  void SearchInStripLineSeg(SegmentsColl& segments, int4* L_) {
    auto last = L_ + L_size;
    if constexpr (has_sentinels<SegmentsColl>) 
      *last= segments.get_sentinel(true);
    decltype(last) cur = last - 1, next;
    if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::needs_SetCurSegCutBE_at_start) == _Coll_flag_state::state_true)
      segments.SetCurSegCutBE(*cur);//we don't need SetCurSegCutBE in all cases, only for Y caching collections 
    while (cur != L_) {// insertion sort from last to first
      next = cur--;
      segments.SetCurSegCutBE(*cur);
      auto no_int = segments.FindCurSegIntUpWith(next, last);
      std::rotate(cur, next, no_int);
    }
  }

};


uint4 get_max_call(uint4 N) {
    uint4 max_call = 24;
    for (; N; N >>= 1) max_call += 2;
    return max_call;
};

//#undef MY_FREE_ARR_MACRO
#endif