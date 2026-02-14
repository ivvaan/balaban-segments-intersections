#pragma once
/*
*
*      Copyright (c)  2011-2026  Ivan Balaban
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
    CSentinel(SegmentsColl& coll, int4& ini, bool is_top = false) {};
  };

  template<class SegmentsColl>
  struct CSentinel <SegmentsColl, true> {
    CSentinel(SegmentsColl& coll, int4& ini, bool is_top = false) : to_restore(ini) {
      prev_val = to_restore;
      to_restore = coll.get_sentinel(is_top);
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

  uint4 GetDivPow(uint4 l) {
    if (l < 64) return 2;
    int4 res = 0;
    for (uint4 target = 2.0*l / MAX(avr_segm_on_vline, 8.0);target; target >>= 1)
      ++res;
    return MAX(res,4);
  };

  //functions for fast algorithm
  template <class SegmentsColl>
  static void FindInt(SegmentsColl& segments, int4* const qb, int4* const qe, int4* l)
  {
    // Use range overloads with sentinel instead of single-element overloads.
    // Bottom sentinel is placed via CSentinel RAII (restores *qb on scope exit).
    // Top sentinel is placed directly at *qe (no restore needed — qe is one past
    // the last stair, so the caller does not rely on its value).
    CSentinel bsentinel(segments, *qb);
    auto c = segments.FindCurSegIntDownWith(l, qb); //first get intersections below

    constexpr const bool line_seg = (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) == _Coll_flag_state::state_true);

    if (line_seg && (c != l))
      return; //if found and segment is line there can't be any more intersections above

    CSentinel tsentinel(segments, *qe, true);
    segments.FindCurSegIntUpWith(l + 1, qe); // get intersections above
  };

  template<class SegmentsColl>
  void SearchInStripNonLineSeg(SegmentsColl& segments, int4* _L, int4* _Q)//simplified version for SearchInStrip
  {
    auto last_L = _L + L_size;
    auto first_L = _L + 1;
    int4 * _Q_pos = _Q + 1;
    auto T_pos = GetQTail();
    do
    {
      auto new_L_pos = _L;
      *_Q_pos = _L[0];
      if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::needs_SetCurSegCutBE_at_start) == _Coll_flag_state::state_true)
        segments.SetCurSegCutBE(_L[0]);
      auto Q_tail = T_pos;
      auto cur_L = first_L;
      for (CSentinel sentinel(segments , *_Q); cur_L < last_L; ++cur_L)
      {
        segments.SetCurSegCutBE(*cur_L);
        
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
      _Q = _Q_pos - 1;
      last_L = new_L_pos;
    } while (last_L > first_L);//L contains >1 segment
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


/*inline uint4 get_max_call(uint4 N) {
    uint4 max_call = 24;
    for (; N; N >>= 1) max_call += 2;
    return max_call;
};
*/

