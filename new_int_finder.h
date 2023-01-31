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
    bool dont_cut_stripe;

    template <class SegmentsColl>
    void prepare_ends(SegmentsColl& segments)
    {
      auto N= segments.GetSegmNumb();
      ENDS = new uint4[2 * N];
      segments.PrepareEndpointsSortedList(ENDS);
      SegL = new uint4[N];
      SegR = new uint4[N];
      uint4 max_segm_on_vline=0,nsegm_on_vline = 0;
      avr_segm_on_vline = 0;
      for (uint4 i = 0; i < 2 * N; i++) {
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
        avr_segm_on_vline += nsegm_on_vline;
      }
      avr_segm_on_vline /= N * 2;
      LR_len = max_segm_on_vline + 3;//+3 for sentinels(guards) to have some space
      nTotSegm = N;
          
    };

    ~CommonImpl() { FreeMem(); };
protected:
  constexpr static uint4 cut_margin=24;

  uint4 LR_len = 0;
  uint4 nTotSegm = 0;
  uint4 len_of_Q = 0;
  uint4 L_size = 1;

  uint4* SegL = nullptr, * SegR = nullptr, * ENDS = nullptr;
  int4* Q = nullptr;
  int4* L = nullptr;
  double avr_segm_on_vline = 0;

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
    MY_FREE_ARR_MACRO(SegL);
    MY_FREE_ARR_MACRO(SegR);
    MY_FREE_ARR_MACRO(ENDS);
  };

  template <class IntersectionFinder, class SegmentsColl>
  static void SISFindR(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos)
  {
    auto ENDS = i_f.ENDS;
    for (auto i = interval_left_index + 1; i != interval_right_index; ++i) {
      auto pt = ENDS[i];
      if (i_f.L_size > 1) {
        segments.SetCurStripe(ENDS[i - 1], pt);
        i_f.SearchInStrip(segments, ladder_start_index);
      }
      i_f.InsDel(segments, pt, stack_pos);
    }
    //the stripe right bound needs to be installed correctly even if i_f.L_size <= 1
    segments.SetCurStripe(ENDS[interval_right_index - 1], ENDS[interval_right_index]);
    if (i_f.L_size > 1)i_f.SearchInStrip(segments, ladder_start_index);
    return;
  }


  template <class IntersectionFinder, class SegmentsColl>
  static void FindRNoChecks(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos)
  {
    auto ENDS = i_f.ENDS;
    do {
      if (i_f.L_size == 0)break;
      auto Q_pos = i_f.Split(segments, i_f.Q + ladder_start_index, interval_right_index);
      if (Q_pos == 0)break;
      Q_pos += ladder_start_index;
      ProgramStackRec stack_rec(Q_pos, interval_right_index, stack_pos);
      if (i_f.dont_cut_stripe) { //if found a lot of intersections repeat FindR
        i_f.dont_cut_stripe = false;
        FindRNoChecks(i_f, segments, Q_pos, interval_left_index, interval_right_index, &stack_rec);
        i_f.Merge(segments, interval_left_index, ladder_start_index, Q_pos);
        return;
      }
      uint4 m = (interval_left_index + interval_right_index) / 2;
      _FindR(i_f, segments, Q_pos, interval_left_index, m, &stack_rec);
      i_f.InsDel(segments, ENDS[m], &stack_rec);
      _FindR(i_f, segments, Q_pos, m, interval_right_index, &stack_rec);
   
//actually works without SetCurStripeLeft, but it simplifies segment collection class protocol
      segments.SetCurStripeLeft(ENDS[interval_left_index]);

      i_f.Merge(segments, interval_left_index, ladder_start_index, Q_pos);
      return;
    } while (false);
    //if L or Q empty cut into 16 parts
    constexpr const uint4 divide_pow = 4;
//    constexpr const uint4 divide_into = 1 << divide_pow;//16
    constexpr const uint4 min_size = 128;//must be >= 1<<(divide_pow+1)
//    constexpr const uint4 min_step = 12;

    auto step = interval_right_index - interval_left_index;

    if (step <= min_size) {
      return SISFindR(i_f, segments, ladder_start_index, interval_left_index, interval_right_index, stack_pos);
    }

    //step = MAX(step, (1 << divide_pow)*min_step);
    auto rb = step + (interval_left_index << divide_pow);//* divide_into
    auto left_bound= interval_left_index;
    auto right_bound = rb >> divide_pow; // / divide_into;
    auto stripe_left = ENDS[left_bound];
    do {
      auto stripe_right = ENDS[right_bound];
      segments.SetCurStripe(stripe_left, stripe_right);
      FindRNoChecks(i_f, segments, ladder_start_index, left_bound, right_bound, stack_pos);
      i_f.InsDel(segments, stripe_right, stack_pos);
      left_bound = right_bound;
      stripe_left = stripe_right;
      right_bound = (rb += step) >> divide_pow; // / divide_into;
    } while (right_bound < interval_right_index);
    //_FindR(i_f, segments, ladder_start_index, left_bound, interval_right_index, stack_pos);
    assert(interval_right_index - left_bound > 1);
    segments.SetCurStripe(stripe_left, ENDS[interval_right_index]);
    FindRNoChecks(i_f, segments, ladder_start_index, left_bound, interval_right_index, stack_pos);
  };

  template <class IntersectionFinder, class SegmentsColl>
  static void _FindR(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos)
  {
    segments.SetCurStripe(i_f.ENDS[interval_left_index], i_f.ENDS[interval_right_index]);
    if (interval_right_index - interval_left_index == 1) {
      if (i_f.L_size > 1) {
        i_f.SearchInStrip(segments, ladder_start_index);
      }
      return;
    }
    FindRNoChecks(i_f, segments, ladder_start_index, interval_left_index, interval_right_index, stack_pos);
  }

  static uint4 GetDivPow(uint4 l, double avr_segm_on_vline) {
    int4 res = 0;
    for (uint4 target = 2.0 * l / MAX(avr_segm_on_vline, 8.0);target; target >>= 1)
      ++res;
    return MAX(res,4);
  };

  template <class IntersectionFinder, class SegmentsColl>
  static void FindR(IntersectionFinder& i_f, SegmentsColl& segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos)
  {
    auto ENDS = i_f.ENDS;
    long long step = interval_right_index - interval_left_index;
    uint4 divide_pow = GetDivPow(step, i_f.avr_segm_on_vline);
    long long rb = step + (((long long)interval_left_index) << divide_pow);
    uint4 left_bound = interval_left_index;
    uint4 right_bound = rb >> divide_pow; 
    auto stripe_left = ENDS[left_bound];
    do {
      auto stripe_right = ENDS[right_bound];
      segments.SetCurStripe(stripe_left, stripe_right);
      FindRNoChecks(i_f, segments, ladder_start_index, left_bound, right_bound, stack_pos);
      i_f.InsDel(segments, stripe_right, stack_pos);
      left_bound = right_bound;
      stripe_left = stripe_right;
      right_bound = (rb += step) >> divide_pow; 
    } while (right_bound < interval_right_index);
    assert(interval_right_index - left_bound > 1);
    segments.SetCurStripe(stripe_left, ENDS[interval_right_index]);
    FindRNoChecks(i_f, segments, ladder_start_index, left_bound, interval_right_index, stack_pos);
  }

  //functions for fast algorithm
  template <class SegmentsColl>
  static void FindInt(SegmentsColl& segments, int4* const qb, int4* const qe, int4* l)
  {
    auto c = l;
    while ((c != qb) && segments.FindCurSegIntDownWith(*c)) //first get intersections below
      --c;
    if ((SegmentsColl::is_line_segments && (c != l)))
      return; //if found and segment is line or no stair above it can't be any more
    do {
      ++l;
    } while ((l != qe) && segments.FindCurSegIntUpWith(*l)); // get intersections above
  };

  template <class SegmentsColl>
  void FindIntI(SegmentsColl& segments, uint4 r_index, ProgramStackRec* stack_pos) const
  {
    while (stack_pos->right_bound <= r_index)
      stack_pos = stack_pos->prev; // go from bottom to top and find staircase to start
    int4  m, len;
    int4 *Qb,  *Ql, *Qe = Q + (stack_pos->Q_pos + 1);
    for (stack_pos = stack_pos->prev; stack_pos != nullptr; stack_pos = stack_pos->prev) {
      Qb = Ql = Q + stack_pos->Q_pos;
      len = Qe - Qb;
      while (len > 1) // binary search
      {
        m = len / 2; //
        if(segments.UnderCurPoint(Ql[m]))
          Ql += m;
        len -= m;
      }
      FindInt(segments, Qb, Qe, Ql);
      Qe = Qb+1; // move staircase bound to parent staircase
    };
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
    auto Size = L_size;
    auto _L = L_ - 1;
    segments.SetCurSegCutBE(L_[0]);
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
    uint4 max_call = 24;
    for (; N; N >>= 1) max_call += 2;
    return max_call;
};

//#undef MY_FREE_ARR_MACRO
#endif