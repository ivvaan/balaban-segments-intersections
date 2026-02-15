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

#include <algorithm>
#include <vector>
#include <thread>
#include <cassert>
#include <type_traits>

class CFastIntFinder 
{
  template<class SegmentsColl, bool HasSentinels = has_sentinels<SegmentsColl> >
  struct CSentinel {
    CSentinel(SegmentsColl& coll, uint4& ini, bool is_top = false) {};
  };

  template<class SegmentsColl>
  struct CSentinel <SegmentsColl, true> {
    CSentinel(SegmentsColl& coll, uint4& ini, bool is_top = false) : to_restore(ini) {
      prev_val = to_restore;
      to_restore = coll.get_sentinel(is_top);
    };
    ~CSentinel() {
      to_restore = prev_val;
    };
  private:
    uint4& to_restore;
    uint4 prev_val;
  };

public:
 
  CFastIntFinder() = default;
  ~CFastIntFinder() { unclone(); FreeMem(); };
  void FreeMem(){};

  template<FastSegColl SegmentsColl>
  void find_intersections(SegmentsColl  &segments,uint4 from=0,uint4 to=0)
  {
    //AllocMem
    len_of_Q = LR_len;
    DECL_RAII_ARR(L, LR_len+1);
    ++L;//to have one cell before L for sentinel
    DECL_RAII_ARR(R, LR_len+1);
    ++R;//to have one cell before R for sentinel
    DECL_RAII_ARR(Q, len_of_Q);

    bool not_parallel = false;
    if (to == 0) {
      from = 0;
      to = nTotX - 1;
      not_parallel = true;
    }
    constexpr uint4 bottom_index = 0;
    ProgramStackRec stack_rec(bottom_index, nTotX); //need to be initialized this way
    if (from == 0) {
      segments.SetCurStripeRight(0);
      InsDel(*this,segments,0,&stack_rec);
    }
    else
      L_size = CalcLAt(segments, from);
    if (not_parallel) {
      if (avr_segm_on_vline < 20)
        SISFindR(segments, bottom_index, from, to, &stack_rec);
      else
        MultipleCutting(segments, bottom_index, from, to, &stack_rec, GetDivPow(to - from));
    }
    else //parallel call
      FindR(segments, bottom_index, from, to, &stack_rec);
    InsDel(*this,segments, to, &stack_rec);
    if(not_parallel)
      segments.IntersectionsFindingDone();
  }

  template<ParallelSegColl SegmentsColl>
  void find_intersections(uint4 n_threads,SegmentsColl& segments)
  {
    using namespace std;
    vector<thread> wrk_threads;
    auto thread_func = [](CFastIntFinder* master, SegmentsColl* segments,uint4 thread_index, uint4 from, uint4 to) {
      SegmentsColl coll(*segments, thread_index);
      CFastIntFinder(master).find_intersections(coll, from, to);
      };
    segments.InitClone(n_threads);
    auto n = segments.GetSegmNumb();
    double part = nTotX / (double)n_threads;
    uint4 start_from = part;
    for (uint4 i = 2, from = start_from, to; i <= n_threads; ++i) {
      to = (i == n_threads) ? nTotX-1 : (uint4)(part * i);
#if defined(DEBUG) || defined(_DEBUG) 
      thread_func(this, &segments,i-1, from, to); // starts intersection finding in a stripe <from,to>
#else
      wrk_threads.emplace_back(thread_func, this, &segments,i-1, from, to); // starts intersection finding in a stripe <from,to>
#endif
      from = to;
    }

    find_intersections(segments, 0, start_from);
#if !(defined(DEBUG) || defined(_DEBUG))
    for (auto& cur_thread : wrk_threads)
      cur_thread.join(); //waiting for calculation of all threads are finished
#endif
    segments.IntersectionsFindingDone();
  }


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


protected:
  CFastIntFinder(CFastIntFinder* master) {
    clone(master);
  };

  struct ProgramStackRec //structure to use program stack to store list of starcases above current call
  {
    uint4 Q_pos; //starting position and 
    uint4 right_bound;//right boind of current staircase
    ProgramStackRec* prev = nullptr; //link to program stack record storing parent staircase information 
    ProgramStackRec(uint4 qp) :Q_pos(qp) {};
    ProgramStackRec(uint4 qp, uint4 rb) : Q_pos(qp), right_bound(rb) {};
    ProgramStackRec(uint4 qp, uint4 rb, ProgramStackRec* pr) : Q_pos(qp), right_bound(rb), prev(pr) {};
    ProgramStackRec* Set(ProgramStackRec* p, uint4 rb) { prev = p, right_bound = rb; return this; };
    bool isnot_top() { return prev != nullptr; };
  };

  constexpr static const uint4 cut_margin = 24;
  constexpr static const uint4 min_strip_width = 6;// at least 2

  uint4* Q = nullptr;
  uint4* L = nullptr;
  uint4* R = nullptr;

  CFastIntFinder* clone_of = nullptr;

  double avr_segm_on_vline = 0;
  uint4 LR_len = 0;
  uint4 nTotSegm = 0;
  uint4 len_of_Q = 0;
  uint4 L_size = 0;
  uint4 nTotX = 0;
  bool dont_cut_stripe;



  uint4* GetQTail() {
    return Q + len_of_Q;
  }

  uint4 GetDivPow(uint4 l) {
    if (l < 64) return 2;
    uint4 res = 0;
    for (uint4 target = 2.0 * l / MAX(avr_segm_on_vline, 8.0); target; target >>= 1)
      ++res;
    return MAX(res, 4);
  };

  template<class IntersectionFinder,class SegmentsColl>
  static void InsDel(IntersectionFinder &self, SegmentsColl& segments, uint4 end_rank, ProgramStackRec* stack_pos) {
    segments.InsDel(self.L_size, self.L, end_rank);
    if (stack_pos->isnot_top()) {
      // Find intersections for all inserted (first endpoints) and vertical segments
      // in staircase levels above (in the stack).
      for (auto p = segments.GetEndsListFirst(end_rank); p; p = segments.GetEndsListNext(p)) {
        auto sn = SegmentsColl::get_segm(*p);
        segments.SetCurSeg4Bubble(sn);
        self.FindIntI(segments, sn, stack_pos); // get internal intersections
      }
    }

  }



  template <class SegmentsColl>
  static void FindInt(SegmentsColl& segments, uint4* const qb, uint4* const qe, uint4* l)
  {
    // Use range overloads with sentinel instead of single-element overloads.
    // Bottom and top sentinel is placed via CSentinel RAII 
    CSentinel bsentinel(segments, *qb);
    auto c = segments.FindCurSegIntDownWith(l, qb); //first get intersections below

    constexpr const bool line_seg = (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) == _Coll_flag_state::state_true);

    if (line_seg && (c != l))
      return; //if found and segment is line there can't be any more intersections above

    CSentinel tsentinel(segments, *qe, true);
    segments.FindCurSegIntUpWith(l + 1, qe); // get intersections above
  };

  template<class SegmentsColl>
  void SearchInStripNonLineSeg(SegmentsColl& segments, uint4* _L, uint4* _Q)//simplified version for SearchInStrip
  {
    auto last_L = _L + L_size;
    auto first_L = _L + 1;
    auto* _Q_pos = _Q + 1;
    auto T_pos = GetQTail();
    do
    {
      auto new_L_pos = _L;
      *_Q_pos = _L[0];
      if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::needs_SetCurSegCutBE_at_start) == _Coll_flag_state::state_true)
        segments.SetCurSegCutBE(_L[0]);
      auto Q_tail = T_pos;
      auto cur_L = first_L;
      for (CSentinel sentinel(segments, *_Q); cur_L < last_L; ++cur_L)
      {
        segments.SetCurSegCutBE(*cur_L);

        if (_Q_pos == segments.FindCurSegIntDownWith(_Q_pos, _Q))
          *++_Q_pos = *cur_L;
        else
        {
          *new_L_pos++ = *cur_L;
          *--Q_tail = _Q_pos - _Q;
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
      for (cur_L = _L; cur_Q != _Q_pos; ++cur_L)
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
  void SearchInStripLineSeg(SegmentsColl& segments, uint4* L_) {
    auto last = L_ + L_size;
    if constexpr (has_sentinels<SegmentsColl>)
      *last = segments.get_sentinel(true);
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

  template <class SegmentsColl>
  void FindIntI(SegmentsColl& segments, uint4 sn, ProgramStackRec* stack_pos) const
  {
    auto r_index = segments.GetSegR(sn);
    while (stack_pos->right_bound <= r_index)
      stack_pos = stack_pos->prev; // go from bottom to top and find staircase to start
    int4  m, len;
    uint4* Qb, * Ql, * Qe = Q + (stack_pos->Q_pos + 1);
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


  template<class SegmentsColl>
  void SearchInStrip(SegmentsColl &segments, uint4 qp)
  {
    if constexpr(SegmentsColl::get_coll_flag(_Coll_flags::line_segments)==_Coll_flag_state::state_true)
    {
      //For line segments we can do more efficient insertion sorting using intersection check as comparison.
      //If s1<s2 at the left bound then s1 intersects s2 inside the stripe means s1>s2 at the right bound of the stripe.
      //All segments are sorted at the left bound by precondition, so intesection is the same as comparison
      // at the right bound. Byproduct of the sorting is intersections detection.
      SearchInStripLineSeg(segments, L);
    }
    else {
      //We use R as temporary Q making sequential simplified Split
      SearchInStripNonLineSeg(segments, L, R - 1);
      //at the and we have all the intersections found
      //and all the segments moved from L to Q
      //so we move them back
      std::swap(L, R);
      //now segments in L once again, but in icorrect order
      //sort it
      std::sort(L, L + L_size, [&](uint4 s1, uint4 s2) {return segments.RBelow(s1, s2); });
    }
  };

  template <class SegmentsColl>
  void SISFindR(SegmentsColl& segments, uint4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    for (auto i = interval_left_rank + 1; i != interval_right_rank; ++i) {
      if (L_size > 1) {
        segments.SetCurStripe(i - 1, i);
        SearchInStrip(segments, ladder_start_index);
      }
      else
        segments.SetCurStripeRight(i);
      InsDel(*this,segments, i, stack_pos);
    }
    //the stripe right bound needs to be installed correctly even if L_size <= 1
    segments.SetCurStripe(interval_right_rank - 1, interval_right_rank);
    if (L_size > 1) SearchInStrip(segments, ladder_start_index);
  }

  template <class SegmentsColl>
  void FindRNoChecks(SegmentsColl& segments, uint4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    uint4 Q_pos;
    //if L or Q empty cut the strip into 2^divide_pow substrips
    constexpr const uint4 divide_pow = 4;
    if ((L_size == 0)||(0 == (Q_pos = Split(segments, Q + ladder_start_index, interval_right_rank))))
      return MultipleCutting(segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos, divide_pow);
    //normal processing otherwise
    Q_pos += ladder_start_index;
    ProgramStackRec stack_rec(Q_pos, interval_right_rank, stack_pos);// creating new stack list entry for "bubbling up" 
    if (dont_cut_stripe ) { //if found a lot of intersections repeat FindR
      FindRNoChecks(segments, Q_pos, interval_left_rank, interval_right_rank, &stack_rec);
    }else{// cut at the middle into two substrips 
      uint4 middle = (interval_left_rank + interval_right_rank) / 2;
      FindR(segments, Q_pos, interval_left_rank, middle, &stack_rec);
      InsDel(*this,segments, middle, &stack_rec);
      FindR(segments, Q_pos, middle, interval_right_rank, &stack_rec);
      //actually works without SetCurStripeLeft, but it simplifies segment collection class protocol
      //segments.SetCurStripeLeft(interval_left_rank);
    }
    Merge(segments, interval_left_rank, ladder_start_index, Q_pos);
  };

  template <class SegmentsColl>
  void FindR(SegmentsColl& segments, uint4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    if (interval_right_rank < interval_left_rank + min_strip_width) {//end of recursion for narrow stripes
      SISFindR(segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos);
    } else {// normal step
      segments.SetCurStripe(interval_left_rank, interval_right_rank);
      FindRNoChecks(segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos);
    }
  }

  //cuts the strip into 2^divide_pow substrips
  template <class SegmentsColl>
  void MultipleCutting(SegmentsColl& segments, uint4 ladder_start_index, uint4 left_bound, uint4 interval_right_rank, ProgramStackRec* stack_pos, uint4 divide_pow)
  {
    if (interval_right_rank < left_bound + (min_strip_width << divide_pow)) // end of recursion for narrow stripes
      return SISFindR(segments, ladder_start_index, left_bound, interval_right_rank, stack_pos);  //if strip narrow just apply SISFindR

    uint8 step = interval_right_rank - left_bound;
    uint8 rb = (uint8)left_bound << divide_pow;
    uint4 right_bound = (rb += step) >> divide_pow;
    segments.SetCurStripe(left_bound, right_bound);
    FindRNoChecks(segments, ladder_start_index, left_bound, right_bound, stack_pos);
    do {
      left_bound = right_bound;
      right_bound = (rb += step) >> divide_pow;
      InsDel(*this,segments, left_bound, stack_pos);
      segments.SetCurStripe(left_bound, right_bound);
      FindRNoChecks(segments, ladder_start_index, left_bound, right_bound, stack_pos);
    } while (right_bound != interval_right_rank);
  }

  template<class SegmentsColl>
  void Merge(SegmentsColl &segments, uint4 LBoundIdx, uint4 QB, uint4 QE)
  {
    const auto Size = L_size;
    uint4 cur_R_pos = 0, new_size = 0;
    auto cur_stair = Q + QB + 1;
    uint4 * const __restrict _R = L;
    uint4 * const __restrict _L = R;
    auto cur_seg = _R[cur_R_pos];
    auto const bot_Q = Q + QB;
    auto const top_Q = Q + QE + 1;
    CSentinel sentinel(segments, *bot_Q);
    if constexpr (has_sentinels<SegmentsColl>)
      *top_Q = segments.get_sentinel(true);//we don't need to restore top_Q just place sentinel
    while ((cur_stair!= top_Q) && (cur_R_pos<Size))
    {
      if (segments.RBelow(cur_seg, *cur_stair))
      {
        if (segments.GetSegL(cur_seg) > LBoundIdx)
        {
          segments.SetCurSegCutEnd(cur_seg);
          {
            auto c = segments.FindCurSegIntUpWith(cur_stair, top_Q); 
            constexpr const bool non_line_seg = (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) != _Coll_flag_state::state_true);
            if (non_line_seg || (c == cur_stair))
              segments.FindCurSegIntDownWith(cur_stair-1, bot_Q);
          };
        }
        _L[new_size++] = cur_seg;
        cur_seg = _R[++cur_R_pos];
      }
      else
        _L[new_size++] = *cur_stair++;
    }
    for (auto last_Q = top_Q - 1; cur_R_pos < Size;)
    {
      if (segments.GetSegL(cur_seg) > LBoundIdx)
      {
        segments.SetCurSegCutEnd(cur_seg);
        segments.FindCurSegIntDownWith(last_Q,bot_Q);
      }
      _L[new_size++] = cur_seg;
      cur_seg = _R[++cur_R_pos];
    }
    while (cur_stair!=top_Q)
      _L[new_size++] = *cur_stair++;
    L = _L; R = _R;
    L_size = new_size;
  };

  template<class SegmentsColl>
  uint4 Split4LineSeg(SegmentsColl& segments, uint4* _Q, uint4 RBoundIdx)
  {
    auto R_pos = R;
    auto new_L_pos = L;
    auto Q_tail = GetQTail();
    assert(Q_tail > _Q + L_size);
    auto const last_L = L + L_size;
 //at first, place all lowest segments not covering current stripe to R (and implicitly L)
    while ((new_L_pos < last_L) && (segments.GetSegR(*new_L_pos) < RBoundIdx)) {
            *(R_pos++) = *(new_L_pos++);//place segment in  R
            //storing segment position in Q_tail 
  //we have the empty ledder so current position is 0 and position above is 1
            *(--Q_tail) = 0;  //it should be position above i.e. 1. We save 
            // one addition per loop by incrementing _Q later.
    }
    if (new_L_pos == last_L){
      return 0;
    }

    long long n_int = 0;
    auto  _Q_pos = _Q;
    //first segment covering current stripe we place to Q 
    //it can't intersect any of the steps(stairs) because there no stairs yet 
    *++_Q_pos = *new_L_pos;
    if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::needs_SetCurSegCutBE_at_start) == _Coll_flag_state::state_true)
      segments.SetCurSegCutBE(*new_L_pos);//we don't need SetCurSegCutBE in all cases, only for Y caching collections 
    auto cur = new_L_pos + 1;
    for (CSentinel sentinel(segments, *_Q); cur < last_L; ++cur) {
      auto cur_seg = *cur;
      segments.SetCurSegCutBE(cur_seg);
      auto k= _Q_pos - segments.FindCurSegIntDownWith(_Q_pos, _Q);
      if (k) {//segment  intersects upper ladder stair
        // add number of intersections with stairs to counter
        n_int += k;//increment found int number
        *(new_L_pos++) = cur_seg;//place segment in L
        continue;
      }

      //segment doesn't intersect the ladder stairs
      if (segments.GetSegR(cur_seg) >= RBoundIdx) {//segment is covering current stripe 
        *++_Q_pos = cur_seg; //place segment in Q
        continue;
      }
      //segment is not covering current stripe 
      *(R_pos++) = *(new_L_pos++) = cur_seg; //place segment in L and R
      //storing segment position in Q_tail
      *(--Q_tail) = _Q_pos - _Q;  //it should be _Q_pos+1. We save 
      // one addition per loop by incrementing _Q later.
    }
    L_size = new_L_pos-L;
    uint4 Q_size=_Q_pos - _Q;
    *--Q_tail = Q_size;//place a guard for the loop below

    // important to start from stair above current segm, meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we incremement _Q and _Q[*Q_tail] become stair above
    ++_Q_pos;

    if constexpr (has_sentinels<SegmentsColl>)//placing sentinel if collection supports
      *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel
    
    Q_tail = GetQTail() - 1;
    cur = R;
    for (auto Q_idx = *Q_tail; Q_idx != Q_size; Q_idx = *--Q_tail, ++cur)
    {
      segments.SetCurSegCutBE(*cur);//SetCurSegCutBeg
//      auto Q_loc = _Q + Q_idx;
//      n_int += segments.FindCurSegIntUpWith(Q_loc, _Q_pos) - Q_loc;
      segments.FindCurSegIntUpWith(_Q + Q_idx, _Q_pos);
    }
    dont_cut_stripe = n_int > L_size + cut_margin;
    return Q_size;
  };

  template<class SegmentsColl>
  uint4 Split4NonLineSeg(SegmentsColl& segments, uint4* _Q, uint4 RBoundIdx)
  {
    auto  _Q_pos = _Q;
    long long n_int = 0;
    auto Q_tail = GetQTail();
    auto new_L_pos = L;
    CSentinel sentinel(segments, *_Q);
    for (auto L_pos = L, last_L = L + L_size; L_pos < last_L; ++L_pos) {
      auto cur_seg = *L_pos;
      segments.SetCurSegCutBE(cur_seg);
      auto cur_Q = segments.FindCurSegIntDownWith(_Q_pos, _Q);
      if ((cur_Q == _Q_pos) && (segments.GetSegR(cur_seg) >= RBoundIdx)) //segment doesn't intersect ladder stairs
        //and is covering current stripe and doesn't intersect ladder stairs
          *++_Q_pos = cur_seg;
      else {   //segment  intersects ladder stairs or is not covering 
        n_int += _Q_pos - cur_Q; //inc int counter
        *(new_L_pos++)= cur_seg; //place segment in L
        //storing segment position in Q_tail
        *(--Q_tail) = _Q_pos -_Q;  //it should be _Q_pos+1. We save 
        // one addition per loop by incrementing _Q later.
      }
    }
    if (_Q_pos == _Q) {
      return 0;
    }
    L_size = new_L_pos - L;
    uint4 Q_size = _Q_pos - _Q;
    *--Q_tail = Q_size;//place a guard for the loop below

    // important to start from stair above current segm, meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we incremement _Q and _Q[*Q_tail] become stair above
    ++_Q_pos;
    if constexpr (has_sentinels<SegmentsColl>)
      *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel
    
    Q_tail = GetQTail() - 1;
    auto L_pos = L;
    for (auto Q_idx =*Q_tail; Q_idx != Q_size; Q_idx = *--Q_tail, ++L_pos)
    {
      segments.SetCurSegCutBE(*L_pos);
      segments.FindCurSegIntUpWith(_Q + Q_idx, _Q_pos);
    }
    dont_cut_stripe = n_int > L_size + cut_margin;
    return Q_size;
  };


  template<class SegmentsColl>
  uint4 Split(SegmentsColl& segments, uint4* _Q, uint4 RBoundIdx)
  {
    if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::line_segments)==_Coll_flag_state::state_true)
      return Split4LineSeg(segments, _Q, RBoundIdx);
    return Split4NonLineSeg(segments, _Q, RBoundIdx);
  };


 
  void clone(CFastIntFinder* master)
  {
      nTotSegm = master->nTotSegm;
      LR_len = master->LR_len;
      len_of_Q = LR_len;
      avr_segm_on_vline = master->avr_segm_on_vline;
      nTotX = master->nTotX;
      clone_of = master;
  };

  void unclone()
  {
      if (clone_of != nullptr) {
          clone_of = nullptr;
      }
  };



 

  template<class SegmentsColl>
  uint4 CalcLAt(SegmentsColl& segments, uint4 end_rank)
  {
      uint4 Size = 0;
      for (uint4 cur_seg = 0; cur_seg < nTotSegm; ++cur_seg)
          if ((segments.GetSegL(cur_seg) <= end_rank) && (segments.GetSegR(cur_seg) > end_rank))
              L[Size++] = cur_seg;
      segments.SortAt(end_rank, Size, L);
      return Size;
  };

};

