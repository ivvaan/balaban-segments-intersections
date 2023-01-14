#pragma once
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

#include "new_int_finder.h"
#include <algorithm>
#include <vector>
#include <thread>
#include <cassert>
#include <type_traits>

class CFastIntFinder : public CommonImpl
{
public:
  using CTHIS = CFastIntFinder;
  using CIMP = CommonImpl;
 
  ~CFastIntFinder() { unclone(); FreeMem(); };

  template<class SegmentsColl>
  void find_intersections(SegmentsColl  &segments,uint4 from=0,uint4 to=0)
  {
    //AllocMem
    len_of_Q = LR_len;
    DECL_RAII_ARR(L, LR_len);
    ++L;//to have one cell before L for sentinel
    DECL_RAII_ARR(R, LR_len);
    ++R;//to have one cell before R for sentinel
    DECL_RAII_ARR(Q, len_of_Q);

    if (to == 0) {
      from = 0;
      to = 2 * nTotSegm - 1;
    }
    if (from == 0) {
      L[0] = SegmentsColl::get_segm(ENDS[0]);
      L_size = 1;
    }
    else
      L_size = CalcLAt(segments, from);

    constexpr int4 bottom_index = 0;
    ProgramStackRec stack_rec(bottom_index, 2 * nTotSegm); //need to be initialized this way
    FindR(*this, segments, bottom_index, from, to, &stack_rec, 0, get_max_call(to-from));
  }

 template<template <class> class SegmentsColl, class CIntRegistrator >
  void find_intersections(SegmentsColl<CIntRegistrator>& segments, uint4 n_threads, CIntRegistrator* regs[])
  {
    using namespace std;
    vector<thread> wrk_threads;
    auto thread_func = [](CTHIS *master, SegmentsColl<CIntRegistrator> *segments, uint4 from, uint4 to, CIntRegistrator* add_reg) {
      SegmentsColl<CIntRegistrator> coll(*segments, add_reg);
      CTHIS(master).find_intersections(coll, from, to);
    };
    auto n = segments.GetSegmNumb();

    double part = 2 * n / (double)n_threads;
    uint4 start_from = part;
    for (uint4 i = 2, from = start_from, to; i <= n_threads; ++i) {
      to = (i == n_threads) ? 2 * n - 1 : (uint4)(part * i);
      wrk_threads.emplace_back(thread_func, this, &segments, from, to, regs[i - 2]); // starts intersection finding in a stripe <from,to>
      from = to;
    }

    find_intersections(segments, 0, start_from);
    for (auto& cur_thread:wrk_threads)
      cur_thread.join(); //waiting for calculation of all threads are finished
  }

  template<class SegmentsColl>
  void SearchInStrip(SegmentsColl &segments, int4 qp)
  {
    if constexpr(SegmentsColl::is_line_segments)
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
      std::sort(L, L + L_size, [&](int4 s1, int4 s2) {return segments.RBelow(s1, s2); });
    }
  };

  template<class SegmentsColl>
  void InsDel(SegmentsColl &segments, uint4 end_index, ProgramStackRec * stack_pos)
  {
    int4 Size = L_size;
    int4 i;
    auto sn = SegmentsColl::get_segm(ENDS[end_index]);
    if (SegmentsColl::is_last(ENDS[end_index])) // if endpoint - delete
    {
      i = --Size;
      auto cur = L[i];
      while (cur != sn)
      {
        --i;
        auto buf = cur;
        cur = L[i];
        L[i] = buf;
      }
    }
    else// if startpoint - insert
    {
      segments.SetCurSegAndPoint(sn);
      for (i = Size - 1; (i > -1) && (!segments.UnderCurPoint(L[i])); --i)
        L[i + 1] = L[i];
      L[i + 1] = sn;
      Size++;
      FindIntI(segments, SegR[sn], stack_pos);// get internal intersections
    }
    L_size = Size;
  }

  template<class SegmentsColl>
  void Merge(SegmentsColl &segments, uint4 LBoundIdx, int4 QB, int4 QE)
  {
    const auto Size = L_size;
    uint4 cur_R_pos = 0, new_size = 0;
    auto cur_stair = Q + QB + 1;
    int4 * const __restrict _R = L;
    int4 * const __restrict _L = R;
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
        if (SegL[cur_seg]>LBoundIdx)
        {
          segments.SetCurSegCutEnd(cur_seg);
          {
            auto l = cur_stair - 1;
            auto c = segments.FindCurSegIntDownWith(l, bot_Q);
            if ((!SegmentsColl::is_line_segments || (c == l)))
              segments.FindCurSegIntUpWith(cur_stair, top_Q);
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
      if (SegL[cur_seg]>LBoundIdx)
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
  int4 Split4LineSeg(SegmentsColl& segments, int4* _Q, uint4 RBoundIdx)
  {
    auto R_pos = R;
    auto new_L_pos = L;
    auto Q_tail = Q + len_of_Q;
    assert(Q_tail > _Q + L_size);
    auto const last_L = L + L_size;
 //at first, place all lowest segments not covering current stripe to R (and implicitly L)
    while ((new_L_pos < last_L) && (SegR[*new_L_pos] < RBoundIdx)) {
            *(R_pos++) = *(new_L_pos++);//place segment in  R
            //storing segment position in Q_tail 
  //we have the empty ledder so current position is 0 and position above is 1
            *(--Q_tail) = 0;  //it should be position above i.e. 1. We save 
            // one addition per loop by incrementing _Q later.
    }
    if (new_L_pos == last_L){
      dont_split_stripe = false;
      return 0;
    }

    long long n_int = 0;
    auto  _Q_pos = _Q;
    //first segment covering current stripe we place to Q 
    //it can't intersect any of the steps(stairs) because there no stairs yet     
    segments.SetCurSegCutBE(*++_Q_pos = *new_L_pos);//we don't need SetCurSegCutBE but it is allow to speedup by y values caching 
    auto cur_L = new_L_pos + 1;
    for (CSentinel sentinel(segments, *_Q); cur_L < last_L; ++cur_L) {
      auto cur_seg = *cur_L;
      segments.SetCurSegCutBE(cur_seg);
      if (segments.FindCurSegIntDownWith(*_Q_pos)) {//segment  intersects upper ladder stair
          // finding another ledder intersections
        n_int += _Q_pos - segments.FindCurSegIntDownWith(_Q_pos - 1, _Q);//increment found int number
        *(new_L_pos++) = cur_seg;//place segment in L
        continue;
      }

      //segment doesn't intersect the ladder stairs
      if (SegR[cur_seg] >= RBoundIdx) {//segment is covering current stripe 
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
    // important to start from stair above current segm, meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we incremement _Q and _Q[*Q_tail] become stair above
    ++_Q_pos;
    if constexpr (has_sentinels<SegmentsColl>)
      *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel
    while (R_pos!=R)
    {
      --R_pos;
      segments.SetCurSegCutBeg(*R_pos);
      segments.FindCurSegIntUpWith(_Q + *Q_tail, _Q_pos);
      ++Q_tail;
    } 
    dont_split_stripe = n_int > L_size;
    return _Q_pos-_Q;
  };

  template<class SegmentsColl>
  int4 Split4NonLineSeg(SegmentsColl& segments, int4* _Q, uint4 RBoundIdx)
  {
    auto  _Q_pos = _Q;
    long long n_int = 0;
    auto Q_tail = Q + len_of_Q;
    auto new_L_pos = L;
    CSentinel sentinel(segments, *_Q);
    for (auto L_pos = L, last_L = L + L_size; L_pos < last_L; ++L_pos) {
      auto cur_seg = *L_pos;
      segments.SetCurSegCutBE(cur_seg);
      auto cur_Q = segments.FindCurSegIntDownWith(_Q_pos, _Q);
      if ((cur_Q == _Q_pos) && (SegR[cur_seg] >= RBoundIdx)) //segment doesn't intersect ladder stairs
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
      dont_split_stripe = false;
      return 0;
    }
    L_size = new_L_pos - L;
    Q_tail = Q + len_of_Q - 1;
    // important to start from stair above current segm, meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we incremement _Q and _Q[*Q_tail] become stair above
    ++_Q_pos;
    if constexpr (has_sentinels<SegmentsColl>)
      *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel
    for (auto L_pos=L; L_pos != new_L_pos ; ++L_pos, --Q_tail)
    {
      segments.SetCurSegCutBE(*L_pos);
      segments.FindCurSegIntUpWith(_Q + *Q_tail, _Q_pos);
    }
    dont_split_stripe = n_int > L_size;
    return _Q_pos-_Q;
  };


  template<class SegmentsColl>
  int4 Split(SegmentsColl& segments, int4* _Q, uint4 RBoundIdx)
  {
    if constexpr (SegmentsColl::is_line_segments)
      return Split4LineSeg(segments, _Q, RBoundIdx);
    return Split4NonLineSeg(segments, _Q, RBoundIdx);
  };

  CFastIntFinder() {};

  protected:
  int4 *R = nullptr;
  CTHIS* clone_of = nullptr;
 
  void clone(CTHIS* master)
  {
      nTotSegm = master->nTotSegm;
      LR_len = master->LR_len;
      len_of_Q = LR_len;
      SegL = master->SegL;
      SegR = master->SegR;
      ENDS = master->ENDS;
      clone_of = master;
  };

  void unclone()
  {
      if (clone_of != nullptr) {
          SegL = nullptr;
          SegR = nullptr;
          ENDS = nullptr;
          clone_of = nullptr;
      }
  };


  CFastIntFinder(CTHIS* master) {
    clone(master);
  };
 

  template<class SegmentsColl>
  uint4 CalcLAt(SegmentsColl& segments, uint4 end_index)
  {
      uint4 i, Size = 0;
      for (i = 0; i < nTotSegm; i++)
          if ((SegL[i] <= end_index) && (SegR[i] > end_index))
              L[Size++] = i;
      segments.SortAt(ENDS[end_index], Size, L);
      return Size;
  };

public:

};

