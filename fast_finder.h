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

class CFastIntFinder : public CommonImpl
{

  constexpr static int4 bottom_step_index = 0;

  template<typename SegmentsColl>
  struct has_get_sentinel_idx
  {
  private:
    template<typename SegColl> static auto test() -> decltype(std::declval<SegColl>().get_sentinel_idx(true) == 1, std::true_type());

    template<typename> static std::false_type test(...);

  public:
    static constexpr bool value = std::is_same<decltype(test<SegmentsColl>()), std::true_type>::value;
  };

  template< class T >
  constexpr static bool has_get_sentinel_idx_member = has_get_sentinel_idx<T>::value;
  
  template<class SegmentsColl, bool HasSentinel = has_get_sentinel_idx_member<SegmentsColl> >
  struct CSentinel {
    CSentinel(SegmentsColl& coll, int4& ini) {};
  };

  template<class SegmentsColl>
  struct CSentinel <SegmentsColl, true> {
    CSentinel(SegmentsColl& coll, int4& ini) : to_restore(ini) {
      prev_val = to_restore;
      to_restore = coll.get_sentinel_idx(false);
    };
    ~CSentinel() {
      to_restore = prev_val;
    };
  private:
    int4& to_restore;
    int4 prev_val;
  };

public:
  using CTHIS = CFastIntFinder;
  using CIMP = CommonImpl;
  using CIMP::dont_split_stripe;
  using CIMP::SegL; using CIMP::SegR; using CIMP::ENDS; using CIMP::Q;
  using CIMP::prepare_ends; using CIMP::FindInt; using CIMP::FindIntI;
  using CIMP::SearchInStripNonLineSeg; using CIMP::SearchInStripLineSeg;

  ~CFastIntFinder() { unclone(); FreeMem(); };

  template<class SegmentsColl>
  void find_intersections(SegmentsColl &segments)
  {
    AllocMem(segments);
    ProgramStackRec stack_rec(bottom_step_index, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl::get_segm(ENDS[0]);
    L_size = 1;
    FindR(*this, segments, bottom_step_index, 0, 2 * nTotSegm - 1, &stack_rec, 0, get_max_call(2 * nTotSegm));
    FreeMem();
  }

/*  template<class SegmentsColl, class CIntRegistrator >
  void find_intersections(SegmentsColl& segments, uint4 n_threads, CIntRegistrator* regs[])*/
  template<template <class> class SegmentsColl, class CIntRegistrator >
  void find_intersections(SegmentsColl<CIntRegistrator>& segments, uint4 n_threads, CIntRegistrator* regs[])
  {
    AllocMem(segments);
    using namespace std;
    vector<thread> wrk_threads;
    auto thread_func = [](CTHIS *master, SegmentsColl<CIntRegistrator> *segments, uint4 from, uint4 to,uint4 max_call, CIntRegistrator* add_reg) {
      CTHIS i_f;
      SegmentsColl<CIntRegistrator> coll;
      i_f.clone(master);
      coll.clone(*segments, add_reg);
      ProgramStackRec psr(bottom_step_index, 2 * i_f.nTotSegm);
      i_f.L_size = i_f.CalcLAt(coll, from);
      FindR(i_f,coll, bottom_step_index, from, to, &psr, 0,max_call);
      coll.unclone();
      i_f.unclone();
    };
    auto n = segments.GetSegmNumb();

    double part = 2 * n / (double)n_threads;
    uint4 max_call = get_max_call(part);
    int4 i = 1;
    uint4 start_from = part;
    uint4 from, to = start_from;
    while (i < n_threads) {
      from = to;
      i++;
      to = (i == n_threads) ? 2 * n - 1 : part * i;
      wrk_threads.emplace_back(thread_func, this, &segments, from, to, max_call,regs[i - 2]); // starts intersection finding in a stripe <from,to>
    }

    ProgramStackRec stack_rec(bottom_step_index, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl<CIntRegistrator>::get_segm(ENDS[0]);
    L_size = 1;
    FindR(*this,segments, bottom_step_index, 0, start_from, &stack_rec,  0,max_call);
    for (auto cur_thread = wrk_threads.begin(); cur_thread != wrk_threads.end(); cur_thread++)
      cur_thread->join(); //waiting for calculation of all threads are finished
    FreeMem();
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
      std::sort(L, L + L_size, [segments](int4 s1, int4 s2) {return segments.RBelow(s1, s2); });
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
      segments.SetCurPointAtBeg(sn);
      for (i = Size - 1; (i > -1) && (!segments.UnderCurPoint(L[i])); --i)
        L[i + 1] = L[i];
      L[i + 1] = sn;
      Size++;
      segments.SetCurSegAE(sn);
      FindIntI(segments, SegR[sn], stack_pos);// get internal intersections
    }
    L_size = Size;
  }


  template<class SegmentsColl>
  void Merge(SegmentsColl &segments, uint4 LBoundIdx, int4 QB, int4 QE)
  {
    auto Size = L_size;
    int4 cur_R_pos = 0, new_size = 0;
    auto _R = L, _L = R;
    int4 cur_seg = _R[cur_R_pos];
    auto bot_Q = Q + QB;
    auto top_Q = Q + QE + 1;
    auto cur_stair = bot_Q + 1;
    CSentinel sentinel(segments, *bot_Q);
    if constexpr (has_get_sentinel_idx_member<SegmentsColl>)
      *top_Q = segments.get_sentinel_idx(true);

    while ((cur_stair<top_Q) && (cur_R_pos<Size))
    {
      if (segments.RBelow(cur_seg, *cur_stair))
      {
        if (SegL[cur_seg]>LBoundIdx)
        {
          segments.SetCurSegCutEnd(cur_seg);
          auto c = cur_stair - 1;
          while ((has_get_sentinel_idx_member<SegmentsColl> || (c != bot_Q)) &&
            (segments.FindCurSegIntDownWith(*c)))
              --c; //first get intersections below
          if (!SegmentsColl::is_line_segments || (c == cur_stair - 1)) {//if found and segment is line it can't be any more
            c = cur_stair;
            while ((has_get_sentinel_idx_member<SegmentsColl>||(c != top_Q)) && 
              (segments.FindCurSegIntUpWith(*c))) // get intersections above
                ++c;
          }
        }
        _L[new_size++] = cur_seg;
        cur_seg = _R[++cur_R_pos];
      }
      else
        _L[new_size++] = *(cur_stair++);
    }
    while (cur_R_pos<Size)
    {
      if (SegL[cur_seg]>LBoundIdx)
      {
        segments.SetCurSegCutEnd(cur_seg);
        for (auto c = top_Q; (has_get_sentinel_idx_member<SegmentsColl> || (c > bot_Q))
          && segments.FindCurSegIntDownWith(*--c);); //get intersections below
      }
      _L[new_size++] = cur_seg;
      cur_seg = _R[++cur_R_pos];
    }
    while (cur_stair<top_Q)
      _L[new_size++] = *(cur_stair++);
    L = _L; R = _R;
    L_size = new_size;
  };

  template<class SegmentsColl>
  int4 Split4LineSeg(SegmentsColl& segments, int4* _Q, uint4 RBoundIdx)
  {
    auto R_pos = R;
    auto new_L_pos = L;
    auto Q_tail = Q + len_of_Q;
    auto last_L = L + L_size;
 //at first, place all lowest segments not covering current stripe to R (and implicitly L)
    while ((new_L_pos < last_L) && (SegR[*new_L_pos] < RBoundIdx)) {
            *(R_pos++) = *(new_L_pos++);//place segment in  R
            //storing segment position in Q_tail 
  //we have the empty ledder so current position is 0 and position above is 1
            *(--Q_tail) = 0;  //it should be position above i.e. 1. We save 
            // one addition per loop by incrementing _Q later.
    }

    long long n_int = 0;
    auto  _Q_pos = _Q;
    if (new_L_pos < last_L) {
        //first segment covering current stripe we place to Q 
        //it can't intersect any of the steps(stairs) because there no stairs yet     
        *++_Q_pos = *new_L_pos;
        CSentinel<SegmentsColl> sentinel(segments, *_Q);
        for (auto cur_L_seg = new_L_pos + 1; cur_L_seg < last_L; ++cur_L_seg) {
            segments.SetCurSegCutBE(*cur_L_seg);
            auto cur_Q = _Q_pos;
            if (segments.FindCurSegIntDownWith(*(cur_Q--))){//segment  intersects upper ladder stair
                // finding another ledder intersections
                while ((has_get_sentinel_idx_member<SegmentsColl>||(cur_Q != _Q))
                  && (segments.FindCurSegIntDownWith(*cur_Q)))
                    --cur_Q;
                n_int += _Q_pos-cur_Q;//increment found int number
 
                *(new_L_pos++) = *cur_L_seg;//place segment in L
              continue;
            }
 
            //segment doesn't intersect the ladder stairs
            if (SegR[*cur_L_seg] >= RBoundIdx) {//segment is covering current stripe 
                *++_Q_pos = *cur_L_seg; //place segment in Q
                continue;
            }
            //segment is not covering current stripe 
            *(R_pos++) = *(new_L_pos++) = *cur_L_seg; //place segment in L and R
            //storing segment position in Q_tail
            *(--Q_tail) = _Q_pos - _Q;  //it should be _Q_pos+1. We save 
            // one addition per loop by incrementing _Q later.
        }
    }
    L_size = new_L_pos-L;
    if (_Q_pos == _Q) {
      dont_split_stripe = false;
      return 0;
    }
    // important to start from stair above current segm, meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we incremement _Q and _Q[*Q_tail] become stair above
    ++_Q_pos;
    if constexpr(has_get_sentinel_idx_member<SegmentsColl>)
      *_Q_pos=segments.get_sentinel_idx(true);
    while (R_pos>R)
    {
      segments.SetCurSegCutBeg(*--R_pos);
      auto cur_Q = _Q + *Q_tail++;          // getting position stored in tail of Q;
      while ((has_get_sentinel_idx_member<SegmentsColl>||(cur_Q != _Q_pos)) &&
        (segments.FindCurSegIntUpWith(*cur_Q)))
          ++cur_Q;
    } 
    dont_split_stripe = n_int > L_size;
    return _Q_pos-_Q;
  };

  template<class SegmentsColl>
  int4 Split4NonLineSeg(SegmentsColl& segments, int4* _Q, uint4 RBoundIdx)
  {
    uint4  _Q_pos = 0;
    long long n_int = 0;
    auto Q_tail = Q + len_of_Q;
    auto new_L_pos = L;
    for (auto L_pos = L, last_L = L + L_size; L_pos < last_L; ++L_pos) {
      auto cur_seg = *L_pos;
      auto step = _Q_pos;
      segments.SetCurSegCutBE(cur_seg);
      for (auto cur_Q = _Q + step; (step != 0) && (segments.FindCurSegIntDownWith(*cur_Q)); --step, --cur_Q);

      if ((_Q_pos == step) && (SegR[cur_seg] >= RBoundIdx)) //segment doesn't intersect ladder stairs
        //and is covering current stripe and doesn't intersect ladder stairs
          _Q[++_Q_pos] = cur_seg;
      else {   //segment  intersects ladder stairs or is not covering 
        n_int += _Q_pos - step; //inc int counter
        *(new_L_pos++)= cur_seg; //place segment in L
        //storing segment position in Q_tail
        *(--Q_tail) = _Q_pos;  //it should be _Q_pos+1. We save 
        // one addition per loop by incrementing _Q later.
      }
    }
    L_size = new_L_pos - L;
    if (_Q_pos == 0){
      dont_split_stripe = false;
      return 0;
    }
    Q_tail = Q + len_of_Q - 1;
    // important to start from stair above current segm, meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we incremement _Q and _Q[*Q_tail] become stair above
    for (auto L_pos=L, last_Q = _Q + _Q_pos; L_pos < new_L_pos ; ++L_pos, --Q_tail)
    {
      segments.SetCurSegCutBE(*L_pos);
      auto cur_Q = _Q + *Q_tail;          // getting position stored in tail of Q;
      while ((cur_Q < last_Q) && (segments.FindCurSegIntUpWith(*cur_Q)))++cur_Q;
    }
    dont_split_stripe = n_int > L_size;
    return _Q_pos;
  };


  template<class SegmentsColl>
  int4 Split(SegmentsColl& segments, int4* _Q, uint4 RBoundIdx)
  {
    if constexpr (SegmentsColl::is_line_segments)
      return Split4LineSeg(segments, _Q, RBoundIdx);
    return Split4NonLineSeg(segments, _Q, RBoundIdx);
  };

 
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
      L = new int4[LR_len];
      Q = new int4[len_of_Q];
      R = new int4[LR_len];
      clone_of = master;
  };

  void unclone()
  {
      if (clone_of != nullptr) {
          SegL = nullptr;
          SegR = nullptr;
          ENDS = nullptr;
          FreeMem();
          clone_of = nullptr;
      }
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

  template<class SegmentsColl>
  void AllocMem(SegmentsColl& segments)
  {
      assert(nTotSegm == segments.GetSegmNumb());
      len_of_Q = LR_len;
      L = new int4[LR_len];
      R = new int4[LR_len];
      Q = new int4[len_of_Q];
  };
public:
  void FreeMem()
  {
      MY_FREE_ARR_MACRO(L);
      MY_FREE_ARR_MACRO(R);
      MY_FREE_ARR_MACRO(Q);
      //CIMP::FreeMem();
  };

};

