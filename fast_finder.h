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


class CFastIntFinder : public CommonImpl
{
public:
  using CTHIS = CFastIntFinder;
  using  CIMP = CommonImpl;
  using CIMP::max_call;using CIMP::dont_split_stripe;
  using CIMP::SegL; using CIMP::SegR; using CIMP::ENDS; using CIMP::Q;
  using CIMP::prepare_ends; using CIMP::FindInt; using CIMP::FindIntI;

  ~CFastIntFinder() { FreeMem(); };

  template<class SegmentsColl>
  void find_intersections(SegmentsColl *segments)
  {
    AllocMem(segments);
    ProgramStackRec stack_rec(-1, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl::get_segm(ENDS[0]);

    FindR(this, segments, -1, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0, get_max_call(2 * nTotSegm));
    //FreeMem();
  }

  template<class SegmentsColl, class CIntRegistrator >
  void find_intersections(SegmentsColl* segments, uint4 n_threads, CIntRegistrator* regs[])
  {
    AllocMem(segments);
    prepare_ends(segments);
    using namespace std;
    vector<thread> wrk_threads;
    auto thread_func = [](CTHIS* master, SegmentsColl* segments, uint4 from, uint4 to,uint4 max_call, CIntRegistrator* add_reg) {
      CTHIS i_f;
      SegmentsColl coll;
      i_f.clone(master);
      coll.clone(segments, add_reg);
      FindR(&i_f,&coll, -1, from, to, &ProgramStackRec(-1, 2 * i_f.nTotSegm), i_f.CalcLAt(&coll, from), 0,max_call);
      coll.unclone();
      i_f.unclone();
    };
    auto n = segments->GetSegmNumb();

    double part = 2 * n / (double)n_threads;
    uint4 max_call = get_max_call(part);
    int4 i = 1;
    uint4 start_from = part;
    uint4 from, to = start_from;
    while (i < n_threads) {
      from = to;
      i++;
      to = (i == n_threads) ? 2 * n - 1 : part * i;
      wrk_threads.emplace_back(thread_func, this, segments, from, to, max_call,regs[i - 2]); // starts intersection finding in a stripe <from,to>
    }

    ProgramStackRec stack_rec(-1, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl::get_segm(ENDS[0]);
    FindR(this,segments, -1, 0, start_from, &stack_rec, 1, 0,max_call);
    //      FindR(segments, -1, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0);
    for (auto cur_thread = wrk_threads.begin(); cur_thread != wrk_threads.end(); cur_thread++)
      cur_thread->join(); //waiting for calculation of all threads are finished
    FreeMem();
  }

  template<class SegmentsColl>
  int4 SearchInStrip(SegmentsColl *segments, int4 qp, int4 Size)
  {
    if (SegmentsColl::is_line_segments)
    {
      //For line segments we can do more efficient insertion sorting using intersection check as comparison.
      //If s1<s2 at the left bound then s1 intersects s2 inside the stripe means s1>s2 at the right bound of the stripe.
      //All segments are sorted at the left bound by precondition, so intesection is the same as comparison
      // at the right bound. Byproduct of the sorting is intersections detection.
      auto L_ = L;
      auto _L = L_ - 1;
      for (uint4 i = 1; i < Size; i++)
      {
        auto sn = L_[i];
        segments->SetCurSegCutBE(sn);
        if (segments->FindCurSegIntWith(_L[i])) {
          L_[i] = _L[i];
          uint4 j = i - 1;
          for (; (j) && (segments->FindCurSegIntWith(_L[j])); --j)
            L_[j] = _L[j];
          L_[j] = sn;
        }
      }
    }
    else
    {
      auto _L = L, q = Q + qp + 1, last = L + Size;
      int4 n_split = 0, size = Size;
      while (size = SplitSIS(segments, qp, size))++n_split; //it change QP must be after q = Q + QP + 1 
                                                            // at this point we can just place Q starting from QP+1 to L and sort it
      for (; _L < last; ++_L, ++q)  *_L = *q;
      if (n_split)std::sort(L, L + Size, [segments](int4 s1, int4 s2) {return segments->RBelow(s1, s2); });
    }
    return Size;
  };
  template<class SegmentsColl>
  int4 InsDel(SegmentsColl *segments, uint4 end_index, ProgramStackRec * stack_pos, int4 Size)
  {
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
      segments->SetCurPointAtBeg(sn);
      for (i = Size - 1; (i > -1) && (!segments->UnderCurPoint(L[i])); --i)
        L[i + 1] = L[i];
      L[i + 1] = sn;
      Size++;
      segments->SetCurSeg(sn);
      FindIntI(segments, SegR[sn], stack_pos);// get internal intersections
    }
    return Size;
  }

  template<class SegmentsColl>
  int4 Merge(SegmentsColl *segments, uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
  {
    int4 cur_R_pos = 0, new_size = 0, cur_stair = QB;
    auto _R = L, _L = R;
    int4 cur_seg = _R[cur_R_pos];
    while ((cur_stair<QE) && (cur_R_pos<Size))
    {

      if (segments->RBelow(cur_seg, Q[cur_stair + 1]))
      {
        if (SegL[cur_seg]>LBoundIdx)
        {
          segments->SetCurSegCutEnd(cur_seg);
          FindInt(segments, QB, QE, cur_stair);
        }
        _L[new_size++] = cur_seg;
        cur_seg = _R[++cur_R_pos];
      }
      else
        _L[new_size++] = Q[++cur_stair];
    }
    while (cur_R_pos<Size)
    {
      if (SegL[cur_seg]>LBoundIdx)
      {
        segments->SetCurSegCutEnd(cur_seg);
        FindInt(segments, QB, QE, QE);
      }
      _L[new_size++] = cur_seg;
      cur_seg = _R[++cur_R_pos];
    }
    while (cur_stair<QE)
      _L[new_size++] = Q[++cur_stair];
    L = _L; R = _R;
    return new_size;
  };
  template<class SegmentsColl>
  int4 Split(SegmentsColl* segments, uint4 RBoundIdx, int4& _step_index, int4 Size)
  {
    auto _Q = Q + _step_index;
    uint4 new_L_size = 0, _Q_pos = 0;
    auto Q_tail = Q + len_of_Q;
    int4 *R_pos;
    long long n_int=0;
    if (SegmentsColl::is_line_segments)R_pos = R;
    for (int4 cur_L_pos = 0; cur_L_pos<Size; cur_L_pos++)
    {
      int4 cur_seg = L[cur_L_pos];
      int4 step = _Q_pos;
      segments->SetCurSegCutBE(cur_seg);
      for (auto cur_Q = _Q + step; (step!=0) && (segments->FindCurSegIntWith(*cur_Q)); --step, --cur_Q);
        

      if (_Q_pos == step)//segment doesn't intersect ladder stairs
      {
        if (SegR[cur_seg] >= RBoundIdx)//segment is covering current stripe and doesn't intersect ladder stairs
          _Q[++_Q_pos] = cur_seg;
        else
        {
          //place segment in L
          L[new_L_size++] = cur_seg;
          //storing segment position in Q_tail
          *(--Q_tail) = _Q_pos;
          //storing segment number in R
          if (SegmentsColl::is_line_segments)*(R_pos++) = cur_seg;
        }
      }
      else//segment  intersects ladder stairs
      {
        n_int+= _Q_pos -step;
        //place segment in L
        L[new_L_size++] = cur_seg;
        //storing segment position in Q_tail, only for nonline segments
        if (!SegmentsColl::is_line_segments) *(--Q_tail) = _Q_pos;
      }
    }
    if(_Q_pos == 0)
    { 
      dont_split_stripe = false;
      return new_L_size;
    }
    Q_tail = Q + len_of_Q;
    if (!SegmentsColl::is_line_segments)
      R_pos = L + new_L_size;
    // important to start from stair above current segm, meanwhile _Q[loc] is stair below
    _Q++;// so we incremement _Q and _Q[loc] become stair above
    for (auto _R = SegmentsColl::is_line_segments ? R : L; _R < R_pos; ++_R)
    {
      segments->SetCurSegCutBE(*_R);
      auto loc = *(--Q_tail);          // getting position from tail of Q
      auto c=loc;
      for (auto cur_Q = _Q + loc; (loc < _Q_pos) && (segments->FindCurSegIntWith(*cur_Q)); ++cur_Q, ++loc);
      n_int+=loc-c;
    }
    dont_split_stripe = n_int > new_L_size;
    _step_index += _Q_pos;
    return new_L_size;
  };

 
 
  protected:
  int4 *L = nullptr, *R = nullptr;
  CTHIS* clone_of = nullptr;

  uint4 nTotSegm, len_of_Q;

  template<class SegmentsColl>
  int4 SplitSIS(SegmentsColl* segments, int4& step_index, int4 Size)//simplified version for SearchInStrip
  {
    int4 father_last_step = step_index, new_L_size = 0;
    auto Q_tail = Q + len_of_Q;
    for (int4 cur_L_pos = 0; cur_L_pos<Size; cur_L_pos++)
    {
      int4 cur_seg = L[cur_L_pos];
      int4 step = step_index;
      segments->SetCurSegCutBE(cur_seg);
      while ((father_last_step<step) && (segments->FindCurSegIntWith(Q[step])))
        step--;

      if (step_index==step)
          Q[++step_index] = cur_seg;
      else
      {
        L[new_L_size++] = cur_seg;
        *(--Q_tail) = step_index + 1;
      }
    }
    Q_tail = Q + len_of_Q;
    for (int4 i = 0; i < new_L_size; ++i)
    {
      auto c = *(--Q_tail);
      segments->SetCurSegCutBE(L[i]);
      while ((c <= step_index) && (segments->FindCurSegIntWith(Q[c++])));
    }

    return new_L_size;
  };



  void clone(CTHIS* master)
  {
      nTotSegm = master->nTotSegm;
      len_of_Q = nTotSegm;
      SegL = master->SegL;
      SegR = master->SegR;
      ENDS = master->ENDS;
      L = new int4[nTotSegm];
      Q = new int4[nTotSegm];
      R = new int4[nTotSegm];
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
  uint4 CalcLAt(SegmentsColl* segments, uint4 end_index)
  {
      uint4 i, Size = 0;
      for (i = 0; i < nTotSegm; i++)
          if ((SegL[i] <= end_index) && (SegR[i] > end_index))
              L[Size++] = i;
      segments->SortAt(ENDS[end_index], Size, L);
      return Size;
  };

  template<class SegmentsColl>
  void AllocMem(SegmentsColl* segments)
  {
      nTotSegm = len_of_Q = segments->GetSegmNumb();
      L = new int4[nTotSegm];
      R = new int4[nTotSegm];
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

