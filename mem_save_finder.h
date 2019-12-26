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


template <class SegmentsColl>
class CMemSaveIntFinder : private CommonImpl<SegmentsColl> {
public:
  using CIMP = CommonImpl<SegmentsColl>;
  using CIMP::max_call;using CIMP::dont_split_stripe;
  using CIMP::prepare_ends; using CIMP::FindInt; using CIMP::FindIntI;
  using CIMP::SegL; using CIMP::SegR; using CIMP::ENDS; using CIMP::Q;

  void find_intersections(SegmentsColl* segments)
  {
    AllocMem(segments);
    prepare_ends(segments);
    ProgramStackRec stack_rec(-1, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl::get_segm(ENDS[0]);
    from_begin = true;
    FindR(this,segments, -1, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0,get_max_call(2*nTotSegm));
    FreeMem();
  }
  int4 SearchInStrip(SegmentsColl* segments, int4 qp, int4 Size)
  {
    if (SegmentsColl::is_line_segments)
    {
      auto L_ = from_begin ? L : L + (nTotSegm - Size);
      auto _L = L_ - 1;
      for (uint4 i = 1; i < Size; i++)
      {
        auto sn = L_[i];
        segments->SetCurSegCutBE(sn);
        uint4 j = i;
        for (; (j) && (segments->FindCurSegIntWith(_L[j])); --j)
          L_[j] = _L[j];
        L_[j] = sn;
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


  int4 InsDel(SegmentsColl* segments, uint4 end_index, ProgramStackRec* stack_pos, int4 Size)
  {
    int4 i;
    auto sn = SegmentsColl::get_segm(ENDS[end_index]);
    if (from_begin)
    {
      if (SegmentsColl::is_last(ENDS[end_index])) // if endpoint - delete
      {
        i = --Size;
        auto cur = L[i];
        while (cur != sn) {
          --i;
          auto buf = cur;
          cur = L[i];
          L[i] = buf;
        }
      }
      else // if startpoint - insert
      {
        segments->SetCurPointAtBeg(sn);
        for (i = Size - 1; (i > -1) && (!segments->UnderCurPoint(L[i])); --i)
          L[i + 1] = L[i];
        L[i + 1] = sn;
        Size++;
        segments->SetCurSeg(sn);
        FindIntI(segments, SegR[sn], stack_pos); // get internal intersections
      }

    }
    else
    {
      if (SegmentsColl::is_last(ENDS[end_index])) // if endpoint - delete
      {
        i = nTotSegm - Size;
        --Size;
        auto cur = L[i];
        while (cur != sn) {
          ++i;
          auto buf = cur;
          cur = L[i];
          L[i] = buf;
        }
      }
      else // if startpoint - insert
      {
        segments->SetCurPointAtBeg(sn);
        for (i = nTotSegm - Size; (i < nTotSegm) && (segments->UnderCurPoint(L[i])); ++i)
          L[i - 1] = L[i];
        L[i - 1] = sn;
        ++Size;

        segments->SetCurSeg(sn);
        FindIntI(segments, SegR[sn], stack_pos); // get internal intersections
      }
    }
    return Size;
  };

  int4 Merge(SegmentsColl* segments, uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
  {
    int4 cur_R_pos = 0, new_size = 0;
    int4 cur_seg;
    if (from_begin)
    {
      auto _R = L + Size - 1;
      auto _L = L + (nTotSegm - 1);
      auto cur_stair = QE, _Size = -Size;
      cur_seg = _R[cur_R_pos];
      while ((cur_stair > QB) && (cur_R_pos > _Size)) {
        if (segments->RBelow(Q[cur_stair], cur_seg)) {
          if (SegL[cur_seg] > LBoundIdx) {
            segments->SetCurSegCutEnd(cur_seg);
            FindInt(segments, QB, QE, cur_stair);
          }
          _L[new_size--] = cur_seg;
          cur_seg = _R[--cur_R_pos];
        }
        else
          _L[new_size--] = Q[cur_stair--];
      }
      while (cur_R_pos > _Size) {
        if (SegL[cur_seg] > LBoundIdx) {
          segments->SetCurSegCutEnd(cur_seg);
          FindInt(segments, QB, QE, QB);
        }
        _L[new_size--] = cur_seg;
        cur_seg = _R[--cur_R_pos];
      }
      while (cur_stair > QB)
        _L[new_size--] = Q[cur_stair--];
      from_begin = false;
      return -new_size;
    }
    auto _R = L + (nTotSegm - Size);
    int4 cur_stair = QB;
    cur_seg = _R[cur_R_pos];
    while ((cur_stair < QE) && (cur_R_pos < Size)) {
      if (segments->RBelow(cur_seg, Q[cur_stair + 1])) {
        if (SegL[cur_seg] > LBoundIdx) {
          segments->SetCurSegCutEnd(cur_seg);
          FindInt(segments, QB, QE, cur_stair);
        }
        L[new_size++] = cur_seg;
        cur_seg = _R[++cur_R_pos];
      }
      else
        L[new_size++] = Q[++cur_stair];
    }
    while (cur_R_pos < Size) {
      if (SegL[cur_seg] > LBoundIdx) {
        segments->SetCurSegCutEnd(cur_seg);
        FindInt(segments, QB, QE, QE);
      }
      L[new_size++] = cur_seg;
      cur_seg = _R[++cur_R_pos];
    }
    while (cur_stair < QE)
      L[new_size++] = Q[++cur_stair];
    from_begin = true;
    return new_size;
  };

  int4 Split(SegmentsColl* segments, uint4 RBoundIdx, int4& _step_index, int4 Size)
  {
    auto _L = from_begin ? L : L + (nTotSegm - Size);
    int4 step_index = _step_index;
    int4 father_last_step = step_index, new_L_size = 0;
    auto Q_tail = Q + len_of_Q;
    long long n_int=0;
    for (int4 cur_L_pos = 0; cur_L_pos < Size; cur_L_pos++)
    {
      int4 cur_seg = _L[cur_L_pos];
      int4 step = step_index;
      segments->SetCurSegCutBE(cur_seg);
      while ((father_last_step < step) && (segments->FindCurSegIntWith(Q[step])))
        step--;

      if (step_index == step)
      {
        if (SegR[cur_seg] >= RBoundIdx) //segment is covering current stripe and doesn't intersect ladder stairs
          Q[++step_index] = cur_seg;
        else {
          //place segment in L
          L[new_L_size++] = cur_seg;
          //storing segment position in Q_tail
          *(--Q_tail) = step_index + 1;
        }
      }
      else
      {
        n_int+=step_index-step;
        //place segment in L
        L[new_L_size++] = cur_seg;
        //storing segment position in Q_tail
        *(--Q_tail) = SegmentsColl::is_line_segments ? INT_MAX : step_index + 1;
      }
    }
    Q_tail = Q + len_of_Q;
    int4 c;
    for (int4 i = 0; i < new_L_size; ++i)
      if ((c = *(--Q_tail)) != INT_MAX) {
        segments->SetCurSegCutBE(L[i]);
        auto loc=c;
        while ((c <= step_index) && (segments->FindCurSegIntWith(Q[c++])));
        n_int+=c-loc;
      }
//    dont_split_stripe = step_index != _step_index;
    dont_split_stripe = n_int>new_L_size;
    _step_index = step_index;
    from_begin = true;
    return new_L_size;
  };


private:
  int4 *L = nullptr;
  uint4 nTotSegm, len_of_Q;
  bool from_begin = true;

  void AllocMem(SegmentsColl* segments)
  {
    nTotSegm = len_of_Q = segments->GetSegmNumb();
    L = new int4[nTotSegm];
    Q = new int4[len_of_Q];
  };

  void FreeMem()
  {
    MY_FREE_ARR_MACRO(L);
    MY_FREE_ARR_MACRO(Q);
    CIMP::FreeMem();
  };

  //functions for fast algorithm
  int4 SplitSIS(SegmentsColl* segments, int4& step_index, int4 Size)//simplified version for SearchInStrip
  {
    auto _L = from_begin ? L : L + (nTotSegm - Size);
    int4 father_last_step = step_index, new_L_size = 0;
    auto Q_tail = Q + len_of_Q;
    for (int4 cur_L_pos = 0; cur_L_pos < Size; cur_L_pos++)
    {
      int4 cur_seg = _L[cur_L_pos];
      int4 step = step_index;
      segments->SetCurSegCutBE(cur_seg);
      while ((father_last_step < step) && (segments->FindCurSegIntWith(Q[step])))
        step--;

      if (step_index == step)
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
        auto c =  *(--Q_tail);
        segments->SetCurSegCutBE(L[i]);
        while ((c <= step_index) && (segments->FindCurSegIntWith(Q[c++])));
      }

    from_begin = true;
    return new_L_size;
  };


};
