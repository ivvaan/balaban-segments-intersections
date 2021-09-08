#ifndef MEM_SAVE_FINDER_FOR_SEGMENT_INTERSECTION
#define MEM_SAVE_FINDER_FOR_SEGMENT_INTERSECTION
 
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



class CMemSaveIntFinder : public CommonImpl 
{
public:
  using CIMP = CommonImpl;
  using CIMP::dont_split_stripe;
  using CIMP::prepare_ends; using CIMP::FindInt; using CIMP::FindIntI;
  using CIMP::SegL; using CIMP::SegR; using CIMP::ENDS; using CIMP::Q;
  using CIMP::SearchInStripNonLineSeg; using CIMP::SearchInStripLineSeg;


  template<class SegmentsColl>
  void find_intersections(SegmentsColl* segments)
  {
    AllocMem(segments);
    ProgramStackRec stack_rec(-1, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl::get_segm(ENDS[0]);
    from_begin = true;
    FindR(this, segments, -1, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0, get_max_call(2 * nTotSegm));
    FreeMem();
  }

  template<class SegmentsColl>
  int4 SearchInStrip(SegmentsColl* segments, int4 qp, int4 Size)
  {
    auto _L = from_begin ? L : L + (LR_len - Size);
    if constexpr (SegmentsColl::is_line_segments)
    {
      //For line segments we can do more efficient insertion sorting using intersection check as comparison.
      //If s1<s2 at the left bound then s1 intersects s2 inside the stripe means s1>s2 at the right bound of the stripe.
      //All segments are sorted at the left bound by precondition, so intesection is the same as comparison
      // at the right bound. Byproduct of the sorting is intersections detection.
      SearchInStripLineSeg(segments, _L, Size);
    }
    else {

      auto _Q = Q + qp;  //prepere Q pointer for SearchInStripNonLineSeg
      //Make sequential simplified Split
      SearchInStripNonLineSeg(segments, _L, _Q++, Size); //_Q pointed to element before the first new (moved from L) 
      //segment in Q and we need to increment Q for correct copy (see below)

      //as a result of SearchInStripNonLineSeg we have all the intersections 
      //found and all the segments moved from L to Q
      //so we need to copy them back and
      //once segments in L again (but in icorrect order)
      //we sort it
      std::sort(_L, std::copy(_Q, _Q + Size, _L),
        [segments](int4 s1, int4 s2) {return segments->RBelow(s1, s2); });
    }

    return Size;

  };

  template<class SegmentsColl>
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
        i = LR_len - Size;
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
        for (i = LR_len - Size; (i < LR_len) && (segments->UnderCurPoint(L[i])); ++i)
          L[i - 1] = L[i];
        L[i - 1] = sn;
        ++Size;

        segments->SetCurSeg(sn);
        FindIntI(segments, SegR[sn], stack_pos); // get internal intersections
      }
    }
    return Size;
  };

  template<class SegmentsColl>
  int4 Merge(SegmentsColl* segments, uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
  {
    int4 cur_R_pos = 0, new_size = 0;
    int4 cur_seg;
    if (from_begin)
    {
      auto _R = L + Size - 1;
      auto _L = L + (LR_len - 1);
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
    auto _R = L + (LR_len - Size);
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

  template<class SegmentsColl>
  int4 Split(SegmentsColl* segments, uint4 RBoundIdx, int4& _step_index, int4 Size)
  {
    auto _Q = Q + _step_index;
    auto _L = from_begin ? L : L + (LR_len - Size);
    int4  new_L_size = 0,_Q_pos=0;
    auto Q_tail = Q + len_of_Q;
    long long n_int=0;
    for (int4 cur_L_pos = 0; cur_L_pos < Size; cur_L_pos++)
    {
      int4 cur_seg = _L[cur_L_pos];
      auto step = _Q_pos;
      segments->SetCurSegCutBE(cur_seg);
      for (auto cur_Q = _Q + step; (step!=0) && (segments->FindCurSegIntWith(*cur_Q)); --step, --cur_Q);

      if (_Q_pos == step)
      {
        if (SegR[cur_seg] >= RBoundIdx) //segment is covering current stripe and doesn't intersect ladder stairs
          _Q[++_Q_pos] = cur_seg;
        else {
          //place segment in L
          L[new_L_size++] = cur_seg;
          //storing segment position in Q_tail
          *(--Q_tail) = _Q_pos;
        }
      }
      else
      {
        n_int+= _Q_pos -step;
        //place segment in L
        L[new_L_size++] = cur_seg;
        //storing segment position in Q_tail
        *(--Q_tail) = SegmentsColl::is_line_segments ? INT_MAX : _Q_pos;
      }
    }
    from_begin = true;
    if (_Q_pos == 0)
    {
      dont_split_stripe = false;
      return new_L_size;
    }
    Q_tail = Q + len_of_Q;
    int4 loc;
    _L = L;
    // important to start from stair above current segm, meanwhile _Q[loc] is stair below
    _Q++;// so we incremement _Q and _Q[loc] become stair above
    auto last_Q = _Q + _Q_pos;
    for (int4 i = 0; i < new_L_size; ++i)
      if ((loc = *(--Q_tail)) != INT_MAX) {
        segments->SetCurSegCutBE(_L[i]);
        auto cur_Q = _Q + loc;
        auto ini_Q = cur_Q;
        while ((cur_Q < last_Q) && (segments->FindCurSegIntWith(*cur_Q)))++cur_Q;
        n_int+= cur_Q - ini_Q;
      }
//    dont_split_stripe = location != _step_index;
    dont_split_stripe = n_int>new_L_size;
    _step_index += _Q_pos;
    
    return new_L_size;
  };

  void FreeMem()
  {
    MY_FREE_ARR_MACRO(L);
    MY_FREE_ARR_MACRO(Q);
    //CIMP::FreeMem();
  };


private:
  bool from_begin = true;

  template<class SegmentsColl>
  void AllocMem(SegmentsColl* segments)
  {
    assert(nTotSegm  == segments->GetSegmNumb());
    len_of_Q = LR_len;
    L = new int4[LR_len];
    Q = new int4[len_of_Q];
  };

};

#endif