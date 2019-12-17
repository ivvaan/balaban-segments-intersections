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
#include "new_int_finder.h"
#include <algorithm>
#include <vector>
#include <thread>

template<class SegmentsColl>
class CFastIntFinder : protected CommonImpl<SegmentsColl>
{
public:
  using CTHIS = CFastIntFinder<SegmentsColl>;
  using  CIMP = CommonImpl<SegmentsColl>;
  using CIMP::SegL; using CIMP::SegR; using CIMP::ENDS; using CIMP::Q;
  using CIMP::prepare_ends; using CIMP::FindInt; using CIMP::FindIntI;

  void find_intersections(SegmentsColl *segments)
  {
    AllocMem(segments);
    prepare_ends(segments);
    ProgramStackRec stack_rec(-1, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl::get_segm(ENDS[0]);
    FindR(segments, -1, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0);
    FreeMem();
  }
  template<class CIntRegistrator>
  void find_intersections(SegmentsColl* segments, uint4 n_threads, CIntRegistrator* regs[])
  {
    AllocMem(segments);
    prepare_ends(segments);
    using namespace std;
    vector<thread> wrk_threads;
    auto thread_func = [](CTHIS* master, SegmentsColl* segments, uint4 from, uint4 to, CIntRegistrator* add_reg) {
      CTHIS i_f;
      SegmentsColl coll;
      i_f.clone(master);
      coll.clone(segments, add_reg);
      i_f.FindR(&coll, -1, from, to, &ProgramStackRec(-1, 2 * i_f.nTotSegm), i_f.CalcLAt(&coll, from), 0);
      coll.unclone();
      i_f.unclone();
    };
    auto n = segments->GetSegmNumb();
    double part = 2 * n / (double)n_threads;
    int4 i = 1;
    uint4 start_from = part;
    uint4 from, to = start_from;
    while (i < n_threads) {
      from = to;
      i++;
      to = (i == n_threads) ? 2 * n - 1 : part * i;
      wrk_threads.emplace_back(thread_func, this, segments, from, to, regs[i - 2]); // starts intersection finding in a stripe <from,to>
    }

    ProgramStackRec stack_rec(-1, 2 * nTotSegm); //need to be initialized this way
    L[0] = SegmentsColl::get_segm(ENDS[0]);
    FindR(segments, -1, 0, start_from, &stack_rec, 1, 0);
    //      FindR(segments, -1, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0);
    for (auto cur_thread = wrk_threads.begin(); cur_thread != wrk_threads.end(); cur_thread++)
      cur_thread->join(); //waiting for calculation of all threads are finished
    FreeMem();
  }


protected:

  int4 *L = nullptr, *R = nullptr;
  CTHIS* clone_of = nullptr;

  uint4 nTotSegm, len_of_Q, int_numb;

  uint4 CalcLAt(SegmentsColl* segments, uint4 end_index)
  {
    uint4 i, Size = 0;
    for (i = 0; i < nTotSegm; i++)
      if ((SegL[i] <= end_index) && (SegR[i] > end_index))
        L[Size++] = i;
    segments->SortAt(ENDS[end_index], Size, L);
    return Size;
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
    if (clone_of != nullptr)
    {
      SegL = nullptr;
      SegR = nullptr;
      ENDS = nullptr;
      FreeMem();
      clone_of = nullptr;
    }
  };

  int4 SearchInStrip(SegmentsColl *segments, int4 qp, int4 Size)
  {
    if (SegmentsColl::is_line_segments)
    {
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
      while (size = Split(segments, 0, qp, size))++n_split; //it change QP must be after q = Q + QP + 1 
                                                            // at this point we can just place Q starting from QP+1 to L and sort it
      for (; _L < last; ++_L, ++q)  *_L = *q;
      if (n_split)std::sort(L, L + Size, [segments](int4 s1, int4 s2) {return segments->Below(s1, s2); });
    }
    return Size;
  };

  void AllocMem(SegmentsColl *segments)
  {
    nTotSegm = len_of_Q = segments->GetSegmNumb();
    L = new int4[nTotSegm];
    R = new int4[nTotSegm];
    Q = new int4[len_of_Q];
  };

  void FreeMem()
  {
    MY_FREE_ARR_MACRO(L);
    MY_FREE_ARR_MACRO(R);
    MY_FREE_ARR_MACRO(Q);
    CIMP::FreeMem();
  };
  //void FindIntL(int4 qe,int4 segm_numb);

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


  int4 Merge(SegmentsColl *segments, uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
  {
    int4 cur_R_pos = 0, new_size = 0, cur_stair = QB;
    auto _R = L, _L = R;
    int4 cur_seg = _R[cur_R_pos];
    //segments->SetCurVLine(stripe_right);
    while ((cur_stair<QE) && (cur_R_pos<Size))
    {

      if (segments->Below(cur_seg, Q[cur_stair + 1]))
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

  int4 Split(SegmentsColl* segments, uint4 RBoundIdx, int4& step_index, int4 Size)
  {
    int4 father_last_step = step_index, new_L_size = 0;
    auto Q_tail = Q + len_of_Q;
    int4 *R_pos;
    if (SegmentsColl::is_line_segments)R_pos = R;
    for (int4 cur_L_pos = 0; cur_L_pos<Size; cur_L_pos++)
    {
      int4 cur_seg = L[cur_L_pos];
      int4 step = step_index;
      segments->SetCurSegCutBE(cur_seg);
      while ((father_last_step<step) && (segments->FindCurSegIntWith(Q[step])))
        step--;
      int4 n_int = step_index - step;

      if (n_int == 0)
      {
        if (SegR[cur_seg] >= RBoundIdx)//segment is covering current stripe and doesn't intersect ladder stairs
          Q[++step_index] = cur_seg;
        else
        {
          //place segment in L
          L[new_L_size++] = cur_seg;
          //storing segment position in Q_tail
          *(--Q_tail) = step_index + 1;
          //storing segment number in R
          if (SegmentsColl::is_line_segments)*(R_pos++) = cur_seg;
        }
      }
      else
      {
        int_numb += n_int;
        //place segment in L
        L[new_L_size++] = cur_seg;
        //storing segment position in Q_tail, only for nonline segments
        //*(--Q_tail) = is_line_seg? INT_MAX :step_index+1;
        if (!SegmentsColl::is_line_segments) *(--Q_tail) = step_index + 1;
      }
    }
    Q_tail = Q + len_of_Q;
    if (!SegmentsColl::is_line_segments)
      R_pos = L + new_L_size;
    for (auto _R = SegmentsColl::is_line_segments ? R : L; _R < R_pos; ++_R)
    {
      segments->SetCurSegCutBE(*_R);
      int4 loc = *(--Q_tail);          // getting ts position from tail of Q
      int4 c = loc;
      while ((c <= step_index) && (segments->FindCurSegIntWith(Q[c++])));
      int_numb += c - loc;
    }

    return new_L_size;
  };

  int4 FindR(SegmentsColl *segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec *stack_pos, int4 Size, int4 call_numb)
  {
    segments->SetCurStripe(ENDS[interval_left_index], ENDS[interval_right_index]);
    if (interval_right_index - interval_left_index == 1)
      return Size>1 ? SearchInStrip(segments, ladder_start_index, Size) : Size;
    ProgramStackRec stack_rec(ladder_start_index);
    int_numb = 0;// variable to count intersections on Split stage
    if (Size>0)
    {
      Size = Split(segments, interval_right_index, stack_rec.Q_pos, Size);
      if ((ladder_start_index<stack_rec.Q_pos))
        stack_pos = stack_rec.Set(stack_pos, interval_right_index);
    };
    if ((int_numb>Size) && (call_numb<max_call)) //if found a lot of intersections repeat FindR
      Size = FindR(segments, stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
    else //cut stripe 
    {

      uint4 m = (interval_left_index + interval_right_index) / 2;
      if (call_numb>1)
      { // if L contains a lot of segments then cut on two parts
        Size = FindR(segments, stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
        Size = InsDel(segments, m, stack_pos, Size);
        Size = FindR(segments, stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
      }
      else
      {// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
        uint4 q = (interval_left_index + m) / 2;
        if (interval_left_index != q) {
          Size = FindR(segments, stack_rec.Q_pos, interval_left_index, q, stack_pos, Size, 0);
          Size = InsDel(segments, q, stack_pos, Size);
        }
        Size = FindR(segments, stack_rec.Q_pos, q, m, stack_pos, Size, 0);
        Size = InsDel(segments, m, stack_pos, Size);
        q = (interval_right_index + m) / 2;
        if (q != m) {
          Size = FindR(segments, stack_rec.Q_pos, m, q, stack_pos, Size, 0);
          Size = InsDel(segments, q, stack_pos, Size);
        }
        Size = FindR(segments, stack_rec.Q_pos, q, interval_right_index, stack_pos, Size, 0);
      }
    }
    if (ladder_start_index >= stack_rec.Q_pos) return Size;
    // segments->SetCurStripe(ENDS[interval_left_index], ENDS[interval_right_index]);
    return Merge(segments, interval_left_index, ladder_start_index, stack_rec.Q_pos, Size);
  };

};

