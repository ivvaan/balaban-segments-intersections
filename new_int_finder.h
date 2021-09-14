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
//#include <array>
//#include <memory>
//#include "test_coll.h"



class CTrivialIntFinder
{
public:
  template<class SegmentsColl>
  void find_intersections(SegmentsColl *segments) 
  {
    uint4 N = segments->GetSegmNumb();
    for (uint4 i = 0; i < N; ++i)
    {
      segments->SetCurSeg(i);
      for (uint4 j = i + 1; j < N; ++j)
        segments->TrivCurSegIntWith(j);

    }
  };

};



class CSimpleSweepIntFinder {
public:
  template <class SegmentsColl>
    void find_intersections(SegmentsColl* segments)
    {
        uint4 n = segments->GetSegmNumb();
        uint4* sgm = new uint4[n * 3];
        uint4 *Ends = sgm+n;
        segments->PrepareEndpointsSortedList(Ends);

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
                segments->SetCurSeg(s);
                for (auto j = sgm; j < pos; ++j)
                    segments->SSCurSegIntWith(*j);
                *pos = s;
                ++pos;
            }
        }
        delete[] sgm;
    };
};


  struct ProgramStackRec //structure to use program stack to store list of starcases above current call
  {
    ProgramStackRec *prev=nullptr; //link to program stack record storing parent staircase information 
    int4 Q_pos; //starting position and 
    uint4 right_bound;//right boind of current staircase
    inline ProgramStackRec(int4 qp) :Q_pos(qp) {};
    inline ProgramStackRec(int4 qp, uint4 rb) : Q_pos(qp), right_bound(rb) {};
    inline ProgramStackRec *Set(ProgramStackRec *p, uint4 rb) { prev = p, right_bound = rb; return this; };
  };


class CommonImpl 
{
public:
    bool dont_split_stripe;

    template <class SegmentsColl>
    void prepare_ends(SegmentsColl* segments)
    {
      auto N= segments->GetSegmNumb();
      AllocMem(N);
      segments->PrepareEndpointsSortedList(ENDS);
      
      uint4 max_segm_on_vline=0,nsegm_on_vline = 0;
      for (uint4 i = 0; i < 2 * N; i++)
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
      LR_len = max_segm_on_vline+1;
      nTotSegm = N;
          
    };

    ~CommonImpl() { FreeMem(); };
protected:

  uint4 LR_len = 0;
  uint4 nTotSegm = 0;
  uint4 len_of_Q = 0;

  uint4* SegL = nullptr, * SegR = nullptr, * ENDS = nullptr;
  int4* Q = nullptr;
  int4* L = nullptr;

  void AllocMem(uint4 N)
  {
    SegL = new uint4[N];
    SegR = new uint4[N];
    ENDS = new uint4[2 * N];
  };

  void FreeMem()
  {
    MY_FREE_ARR_MACRO(SegL);
    MY_FREE_ARR_MACRO(SegR);
    MY_FREE_ARR_MACRO(ENDS);
  };

  template <class IntersectionFinder, class SegmentsColl>
  static  int4 FindR(IntersectionFinder* i_f, SegmentsColl* segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos, uint4 Size, uint4 call_numb, uint4 max_call = 30)
  {
    segments->SetCurStripe(i_f->ENDS[interval_left_index], i_f->ENDS[interval_right_index]);
    if (interval_right_index - interval_left_index == 1)
      return Size > 1 ? i_f->SearchInStrip(segments, ladder_start_index, Size) : Size;

    ProgramStackRec stack_rec(ladder_start_index);
    if (Size > 0)
    {
      Size = i_f->Split(segments, interval_right_index, stack_rec.Q_pos, Size);
      if (ladder_start_index < stack_rec.Q_pos)
        stack_pos = stack_rec.Set(stack_pos, interval_right_index);
    };
    if (i_f->dont_split_stripe && (call_numb < max_call)) //if found a lot of intersections repeat FindR
      Size = FindR(i_f, segments, stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1, max_call);
    else //cut stripe 
    {
      uint4 m = (interval_left_index + interval_right_index) >> 1;
      if (call_numb > 1)
      { // if L contains a lot of segments then cut on two parts
        max_call -= 2;
        Size = FindR(i_f, segments, stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0, max_call);
        Size = i_f->InsDel(segments, m, stack_pos, Size);
        Size = FindR(i_f, segments, stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0, max_call);
      }
      else
      {// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
        max_call -= 4;
        uint4 q = (interval_left_index + m) >> 1;
        if (interval_left_index != q) {
          Size = FindR(i_f, segments, stack_rec.Q_pos, interval_left_index, q, stack_pos, Size, 0, max_call);
          Size = i_f->InsDel(segments, q, stack_pos, Size);
        }
        Size = FindR(i_f, segments, stack_rec.Q_pos, q, m, stack_pos, Size, 0, max_call);
        Size = i_f->InsDel(segments, m, stack_pos, Size);
        q = (interval_right_index + m) >> 1;
        if (q != m) {
          Size = FindR(i_f, segments, stack_rec.Q_pos, m, q, stack_pos, Size, 0, max_call);
          Size = i_f->InsDel(segments, q, stack_pos, Size);
        }
        Size = FindR(i_f, segments, stack_rec.Q_pos, q, interval_right_index, stack_pos, Size, 0, max_call);
      }
    }
    if (ladder_start_index >= stack_rec.Q_pos) return Size;
    return i_f->Merge(segments, interval_left_index, ladder_start_index, stack_rec.Q_pos, Size);
  };


  //functions for fast algorithm
  template <class SegmentsColl>
  void FindInt(SegmentsColl* segments, int4 qb, int4 qe, int4 l) const
  {
    int4 c = l;
    segments->SetSearchDirDown(true);
    while ((c > qb) && (segments->FindCurSegIntWith(Q[c]))) //first get intersections below
      --c;
    if (SegmentsColl::is_line_segments && (c != l))
      return; //if found and segment is line it can't be any more
    c = l + 1;
    segments->SetSearchDirDown(false);
    while ((c <= qe) && (segments->FindCurSegIntWith(Q[c]))) // get intersections above
      ++c;
  };

  template <class SegmentsColl>
  void FindIntI(SegmentsColl* segments, uint4 r_index, ProgramStackRec* stack_pos) const
  {
    while (stack_pos->right_bound <= r_index)
      stack_pos = stack_pos->prev; // go from bottom to top and find staircase to start
    int4 l, r, m, QB, QE = stack_pos->Q_pos;
    for (stack_pos = stack_pos->prev; stack_pos != nullptr; stack_pos = stack_pos->prev) {
      l = QB = stack_pos->Q_pos;
      r = QE + 1;
      while ((r - l) > 1) // binary search
      {
        m = (r + l) >> 1; //        m=(r+l)/2;
        if (segments->UnderCurPoint(Q[m]))
          l = m;
        else
          r = m;
      }
      FindInt(segments, QB, QE, l);
      QE = QB; // move staircase bound to parent staircase
    };
  };

  template<class SegmentsColl>
  void SearchInStripNonLineSeg(SegmentsColl* segments, int4* _L, int4* _Q, int4 size)//simplified version for SearchInStrip
  {
    auto Size = size;
    do
    {
      auto Q_tail = Q + len_of_Q;
      int4 new_L_size = 0;
      int4 step_index = 1;
      _Q[1] = _L[0];
      for (int4 cur_L_pos = 1; cur_L_pos < Size; cur_L_pos++)
      {
        auto step = step_index;
        auto cur_seg = _L[cur_L_pos];
        segments->SetCurSegCutBE(cur_seg);
        while ((step) && (segments->FindCurSegIntWith(_Q[step])))
          --step;

        if (step_index == step)
          _Q[++step_index] = cur_seg;
        else
        {
          _L[new_L_size++] = cur_seg;
          *(--Q_tail) = step_index;
        }
      }
      Q_tail = Q + len_of_Q;
      for (int4 i = 0; i < new_L_size; ++i)
      {
        auto c = *(--Q_tail);
        if (c >= step_index) break;
        segments->SetCurSegCutBE(_L[i]);
        while ((c < step_index) && (segments->FindCurSegIntWith(_Q[++c])));
      }
      _Q += step_index;
      Size = new_L_size;
    } while (Size);

  };

  template<class SegmentsColl>
  void SearchInStripLineSeg(SegmentsColl* segments, int4* L_, int4 Size) {
    auto _L = L_ - 1;
    for (uint4 i = 1; i < Size; i++)
    {
      auto sn = L_[i];
      segments->SetCurSegCutBE(sn);
      segments->SetSearchDirDown(true);
      if (segments->FindCurSegIntWith(_L[i])) {
        L_[i] = _L[i];
        uint4 j = i - 1;
        for (; (j) && (segments->FindCurSegIntWith(_L[j])); --j)
          L_[j] = _L[j];
        L_[j] = sn;
      }
    }
  }

};


uint4 get_max_call(uint4 N) {
    uint4 max_call = 8;
    for (; N; N >>= 1) max_call += 2;
    return max_call;
};

//#undef MY_FREE_ARR_MACRO
#endif