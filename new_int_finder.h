#pragma once
/*
*
*      Copyright (c)  2011-2017  Ivan Balaban
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

#include "utils.h"
#include <math.h>
#include <memory.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include <thread>
//#include <array>
//#include <memory>
//#include "test_coll.h"

#define MY_FREE_ARR_MACRO(a) \
    if (a!=nullptr) {                 \
        delete[] a;          \
        a = nullptr;            \
    }

template<class SegmentsColl>
class CTrivialIntFinder
{
public:

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


template <class SegmentsColl>
class CSimpleSweepIntFinder {
public:
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
    ProgramStackRec *prev; //link to program stack record storing parent staircase information 
    int4 Q_pos; //starting position and 
    uint4 right_bound;//right boind of current staircase
    inline ProgramStackRec(int4 qp) :Q_pos(qp) {};
    inline ProgramStackRec(int4 qp, uint4 rb) :prev(NULL), Q_pos(qp), right_bound(rb) {};
    inline ProgramStackRec *Set(ProgramStackRec *p, uint4 rb) { prev = p, right_bound = rb; return this; };
  };

template <class SegmentsColl>
class CommonImpl 
{
protected:
    uint4 *SegL = nullptr, *SegR = nullptr, *ENDS = nullptr;
    int4 *Q = nullptr;
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

    void prepare_ends(SegmentsColl* segments)
    {
        auto N = segments->GetSegmNumb();
        AllocMem(N);
        segments->PrepareEndpointsSortedList(ENDS);
        for (uint4 i = 0; i < 2 * N; i++)
            if (SegmentsColl::is_last(ENDS[i]))
                SegR[SegmentsColl::get_segm(ENDS[i])] = i;
            else
                SegL[SegmentsColl::get_segm(ENDS[i])] = i;
    };

    //functions for fast algorithm
    void FindInt(SegmentsColl* segments, int4 qb, int4 qe, int4 l)
    {
        int4 c = l;
        while ((c > qb) && (segments->FindCurSegIntWith(Q[c]))) //first get intersections below
            --c;
        if (SegmentsColl::is_line_segments && (c != l))
            return; //if found and segment is line it can't be any more
        c = l + 1;
        while ((c <= qe) && (segments->FindCurSegIntWith(Q[c]))) // get intersections above
            ++c;
    };
    ;
    void FindIntI(SegmentsColl* segments, uint4 r_index, ProgramStackRec* stack_pos /*, uint4 seg*/)
    {
        while (stack_pos->right_bound <= r_index)
            stack_pos = stack_pos->prev; // go from bottom to top and find staircase to start
        //    segments->SetCurPointAtBeg(seg);
        //    segments->SetCurSeg(seg);
        int4 l, r, m, QB, QE = stack_pos->Q_pos;
        for (stack_pos = stack_pos->prev; stack_pos; stack_pos = stack_pos->prev) {
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
    } //*/
};


template<class SegmentsColl>
class CFastIntFinder : protected CommonImpl<SegmentsColl> 
{
public:
using CTHIS=CFastIntFinder<SegmentsColl>;
using  CIMP=CommonImpl<SegmentsColl>;
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
          coll.clone(segments,add_reg);
          i_f.FindR(&coll,-1, from, to, &ProgramStackRec(-1, 2 * i_f.nTotSegm), i_f.CalcLAt(&coll,from), 0);
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
  CTHIS* clone_of=nullptr;

  uint4 nTotSegm,len_of_Q,int_numb;
  
  uint4 CalcLAt(SegmentsColl* segments,uint4 end_index)
  {
      uint4 i, Size = 0;
      for (i = 0; i < nTotSegm; i++) 
          if ((SegL[i] <= end_index) && (SegR[i] > end_index))
              L[Size++] = i;
      segments->SortAt(ENDS[end_index],Size,L);
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
    if (clone_of!=nullptr) 
      {
          SegL = nullptr;
          SegR = nullptr;
          ENDS = nullptr;
          FreeMem();
          clone_of = nullptr;
      }
  };

  int4 SearchInStrip(SegmentsColl *segments,int4 qp, int4 Size)
  {
    if (SegmentsColl::is_line_segments)
      {
        auto _L = L-1;
        for (uint4 i = Size; i > 1; )
        {
          auto sn = _L[i];
          segments->SetCurSegCutBE(sn);
          uint4 j=--i;
          while ((j)&&(segments->FindCurSegIntWith(_L[j])))
          {
            _L[j + 1] = _L[j];
            --j;
          }
          _L[j + 1] = sn;
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
      FindIntI(segments,SegR[sn], stack_pos);// get internal intersections
    }
    return Size;
  }


  int4 Merge(SegmentsColl *segments, uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
  {
    int4 cur_R_pos = 0, new_size = 0, cur_stair = QB;
    auto _R = L,_L=R;
   int4 cur_seg = _R[cur_R_pos];
    //segments->SetCurVLine(stripe_right);
    while ((cur_stair<QE) && (cur_R_pos<Size))
    {

      if (segments->Below(cur_seg, Q[cur_stair + 1]))
      {
        if (SegL[cur_seg]>LBoundIdx)
        { 
          segments->SetCurSegCutEnd(cur_seg);
          FindInt(segments,QB, QE, cur_stair);
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
    L=_L;R=_R;
    return new_size;
  };

  int4 Split(SegmentsColl* segments,uint4 RBoundIdx,int4& step_index, int4 Size)
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
      return Size>1 ? SearchInStrip(segments,ladder_start_index, Size) : Size;
    ProgramStackRec stack_rec(ladder_start_index);
    int_numb = 0;// variable to count intersections on Split stage
    if (Size>0)
    {
        Size = Split(segments, interval_right_index,stack_rec.Q_pos, Size);
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

template <class SegmentsColl>
class CMemSaveIntFinder : private CommonImpl<SegmentsColl> {
public:
using CIMP = CommonImpl<SegmentsColl>;
using CIMP::prepare_ends; using CIMP::FindInt; using CIMP::FindIntI;
using CIMP::SegL; using CIMP::SegR; using CIMP::ENDS; using CIMP::Q;

    void find_intersections(SegmentsColl* segments)
    {
        AllocMem(segments);
        prepare_ends(segments);
        ProgramStackRec stack_rec(-1, 2 * nTotSegm); //need to be initialized this way
        L[0] = SegmentsColl::get_segm(ENDS[0]);
        from_begin=true;
        FindR(segments, -1, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0);
        FreeMem();
    }

private:
    int4 *L = nullptr;
    uint4 nTotSegm, len_of_Q, int_numb;
    bool from_begin = true;

    int4 SearchInStrip(SegmentsColl* segments, int4 qp, int4 Size)
    {
      if (SegmentsColl::is_line_segments)
      {
        auto _L = L - 1;
        for (uint4 i = Size; i > 1; )
        {
          auto sn = _L[i];
          segments->SetCurSegCutBE(sn);
          uint4 j = --i;
          while ((j) && (segments->FindCurSegIntWith(_L[j])))
          {
            _L[j + 1] = _L[j];
            --j;
          }
          _L[j + 1] = sn;
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
                if (segments->Below(Q[cur_stair], cur_seg)) {
                    if (SegL[cur_seg] > LBoundIdx) {
                        segments->SetCurSegCutEnd(cur_seg);
                        FindInt(segments, QB, QE, cur_stair);
                    }
                    _L[new_size--] = cur_seg;
                    cur_seg = _R[--cur_R_pos];
                } else
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
            if (segments->Below(cur_seg, Q[cur_stair + 1])) {
                if (SegL[cur_seg] > LBoundIdx) {
                    segments->SetCurSegCutEnd(cur_seg);
                    FindInt(segments, QB, QE, cur_stair);
                }
                L[new_size++] = cur_seg;
                cur_seg = _R[++cur_R_pos];
            } else
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

    int4 Split(SegmentsColl* segments, uint4 RBoundIdx, int4& step_index, int4 Size)
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
            int4 n_int = step_index - step;

            if (n_int == 0) 
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
                int_numb += n_int;
                //place segment in L
                L[new_L_size++] = cur_seg;
                //storing segment position in Q_tail
                *(--Q_tail) = SegmentsColl::is_line_segments ? INT_MAX : step_index + 1;
            }
        }
        Q_tail = Q + len_of_Q;
        int4 loc;
        for (int4 i = 0; i < new_L_size; ++i)
            if ((loc = *(--Q_tail)) != INT_MAX) {
                auto c = loc;
                segments->SetCurSegCutBE(L[i]);
                while ((c <= step_index) && (segments->FindCurSegIntWith(Q[c++])));
                int_numb += c - loc;
            }

        from_begin = true;
        return new_L_size;
    };

    int4 FindR(SegmentsColl* segments, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec* stack_pos, int4 Size, int4 call_numb)
    {
        segments->SetCurStripe(ENDS[interval_left_index], ENDS[interval_right_index]);
        if (interval_right_index - interval_left_index == 1)
            return Size > 1 ? SearchInStrip(segments, ladder_start_index, Size) : Size;
        ProgramStackRec stack_rec(ladder_start_index);
        int_numb = 0; // variable to count intersections on Split stage
        if (Size > 0) {
            Size = Split(segments, interval_right_index, stack_rec.Q_pos, Size);
            if ((ladder_start_index < stack_rec.Q_pos))
                stack_pos = stack_rec.Set(stack_pos, interval_right_index);
        };
        if ((int_numb > Size) && (call_numb < max_call)) //if found a lot of intersections repeat FindR
            Size = FindR(segments, stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
        else //cut stripe
        {

            uint4 m = (interval_left_index + interval_right_index) / 2;
            if (call_numb > 1) 
            { // if L contains a lot of segments then cut on two parts
                Size = FindR(segments, stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
                Size = InsDel(segments, m, stack_pos, Size);
                Size = FindR(segments, stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
            } 
            else 
            { // if L contains not so many segments than cut on four parts (works faster for some segment distributions)
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



#undef MY_FREE_ARR_MACRO