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

#include "fast_finder.h"
#include <algorithm>
#include <vector>
#include <thread>

class COptimalIntFinder : public CFastIntFinder
{
public:
  using  CIMP = CommonImpl;
  using CFAST = CFastIntFinder;
  using CIMP::SegL; using CIMP::SegR; using CIMP::ENDS; using CIMP::Q;
  using CIMP::prepare_ends; CIMP::max_call;
  using CFAST::L; using CFAST::R; CFAST::nTotSegm; using CFAST::len_of_Q; 

  ~COptimalIntFinder() { FreeMem(); };
  template<class SegmentsColl>
  void find_intersections(SegmentsColl *segments)
  {
    AllocMem(segments);
    ProgramStackRec stack_rec(inherit_each, 2 * nTotSegm);  //need to be initialized this way 
    L[0] = SegmentsColl::get_segm(ENDS[0]);
    FindR(segments, inherit_each + 1, inherit_each, 0, 2 * nTotSegm - 1, &stack_rec, 1, 0);
    //FreeMem();
  }

  void FreeMem()
  {
    MY_FREE_ARR_MACRO(father_loc);
  };


private:
  static const int4 undef_loc = 0;

  static const int4 inherit_each = 16; // defines how often in optimal algorithm stairs are inherited,  one inherited per inherit_each
  static const int4 inherit_offset = inherit_each / 2; // first stair to inherit; must be in [1..inherit_each-1]
                                                       // so inherited stair positions are: inherit_offset,inherit_offset+inherit_each,inherit_offset+2*inherit_each,...
  static const int4 big_staircase_threshold = inherit_each*8;// 1024;//used in optFindR

  int4 *father_loc = nullptr;

  template<class SegmentsColl>
  void AllocMem(SegmentsColl *segments)
  {
    nTotSegm =  segments->GetSegmNumb();
    L = new int4[nTotSegm];
    R = new int4[nTotSegm];
    len_of_Q = nTotSegm + nTotSegm / (inherit_each - 1) + inherit_each + 1;
    Q = new int4[len_of_Q];
    father_loc = new int4[len_of_Q];
    for (int4 i = 0; i < len_of_Q; i++)father_loc[i] = undef_loc;
    //father_loc[undef_loc] = undef_loc;//nesessary precondition for optFindIntI
  };


  static inline bool is_original(int4 v) { return v < 1; };

  template<class SegmentsColl>
  void FindInt(SegmentsColl* segments, int4 qb, int4 qe, int4 l)
  {
    
    int4 c = l;
    while (c>qb)
    {
      if (is_original(father_loc[c]))
      { //if original stair try to find and register intersection
        if (!segments->FindCurSegIntWith(Q[c])) break;
      }
      else //if inherited stair just check if intersects
        if (!segments->IsIntersectsCurSeg(Q[c])) break;
      c--;
    }
    if (SegmentsColl::is_line_segments && (c != l))return;
    c = ++l;
    while (c <= qe)
    {
      if (is_original(father_loc[c]))
      {
        if (!segments->FindCurSegIntWith(Q[c])) break;
      }
      else
        if (!segments->IsIntersectsCurSeg(Q[c])) break;
      c++;
    }
  };

  template<class SegmentsColl>
  void FindIntI(SegmentsColl* segments, uint4 r_index, ProgramStackRec* stack_pos)
  {
    while (stack_pos->right_bound <= r_index)stack_pos = stack_pos->prev;// go from bottom to top and find staircase to start
    if (stack_pos->prev==nullptr)return;
    int4  l = undef_loc, r, m,  QE = stack_pos->Q_pos + 1;
    int4 QB = undef_loc+1;//just to differ from l at first loop
    for (stack_pos = stack_pos->prev; stack_pos!=nullptr; stack_pos = stack_pos->prev)// for all staircases above
    {
        if ((big_staircase_threshold>inherit_each) && (l == QB)) //if location is below first stair we don't have data in father_loc for that location
            // because the number of locations is greater (by one) than number of stairs 
            // (for example one stair has two locations (below and above), but we keep in father_loc data for only one (above))
            // so we use data from next location QE==l+1
        {
            l = QB = stack_pos->Q_pos;// set lower location bound to the location below first step of current staircase
            r = father_loc[QE] ? QB + inherit_each : QE ;// father_loc[QE] !=undef_loc (i.e. 0) means staircase was created by optSplit 
            //and we use range [QB;QB+inherit_each] to search, otherwise it was created by Split and full range search [QB;QE] is used
        }
        else
        {
            QB = stack_pos->Q_pos;
            m=abs(father_loc[l]);//using father_loc get approximate location in parent staircase; 
            // line above sometimes executed when l == undef_loc, so we must initialize father_loc[undef_loc]=undef_loc in AllocMem()
            //to keep l unchanged in this case; undef_loc was chosen to be zero so abs(undef_loc)==undef_loc
            //otherwise we need additional checks
            l=max(QB,m);//here we use a fact that always undef_loc < QB and if m==undef_loc or m<QB l should be QB
            if (m) { r = l + inherit_each; if (r > QE)r = QE; }
            else r = QE;//if m==0 i.e. undef_loc we use [QB;QE] as a range
        }
      while ((r - l) > 1) // binary search
      {
        m = (r + l) >> 1; //        m=(r+l)/2;
        if (segments->UnderCurPoint(Q[m]))
          l = m;
        else
          r = m;
      }
      FindInt(segments, QB, QE - 1, l);
      QE = QB + 1; //set upper location bound to the location above last step of prev (which will be current) staircase
    };
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
        if (is_original(father_loc[++cur_stair]))
          _L[new_size++] = Q[cur_stair];
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
      if (is_original(father_loc[++cur_stair]))
        _L[new_size++] = Q[cur_stair];
    L = _L; R = _R;
    return new_size;
  };

  template<class SegmentsColl>
  int4 Split(SegmentsColl* segments, uint4 RBoundIdx, int4 cur_father_pos, int4& _step_index, int4 Size)
  {
    auto step_index = _step_index;
    int4 father_last_step = step_index, new_L_size = 0, cur_L_pos = 0;
    int4 cur_seg;
    auto location = R;
    while ((cur_L_pos<Size) && (cur_father_pos <= father_last_step))
    {
      cur_seg = L[cur_L_pos];
      segments->SetCurSegCutBE(cur_seg);
      if (segments->LBelow( cur_seg, Q[cur_father_pos])) // current segment below candidate to inherit
      {
        if ((SegR[cur_seg] < RBoundIdx)// current segment not covering strip
          || (step_index>father_last_step) && segments->IsIntersectsCurSeg(Q[step_index])// or intersects last stair
          || segments->IsIntersectsCurSeg(Q[cur_father_pos])) //or intersects candidate to inherit 
        {
          //place segment in L
          L[new_L_size] = cur_seg;
          //storing segment position in location
          location[new_L_size++] = step_index;
        }
        else
        { //add it to Q as original
          step_index++;
          // if segment is over this stair it is over previous inherited segment in father staircase
          // negative value of  father_loc[step_index] means it is original segment ( always cur_father_pos  >  inherit_each )              
          father_loc[step_index] = inherit_each - cur_father_pos;
          Q[step_index] = cur_seg;
        }
        cur_L_pos++; // move to next in L
      }
      else
      {  // add candidate to inherit to Q
        step_index++;
        // if segment is over this stair it is over the same inherited segment in father staircase  
        // positive value of  father_loc[step_index] means it is inherited segment                  
        father_loc[step_index] = cur_father_pos;
        Q[step_index] = Q[cur_father_pos];
        cur_father_pos += inherit_each;
      };
    }
    for (; cur_L_pos<Size; cur_L_pos++)  // father staircase is over. we are adding rest segments the same manner as in suboptimal case
    {
      cur_seg = L[cur_L_pos];
      segments->SetCurSegCutBE(cur_seg);
      if ((SegR[cur_seg] < RBoundIdx) ||//segment is not covering current stripe
        (step_index>father_last_step) && segments->IsIntersectsCurSeg(Q[step_index]))// if it intesects last stair
      {
        //place segment in L
        L[new_L_size] = cur_seg;
        //storing segment position in location
        location[new_L_size++] = step_index;
      }
      else
      {  // add it to Q
        step_index++;
        // negative value of father_loc[step_index] means it is original segment ( always cur_father_pos  >  inherit_each )              
        father_loc[step_index] = inherit_each - cur_father_pos;
        Q[step_index] = cur_seg;
      }
    };
    for (; cur_father_pos <= father_last_step; cur_father_pos += inherit_each) // L segments is over. We are adding rest inherited father segments to Q
    {
      step_index++;
      // positive value of  father_loc[step_index] means it is inherited segment                  
      father_loc[step_index] = cur_father_pos;
      Q[step_index] = Q[cur_father_pos];
    };

    for (int4 i = 0; i < new_L_size; i++)
    {
      segments->SetCurSegCutBE(L[i]);
      FindInt(segments,father_last_step, step_index, location[i]);//location[i] should contain the location of the segment L[i]
    }
    _step_index = step_index;
    return new_L_size;

  };

  template<class SegmentsColl>
  int4 FindR(SegmentsColl *segments, int4 father_first_step, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index, ProgramStackRec *stack_pos, int4 Size, int4 call_numb)
  {
    segments->SetCurStripe(ENDS[interval_left_index], ENDS[interval_right_index]);
    if (interval_right_index - interval_left_index == 1)
        return Size>1 ? CFAST::SearchInStrip(segments, ladder_start_index, Size) : Size;
     
    ProgramStackRec stack_rec(ladder_start_index);
    bool use_opt = (ladder_start_index - father_first_step > big_staircase_threshold);
    if (use_opt)
    {// use optimal variant if father staircase is big
//      segments->SetCurVLine(ENDS[interval_left_index]);
      Size = Split(segments, interval_right_index, father_first_step, stack_rec.Q_pos, Size);
      stack_pos = stack_rec.Set(stack_pos, interval_right_index);
      father_first_step = ladder_start_index + inherit_offset;
    }
    else
    {// use fast variant if new staircase is small 
      Size = CFAST::Split(segments,interval_right_index,stack_rec.Q_pos, Size);
      if ((ladder_start_index < stack_rec.Q_pos))
      {
        for (int4 i = ladder_start_index + 1; i <= stack_rec.Q_pos; ++i)
          father_loc[i] = undef_loc;
        stack_pos = stack_rec.Set(stack_pos, interval_right_index);
        father_first_step = ladder_start_index + inherit_offset;
      }
    }
    if ((ladder_start_index < stack_rec.Q_pos) && (call_numb < max_call)) //if found a lot of intersections repeat optFindR
      Size = FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
    else //cut stripe 
    {
      uint4 m = (interval_left_index + interval_right_index) >>1;
      if (call_numb > 1)
      {// if L contains a lot of segments then cut on two parts
        Size = FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
        Size = InsDel(segments,m, stack_pos, Size);
        Size = FindR(segments,father_first_step, stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
      }
      else
      {// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
        uint4 q = (interval_left_index + m) >>1;
        if (interval_left_index != q) {
          Size = FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_index, q, stack_pos, Size, 0);
          Size = InsDel(segments, q, stack_pos, Size);
        }
        if (q != m) {
          Size = FindR(segments, father_first_step, stack_rec.Q_pos, q, m, stack_pos, Size, 0);
          Size = InsDel(segments, m, stack_pos, Size);
        }
        q = (interval_right_index + m) / 2;
        if (q != m) {
          Size = FindR(segments, father_first_step, stack_rec.Q_pos, m, q, stack_pos, Size, 0);
          Size = InsDel(segments, q, stack_pos, Size);
        }
        Size = FindR(segments, father_first_step, stack_rec.Q_pos, q, interval_right_index, stack_pos, Size, 0);
      }

    }
    if (ladder_start_index >= stack_rec.Q_pos) return Size;
//    segments->SetCurVLine(ENDS[interval_right_index]);
    if (use_opt)
      return Merge(segments, interval_left_index, ladder_start_index, stack_rec.Q_pos, Size);
    return CFAST::Merge(segments, interval_left_index, ladder_start_index, stack_rec.Q_pos, Size);


  };

};

#pragma once
