#ifndef OPTIMAL_FINDER_FOR_SEGMENT_INTERSECTION
#define OPTIMAL_FINDER_FOR_SEGMENT_INTERSECTION
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

#include "fast_finder.h"
#include <algorithm>
#include <vector>
#include <thread>
#include <cassert>


class COptimalIntFinder : public CFastIntFinder
{
public:
  using  imp_T = CommonImpl;
  using CFAST = CFastIntFinder;
 
  ~COptimalIntFinder() { FreeMem(); };
  template<class SegmentsColl>
  void find_intersections(SegmentsColl &segments)
  {
    //AllocMem
    len_of_Q = LR_len + LR_len / (inherit_each - 1) + inherit_each + 2;
    DECL_RAII_ARR(L, LR_len+1); 
    ++L;//to have one cell before L for sentinel
    DECL_RAII_ARR(R, LR_len+1);
    ++R;//to have one cell before R for sentinel
    DECL_RAII_ARR(Q, len_of_Q);
    DECL_RAII_ARR(father_loc, len_of_Q);
    std::fill_n(father_loc, len_of_Q, undef_loc);

    constexpr int4 bottom_index = inherit_each + 1;
    ProgramStackRec stack_rec(bottom_index, nTotX);  //need to be initialized this way
    segments.SetCurStripeRight(0);
    InsDel(segments, 0, &stack_rec);
    FindR(segments, bottom_index + 1, bottom_index, 0, nTotX - 1, &stack_rec, 0);
    InsDel(segments, nTotX - 1, &stack_rec);
    segments.IntersectionsFindingDone();
 }

  template<class SegmentsColl>
  void FindIntI(SegmentsColl& segments, int4 cur_seg, ProgramStackRec* stack_pos) const
  {
    auto r_index = segments.GetSegR(cur_seg);
    while (stack_pos->right_bound <= r_index)stack_pos = stack_pos->prev;// go from bottom to top and find staircase to start
    if (stack_pos->prev == nullptr)return;
    int4  l = undef_loc, r, m, QE = stack_pos->Q_pos + 1;
    int4 QB = undef_loc + 1;//just to differ from l at first loop
    for (stack_pos = stack_pos->prev; stack_pos != nullptr; stack_pos = stack_pos->prev)// for all staircases above
    {
      if ((big_staircase_threshold > inherit_each) && (l == QB)) //if location is below first stair we don't have data in father_loc for that location
        // because the number of locations is greater (by one) than number of stairs 
        // (for example one stair has two locations (below and above), but we keep in father_loc data for only one (above))
        // so we use data from next location QE==l+1
      {
        l = QB = stack_pos->Q_pos;// set lower location bound to the location below first step of current staircase
        r = father_loc[QE] ? QB + inherit_each : QE;// father_loc[QE] !=undef_loc (i.e. 0) means staircase was created by optSplit 
        //and we use range [QB;QB+inherit_each] to search, otherwise it was created by Split and full range search [QB;QE] is used
      }
      else
      {
        QB = stack_pos->Q_pos;
        m = abs(father_loc[l]);//using father_loc get approximate location in parent staircase; 
        // line above sometimes executed when l == undef_loc, so we must initialize father_loc[undef_loc]=undef_loc in AllocMem()
        //to keep l unchanged in this case; undef_loc was chosen to be zero so abs(undef_loc)==undef_loc
        //otherwise we need additional checks
        l = MAX(QB, m);//here we use a fact that always undef_loc < QB and if m==undef_loc or m<QB l should be QB
        r = l + inherit_each;
        if ((m == undef_loc) || (r > QE))r = QE;
        //if m == undef_loc we use [QB;QE] as a range
      }
      auto len = r - l;
      auto Ql = Q + l;
      while (len > 1) // binary search
      {
        m = len / 2; // 
        if (segments.UnderCurPoint(Ql[m]))
          Ql += m;
        len -= m;
      }
      l = Ql - Q;
      FindInt(segments, QB, QE - 1, l);
      QE = QB + 1; //set upper location bound to the location above last step of prev (which will be current) staircase
    };
  };


private:
  static constexpr int4 undef_loc = 0;

  static constexpr int4 max_call = 32; //max number of sequential recursive call (opt)FindR before dividing current strip

  static constexpr int4 inherit_each = 15; // defines how often in optimal algorithm stairs are inherited,  one inherited per inherit_each
  static constexpr int4 inherit_offset = inherit_each / 2; // first stair to inherit; must be in [1..inherit_each-1]
                                                       // so inherited stair positions are: inherit_offset,inherit_offset+inherit_each,inherit_offset+2*inherit_each,...
  static constexpr int4 big_staircase_threshold = 64;// 1024;//used in optFindR

  int4 *father_loc = nullptr;


  static bool is_original(int4 v) { return v < 1; };

  template<class SegmentsColl>
  int4 FindInt(SegmentsColl& segments, int4 qb, int4 qe, int4 l) const
  {
    
    int4 c = l;
    while (c>qb)
    {
      if (is_original(father_loc[c]))
      { //if original stair try to find and register intersection
        if (!segments.FindCurSegIntDownWith(Q[c])) break;
      }
      else //if inherited stair just check if intersects
        if (!segments.IsIntersectsCurSegDown(Q[c])) break;
      c--;
    }

    constexpr const bool line_seg = (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) == _Coll_flag_state::state_true);
    
    if (line_seg && (c != l))
      return l-c;
    ++c; ++l;
    while (l <= qe)
    {
      if (is_original(father_loc[l]))
      {
        if (!segments.FindCurSegIntUpWith(Q[l])) break;
      }
      else
        if (!segments.IsIntersectsCurSegUp(Q[l])) break;
      l++;
    }
    return l-c;
  };


  template<class SegmentsColl>
  void InsDel(SegmentsColl& segments, uint4 end_rank, ProgramStackRec* stack_pos)  {
    segments.InsDel(*this, L_size, L, end_rank, stack_pos);
  }

  template<class SegmentsColl>
  void Merge(SegmentsColl &segments, uint4 LBoundIdx, int4 QB, int4 QE)
  {
    auto Size = L_size;
    uint4 cur_R_pos = 0, new_size = 0;
    int4 cur_stair = QB;
    auto _R = L, _L = R;
    int4 cur_seg = _R[cur_R_pos];
    while ((cur_stair<QE) && (cur_R_pos<Size))
    {
      if (segments.RBelow(cur_seg, Q[cur_stair + 1]))
      {
        if (segments.GetSegL(cur_seg) >LBoundIdx)
        {
          segments.SetCurSegCutEnd(cur_seg);
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
      if (segments.GetSegL(cur_seg) >LBoundIdx)
      {
        segments.SetCurSegCutEnd(cur_seg);
        FindInt(segments, QB, QE, QE);
      }
      _L[new_size++] = cur_seg;
      cur_seg = _R[++cur_R_pos];
    }
    while (cur_stair<QE)
      if (is_original(father_loc[++cur_stair]))
        _L[new_size++] = Q[cur_stair];
    L = _L; R = _R;
    L_size = new_size;
  };

  template<class SegmentsColl>
  void Split(SegmentsColl& segments, uint4 RBoundIdx, int4 cur_father_pos, int4& _step_index)
  {
    auto Size = L_size;
    auto step_index = _step_index;
    int4 father_last_step = step_index;
    uint4 new_L_size = 0, cur_L_pos = 0;
    int4 cur_seg;
    auto location = R;
    while ((cur_L_pos<Size) && (cur_father_pos <= father_last_step))
    {
      cur_seg = L[cur_L_pos];
      segments.SetCurSegCutBE(cur_seg);
      if (segments.LBelow( cur_seg, Q[cur_father_pos])) // current segment below candidate to inherit
      {
        if ((segments.GetSegR(cur_seg) < RBoundIdx)// current segment not covering strip
          || (step_index>father_last_step) && segments.IsIntersectsCurSegDown(Q[step_index])// or intersects last stair
          || segments.IsIntersectsCurSegUp(Q[cur_father_pos])) //or intersects candidate to inherit 
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
      segments.SetCurSegCutBE(cur_seg);
      if ((segments.GetSegR(cur_seg) < RBoundIdx) ||//segment is not covering current stripe
        (step_index>father_last_step) && segments.IsIntersectsCurSegDown(Q[step_index]))// if it intesects last stair
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
    long long n_int = 0;
    for (decltype(new_L_size) i = 0; i != new_L_size; i++)
    {
      segments.SetCurSegCutBE(L[i]);
      n_int +=FindInt(segments,father_last_step, step_index, location[i]);//location[i] should contain the location of the segment L[i]
    }
    dont_cut_stripe = n_int > new_L_size + cut_margin;
    _step_index = step_index;
    L_size = new_L_size;

  };

  template<class SegmentsColl>
  void FindR(SegmentsColl &segments, int4 father_first_step, int4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec *stack_pos, uint4 call_numb,uint4 _max_call=max_call)
  {
    segments.SetCurStripe(interval_left_rank, interval_right_rank);
    if (interval_right_rank - interval_left_rank == 1) {
      if (L_size > 1) CFAST::SearchInStrip(segments, ladder_start_index);
      return;
    }
 
     
    ProgramStackRec stack_rec(ladder_start_index);
    bool use_opt =  (ladder_start_index - father_first_step > big_staircase_threshold);
    if (use_opt)
    {// use optimal variant if father staircase is big
//      segments.SetCurVLine(ENDS[interval_left_rank]);
      Split(segments, interval_right_rank, father_first_step, stack_rec.Q_pos);
      stack_pos = stack_rec.Set(stack_pos, interval_right_rank);
      father_first_step = ladder_start_index + inherit_offset;
    }
    else
    {// use fast variant if new staircase is small 
      stack_rec.Q_pos+=CFAST::Split(segments,Q+ stack_rec.Q_pos,interval_right_rank);
      if ((ladder_start_index < stack_rec.Q_pos))
      {
        for (int4 i = ladder_start_index + 1; i <= stack_rec.Q_pos; ++i)
          father_loc[i] = undef_loc;
        stack_pos = stack_rec.Set(stack_pos, interval_right_rank);
        father_first_step = ladder_start_index + inherit_offset;
      }
    }
    if (dont_cut_stripe && (call_numb < _max_call)) //if found a lot of intersections repeat optFindR
      FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_rank, interval_right_rank, stack_pos, call_numb + 1,_max_call);
    else //cut stripe 
    {
      uint4 m = (interval_left_rank + interval_right_rank) >>1;
      if (call_numb > 1)
      {
        _max_call -=  (_max_call > 4) * 2;
        // if L contains a lot of segments then cut on two parts
        FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_rank, m, stack_pos, 0, _max_call);
        InsDel(segments,m, stack_pos);
        FindR(segments,father_first_step, stack_rec.Q_pos, m, interval_right_rank, stack_pos, 0, _max_call);
      }
      else
      {// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
        _max_call -=  (_max_call > 8) * 4;

        uint4 q = (interval_left_rank + m) >>1;
        if (interval_left_rank != q) {
          FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_rank, q, stack_pos, 0, _max_call);
          InsDel(segments, q, stack_pos);
        }
        if (q != m) {
          FindR(segments, father_first_step, stack_rec.Q_pos, q, m, stack_pos, 0, _max_call);
          InsDel(segments, m, stack_pos);
        }
        q = (interval_right_rank + m) / 2;
        if (q != m) {
          FindR(segments, father_first_step, stack_rec.Q_pos, m, q, stack_pos, 0, _max_call);
          InsDel(segments, q, stack_pos);
        }
        FindR(segments, father_first_step, stack_rec.Q_pos, q, interval_right_rank, stack_pos, 0, _max_call);
      }

    }
    if (ladder_start_index >= stack_rec.Q_pos) return;
    
    segments.SetCurStripeLeft(interval_left_rank);

    if(use_opt)
      Merge(segments, interval_left_rank, ladder_start_index, stack_rec.Q_pos);
    else
      CFAST::Merge(segments, interval_left_rank, ladder_start_index, stack_rec.Q_pos);
  };

};

#endif
