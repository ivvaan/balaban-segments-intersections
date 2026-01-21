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

#include "segments.h"
#include "registrator.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <numeric> 

namespace {
  template<class T>
  T& unmove(T&& t) { return t; };

  enum _Stages {
    stage_split=0, split_up, stage_bubble, stage_merge
  };

}

template<typename Real>
Real get_max_gap_middle(uint4 N, Real arr[]) {
  // function to find a middle of maximum gap between arr elements, Real is like double 

  assert(N > 1);
  std::sort(arr, arr + N);

  // Find the first maximum gap and return its midpoint.
  uint4 best_i = 1;
  auto _arr = arr - 1;
  Real max_gap = arr[best_i] - _arr[best_i];

  for (uint4 i = 2; i < N; ++i) {
    if (Real g = arr[i] - _arr[i]; g > max_gap) {
      max_gap = g;
      best_i = i;
    }
  }

  return (arr[best_i] + _arr[best_i]) / Real(2);
}


template<class IntersectionRegistrator>
class CIntegerSegmentCollection
{
  using CTHIS = CIntegerSegmentCollection< IntersectionRegistrator>;
public:
  static constexpr _Coll_flag_state get_coll_flag(_Coll_flags flag) {
    if (flag == _Coll_flags::line_segments)
      return _Coll_flag_state::state_false;// degenerate segments should be processed as nonlinear

    return _Coll_flag_state::state_unimplemented;
  }

  static  bool is_last(uint4 pt)
  {
    return pt & 1;
  };
  static  bool is_first(uint4 pt)
  {
    return (pt & 1) == 0;
  };
  static  uint4 get_segm(uint4 pt)
  {
    return pt >> 1;
  };
  static  uint4 first_point(uint4 s)
  {
    return s << 1;
  };
  static  uint4 last_point(uint4 s)
  {
    return (s << 1) + 1;
  };
  static  uint4 other_point(uint4 pt)
  {
    return pt^1;
  };


private:
  struct ListBounds {
    uint4 beg = 0;
    uint4 end = 0;
  };

  //===========================================================================
#include "remaper.h"
  //===========================================================================

public:

  ListBounds* x_collide_points=nullptr;
  int4* tmp = nullptr;

  std::pair<uint4, uint4> get_dublicate_stat(uint4 len, uint4* arr) const {
    uint4 unique_count = 0;
    uint4 duplicated_elements_count = 0;
    
    for (uint4 i = 0,j; i < len; i = j) {
      j = i + 1;
      while (j < len && pts[arr[j]].x == pts[arr[i]].x) 
        ++j;
      ++unique_count;
      if (j >  i + 1) {
        ++duplicated_elements_count;
      }
    }
    return { unique_count , duplicated_elements_count };
  }


  void PrepareEndpointsSortedList(uint4* epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = GetSegmNumb() << 1;
    std::iota(epoints, epoints+NN, 0);
    auto is_below = [pts = pts, coll = collection](uint4 pt1, uint4 pt2) {
      if (pts[pt1].x != pts[pt2].x)
        return pts[pt1].x < pts[pt2].x;
      if (pts[pt1].y != pts[pt2].y)
        return pts[pt1].y < pts[pt2].y;
      auto le1 = is_last(pt1), le2 = is_last(pt2);
      if (le1 != le2)
        return le1 < le2;
      auto& s1 = coll[get_segm(pt1)];
      auto& s2 = coll[get_segm(pt2)];
      auto angle_like = s1.shift % s2.shift;
      if (angle_like != 0)
        return le1 != (0 < angle_like);
      return (pt1 < pt2);
      };
      //    std::sort(epoints, epoints + NN, std::ref(*this));
    std::sort(epoints, epoints + NN, is_below);
  };

  bool is_multiple(uint4 pt) const {
    return pt >= (nSegments << 1);
  };
  auto mark_multiple(uint4 pt) const {
    return pt + (nSegments << 1);
  };
  auto unmark_multiple(uint4 pt) const {
    return pt - (nSegments << 1);
  };

  PrepareResult Prepare()
  {
    Reset();
    uint4 Nn = GetSegmNumb();
    if (Nn == 0) return {};
    uint4 NN = Nn << 1;
    //RAII_ARR(uint4,tmpENDS,NN);
    auto __tmpENDS__ = std::make_unique<uint4[]>(NN); 
      uint4* tmpENDS = __tmpENDS__.get();
    tmp =new int4[NN];
    PrepareEndpointsSortedList(tmpENDS);

    auto [nX,nV]=get_dublicate_stat(NN, tmpENDS);
    auto nAll = NN+nV;
    //auto nD = nAll  - nX;
    //nV - number of unique non unique X,e.g. for X={0,1,3,3,3,4,4} 3 and 4 are non unique, the number is 2
    //nD - number of non unique non unique X, i.e. 5: 3,3,3,4,4

    if (nV != 0) {
      ENDS = new uint4[nAll];
      x_collide_points = new ListBounds[nV];
      uint4 D_last = nX;
      for (uint4 i = 0, xcp_pos = 0, ne_pos = 0, j; i < NN; i = j) {
        j = i + 1;
        while (j < NN && pts[tmpENDS[j]].x == pts[tmpENDS[i]].x)
          ++j;
        auto len = j - i;
        if (len > 1) {
          x_collide_points[xcp_pos] = { D_last,D_last + len };
          std::copy(tmpENDS + i, tmpENDS + j, ENDS + D_last);
          ENDS[ne_pos] = mark_multiple(xcp_pos);
          std::fill(tmp + i, tmp + j, ne_pos);
          ++xcp_pos;
          D_last += len;
        }
        else {
          ENDS[ne_pos] = tmpENDS[i];
          tmp[i] = ne_pos;
        }
        ++ne_pos;
      }
      assert(D_last == nAll);
    }
    else {
      std::iota(tmp, tmp + NN, 0);
      ENDS = tmpENDS;
      __tmpENDS__.release();
    };

    seg_L_rank = new uint4[Nn];
    seg_R_rank = new uint4[Nn];
    uint4 max_segm_on_vline = 0, nsegm_on_vline = 0;
    double avr = 0;
    for (uint4 i = 0; i < NN; ++i) {
      if (is_last(tmpENDS[i])) {
        seg_R_rank[get_segm(tmpENDS[i])] = tmp[i];
        --nsegm_on_vline;
      }
      else {
        seg_L_rank[get_segm(tmpENDS[i])] = tmp[i];
        ++nsegm_on_vline;
        if (nsegm_on_vline > max_segm_on_vline) max_segm_on_vline = nsegm_on_vline;
      }
      avr += nsegm_on_vline;
    }
    avr /= (double)nX;
    nTotX = nX; nCollideX = nV;
    return { nX, max_segm_on_vline, avr };
  }

  bool is_remapped() const {
    return is_collection_remapped;
  }

  int4 GetMaxCollideXIdx() const {
    if (nCollideX < 1)
      return -1;
    auto less=[](ListBounds a, ListBounds b) {
      return (a.end - a.beg) < (b.end - b.beg);
      };
    return std::max_element(x_collide_points, x_collide_points + nCollideX,less)
      - x_collide_points;
  }

  template<class IntersectionFinder, class ProgramStackRec>
  void InsDelOrdinary(IntersectionFinder& i_f, uint4& L_size, int4* L, uint4 pt, ProgramStackRec* stack_pos)
  {
    auto sn = get_segm(pt);
    if (is_last(pt)) // if endpoint - remove
    {
      auto last = std::remove(L, L + L_size, sn);
      assert(last + 1 == L + L_size);
      --L_size;
    }
    else// if startpoint - insert
    {
      if (stack_pos->isnot_top()) {
        SetCurSegAndPoint(sn);
        i_f.FindIntI(*this, sn, stack_pos);// get internal intersections
      }
      else {
        SetCurPointAtBeg(sn);//sets collection current point at the begin of the segment sn
      }
      int4 i = L_size;
      for (auto _L = L - 1; (i != 0) && (!UnderCurPoint(_L[i])); --i)
        L[i] = _L[i];
      L[i] = sn; ++L_size;
    }
  }



  template<class IntersectionFinder, class ProgramStackRec>
  void InsDel(IntersectionFinder& i_f, uint4& L_size, int4* L, uint4 end_rank, ProgramStackRec* stack_pos) {
    auto pt = PointAtRank(end_rank);
    if (is_multiple(pt)) {
      DelStep(end_rank, L_size, L);

      auto actual_pt = unmark_multiple(pt);
      auto [f, l] = get_pt_list_bounds(actual_pt);
      ReorderStep(pts[ENDS[f]].x, L_size, L);
      AllIntCurLine(f, l, L_size, L);
      InsStep(f, l, L_size, L);
      if (stack_pos->isnot_top())// find intersections for all inserted and vertical  
        //segments in the staircases above (in the stack)
        for (auto i = f; i < l; ++i) {
          auto cur_pt = ENDS[i];
          if (is_first(cur_pt)) {
            auto sn = get_segm(cur_pt);
            SetCurSegAndPoint(sn);
            i_f.FindIntI(*this, sn, stack_pos);
          }
        }
    }
    else
      InsDelOrdinary(i_f, L_size, L, pt, stack_pos);
  }

  ListBounds get_pt_list_bounds(uint4 pt) const {
    return x_collide_points[pt];
  }

  uint4 rank_to_rank(uint4 multiple_rank) {
    auto pt = PointAtRank(multiple_rank);
    auto seg = get_segm(pt);
    return is_first(pt)?seg_L_rank[seg]:seg_R_rank[seg];
  }
  void DelStep(uint4 end_rank, uint4& L_size, int4* L)
  {
    uint4 new_size = 0;
    for (uint4 i = 0; i < L_size; ++i)
      if (seg_R_rank[L[i]]!= end_rank) {
        L[new_size++] = L[i];
      }
    L_size = new_size;
  }

  template <class Comp, class T>
  static auto merge_two_ranges(
    const T* it1, const T* last1,
    const T* it2, const T* last2,
    T* out,  Comp comp   // comp(x_from_1, y_from_2)
  ) {
    while (it1 != last1 && it2 != last2) 
      *out++ = comp(*it1, *it2) ? *it1++ : *it2++;
    while (it1 != last1)
      *out++ = *it1++;
    while (it2 != last2)
      *out++ = *it2++;
    return out; 
  }

  void InsStep(uint4 f, uint4 l, uint4& L_size, int4* L)
  {
    auto first_to_insert = tmp;
    auto last_to_insert = first_to_insert;
    auto first_L = tmp + GetSegmNumb();
    auto last_L = std::copy(L, L + L_size, first_L);
    for (uint4 i = f; i < l; ++i) {
      auto cur_pt = ENDS[i];
      if (is_first(cur_pt)) {
        auto s = get_segm(cur_pt);
        if (!collection[s].is_vertical())
          *(last_to_insert++)  = s;
      }
    }
    auto is_below = [coll = this->collection, pts=this->pts](int4 sL, int4 sI) {
      auto& s = coll[sL];
      auto cmp = s.under(pts[first_point(sI)]);
      if (cmp != std::strong_ordering::equal)
        return cmp == std::strong_ordering::less;
      return 0 < s.shift % coll[sI].shift;
      };
    
    auto last_res = merge_two_ranges(first_L, last_L, first_to_insert, last_to_insert, L, is_below);
    L_size = last_res - L;
  }


  void ReorderStep(int4 x, uint4 L_size, int4* L) {
    for (uint4 i = 0, j; i < L_size; i = j) {
      j = i + 1;
      while (j < L_size && 0==(collection[L[i]].YAtX_frac(x) <=> collection[L[j]].YAtX_frac(x)))
        ++j;
      if (j > i + 1) {//reverse order for sublists of equal Y segments (passing segments should be in reverse order)
        std::reverse(L + i, L + j);
      };
    };

  };

  void AllIntCurLine(uint4 f, uint4 l, uint4& L_size, int4* L)//all intersections with vertical segments in the current line are registered
  {
    auto vertical_segments = tmp;//temporary storage for vertical segment indices
    uint4 v_size=0;
    auto non_vertical_pts = tmp + GetSegmNumb();//temporary storage for non vertical segment endpoints
    uint4 nv_size = 0;

    for (uint4 i = f; i < l; ++i) {
      auto cur_pt = ENDS[i];
      auto s = get_segm(cur_pt);
      if (collection[s].is_vertical()) {
        if (is_first(cur_pt))
            vertical_segments[v_size++] = s;
      }
      else
        non_vertical_pts[nv_size++] = cur_pt;
    };
 
    auto prev_vert = vertical_segments[0];
    for (uint4 i = 1; i < v_size; ++i) {//vertical - vertical intersections
      auto cur_vert = vertical_segments[i];
      if (pts[first_point(cur_vert)].y == pts[last_point(prev_vert)].y)
        register_pair(this, prev_vert, cur_vert);
      prev_vert = cur_vert;
    }
     
    for (uint4 i = 0, j = 0; j < L_size && i < v_size; ++i) {// vertical - passing intersections
      auto s = vertical_segments[i];
      //skip segments below the vertical segment
      while (j < L_size && collection[L[j]].exact_under(pts[first_point(s)]))
        j++;
      //register segments intersecting the vertical segment except those exactly ending at its upper end
      for (; j < L_size && collection[L[j]].exact_under(pts[last_point(s)]);++j)
        register_pair(this, L[j], s); 
      //register segments passing exactly through the upper end of the vertical segment
      //keep j unchanged for the next vertical segment beginning at the same point
      for (auto k = j; k < L_size && collection[L[k]].exact_on(pts[last_point(s)]);++k) 
        register_pair(this, L[k], s);
    };


    for (uint4 i = 0, j = 0; j < nv_size && i < v_size; ++i) {// vertical - endpoint intersections
      auto s = vertical_segments[i];
      //skip segments below the vertical segment
      while (j < nv_size && pts[non_vertical_pts[j]].y<pts[first_point(s)].y)
        j++;
      //register segments intersecting the vertical segment except those exactly ending at its upper end
      for (; j < nv_size && pts[non_vertical_pts[j]].y < pts[last_point(s)].y; ++j)
        register_pair(this, get_segm(non_vertical_pts[j]), s);
      //register segments exactly ending at the upper end of the vertical segment
      //keep j unchanged for the next vertical segment beginning at the same point
      for (auto k=j; k < nv_size && pts[non_vertical_pts[k]].y == pts[last_point(s)].y; ++k)
        register_pair(this, get_segm(non_vertical_pts[k]), s);
    };
    LastFirstInt(nv_size, non_vertical_pts);
  }


  void LastFirstInt(uint4 len, int4 * nv_pts)//all intersections with vertical segments in the current line are registered
  {
    for (uint4 i = 0, j; i < len; i = j) {
      j = i + 1;
      while (j < len && pts[nv_pts[j]].y == pts[nv_pts[i]].y)
        ++j;
      if (j > i + 1) {
        auto first_beg = nv_pts + i;
        auto all_ends = nv_pts + j;
        auto first_end = std::find_if(first_beg, all_ends, +[](uint4 pt) { return is_last(pt); });
        assert(all_ends == std::find_if(first_end, all_ends, +[](uint4 pt) { return is_first(pt); }));
        if (first_end != all_ends)
          for (auto beg_pt = first_beg; beg_pt < first_end; ++beg_pt)
            for (auto end_pt = first_end; end_pt < all_ends; ++end_pt)
              register_pair(this, get_segm(*beg_pt), get_segm(*end_pt));
      };
    };
  };

  void Reset()
  {
    MY_FREE_ARR_MACRO(ENDS);
    MY_FREE_ARR_MACRO(seg_L_rank);
    MY_FREE_ARR_MACRO(seg_R_rank);
  }



  static void reg_pair(CIntegerSegmentCollection* c, uint4 s1, uint4 s2) {
    c->remaper.registrator->register_pair(s1, s2);
  };
  static void reg_remapped_pair(CIntegerSegmentCollection* c, uint4 s1, uint4 s2) {
    c->remaper.register_pair(s1, s2);
  };
  uint4  GetSegmNumb() const {
    return nSegments;
  };
  uint4 GetSegR(uint4 sn) const {
    return seg_R_rank[sn];
  };
  uint4 GetSegL(uint4 sn) const {
    return seg_L_rank[sn];
  };
  uint4 PointAtRank(uint4 rank) const {
    return ENDS[rank];
  }

  bool is_right_pt_on_bound(uint4 s) const {
    return stripe_right_rank == seg_R_rank[s];
  }
  bool is_left_pt_on_bound(uint4 s) const {
    return stripe_left_rank == seg_L_rank[s];
  }

  void SetCurStripe(uint4 left_rank, uint4 right_rank) {
    stripe_left_rank = left_rank;
    stripe_right_rank = right_rank;
    B = XAtRank(left_rank);
    E = XAtRank(right_rank);
  };

  void SetCurStripeLeft(uint4 left_rank) {
    stripe_left_rank = left_rank;
    B = XAtRank(left_rank);
  };

  void SetCurStripeRight(uint4 right_rank) {
    stripe_right_rank = right_rank;
    E = XAtRank(right_rank);
  };


  void SetCurPointAtBeg(uint4 s) {
    cur_point_seg = s;
    cur_point = collection[s].BegPoint();
  };

  void SetCurSeg(uint4 s){
    cur_seg_idx = s;
    cur_seg = collection[s];
    is_rstump=false;
  };

  void SetCurSegCutBE(uint4 s){
    SetCurSeg(s);
    stage = _Stages::stage_split;
    is_rstump = E < cur_seg.ex();
    active_end = cur_seg.EndPoint();
    //cur_seg_pt_on_right_bound = is_right_pt_on_bound(s);
  };

  void SetCurSegCutBeg(uint4 s) {
    SetCurSeg(s);
    stage = _Stages::stage_split;
    active_end = cur_seg.EndPoint();
  };

  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    stage = _Stages::stage_merge;
    active_end = cur_seg.BegPoint();
    cur_seg_pt_on_right_bound = is_right_pt_on_bound(s);
  };

  void SetCurSegAndPoint(uint4 s) {
    SetCurPointAtBeg(s);

    SetCurSeg(s);
    stage = _Stages::stage_bubble;
    active_end = pts[last_point(s)];
  };

  enum _Int_reg {
    // segments are not intersecting
    no_int = 0,
    // segments are intersecting, but no registration
    int_no_reg = 1,
    //Segments are intersecting and must be registered
    int_reg = int_no_reg + true //i.e. int_reg == 2 
  };

  constexpr bool static way_up = true;
  constexpr bool static way_down = false;
  // RStump - right stump of cur_seg
  int4 WayDownRStumpIntersect(int4 s_) const { // s_  is_upper than cur_seg right stump
    return RStumpIntersect<way_down>(s_);
  };
  int4 WayUpRStumpIntersect(int4 s_) const { // s_  is_under  cur_seg right stump
    return RStumpIntersect<way_up>(s_);
  };

  int4 BoundIntersect(int4 s_) const {
    assert(stage == _Stages::stage_split || stage == _Stages::stage_merge);
    bool is_left_bound = (stage == _Stages::stage_split);
    auto X = is_left_bound ? B : E;
    auto cmp = collection[s_].YAtX_frac(X) <=> cur_seg.YAtX_frac(X);
    //  on the left bound always must be registered
    // on the right bound no registration except when the bound is defined by the right point of s_
    bool need_reg = is_left_bound || cur_seg_pt_on_right_bound || is_right_pt_on_bound(s_);
    if (cmp == std::strong_ordering::equal)
      return  _Int_reg::int_no_reg + need_reg;// and must be registered
    return _Int_reg::no_int;
  };

  template<bool check_way_up_s_is_under_cur_seg>// s_  is_under or is_upper than cur_seg right stump
  int4 RStumpIntersect(int4 s_) const {
    assert(stage == _Stages::stage_split);
    auto& s = collection[s_];

    auto cmpE = s.YAtX_frac(E) <=> cur_seg.YAtX_frac(E);// cmpE==1 is s.YAtX_frac(E) > cur_seg.YAtX_frac(E)
    if (cmpE == std::strong_ordering::equal) {// inresection on the right bound
      //but no registration except when right bound is defined by the right point of s_
      return _Int_reg::int_no_reg + is_right_pt_on_bound(s_);
    }

    if constexpr (check_way_up_s_is_under_cur_seg)
      return (cmpE == std::strong_ordering::less) * _Int_reg::int_reg; //s < cur_seg
    else //check_way_down_s_is_upper_cur_seg
      return (cmpE == std::strong_ordering::greater) * _Int_reg::int_reg; //s > cur_seg
  };

  constexpr bool static left_bound = false;
  constexpr bool static right_bound = true;

  template<bool is_right>
  bool XBelow(int4 s_1, int4 s_2, int4 X) const {
    auto& s1 = collection[s_1];
    auto& s2 = collection[s_2];//for right bound s2 always a stair
    auto cmp = s1.YAtX_frac(X) <=> s2.YAtX_frac(X);
    if (cmp != std::strong_ordering::equal)
      return cmp == std::strong_ordering::less;
    auto prod = s1.shift % s2.shift;
    assert(prod != 0);
    return is_right == prod < 0;//intersections is on the right of X (prod>0 would be on the left)
  }

  bool LBelow(int4 s_1, int4 s_2)const{  //retuns if s1 below s2 at current vertical line
    return XBelow<left_bound>(s_1, s_2, B);
  };

  bool RBelow(int4 s_1, int4 s_2)const{  //retuns if s1 below s2 at current vertical line
    return XBelow<right_bound>(s_1, s_2, E);
  };

  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    if ((s1.ex() < s2.bx()) || (s2.ex() < s1.bx()))
      return false;//segments are not on vertical line and can't intersect
    return SSCurSegIntWith(s_);
  };


  template <bool reg = true>
  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (both are on some vertical line) and register them
  {
    //!!! can't process zero segments
    //!precondition segments share some x coords
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    if (s1.no_common_y(s2))
      return false;
    //at this moment, we can be sure that segments' x,y-projections intersect each other (1)
    // let l1 is the line through s1 and l2 is the line through s2
    auto delt = s2.BegPoint() - s1.BegPoint();
    auto prod = s1.shift % s2.shift, beg = s1.shift % delt;
    if (prod == 0) { // segments parallel
      if ( (beg == 0) && (s2.shift % delt == 0)) {// they are on the same line
        if constexpr (reg)register_pair(this, cur_seg_idx, s_);
        return true;
      }
      return false;//segments are not on one line and can't intersect
    }
    auto end = beg + prod;
    if (((beg > 0) != (end < 0)) && (beg != 0) && (end != 0))//begin and end of s2 are on the same side of l1 
        return false;

    //auto int_type = _IntType::s2_beg_int * (beg == 0) + _IntType::s2_end_int * (end == 0);

    //at this moment, we can be sure that s2 intesects l1 (2)
    beg = delt % s2.shift;
    end = beg - prod;// if beg==0 && end == 0 -> prod == 0 impossible here  
    if ((beg == 0) ||//begin of s1 lies on l2, taking (1)(2) to account 
        (end == 0) ||//end of s1 lies on l2, taking (2) to account
        ((beg > 0) != (end > 0))) {//begin and end of s1 are on diffenent sides l2

      //int_type += _IntType::s1_beg_int * (beg == 0) + _IntType::s1_end_int * (end == 0);

      if constexpr (reg)register_pair(this, cur_seg_idx, s_);
      return true;
    };
    return false;
  };


  template<bool is_way_up>
  bool IsCurSegIntWith(int4 s_) {//checks if cur_seg and s intersects
    auto cs = cur_seg_idx;
    int4 is_intersect;
    if (is_rstump) // can be in stage_split only
      return BoundIntersect(s_) || RStumpIntersect<is_way_up>(s_);
    if (stage != _Stages::stage_bubble) 
      return BoundIntersect(s_) || ActiveEndIntersect<is_way_up>(s_);
    // no need to check BoundIntersect in stage_bubble but check cur_point intersection
    return CurPointIntersect(s_)||ActiveEndIntersect<is_way_up>(s_);
  };

  bool IsIntersectsCurSegDown(int4 s_) //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return IsCurSegIntWith<way_down>(s_);
  };

  bool IsIntersectsCurSegUp(int4 s_) //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return IsCurSegIntWith<way_up>(s_);
  };


  template<bool is_way_up>
  auto FindCurSegIntWith(int4 s_){//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    auto cs = cur_seg_idx;
    int4 is_intersect;
    if (is_rstump) {//can be in stage_split only
      if ((is_intersect = BoundIntersect(s_)) || (is_intersect = RStumpIntersect<is_way_up>(s_))) {
        if (is_intersect == _Int_reg::int_reg)
          register_pair(this, cs, s_);
        return true;
      };
      return false;
    };
    if (stage != _Stages::stage_bubble) {
      if ((is_intersect = BoundIntersect(s_)) || (is_intersect = ActiveEndIntersect<is_way_up>(s_))) {
        if (is_intersect == _Int_reg::int_reg)
          register_pair(this, cs, s_);
        return true;
      };
      return false;
    };
    // no need to check BoundIntersect in stage_bubble but check cur_point intersection
    if ((is_intersect = CurPointIntersect(s_))||(is_intersect = ActiveEndIntersect<is_way_up>(s_))) {
      if (is_intersect == _Int_reg::int_reg)
        register_pair(this, cs, s_);
      return true;
    };
    return false;
  };

  template<bool is_way_up>
  auto FindCurSegIntWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    auto cs = cur_seg_idx;
    int4 is_intersect;
    int4 increment = is_way_up ? 1 : -1;
    if (is_rstump) {//can be in stage_split only
      while ((s_ != last) && ((is_intersect = BoundIntersect(*s_)) ||
        (is_intersect = RStumpIntersect<is_way_up>(*s_)))) {
        if (is_intersect == _Int_reg::int_reg)
          register_pair(this, cs, *s_);
        s_ += increment;
      }
      return s_;
    }
    if (stage != _Stages::stage_bubble) {
      while ((s_ != last) && ((is_intersect = BoundIntersect(*s_)) ||
        (is_intersect = ActiveEndIntersect<is_way_up>(*s_)))) {
        if (is_intersect == _Int_reg::int_reg)
          register_pair(this, cs, *s_);
        s_ += increment;
      }
      return s_;
    }
    // no need to check BoundIntersect in stage_bubble
    while ((s_ != last) && ((is_intersect = CurPointIntersect(*s_)) || 
      (is_intersect = ActiveEndIntersect<is_way_up>(*s_)))) {
      if (is_intersect == _Int_reg::int_reg)
        register_pair(this, cs, *s_);
      s_ += increment;
    }
    return s_;
  };

  auto FindCurSegIntDownWith(int4 s_) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_down>(s_);
  };

  auto FindCurSegIntUpWith(int4 s_) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_up>(s_);
  };

  auto FindCurSegIntDownWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_down>(s_, last);
  };

  auto FindCurSegIntUpWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_up>(s_, last);
  };


  bool UnderCurPoint(int4 s_) const { //returns true if s_ is under current point 
    //write as ActiveEndIntersect for consistency
    auto& s = collection[s_]; 
    auto pp = s.point_pos(cur_point);
    if(pp != 0)
      return pp < 0;//s_ placed under current point
    auto prod = collection[cur_point_seg].shift % s.shift;
    assert(prod != 0);
    return prod < 0;//if prod<0 then s_ is clockwise from cur_seg at cur_point, 
    // taking into account that cur_seg must be on the right of cur_point 
    // (current point is always the begin point of cur_seg),
    //s_ is under cur_seg
  };

  auto WayUpIntersect(int4 s_) const { //returns true if s is under active end
    return ActiveEndIntersect<way_up>(s_);
  };

  auto WayDownIntersect(int4 s_) const { //returns true if s is upper active end
    return ActiveEndIntersect<way_down>(s_);
  };

  template<bool is_way_up>
  auto ActiveEndIntersect(int4 s_) const { //
    auto pp = collection[s_].point_pos(active_end);
    if constexpr(is_way_up)
      return (pp <= 0)*_Int_reg::int_reg;//s_ placed under active_end
    else
      return (pp >= 0)*_Int_reg::int_reg;//s_ placed over active_end
  };

  auto CurPointIntersect(int4 s_) const { //
    return (collection[s_].point_pos(cur_point) == 0) * _Int_reg::int_reg;//current point placed on  s_
  };


  void clone(CIntegerSegmentCollection &c, IntersectionRegistrator *r)
  {
      clone_of = &c;
      nSegments = c.nSegments;
      ENDS = c.ENDS;
      seg_L_rank = c.seg_L_rank;
      seg_R_rank = c.seg_R_rank;
      collection = c.collection;
      //!!!!!clone remaper!!!!!!!!!!!!
      SetRegistrator(r);
  };

  void unclone() { 
    if (clone_of == nullptr)return; 
    collection = nullptr; 
    clone_of = nullptr; 
    ENDS = nullptr;
    seg_L_rank = nullptr;
    seg_R_rank = nullptr;
    //unclone remaper !!!!!!!!

  };

  void SortAt(uint4 rank, uint4 n, int4 *L)
  {
    SetCurStripeLeft(rank);
    std::sort(L, L + n, [this](int4 s1, int4 s2) {return LBelow(s1, s2); });
  };

  void SetRegistrator(IntersectionRegistrator* r)
  {
    assert(("integer collection registers only intersecting pairs:\n\
       intersection points don't have integer coords ",
      (IntersectionRegistrator::reg_type & _RegistrationType::point) == 0));
    remaper.registrator = r;
  };

  IntersectionRegistrator* GetRegistrator() { return  remaper.registrator; };

  CIntegerSegmentCollection(uint4 n, TLineSegment1* c, IntersectionRegistrator *r, int4 range)
  {
    assert(("integer collection registers only intersecting pairs:\n\
       intersection points don't have integer coords ",
      (IntersectionRegistrator::reg_type & _RegistrationType::point) == 0));
    remaper.NoRemapInit(n, c, r, range, *this);
    TurnRemapOn();
   /* is_collection_remapped = false;
    register_pair = reg_pair;
    N = remaper.get_N();
    collection = segments.data();
    pts = points.data();*/

  }

  void FromIntSegVect(std::vector<TIntegerSegment> &v, IntersectionRegistrator* r)
  {
    assert(("integer collection registers only intersecting pairs:\n\
       intersection points don't have integer coords ",
      (IntersectionRegistrator::reg_type & _RegistrationType::point) == 0));
    remaper.FromIntSegVect(v, r, *this);
    is_collection_remapped = false;
    register_pair = reg_pair;
    nSegments = remaper.get_N();
    collection = segments.data();
    pts = points.data();

  }


  bool TurnRemapOn() {
    is_collection_remapped = remaper.TurnRemapOn(*this);

#ifdef DEBUG_INTERSECTION_SET
    register_pair = reg_pair;
    is_collection_remapped = false;
#else
    register_pair = is_collection_remapped ? reg_remapped_pair : reg_pair;
#endif

    nSegments = remaper.get_N();
    collection = segments.data();
    pts = points.data();
    return is_collection_remapped;
  }

  TIntegerSegment& get_original_segment(uint4 s) {
    if (is_collection_remapped)
      return remaper.get_original_segment(s);
    return collection[s];
  }

  // Print difference vectors produced by multiset_difference.
  void print_difference_vectors(const std::vector<int_rec>& a_minus_b,
    const std::vector<int_rec>& b_minus_a,
    const char* nameA = "reg1",
    const char* nameB = "reg2")
  {
    printf("Pairs present in %s but not in %s (with multiplicity): %zu\n", nameA, nameB, a_minus_b.size());
    for (const auto& t : a_minus_b) {
      auto n1 = std::get<0>(t);
      auto n2 = std::get<1>(t);
      auto& s1 = get_original_segment(n1);
      auto& s2 = get_original_segment(n2);
      printf("%u %u:[(%i,%i)(%i,%i){%i,%i}] [(%i,%i)(%i,%i){%i,%i}]\n", n1, n2,
        s1.org.x, s1.org.y, s1.org.x + s1.shift.x, s1.org.y + s1.shift.y, s1.shift.x, s1.shift.y,
        s2.org.x, s2.org.y, s2.org.x + s2.shift.x, s2.org.y + s2.shift.y, s2.shift.x, s2.shift.y);
    }

    printf("Pairs present in %s but not in %s (with multiplicity): %zu\n", nameB, nameA, b_minus_a.size());
    for (const auto& t : b_minus_a) {
      auto n1 = std::get<0>(t);
      auto n2 = std::get<1>(t);
      auto& s1 = get_original_segment(n1);
      auto& s2 = get_original_segment(n2);
      printf("%u %u:[(%i,%i)(%i,%i){%i,%i}] [(%i,%i)(%i,%i){%i,%i}]\n", n1, n2,
        s1.org.x, s1.org.y, s1.org.x + s1.shift.x, s1.org.y + s1.shift.y, s1.shift.x, s1.shift.y,
        s2.org.x, s2.org.y, s2.org.x + s2.shift.x, s2.org.y + s2.shift.y, s2.shift.x, s2.shift.y);
    }
  }

  CIntegerSegmentCollection() {};
 
  CIntegerSegmentCollection(CIntegerSegmentCollection& coll, IntersectionRegistrator* r)
  {
    //clone(coll, r);
  }

  ~CIntegerSegmentCollection()
  {
    unclone();
    Reset();
   // MY_FREE_ARR_MACRO(SegR);
  }


 
  void coll_to_SVG(chostream* SVG_stream) {
    if (!SVG_stream)return;
    int4 n = MIN(max_SVG_items, nSegments);
    auto mmr = get_mmrect01(collection, n);
    *SVG_stream << "<svg height='100%' width='100%' viewBox='";
    *SVG_stream << mmr.ld.x << " " << mmr.ld.y << " "
      << mmr.rt.x - mmr.ld.x << " " << mmr.rt.y - mmr.ld.y << "' transform='scale(1, -1)'>\n";
    for (int4 i = 0; i < n; ++i)collection[i].write_SVG(i, SVG_stream);
  };
  auto &get_segments() {
    return segments;
  };

  auto XAtRank(uint4 rank) const {
    auto pt = PointAtRank(rank);
    return GetXEx(pt);
  };

private:
  auto GetXEx(uint4 pt) const {
    if (is_multiple(pt)){
      auto [f, l] = get_pt_list_bounds(unmark_multiple(pt));
      return pts[ENDS[f]].x;
    };
    return GetX(pt);
  };

  auto GetX(uint4 pt) const {
    assert(pt < 2 * nSegments);
    return pts[pt].x;
  };

  TIntegerSegment* collection = nullptr;
  TIntegerVect* pts = nullptr;
  //uint4* SegL = nullptr; 
  uint4* seg_L_rank = nullptr, * seg_R_rank = nullptr, * ENDS = nullptr;

  TIntegerSegment cur_seg;
  TIntegerVect cur_point, active_end;
  std::vector<TIntegerSegment> segments;
  std::vector<TIntegerVect> points;
  CRemaper remaper;
  decltype(reg_pair) *register_pair = reg_pair;
  CIntegerSegmentCollection *clone_of = nullptr;
  uint4 nSegments=0,nTotX=0,nCollideX=0;
  int4 B, E,  stripe_left_rank, stripe_right_rank;
  uint4 stage = _Stages::stage_split;

  uint4 cur_seg_idx = 0xFFFFFFFF;// , cur_point_idx = 0xFFFFFFFF;
  uint4 cur_point_seg = 0xFFFFFFFF;// , active_end_idx = 0xFFFFFFFF;
  //TPlaneVect x_transf, y_transf;
  bool is_rstump = false;
  bool is_collection_remapped = false;
  bool cur_seg_pt_on_right_bound = false;
//  bool rbelong_to_stripe=true;
};

// some string
