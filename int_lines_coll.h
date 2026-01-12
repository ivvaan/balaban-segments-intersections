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
  //static constexpr bool is_line_segments = false;//they can have more than one point of intersection with each other!!
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


private:

  //===========================================================================
  class CRemaper {

    struct RemappedRecord {
      uint4 beg = 0;
      uint4 end = 0;
    };

    static REAL get_rot_angle(uint4 n, TLineSegment1 sc[]) {
      REAL* arr;
      DECL_RAII_ARR(arr, n + 2);
      arr[0] = - M_PI / 2.;
      std::transform(sc, sc + n, arr,
        [](const TLineSegment1& s) {
          return std::atan2(s.shift.y, s.shift.x);
        });
      arr[n+1] = M_PI / 2.;
      auto max_gap_mid = get_max_gap_middle(n + 2, arr);
      return max_gap_mid + (max_gap_mid > 0 ? -M_PI / 2. : M_PI / 2.);
    }

  public:

    template<bool remove_zero_seg>
    uint4 int_seg_from_real(uint4 n, TLineSegment1 sc[], std::vector<TIntegerVect>& points) {
      //REAL angle = CRandomValueGen().GetRandomDouble();
      //REAL angle = get_rot_angle(n,sc);//try to rotate the scene to remove vertical segments
      //REAL si = std::sin(angle), co = std::cos(angle);

      REAL si = 0., co = 1.;
      auto mm_rect = get_rot_minmax(n, sc, si, co);
      auto x_transf = transf1D{ .shift = -mm_rect.ld.x - mm_rect.get_width() / 2.,
        .scale = 2.*range / mm_rect.get_width() };
      auto y_transf = transf1D{ .shift = -mm_rect.ld.y - mm_rect.get_height() / 2.,
        .scale = 2.*range / mm_rect.get_height() };
      auto coll_begin = set_size(seg_v, n);
#ifdef PRINT_SEG_AND_INT
      TIntegerSegment::coll_begin = coll_begin;
#endif
      TIntegerSegment s;
      uint4 i = 0;
      for (uint4 m = 0; m != n; ++m) {
        s.Init(sc[m], si, co, x_transf, y_transf);
        if constexpr (remove_zero_seg) {
          if (s.shift.is_zero()) continue;
        }
        seg_v[i] = s;
        points[first_point(i)] = seg_v[i].BegPoint();
        points[last_point(i)] = seg_v[i].EndPoint();
        ++i;
      }
      return i;
    }

    bool prepare_remap(std::vector<uint4>& indexes, TIntegerVect points[],
      std::vector<TIntegerSegment>& res_coll,
      std::vector<TIntegerVect>& res_pts) {

      auto N = indexes.size();
      uint4 remap_size = 0;
      int seg_numb = 0;
      remapped_SN = 0;
      for (auto pt : indexes) {
        seg_numb += is_last(pt) ? -1 : 1;
        assert(seg_numb >= 0);
        remap_size += seg_numb;
        remapped_SN += (seg_numb) ? 1 : 0;
      }
      if ((initial_SN == remapped_SN) && (initial_SN == nonzero_N))
        return true;//not remapped
      auto rrec = set_size(rrec_v, remapped_SN);
      auto remap = set_size(remap_v, remap_size);

      set_size(res_coll, remapped_SN);
      set_size(res_pts, 2 * remapped_SN);
      uint4 remaper_pos = 0, new_seg_num = 0;
      std::vector<uint4> stack;
      stack.reserve(nonzero_N);
      bool not_remapped = initial_SN == nonzero_N;
      for (size_t i = 0; i < N; ++i) {
        auto pt = indexes[i];
        if (is_last(pt)) {
          remove_by_val(stack, get_segm(pt));
        }
        else {
          auto s = get_segm(pt);
          for (auto p : stack) {
            registrator->register_pair(s, p);
          }
          stack.push_back(s);
        };
        auto size = stack.size();

        if (size != 0) {
          assert(i + 1 < N);

          auto beg = points[pt];
          auto end = points[indexes[i + 1]];

          if (beg != end) {
            std::copy(stack.begin(), stack.end(), remap + remaper_pos);
            rrec[new_seg_num] = { remaper_pos,remaper_pos + (uint4)size };
            remaper_pos += size;
            not_remapped &= size < 2;
            res_coll[new_seg_num] = { beg,end };
            res_pts[first_point(new_seg_num)] = beg;
            res_pts[last_point(new_seg_num)] = end;
            ++new_seg_num;

          }
        };
      }
      remapped_SN = new_seg_num;
      remap_size = remaper_pos;
      return not_remapped;
    }

    void register_pair(uint4 s1, uint4 s2, uint4 int_type) {
      auto orig_int_type = seg_v[s1].get_int_type_beg(seg_v[s2]);
      if (int_type == orig_int_type)
        registrator->register_pair(s1, s2);
    }

    void register_pair(uint4 s1, uint4 s2) {
#ifdef PRINT_SEG_AND_INT
      if (s1 > s2)
        std::swap(s1, s2);
      printf("remapped int pair %i %i\n", s1, s2);
#endif 
      auto rr1 = rrec_v[s1];
      auto rr2 = rrec_v[s2];
      auto int_type = remapped_segs[s1].get_int_type_beg(remapped_segs[s2]);
      //trivial remap (should be most of the time)
      if (rr1.end + rr2.end - rr1.beg - rr2.beg < 3) {
        register_pair(remap_v[rr1.beg], remap_v[rr2.beg], int_type);
        return;
      }
      //marking all nonunique segments
      constexpr const uint4 mark = 1 << (sizeof(uint4) * 8 - 1);
      for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
        for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
          if (remap_v[seg1] == remap_v[seg2]) {//same segment found, mark it 
            remap_v[seg1] |= mark;
            remap_v[seg2] |= mark;
            break;
          }

      //register only intersections with nonmarked segments
      for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
        if ((remap_v[seg1] & mark) == 0) {
          for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
            if ((remap_v[seg2] & mark) == 0)
              register_pair(remap_v[seg1], remap_v[seg2], int_type);
        }

      // erase marks, back to original values
      for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
        remap_v[seg1] &= ~mark;
      for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
        remap_v[seg2] &= ~mark;

    };

    template<bool remove_zero_seg = true>
    auto NoRemapInit(uint4 n, void* c, IntersectionRegistrator* r, int4 range, CTHIS& collection) {
      initial_SN = n;
      auto sc = reinterpret_cast<TLineSegment1*>(c);
      std::vector<TIntegerVect> points(2 * n);

      this->range = range;
      remapped_SN = nonzero_N = int_seg_from_real<remove_zero_seg>(n, sc, points);
      if constexpr (remove_zero_seg)
        initial_SN = nonzero_N;
      registrator = r;

      collection.segments = std::move(seg_v);
      collection.points = std::move(points);
      return;
    };

    auto TurnRemapOn(CTHIS& collection) {
      auto points = std::move(collection.points);
      seg_v = std::move(collection.segments);
      std::vector<uint4> indexes;
      auto n = initial_SN;
      indexes.reserve(2 * n);

      for (uint4 i = 0; i < 2 * n; ++i)
        if (seg_v[get_segm(i)].shift.is_non_zero())
          indexes.push_back(i);

      nonzero_N = indexes.size() / 2;

      auto comparator = [segments = seg_v.data(), pts = points.data()](uint4 pt1, uint4 pt2) {
        auto i1 = get_segm(pt1);
        auto i2 = get_segm(pt2);
        if (i1 == i2) // same segment
          return pt1 < pt2; //beging first
        auto& s1 = segments[i1];
        auto& s2 = segments[i2];
        auto S = s1.shift % s2.shift;
        if (S != 0)//segments non parallel
          return S < 0;
        //segments parallel
        auto shift = s1.shift + s2.shift;
        auto oo = pts[pt2] - pts[pt1];
        S = shift % oo;
        if (S != 0)// segments on different lines
          return  S < 0;
        //segments parallel and lies on one line
        assert((shift.x > 0) || ((shift.x == 0) && (shift.y > 0)));
        auto prod = oo * shift;
        if (prod != 0)// points not coinside
          return  0 < prod;
        // points coinside
        if (is_last(pt1) != is_last(pt2))
          return is_last(pt1) > is_last(pt2); //different segments end first
        return pt1 < pt2;
        };

      std::sort(indexes.begin(), indexes.end(), comparator);

      bool not_remapped = prepare_remap(indexes, points.data(), collection.segments, collection.points);
      if (not_remapped) {
        collection.segments = std::move(seg_v);
        collection.points = std::move(points);
        return false;
      }
      remapped_segs = collection.segments.data();
      remapped_ends = collection.points.data();

      return true;
    };

    TIntegerSegment& get_original_segment(uint4 s) {
      return seg_v[s];
    }


    auto get_N() {
      return remapped_SN;
    }

    IntersectionRegistrator* registrator = nullptr;

  private:
    uint4 initial_SN = 0;
    uint4 nonzero_N = 0;
    uint4 remapped_SN = 0;
    uint4 remap_size = 0;
    int4 range = 1;
    std::vector<RemappedRecord> rrec_v;
    std::vector<uint4> remap_v;
    std::vector<TIntegerSegment> seg_v;//if remapped original int segments are stored here
    TIntegerSegment* remapped_segs = nullptr;
    TIntegerVect* remapped_ends = nullptr;

  };
  //===========================================================================

public:


  static void reg_pair(CIntegerSegmentCollection* c, uint4 s1, uint4 s2) {
    c->remaper.registrator->register_pair(s1, s2);
  };
  static void reg_remapped_pair(CIntegerSegmentCollection* c, uint4 s1, uint4 s2) {
    c->remaper.register_pair(s1, s2);
  };
  //TPlaneVect
  uint4  GetSegmNumb() { return N; };
  uint4 get_right_pt_idx(uint4 right_pt) {
    return SegR[get_segm(right_pt)];
  }
  uint4 get_left_pt_idx(uint4 left_pt) {
    return SegL[get_segm(left_pt)];
  }

  bool is_right_pt_in_stripe(uint4 s) const {
    return right_bound_idx <= SegR[s];
  }
  bool is_left_pt_in_stripe(uint4 s) const {
    return left_bound_idx >= SegL[s];
  }
  bool is_right_pt_inside(uint4 s) const {
    return right_bound_idx < SegR[s];
  }
  bool is_left_pt_inside(uint4 s) const {
    return left_bound_idx > SegL[s];
  }
 
  bool is_right_pt_on_bound(uint4 s) const {
    return right_bound_idx == SegR[s];
  }
  bool is_left_pt_on_bound(uint4 s) const {
    return left_bound_idx >= SegL[s];
  }

  void SetCurStripe(uint4 left, uint4 right){//left,right: #points - bounds of the stripe 
    B = GetX(stripe_left = left);
    E = GetX(stripe_right = right);
    left_bound_idx = get_left_pt_idx(left);
    right_bound_idx = get_right_pt_idx(right);
  };

  void SetCurStripeRight(uint4 right){ 
    E = GetX(stripe_right=right); 
  };

  void SetCurStripeLeft(uint4 left) { 
    B = GetX(stripe_left = left);
  };

  /*void SetCurPoint(uint4 pt){
    cur_point = pts[pt];
  };
  };*/

  void SetCurPointAtBeg(uint4 s) {
    cur_point_idx = first_point(s);
    cur_point_seg = s;
    cur_point = collection[s].BegPoint();
  };

  void SetCurSeg(uint4 s){
    cur_seg_idx = s;
    cur_seg = collection[s];
    curB = cur_seg.bx();
    curE = cur_seg.ex();
    is_rstump=false;
  };

  void SetCurSegCutBE(uint4 s){
    SetCurSeg(s);
    stage = _Stages::stage_split;
    curB = MAX(B, curB);// possibly curB=B; !!!!!
    if (E < curE) {//means cur_seg.is_vertical()==false otherwise curE<=E
        curE = E;
        is_rstump = true;
    }
    active_end = cur_seg.EndPoint();
  };

  void SetCurSegCutBeg(uint4 s) {
    SetCurSeg(s);
    stage = _Stages::split_up;
    curB = MAX(B, curB);// possibly curB=B; !!!!!
    active_end = cur_seg.EndPoint();
  };

  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    stage = _Stages::stage_merge;
    curE = MIN(E, curE);// possibly curE=E; !!!!!
    active_end = cur_seg.BegPoint();
  };

  void SetCurSegAndPoint(uint4 s) {
    SetCurPointAtBeg(s);

    SetCurSeg(s);
    stage = _Stages::stage_bubble;
    active_end = pts[last_point(s)];
  };

  bool is_cur_strip_zero() const {
    return B == E;
  }

  bool is_cur_strip_non_zero() const {
    return !is_cur_strip_zero();
  }

  enum _Int_reg {
    // segments are not intersecting
    no_int = 0,
    // segments are intersecting, but no registration
    int_no_reg = 1,
    //Segments are intersecting and must be registered
    int_reg = int_no_reg + true //i.e. int_reg == 2 
  };

  int4 UpperThanRStump(int4 s_) const { // s_  is_upper than cur_seg right stump
    return RStump2CurSeg<false>(s_);
  };
  int4 UnderRStump(int4 s_) const { // s_  is_upper than cur_seg right stump
    return RStump2CurSeg<true>(s_);
  };

  template<bool check_way_up_s_is_under_cur_seg>// s_  is_under or is_upper than cur_seg right stump
  int4 RStump2CurSeg(int4 s_) const { 
    assert(stage == _Stages::stage_split);
    assert(!cur_seg.is_vertical());
    auto &s = collection[s_];
    if (is_cur_strip_non_zero()) {
      auto cmpB = s.YAtX_frac(B) <=> cur_seg.YAtX_frac(B);// 1 is >
      if (cmpB == std::strong_ordering::equal)// inresection on left bound
        return _Int_reg::int_reg;// and must be registered

      auto cmpE = s.YAtX_frac(E) <=> cur_seg.YAtX_frac(E);// cmpE==1 is s.YAtX_frac(E) > cur_seg.YAtX_frac(E)
      if (cmpE == std::strong_ordering::equal) {// inresection on right bound
        //  but no registration except right bound is defined last point of s
        return _Int_reg::int_no_reg + (stripe_right == last_point(s_));
      }

      if constexpr (check_way_up_s_is_under_cur_seg)
        return (cmpE == std::strong_ordering::less) * _Int_reg::int_reg; //s < cur_seg
      else //check_way_down_s_is_upper_cur_seg
        return (cmpE == std::strong_ordering::greater) * _Int_reg::int_reg; //s > cur_seg
    }
    else {// cur strip is zero
      if (stripe_right == last_point(s_)) {// the right bound defined by s last point
        auto cmp_last = cur_seg.upper(pts[stripe_right]);// 1 == false: last p of s is above cur_seg (s>cur_seg), remember last_point(s_)==stripe_right
        if (s.is_vertical()) {
          auto cmp_first = cur_seg.upper(pts[first_point(s_)]);// 1 == false:  beg. of s is above cur_seg
          assert( check_way_up_s_is_under_cur_seg ? cmp_first <= 0 : cmp_last >= 0);
          assert((cmp_first != 0) || (cmp_last != 0));
          if (cmp_first != cmp_last) {
            return _Int_reg::int_reg;
          }
        }
        if constexpr (check_way_up_s_is_under_cur_seg) 
          return (cmp_last <= 0) * _Int_reg::int_reg; //s <= cur_seg
        else //check_way_down_s_is_upper_cur_seg
          return (cmp_last >= 0) * _Int_reg::int_reg; //s >= cur_seg -> int_reg
      }
      else {// cur strip is zero and its right bound IS NOT defined by s last point
        if (s.is_vertical()) {
          auto cmp_last = cur_seg.upper(pts[last_point(s_)]);// -1 if end of s is bellow cur_seg
          if constexpr (check_way_up_s_is_under_cur_seg)
            return (cmp_last <= 0) * _Int_reg::int_reg; //s <= cur_seg
          else //check_way_down_s_is_upper_cur_seg
            return (cmp_last >= 0) * _Int_reg::int_reg; //s >= cur_seg -> int_reg
        }
        // int_no_reg
        auto cmpE = s.YAtX_frac(E) <=> cur_seg.YAtX_frac(E);
        assert(check_way_up_s_is_under_cur_seg ? cmpE >= 0 : cmpE <= 0);
        return cmpE == std::strong_ordering::equal; 
       /* if constexpr (check_way_up_s_is_under_cur_seg)
          return (cmpE <= 0); //s <= cur_seg
        else //check_way_down_s_is_upper_cur_seg
          return (cmpE >= 0); //s >= cur_seg -> int_no_reg*/
      }
    }
  };

  constexpr bool static left_bound = false;
  constexpr bool static right_bound = true;

  template <bool is_right>
  auto StumpPos(int4 s_) const {
    //returns cur_seg Y < s_ Y ? 1:-1;
    // probable error here!!!!!!
    auto X = is_right ? E : B;
    return cur_seg.YAtX_frac(X) <=> collection[s_].YAtX_frac(X);
  };

  bool FindCSIntWith(int4 s_) {
    std::strong_ordering lp, rp;
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    if (stage == _Stages::stage_bubble) {
      if (s2.is_vertical() || s1.is_vertical())
        return SSCurSegIntWith(s_);
      lp = collection[s_].under(cur_point);
      rp = collection[s_].under(active_end);
      // now cur_seg is not vertical, stripe is not degenerate (means s_ is no vertical too)
      if ((lp == 0) || (rp == 0)){
        register_pair(this,cur_seg_idx, s_);
        return true;
      }
      if ((lp < 0) ^ (rp < 0)) {
        register_pair(this,cur_seg_idx, s_);
        return true;
      };

    }

    bool one_s_pt_in_stripe = is_right_pt_in_stripe(cur_seg_idx) || is_right_pt_on_bound(s_);

    if (stage <= _Stages::split_up) {
      if (s1.is_vertical() || s2.is_vertical())
          return one_s_pt_in_stripe?SSCurSegIntWith(s_): SSCurSegIntWith<false>(s_);
      lp = StumpPos<left_bound>(s_);// if negative s_ is under cur left end 
      if (E == B)//degenerate stripe
        if (lp == 0) {
          if (one_s_pt_in_stripe)register_pair(this,cur_seg_idx, s_);
          return true;
        }
      rp = is_rstump ?
        StumpPos<right_bound>(s_) :// if negative s_ is under cur right end 
        collection[s_].under(active_end);// if negative s_ is under cur right end (active_end)
      if (lp == 0 && rp == 0){//segments parallel !!! impossible for nondegenerate strip
        if (one_s_pt_in_stripe)register_pair(this,cur_seg_idx, s_);
        return true;
      }      
      // now cur_seg is not vertical, stripe is not degenerate (means s_ is no vertical too)
      // and segments not parallel
      if (lp == 0)
      {
        register_pair(this,cur_seg_idx, s_);
        return true;
      } 
      if(rp == 0)  {
        if (one_s_pt_in_stripe)
          register_pair(this,cur_seg_idx, s_);
        return true;
      }
      if((lp < 0) ^ (rp < 0)) {
        register_pair(this,cur_seg_idx, s_);
        return true;
      };
    }

    if (stage == _Stages::stage_merge) {
      if (s1.is_vertical() || s2.is_vertical())
        return one_s_pt_in_stripe ? SSCurSegIntWith(s_) : SSCurSegIntWith<false>(s_);
      rp = StumpPos<right_bound>(s_);// if negative s_ is under cur right end 
      if (E == B)//degenerate stripe
        if (rp == 0) {
          if (one_s_pt_in_stripe)register_pair(this,cur_seg_idx, s_);
          return true;
        }
      lp =  collection[s_].point_pos(active_end);// if negative s_ is under cur left end (active_end)
      if (lp == 0 && rp == 0) {//segments parallel
        if (one_s_pt_in_stripe)register_pair(this,cur_seg_idx, s_);
        return true;
      }
      // now cur_seg is not vertical, stripe is not degenerate (means s_ is no vertical too)
      // and segments not parallel
      if (lp == 0)
      {
        register_pair(this,cur_seg_idx, s_);
        return true;
      }
      if (rp == 0) {
        if (one_s_pt_in_stripe)
          register_pair(this,cur_seg_idx, s_);
        return true;
      }
      if ((lp < 0) ^ (rp < 0)) {
        register_pair(this,cur_seg_idx, s_);
        return true;
      };

    }

  };


  template<bool is_right>
  bool XBelow(int4 s_1, int4 s_2, int4 X) const {
    auto& s1 = collection[s_1];
    auto& s2 = collection[s_2];//for right bound s2 always a stair
    if ((s1.is_vertical())&&(s2.is_vertical())) {
      //y1 = is_right_pt_on_bound(s_1) ? s1.ey() : s1.by();
      //y2 = is_right_pt_on_bound(s_2) ? s2.ey() : s2.by();
      //return y1 < y2;
      return s1.ey() < s2.ey();
    };
    {
      auto cmp = s1.YAtX_frac(X) <=> s2.YAtX_frac(X);
      if (cmp != std::strong_ordering::equal)
        return cmp < 0;
    };
    auto prod = s1.shift % s2.shift;
 
    if (prod != 0)
      return prod < 0;//intersections is on the right of X (prod>0 would be on the left)
    return s_1 < s_2;
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



  int4 FindCurSegIntDownWith(int4 s_){//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    if (IsIntersectsCurSegDown(s_)) {
      register_pair(this,cur_seg_idx, s_);
      return true;
    }
    return false;
  };


  bool FindCurSegIntUpWith(int4 s_){//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    if (IsIntersectsCurSegUp(s_)) {
      register_pair(this,cur_seg_idx, s_);
      return true;
    }
    return false;
  };


  bool IsIntersectsCurSegDown(int4 s_) //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    if (is_rstump)
      return UpperThanRStump(s_); // s_is_upper_cur_seg==true
    return UpperActiveEnd(s_); // s_ is upper when active end
  };

  auto FindCurSegIntDownWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  //Caller ensures that "last" points to the allocated memory address accessible for reading and writing.
    auto cs = cur_seg_idx;
    if (is_rstump) {
      while ((THIS_HAS_SENTINELS || (s_ != last)) && UpperThanRStump(*s_)) {
        register_pair(this,cs, *s_);
        --s_;
      }
      return s_;
    }

    while ((THIS_HAS_SENTINELS || (s_ != last)) && UpperActiveEnd(*s_)) {
      register_pair(this,cs, *s_);
      --s_;
    }
    return s_;
  };

  auto FindCurSegIntUpWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  //Caller ensures that "last" points to the allocated memory address accessible for reading and writing.
    auto cs = cur_seg_idx;
    /*if (is_rstump) {
      while ((THIS_HAS_SENTINELS || (s_ != last)) && RStump2CurSeg<true>(*s_)) {
        register_pair(this,cs, *s_);
        ++s_;
      }
      return s_;
    }*/
    while ((THIS_HAS_SENTINELS || (s_ != last)) && UnderActiveEnd(*s_)) {
      register_pair(this,cs, *s_);
      ++s_;
    }
    return s_;
  };


  bool IsIntersectsCurSegUp(int4 s_){ //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  
    //if (is_rstump)return RStump2CurSeg<true>(s_);
    return UnderActiveEnd(s_);
  };

  bool UnderCurPoint(int4 s_) const { //returns true if s_ is under current point 
    auto& s2 = collection[s_];
    if (s2.is_vertical()) {
      auto y = is_right_pt_in_stripe(s_) ? s2.org.y + s2.shift.y : s2.org.y;
      if (y != cur_point.y)
        return y < cur_point.y;
      goto fin;
    }
    {
    auto pp = s2.point_pos(cur_point);
    if (pp != 0)
      return pp < 0;
    }
    fin:
    auto& s1 = collection[cur_point_seg];
    auto prod = s1.shift % s2.shift;
    if (prod != 0)
      return prod > 0;
    return s_ < cur_point_seg;
  };

  bool UnderActiveEnd(int4 s_) const { //returns true if s is under active end
    auto pp = collection[s_].point_pos(active_end);
    if (pp != 0)
      return pp < 0;
    return is_right_pt_in_stripe(cur_seg_idx);
  }; 

  bool UpperActiveEnd(int4 s_) const { //returns true if s is upper current point
    auto pp = collection[s_].point_pos(active_end);
    if (pp != 0)
      return pp > 0;
    return is_right_pt_in_stripe(cur_seg_idx);
  };

  bool operator ()(uint4 pt1, uint4 pt2) {
    if (pts[pt1].x != pts[pt2].x)
      return pts[pt1].x < pts[pt2].x;
    if (pts[pt1].y != pts[pt2].y)
      return pts[pt1].y < pts[pt2].y;
    auto le1 = is_last (pt1), le2 = is_last(pt2);
    if (le1 != le2) 
      return le1 < le2;
    return (pt1 < pt2);
  }



  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N * 2;
    for (uint4 i = 0; i < NN; ++i)  epoints[i] = i;


/*
    static_assert(sizeof(TIntegerSegment) == 4 * sizeof(int4), "one segment should occupy 4 reals space");
    static_assert(offsetof(TIntegerSegment, org.x) == 0, "TLineSegment1.org.x should have 0 offset");
    static_assert(offsetof(TIntegerSegment, shift.x) == 2 * sizeof(int4), "TLineSegment1.shift.x should have 2 reals offset");
*/
    std::sort(epoints, epoints + NN, *this);

    ENDS = epoints;


    /*SegR = new uint4[N];
    for (uint4 i = 0, cur_end; i < NN; ++i) {
      cur_end = ENDS[i] /= 2;
      if (is_last(cur_end))
        SegR[get_segm(cur_end)] = i;
    };*/  

  };

  void set_seg2end_arr(uint4* _SegL, uint4* _SegR) {
    SegL = _SegL;
    SegR = _SegR;
  };

  void clone(CIntegerSegmentCollection &c, IntersectionRegistrator *r)
  {
      clone_of = &c;
      N = c.N;
      ENDS = c.ENDS;
      SegR = c.SegR;
      collection = c.collection;
      //!!!!!clone remaper!!!!!!!!!!!!
      SetRegistrator(r);
  };

  void unclone() { 
    if (clone_of == nullptr)return; 
    collection = nullptr; 
    clone_of = nullptr; 
    ENDS = nullptr;
    SegR = nullptr;
    //unclone remaper !!!!!!!!

  };

  void SortAt(uint4 pt, uint4 n, int4 *L)
  {
    SetCurStripeLeft(pt);
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

  CIntegerSegmentCollection(uint4 n, void* c, IntersectionRegistrator *r, int4 range)
  {
    assert(("integer collection registers only intersecting pairs:\n\
       intersection points don't have integer coords ",
      (IntersectionRegistrator::reg_type & _RegistrationType::point) == 0));
    remaper.NoRemapInit(n, c, r, range, *this);
    is_collection_remapped = false;
    register_pair = reg_pair;
    N = remaper.get_N();
    collection = segments.data();
    pts = points.data();

  }

  bool TurnRemapOn() {
    is_collection_remapped = remaper.TurnRemapOn(*this);
    register_pair = is_collection_remapped?reg_remapped_pair: reg_pair;
    N = remaper.get_N();
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
   // MY_FREE_ARR_MACRO(SegR);
  }


 
  void coll_to_SVG(chostream* SVG_stream) {
    if (!SVG_stream)return;
    int4 n = MIN(max_SVG_items, N);
    auto mmr = get_mmrect01(collection, n);
    *SVG_stream << "<svg height='100%' width='100%' viewBox='";
    *SVG_stream << mmr.ld.x << " " << mmr.ld.y << " "
      << mmr.rt.x - mmr.ld.x << " " << mmr.rt.y - mmr.ld.y << "'>\n";
    for (int4 i = 0; i < n; ++i)collection[i].write_SVG(i, SVG_stream);
  };

private:
  auto GetX(uint4 pt) const {
    return pts[pt].x;
  };

  TIntegerSegment* collection = nullptr;
  TIntegerVect* pts = nullptr;
  //uint4* SegL = nullptr; 
  uint4* SegL = nullptr, * SegR = nullptr, * ENDS = nullptr;

  TIntegerSegment cur_seg;
  TIntegerVect cur_point, active_end;
  //int8 cs_r_numerator = 0, cs_l_numerator = 0;
  std::vector<TIntegerSegment> segments;
  std::vector<TIntegerVect> points;
  CRemaper remaper;
  decltype(reg_pair) *register_pair = reg_pair;
  CIntegerSegmentCollection *clone_of = nullptr;
  uint4 N=0;
  int4 B, E, curB, curE, stripe_left, stripe_right, right_bound_idx, left_bound_idx;
  uint4 stage = _Stages::stage_split;

  uint4 cur_seg_idx = 0xFFFFFFFF, cur_point_idx= 0xFFFFFFFF;
  uint4 cur_point_seg = 0xFFFFFFFF, active_end_idx = 0xFFFFFFFF;
  //TPlaneVect x_transf, y_transf;
  bool is_rstump = false;
  bool is_collection_remapped = false;
//  bool rbelong_to_stripe=true;
};

// some string
