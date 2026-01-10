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
#include <cassert>

namespace {
  template<class T>
  T& unmove(T&& t) { return t; };

  enum _Stages {
    split=0, bubble, merge
  };
}



template<class IntersectionRegistrator>
class CIntegerSegmentCollection
{
public:
  static constexpr bool is_line_segments = false;//they can have more than one point of intersection with each other!!

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
    return (s << 1)+1;
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

   void SetCurStripe(uint4 left, uint4 right)
  {
    B = GetX(left);
    E = GetX(stripe_right = right);
    left_bound_idx = get_left_pt_idx(left);
    right_bound_idx = get_right_pt_idx(right);
   };
   void SetCurStripeRight(uint4 right) { E = GetX(stripe_right=right); };
   void SetCurStripeLeft(uint4 left) { B = GetX(left); };
  void SetCurPoint(uint4 pt)
  {
    cur_point = is_last(pt) ?
      collection[get_segm(pt)].EndPoint() :
      collection[get_segm(pt)].BegPoint();
  };
  void SetCurSegAndPoint(uint4 s)
  {
    cur_point_idx = first_point(s);
    cur_point_seg = s;
    cur_point=collection[s].BegPoint();
  };
  void SetCurPointAtBeg(uint4 s) {
    cur_point_idx = first_point(s);
    cur_point_seg = s;
    cur_point = collection[s].BegPoint();
  };
  void SetCurSeg(uint4 s)
  {
    cur_seg_idx = s;
    cur_seg = collection[s];
    curB = cur_seg.org.x;
    curE = cur_seg.org.x + cur_seg.shift.x;
    is_rstump=false;
  };

  void SetCurSegCutBE(uint4 s)
  {
    SetCurSeg(s);
    stage = _Stages::split;
    curB = MAX(B, curB);
    //cs_l_numerator = cur_seg.YAtX_Numerator(E);
    if (E < curE) {
        //cs_r_numerator = cur_seg.YAtX_Numerator(E);
        curE = E;
        is_rstump = true;
    }
    active_end = cur_seg.EndPoint();
  };

  /*void SetCurSegCutBeg(uint4 s)
  {
    SetCurSeg(s);
    curB = MAX(B, curB);
    active_end = cur_seg.EndPoint();
  };*/

  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    stage = _Stages::merge;
    curE = MIN(E, curE);
    active_end = cur_seg.BegPoint();
  };

  /*void SetCurSegAE(uint4 s)
  {
    SetCurSeg(s);
    stage = _Stages::bubble;
    active_end = cur_seg.EndPoint();
  };*/

#define FLAG_BASE_COMP_RET(v1,v2,f)  if (v1 != v2) return (v1 < v2) && !f || f && (v1 > v2);


  template<bool check_s_under_cur_seg>// s_ is_under or is_upper of cur_seg
  bool RStump2CurSeg(int4 s_)  {
    auto& s = collection[s_];
    auto prod = cur_seg.shift % s.shift;
    if (prod  == 0)return false;
    unsigned __int64 prod_lo1, prod_hi1;//cur_seg Y
    unsigned __int64 prod_lo2, prod_hi2;//s_ Y
    prod_lo1 = _umul128(s.shift.x, cur_seg.YAtX_Numerator(E), &prod_hi1);
    prod_lo2 = _umul128(cur_seg.shift.x, s.YAtX_Numerator(E), &prod_hi2);
    FLAG_BASE_COMP_RET(prod_hi1, prod_hi2, check_s_under_cur_seg);
    FLAG_BASE_COMP_RET(prod_lo1, prod_lo2, check_s_under_cur_seg);
 
    return is_right_pt_in_stripe(s_);
  };

  constexpr bool static left_bound = false;
  constexpr bool static right_bound = true;

  template <bool is_right>
  int8 StumpPos(int4 s_) const {
    auto& s = collection[s_];
    unsigned __int64 prod_lo1, prod_hi1;//cur_seg Y
    unsigned __int64 prod_lo2, prod_hi2;//s_ Y
    prod_lo1 = _umul128(s.shift.x, cur_seg.YAtX_Numerator(is_right ? E : B), &prod_hi1);
    prod_lo2 = _umul128(cur_seg.shift.x, s.YAtX_Numerator(is_right ? E : B), &prod_hi2);
    if (prod_hi1 != prod_hi2)
      return (prod_hi1 < prod_hi2) ? 1 : -1;
    if (prod_lo1 != prod_lo2)
      return (prod_lo1 < prod_lo2) ? 1 : -1;
    return 0;
  };

  bool FindCSIntWith(int4 s_) {
    int8 lp, rp;
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    if (stage == _Stages::bubble) {
      if (s2.shift.x == 0 || s1.shift.x == 0)
        return SSCurSegIntWith(s_);
      lp = collection[s_].point_pos(cur_point);
      rp = collection[s_].point_pos(active_end);
      // now cur_seg is not vertical, stripe is not degenerate (means s_ is no vertical too)
      if ((lp == 0) || (rp == 0)){
        registrator->register_pair(cur_seg_idx, s_);
        return true;
      }
      if ((lp < 0) ^ (rp < 0)) {
        registrator->register_pair(cur_seg_idx, s_);
        return true;
      };

    }

    bool one_s_pt_in_stripe = is_right_pt_in_stripe(cur_seg_idx) || is_right_pt_on_bound(s_);

    if (stage == _Stages::split) {
      if (s1.shift.x == 0 || s2.shift.x==0)
          return one_s_pt_in_stripe?SSCurSegIntWith(s_): SSCurSegIntWith<false>(s_);
      lp = StumpPos<left_bound>(s_);// if negative s_ is under cur left end 
      if (E == B)//degenerate stripe
        if (lp == 0) {
          if (one_s_pt_in_stripe)registrator->register_pair(cur_seg_idx, s_);
          return true;
        }
      rp = is_rstump ?
        StumpPos<right_bound>(s_) :// if negative s_ is under cur right end 
        collection[s_].point_pos(active_end);// if negative s_ is under cur right end (active_end)
      if (lp == 0 && rp == 0){//segments parallel
        if (one_s_pt_in_stripe)registrator->register_pair(cur_seg_idx, s_);
        return true;
      }      
      // now cur_seg is not vertical, stripe is not degenerate (means s_ is no vertical too)
      // and segments not parallel
      if (lp == 0)
      {
        registrator->register_pair(cur_seg_idx, s_);
        return true;
      } 
      if(rp == 0)  {
        if (one_s_pt_in_stripe)
          registrator->register_pair(cur_seg_idx, s_);
        return true;
      }
      if((lp < 0) ^ (rp < 0)) {
        registrator->register_pair(cur_seg_idx, s_);
        return true;
      };
    }

    if (stage == _Stages::merge) {
      if (s1.shift.x == 0 || s2.shift.x == 0)
        return one_s_pt_in_stripe ? SSCurSegIntWith(s_) : SSCurSegIntWith<false>(s_);
      rp = StumpPos<right_bound>(s_);// if negative s_ is under cur right end 
      if (E == B)//degenerate stripe
        if (rp == 0) {
          if (one_s_pt_in_stripe)registrator->register_pair(cur_seg_idx, s_);
          return true;
        }
      lp =  collection[s_].point_pos(active_end);// if negative s_ is under cur left end (active_end)
      if (lp == 0 && rp == 0) {//segments parallel
        if (one_s_pt_in_stripe)registrator->register_pair(cur_seg_idx, s_);
        return true;
      }
      // now cur_seg is not vertical, stripe is not degenerate (means s_ is no vertical too)
      // and segments not parallel
      if (lp == 0)
      {
        registrator->register_pair(cur_seg_idx, s_);
        return true;
      }
      if (rp == 0) {
        if (one_s_pt_in_stripe)
          registrator->register_pair(cur_seg_idx, s_);
        return true;
      }
      if ((lp < 0) ^ (rp < 0)) {
        registrator->register_pair(cur_seg_idx, s_);
        return true;
      };

    }

  };


  template<bool is_right>
  bool XBelow(int4 s_1, int4 s_2, int4 X) const {
    auto& s1 = collection[s_1];
    auto& s2 = collection[s_2];//for right bound s2 always a stair
    if ((s1.shift.x == 0)&&(s2.shift.x==0)) {
      //y1 = is_right_pt_on_bound(s_1) ? s1.org.y + s1.shift.y : s1.org.y;
      //y2 = is_right_pt_on_bound(s_2) ? s2.org.y + s2.shift.y : s2.org.y;
      //return y1 < y2;
      return s1.org.y + s1.shift.y < s2.org.y + s2.shift.y;
    };
    {
      unsigned __int64 prod_lo1, prod_hi1, prod_lo2, prod_hi2;
      prod_lo1 = _umul128(s2.shift.x, s1.YAtX_Numerator(X), &prod_hi1);
      prod_lo2 = _umul128(s1.shift.x, s2.YAtX_Numerator(X), &prod_hi2);
      if (prod_hi1 != prod_hi2)
        return (prod_hi1 < prod_hi2);
      if (prod_lo1 != prod_lo2)
        return (prod_lo1 < prod_lo2);
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
    if ((s1.org.x+s1.shift.x < s2.org.x) || (s2.org.x + s2.shift.x < s1.org.x))
      return false;//segments are not on vertical line and can't intersect
    return SSCurSegIntWith(s_);
  };

  template <bool reg=true>
  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (both are on some vertical line) and register them
  {
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    auto delt = s2.org - s1.org;
    //at this moment, we can be sure that segments' x-projections intersect each other (1)
    // let l1 is the line through s1 and l2 is the line through s2
    auto prod = s1.shift % s2.shift, beg = s1.shift % delt;
    decltype(prod) end;
    if (prod == 0) { // segments parallel
      if (beg == 0) {// they are on the same line
        if ((s1.shift.x != 0) || // if they not vertical
          (s1.shift.y >= delt.y) && (s2.shift.y >= -delt.y)) {// or y- intersects
          if constexpr (reg)registrator->register_pair(cur_seg_idx, s_);
          return true;//segments are not on vertical line and can't intersect
          }
      }
      return false;
    }
 
    do {
      if (beg == 0) {//begin of s2 lies on l1, and taking (1) to account
        if ((delt.x > 0)||((delt.x==0)&&(delt.y>=0)))//if begin of s2 is on the right of begin of s1
        {
          if constexpr(reg)registrator->register_pair(cur_seg_idx, s_);
          return true;
        };
        break;
      }
      end = beg + prod;
      if (end == 0) {//end of s2 lies on l1, and taking (1) to account
        if (delt.x + s2.shift.x <= s1.shift.x)//if end of s2 is on the left of s1 end
        {
          if constexpr (reg)registrator->register_pair(cur_seg_idx, s_);
          return true;
        };
        break;
      }
      if ((beg > 0) ^ (end < 0))////begin and end of s2 are on the same side of l1 
        return false;
    } while (false);
    //at this moment, we can be sure that s2 intesects l1 (2)
    beg = delt % s2.shift;
    end = beg - prod;
    /*if ((beg == 0) || (end == 0)) {
      registrator->register_pair(cur_seg_idx, s_);
      return true;
    };
    if ((beg > 0) ^ (end > 0)) {//begin and end of s1 are on diffenent sides l2
      registrator->register_pair(cur_seg_idx, s_);
      return true;
    };*/

    if ((beg == 0) ||//begin of s1 lies on l2, taking (2) to account
      (end == 0) ||//end of s1 lies on l2, taking (2) to account
      ((beg > 0) ^ (end > 0))) {//begin and end of s1 are on diffenent sides l2 
      if constexpr (reg)registrator->register_pair(cur_seg_idx, s_);
      return true;
    };
    return false;
  };

  int4 FindCurSegIntDownWith(int4 s_){//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    if (IsIntersectsCurSegDown(s_)) {
      registrator->register_pair(cur_seg_idx, s_);
      return true;
    }
    return false;
  };


  bool FindCurSegIntUpWith(int4 s_){//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    if (IsIntersectsCurSegUp(s_)) {
      registrator->register_pair(cur_seg_idx, s_);
      return true;
    }
    return false;
  };


  bool IsIntersectsCurSegDown(int4 s_) //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    if (is_rstump)
      return RStump2CurSeg<false>(s_);
    return UpperActiveEnd(s_);
  };

  auto FindCurSegIntDownWith(int4* s_,int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  //Caller ensures that "last" points to the allocated memory address accessible for reading and writing.
    auto r = registrator;
    auto cs = cur_seg_idx;
    if (is_rstump) {
      while ((THIS_HAS_SENTINELS || (s_ != last)) && RStump2CurSeg<false>(*s_)) {
        r->register_pair(cs, *s_);
        --s_;
      }
      return s_;
    }

    while ((THIS_HAS_SENTINELS || (s_ != last)) && UpperActiveEnd(*s_)) {
      r->register_pair(cs, *s_);
      --s_;
    }
    return s_;
  };

  auto FindCurSegIntUpWith(int4* s_,int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  //Caller ensures that "last" points to the allocated memory address accessible for reading and writing.
    auto r = registrator;
    auto cs = cur_seg_idx;
    if (is_rstump) {
      while ((THIS_HAS_SENTINELS || (s_ != last)) && RStump2CurSeg<true>(*s_)) {
        r->register_pair(cs, *s_);
        ++s_;
      }
      return s_;
    }
    while ((THIS_HAS_SENTINELS || (s_ != last)) && UnderActiveEnd(*s_)) {
      r->register_pair(cs, *s_);
      ++s_;
    }
    return s_;
  };


  bool IsIntersectsCurSegUp(int4 s_){ //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  
    if (is_rstump)
      return RStump2CurSeg<true>(s_);
    return UnderActiveEnd(s_);
  };

  bool UnderCurPoint(int4 s_) const { //returns true if s is under current point 
    auto& s2 = collection[s_];
    if (s2.shift.x == 0) {
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

  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N * 2;
    for (uint4 i = 0; i < NN; ++i)  epoints[i] = i*2;


    for (uint4 i = 0; i < N; ++i) {
      collection[i].shift.x += collection[i].org.x;
    };

    static_assert(sizeof(TIntegerSegment) == 4 * sizeof(int4), "one segment should occupy 4 reals space");
    static_assert(offsetof(TIntegerSegment, org.x) == 0, "TLineSegment1.org.x should have 0 offset");
    static_assert(offsetof(TIntegerSegment, shift.x) == 2 * sizeof(int4), "TLineSegment1.shift.x should have 2 reals offset");

    std::sort(epoints, epoints + NN,
      [x = reinterpret_cast<int4*>(collection)](uint4 pt1, uint4 pt2) {
        if (x[pt1] != x[pt2])
          return x[pt1] < x[pt2];
        auto le1 = (pt1 & 2), le2 = (pt2 & 2);
        if (le1 != le2) //last end (last ending on 0b10, first on 0b00) alway greater than the first
          return le1 < le2;
        return (pt1 < pt2);
      }
    );
    ENDS = epoints;

    for (uint4 i = 0; i < NN; ++i)  epoints[i] /= 2;

    /*SegR = new uint4[N];
    for (uint4 i = 0, cur_end; i < NN; ++i) {
      cur_end = ENDS[i] /= 2;
      if (is_last(cur_end))
        SegR[get_segm(cur_end)] = i;
    };*/  

    for (uint4 i = 0; i < N; ++i) {
      collection[i].shift.x -= collection[i].org.x;
    };
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
      SetRegistrator(r);
  };

  void unclone() { 
    if (clone_of == nullptr)return; 
    collection = nullptr; 
    clone_of = nullptr; 
    ENDS = nullptr;
    SegR = nullptr;
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
    registrator = r;
  };

  IntersectionRegistrator* GetRegistrator() { return registrator; };


  CIntegerSegmentCollection(uint4 n, void* c, IntersectionRegistrator *r, int4 range)
  {
    N = n;
    auto coll = reinterpret_cast<TLineSegment1*>(c);
    {
      auto ck = std::make_unique<TIntegerSegment[]>(N);
      collection_keeper.swap(ck);
    }
    //collection_keeper.swap(unmove(std::make_unique<TIntegerSegment[]>(N)));
    collection = collection_keeper.get();
    auto mm_rect = get_minmax(N,coll);
    auto x_transf = transf1D{.shift = -mm_rect.ld.x,
      .scale = range / (mm_rect.rt.x - mm_rect.ld.x)};
    auto y_transf = transf1D{.shift = -mm_rect.ld.y,
      .scale = range / (mm_rect.rt.y - mm_rect.ld.y)};
    for (uint4 i = 0; i != n; ++i)
      collection[i].Init(coll[i], x_transf, y_transf);
    SetRegistrator(r);

  }

  CIntegerSegmentCollection() {};
 
  CIntegerSegmentCollection(CIntegerSegmentCollection& coll, IntersectionRegistrator* r)
  {
    clone(coll, r);
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
    return collection[get_segm(pt)].PointX(is_last(pt));
  };

  TIntegerSegment* collection = nullptr;
  //uint4* SegL = nullptr; 
  uint4* SegL = nullptr, * SegR = nullptr, * ENDS = nullptr;

  TIntegerSegment cur_seg;
  TIntegerVect cur_point, active_end;
  //int8 cs_r_numerator = 0, cs_l_numerator = 0;
  std::unique_ptr<TIntegerSegment[]> collection_keeper;

  IntersectionRegistrator *registrator = nullptr;
  CIntegerSegmentCollection *clone_of = nullptr;
  uint4 N=0;
  int4 B, E, curB, curE, stripe_right, right_bound_idx, left_bound_idx;
  uint4 stage = _Stages::split;

  uint4 cur_seg_idx = 0xFFFFFFFF, cur_point_idx= 0xFFFFFFFF;
  uint4 cur_point_seg = 0xFFFFFFFF, active_end_idx = 0xFFFFFFFF;
  //TPlaneVect x_transf, y_transf;
  bool is_rstump = false;
//  bool rbelong_to_stripe=true;
};

