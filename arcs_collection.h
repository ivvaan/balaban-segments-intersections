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

template<class IntersectionRegistrator>
class CArcSegmentCollection
{
public:
  using this_T = CArcSegmentCollection;
  static constexpr _Coll_flag_state get_coll_flag(_Coll_flags flag) {
    if (flag == _Coll_flags::line_segments)
      return _Coll_flag_state::state_false;
   
    return _Coll_flag_state::state_unimplemented;
  }

  //static constexpr bool is_line_segments = false;

  static bool is_last(uint4 pt)
  {
    return pt & 1;
  };
  static uint4 get_segm(uint4 pt)
  {
    return pt / 6;
  };

  //TPlaneVect
  uint4  GetSegmNumb() const { return N; };
  void SetCurStripe(uint4 left, uint4 right)
  {
    B = GetX(left);
    E = GetX(right);
  };
  void SetCurStripeRight(uint4 right) { E = GetX(right); };
  void SetCurStripeLeft(uint4 left) { B = GetX(left); };

  void SetCurPointAtBeg(uint4 s) {
    collection[s].BegPoint(cur_point.x, cur_point.y);
  };

  void SetCurSegAndPoint(uint4 s)  {
    collection[s].BegPoint(cur_point.x, cur_point.y);
    SetCurSeg(s);
  };

  void SetCurSeg(uint4 s)
  {
    cur_seg_idx = s;
    cur_seg = collection[s];
  };

  void SetCurSegCutBE(uint4 s)
  {
    SetCurSeg(s);
    cur_seg.x1 = MAX(cur_seg.x1, B);
    cur_seg.x2 = MIN(cur_seg.x2, E);
  };
  void SetCurSegCutBeg(uint4 s)
  {
    SetCurSeg(s);
    cur_seg.x1 = MAX(cur_seg.x1, B);
  };
  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    cur_seg.x2 = MIN(cur_seg.x2, E);
  };

  bool LBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    return collection[s_1].under(collection[s_2].PointAtX(B));
  };

  bool RBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    return collection[s_1].under(collection[s_2].PointAtX(E));
  };

#define REGISTER_INTERSECTION(s1,s2,p)\
  if constexpr((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0)\
    registrator->register_pair_and_point(s1, s2-collection, p);\
  else\
    registrator->register_pair(s1, s2-collection);




  template <bool not_register>
  bool IntPointsInStripe(REAL x1, REAL x2, TArcSegment* s2) {
    auto s1 = &cur_seg;
    auto oo = s2->org - s1->org;
    REAL org_dist2 = oo.get_norm();
    REAL delta = 0.5 * (s1->r2 - s2->r2 + org_dist2);
    REAL l2 = s1->r2 * org_dist2 - delta * delta;
    if (l2 < 0)return false;
    oo /= org_dist2;
    auto m = s1->org + delta * oo;
    oo.rotate_prod(sqrt(l2));
    auto res = m + oo;
    bool ret = false;
    if ((res.x >= x1) && (res.x <= x2) && (s1->IsTheSamePart(res)) && (s2->IsTheSamePart(res)))
    {
      if constexpr (not_register) return true;
      ret = true;
      REGISTER_INTERSECTION(cur_seg_idx, s2, res);
    }
    if (l2 == 0)return ret;
    res = m - oo;
    if ((res.x >= x1) && (res.x <= x2) && (s1->IsTheSamePart(res)) && (s2->IsTheSamePart(res)))
    {
      if constexpr (not_register) return true;
      ret = true;
      REGISTER_INTERSECTION(cur_seg_idx, s2, res);
    }
    return ret;
  };



  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = MAX(cur_seg.x1, s->x1);
    auto x2 = MIN(cur_seg.x2, s->x2);
    if (x1 >= x2)return false;
    return IntPointsInStripe<false>(x1, x2, s);

  };

  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = MAX(cur_seg.x1, s->x1);
    auto x2 = MIN(cur_seg.x2, s->x2);
    return IntPointsInStripe<false>(x1, x2, s);
  };

  bool FindCurSegIntDownWith(int4 s_) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return IntPointsInStripe<false>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  auto FindCurSegIntDownWith(int4 *s_,int4 *last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    while((THIS_HAS_SENTINELS || (s_ != last)) && IntPointsInStripe<false>(cur_seg.x1, cur_seg.x2, collection + *s_))
      --s_;
    return s_;
  };

  auto FindCurSegIntUpWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    while ((THIS_HAS_SENTINELS || (s_ != last)) && IntPointsInStripe<false>(cur_seg.x1, cur_seg.x2, collection + *s_))
      ++s_;
    return s_;
  };

  bool FindCurSegIntUpWith(int4 s_) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return IntPointsInStripe<false>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  bool IsIntersectsCurSegDown(int4 s_){//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
    return IntPointsInStripe<true>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  bool IsIntersectsCurSegUp(int4 s_){//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
    return IntPointsInStripe<true>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  bool UnderCurPoint(int4 s_) const { return collection[s_].under(cur_point); };//returns true if s is under current point 
  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N << 1;
    for (uint4 i = 0; i < NN; ++i)     epoints[i] = i * 3;
    std::sort(epoints, epoints + NN,
      [x = ends](uint4 pt1, uint4 pt2) {
      return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
    }
    );
  };

  void SetRegistrator(IntersectionRegistrator *r){
    registrator = r;
  };

  IntersectionRegistrator *GetRegistrator() { return registrator; };

  void set_seg2end_arr(uint4* SegL, uint4* SegR) {};

  void Init(uint4 n, void* c, IntersectionRegistrator *r)
  {
    N = n;
    collection = reinterpret_cast<TArcSegment*>(c);

    static_assert(sizeof(TArcSegment) == 6 * sizeof(REAL), "one segment should occupy 6 reals space");
    static_assert(offsetof(TArcSegment, x1) == 0, "TArcSegment.x1 should have 0 offset");
    static_assert(offsetof(TArcSegment, x2) == 3 * sizeof(REAL), "TArcSegment.x2 should have 3 reals offset");
    ends = reinterpret_cast<REAL*>(collection);
    SetRegistrator(r);
  };

  void SortAt(uint4 pt, uint4 n, int4* L)
  {
    SetCurStripeLeft(pt);
    std::sort(L, L + n, [this](int4 s1, int4 s2) {return LBelow(s1, s2); });
  };

  void clone(CArcSegmentCollection& c, IntersectionRegistrator* r)
  {
    clone_of = &c;
    Init(c.N, c.collection, r);
  };

  void unclone() {
    if (clone_of == nullptr)return; 
    collection = nullptr; 
    clone_of = nullptr;
  };

  CArcSegmentCollection() {};

  CArcSegmentCollection(uint4 n, void* c, IntersectionRegistrator* r)
  {
    Init(n, c, r);
  }

  CArcSegmentCollection(CArcSegmentCollection& coll, IntersectionRegistrator* r)
  {
    clone(coll,r);
  }

  ~CArcSegmentCollection()
  {
    unclone();
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
  auto GetX(uint4 pt) const { return ends[pt]; };
  //{ return is_last(pt) ? collection[get_segm(pt)].x2 :collection[get_segm(pt)].x1; };


  IntersectionRegistrator *registrator = nullptr;
  CArcSegmentCollection *clone_of = nullptr;
  uint4 N=0;
  REAL B, E, X;
  TPlaneVect cur_point;

  uint4 cur_seg_idx = 0xFFFFFFFF;
  TArcSegment cur_seg;
  TArcSegment *collection= nullptr;
  REAL *ends= nullptr;
};


