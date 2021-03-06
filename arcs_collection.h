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
  static const bool is_line_segments = false;

  static inline bool is_last(uint4 pt)
  {
    return pt & 1;
  };
  static inline uint4 get_segm(uint4 pt)
  {
    return pt / 6;
  };

  //TPlaneVect
  uint4  GetSegmNumb() { return N; };
  inline void SetCurStripe(uint4 left, uint4 right)
  {
    B = GetX(left);
    E = GetX(right);
  };
  inline void SetCurStripeRight(uint4 right) { E = GetX(right); };
  inline void SetCurStripeLeft(uint4 left) { B = GetX(left); };
  void SetCurPoint(uint4 pt)
  {
    if (is_last(pt))
      collection[get_segm(pt)].EndPoint(cur_point.x, cur_point.y);
    else
      collection[get_segm(pt)].BegPoint(cur_point.x, cur_point.y);
  };
  void SetCurPointAtBeg(uint4 s)
  {
    collection[s].BegPoint(cur_point.x, cur_point.y);
  };
  void SetCurPointAtEnd(uint4 s)
  {
    collection[s].EndPoint(cur_point.x, cur_point.y);
  };
  void SetCurSegCutBE(uint4 s)
  {
    SetCurSeg(s);
    cur_seg.x1 = max(cur_seg.x1, B);
    cur_seg.x2 = min(cur_seg.x2, E);
  };
  void SetCurSegCutBeg(uint4 s)
  {
    SetCurSeg(s);
    cur_seg.x1 = max(cur_seg.x1, B);
  };
  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    cur_seg.x2 = min(cur_seg.x2, E);
  };
  inline void SetCurSeg(uint4 s)
  {
    cur_seg_idx = s;
    cur_seg = collection[s];
  };

  bool LBelow(int4 s_1, int4 s_2) //retuns if s1 below s2 at current vertical line
  {
    return collection[s_1].under(collection[s_2].PointAtX(B));
  };

  bool RBelow(int4 s_1, int4 s_2) //retuns if s1 below s2 at current vertical line
  {
    return collection[s_1].under(collection[s_2].PointAtX(E));
  };

#define REGISTER_INTERSECTION(s1,s2,p)\
  if (_RegistrationType::point&IntersectionRegistrator::reg_type)\
  {\
    registrator->begin_registration(1)->\
    register_segments(s1, s2-collection)->\
    register_point(p)->\
    end_registration();\
  }\
  else\
  {\
    registrator->begin_registration(1)->\
    register_segments(s1, s2-collection)->\
    end_registration();\
  }




  template <bool not_register>
  bool IntPointsInStripe(REAL x1, REAL x2, TArcSegment *s2)
  {
    auto s1 = &cur_seg;
    TPlaneVect oo = s2->org - s1->org;
    REAL org_dist2 = oo.get_norm();
    REAL delta = 0.5*(s1->r2 - s2->r2 + org_dist2);
    REAL t = delta / org_dist2;
    REAL l2 = s1->r2 - delta*t;
    if (l2<0)return false;
    TPlaneVect m = s1->org + t*oo;
    if (l2 == 0)
    {
      if ((m.x >= x1) && (m.x <= x2) && (s1->IsTheSamePart(m)) && (s2->IsTheSamePart(m)))
      {
        if (not_register) return true;
        REGISTER_INTERSECTION(cur_seg_idx, s2, m);
        return true;
      }
      return false;
    }
    t = sqrt(l2 / org_dist2);
    oo = TPlaneVect(oo.y*t, -oo.x*t);
    int4 npoints = 0;
    TPlaneVect res = m + oo;
    bool ret = false;
    if ((res.x >= x1) && (res.x <= x2) && (s1->IsTheSamePart(res)) && (s2->IsTheSamePart(res)))
    {
      if (not_register) return true;
      ret = true;
      REGISTER_INTERSECTION(cur_seg_idx, s2, res);
    }
    res = m - oo;
    if ((res.x >= x1) && (res.x <= x2) && (s1->IsTheSamePart(res)) && (s2->IsTheSamePart(res)))
    {
      if (not_register) return true;
      ret = true;
      REGISTER_INTERSECTION(cur_seg_idx, s2, res);
    }
    return ret;
  };



  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = max(cur_seg.x1, s->x1);
    auto x2 = min(cur_seg.x2, s->x2);
    if (x1 >= x2)return false;
    return IntPointsInStripe<false>(x1, x2, s);

  };

  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = max(cur_seg.x1, s->x1);
    auto x2 = min(cur_seg.x2, s->x2);
    return IntPointsInStripe<false>(x1, x2, s);
  };

  bool FindCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    return IntPointsInStripe<false>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  bool IsIntersectsCurSeg(int4 s_)//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return IntPointsInStripe<true>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  bool UnderCurPoint(int4 s_) { return collection[s_].under(cur_point); };//returns true if s is under current point 
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
  void clone(CArcSegmentCollection * c, IntersectionRegistrator *r)
  {
    clone_of = c;
    Init(c->N, c->collection,r);
  };


  void unclone() { if (clone_of == nullptr)return; collection = nullptr; clone_of = nullptr; };
  void SortAt(uint4 pt, uint4 n, int4 *L)
  {
    SetCurStripeLeft(pt);
    std::sort(L, L + n, [this](int4 s1, int4 s2) {return LBelow(s1, s2); });
  };

  void SetRegistrator(IntersectionRegistrator *r)
  {
    registrator = r;
    r->Alloc(GetSegmNumb());
  };
  //IntersectionRegistrator *GetRegistrator() { return registrator; };


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

  CArcSegmentCollection(uint4 n, void* c, IntersectionRegistrator *r)
  {
    Init(n, c, r);
  }
  CArcSegmentCollection() {};
private:
  inline auto GetX(uint4 pt) { return ends[pt]; };
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


