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
class CLine2SegmentCollection
{
public:
  static constexpr bool is_line_segments = true;

  static inline bool is_last(uint4 pt)
  {
    return pt & 2;
  };
  static inline uint4 get_segm(uint4 pt)
  {
    return pt >> 2;
  };

  //TPlaneVect
  uint4  GetSegmNumb() const { return N; };
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
  inline void SetCurSeg(uint4 s)
  {
    cur_seg_idx = s;
    cur_seg = collection[s];
  };
  bool LBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    auto s1 = collection + s_1;
    auto s2 = collection + s_2;
    return ((B * (s2->a - s1->a) + s2->b - s1->b) > 0);
  };
  bool RBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    auto s1 = collection + s_1;
    auto s2 = collection + s_2;
    return ((E * (s2->a - s1->a) + s2->b - s1->b) > 0);
  };


  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = MAX(cur_seg.x1, s->x1);
    auto x2 = MIN(cur_seg.x2, s->x2);
    if (x1 >= x2)return false;
    auto da = cur_seg.a - s->a;
    if (da == 0)return false;
    if (_RegistrationType::point&IntersectionRegistrator::reg_type)
    {
      TPlaneVect p;
      p.x = (s->b - cur_seg.b) / da;
      if ((p.x >= x1) && (p.x <= x2))
      {
        p.y = p.x*cur_seg.a + cur_seg.b;
        
        registrator->register_pair_and_point(cur_seg_idx, s_,p);
        
        return true;
      }
    }
    else
    {
      auto x = (s->b - cur_seg.b) / da;
      if ((x >= x1) && (x <= x2))
      {
        
        registrator->register_pair(cur_seg_idx, s_);
        
        return true;
      }
    }
    return false;

  };

  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto da = cur_seg.a - s->a;
    if (da == 0)return false;
    auto x1 = MAX(cur_seg.x1, s->x1);
    auto x2 = MIN(cur_seg.x2, s->x2);
    if (_RegistrationType::point&IntersectionRegistrator::reg_type)
    {
      TPlaneVect p;
      p.x = (s->b - cur_seg.b) / da;
      if ((p.x >= x1) && (p.x <= x2))
      {
        p.y = p.x*cur_seg.a + cur_seg.b;
        
        registrator->register_pair_and_point(cur_seg_idx, s_,p);
        
        return true;
      }
    }
    else
    {
      auto x = (s->b - cur_seg.b) / da;
      if ((x >= x1) && (x <= x2))
      {
        
        registrator->register_pair(cur_seg_idx, s_);
        
        return true;
      }
    }
    return false;

  };

  bool FindCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto da = cur_seg.a - s->a;
    if (da == 0)return false;
    if (_RegistrationType::point&IntersectionRegistrator::reg_type)
    {
      TPlaneVect p;
      p.x = (s->b - cur_seg.b) / da;
      if ((p.x >= cur_seg.x1) && (p.x <= cur_seg.x2))
      {
        p.y = p.x*cur_seg.a + cur_seg.b;
        
        registrator->register_pair_and_point(cur_seg_idx, s_,p);
        
        return true;
      }
    }
    else
    {
      auto x = (s->b - cur_seg.b) / da;
      if ((x >= cur_seg.x1) && (x <= cur_seg.x2))
      {
        
        registrator->register_pair(cur_seg_idx, s_);
        
        return true;
      }
    }
    return false;

  };

  bool IsIntersectsCurSeg(int4 s_)//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    auto s = collection + s_;
    auto da = cur_seg.a - s->a;
    if (da == 0)return false;
    auto x = (s->b - cur_seg.b) / da;
    return ((x >= cur_seg.x1) && (x <= cur_seg.x2));
  };

  bool UnderCurPoint(int4 s_) const { auto s = collection + s_; return s->a*cur_point.x + s->b < cur_point.y; };//returns true if s is under current point 
  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N << 1;
    for (uint4 i = 0; i < NN; ++i)     epoints[i] = i << 1;
    std::sort(epoints, epoints + NN,
      [x = ends](uint4 pt1, uint4 pt2) {
      return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
    }
    );
  };
  void clone(CLine2SegmentCollection * c, IntersectionRegistrator *r)
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
  void Init(uint4 n, void * c, IntersectionRegistrator *r)
  {
    N = n;
    collection = reinterpret_cast<TLineSegment2*>(c);
    static_assert(sizeof(TLineSegment2) == 4 * sizeof(REAL), "one segment should occupy 4 reals space");
    static_assert(offsetof(TLineSegment2, x1) == 0, "TLineSegment2.x1 should have 0 offset");
    static_assert(offsetof(TLineSegment2, x2) == 2 * sizeof(REAL), "TLineSegment2.x2 should have 2 reals offset");

    ends = reinterpret_cast<REAL*>(collection);
    SetRegistrator(r);
  };

  CLine2SegmentCollection(uint4 n, void* c, IntersectionRegistrator *r)
  {
    Init(n, c,r);
  }

  CLine2SegmentCollection() {};

private:
  inline auto GetX(uint4 pt) const { return ends[pt]; };
  //{ return is_last(pt) ? collection[get_segm(pt)].x2 :collection[get_segm(pt)].x1; };


  IntersectionRegistrator *registrator = nullptr;
  CLine2SegmentCollection *clone_of = nullptr;
  uint4 N=0;
  REAL B, E;
  TPlaneVect cur_point;

  uint4 cur_seg_idx = 0xFFFFFFFF;
  TLineSegment2 cur_seg;
  TLineSegment2 *collection = nullptr;;
  REAL *ends = nullptr;;
};
