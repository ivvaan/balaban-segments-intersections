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
class CLine1SegmentCollection
{
public:
  static constexpr bool is_line_segments = true;
  static constexpr bool has_sentinels = true;

  static inline bool is_last(uint4 pt)
  {
    return pt & 1;
  };
  static inline uint4 get_segm(uint4 pt)
  {
    return pt >> 1;
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
    cur_point = is_last(pt)?
      collection[get_segm(pt)].EndPoint():
      collection[get_segm(pt)].BegPoint();
  };
  void SetCurPointAtBeg(uint4 s)
  {
    cur_point=collection[s].BegPoint();
  };
  void SetCurPointAtEnd(uint4 s)
  {
    cur_point = collection[s].EndPoint();
  };
  void SetCurSeg(uint4 s)
  {
    cur_seg_idx = s;
    cur_seg = collection[s];
    curB = cur_seg.org.x;
    curE = cur_seg.org.x + cur_seg.shift.x;
  };

  void SetCurSegCutBE(uint4 s)
  {
    SetCurSeg(s);
    curB = MAX(B, curB);
    curE = MIN(E, curE);
    active_end.x = curE;
    active_end.y = cur_seg.YAtX(curE);
  };

  void SetCurSegCutBeg(uint4 s)
  {
    SetCurSeg(s);
    curB = MAX(B, curB);
    active_end = cur_seg.EndPoint();
  };

  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    curE = MIN(E, curE);
    active_end = cur_seg.BegPoint();
  };

  void SetCurSegAE(uint4 s)
  {
    SetCurSeg(s);
    active_end = cur_seg.EndPoint();
  };

  bool LBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    auto &s1 = collection[s_1];
    auto &s2 = collection[s_2];
    auto dx1 = s1.shift.x, dx2 = s2.shift.x;
    auto y1 = (s1.org.y*dx1 + (B - s1.org.x)*s1.shift.y)*dx2;
    auto y2 = (s2.org.y*dx2 + (B - s2.org.x)*s2.shift.y)*dx1;
    return y1<y2;
  };

  bool RBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    auto &s1 = collection[s_1];
    auto &s2 = collection[s_2];
    auto dx1 = s1.shift.x, dx2 = s2.shift.x;
    auto y1 = (s1.org.y*dx1 + (E - s1.org.x)*s1.shift.y)*dx2;
    auto y2 = (s2.org.y*dx2 + (E - s2.org.x)*s2.shift.y)*dx1;
    return y1<y2;
  };

  bool FindIntWith(REAL x1, REAL x2, uint4 s_)
  {
    auto &s2 = cur_seg;
    auto &s1 = collection[s_];
    TPlaneVect delt = s2.org - s1.org;
    REAL prod = s1.shift % s2.shift;
    auto mul = s1.shift%delt;
    if ((mul>0) ^ (mul + prod>0))
    {
      mul /= prod;
      auto xc = s2.org.x - mul*s2.shift.x;
      if (((xc <= x1) || (xc>x2))) return false;
      if constexpr(_RegistrationType::point&IntersectionRegistrator::reg_type)
        registrator->register_pair_and_point(cur_seg_idx, s_, TPlaneVect(xc, s2.org.y - mul * s2.shift.y));
      else
        registrator->register_pair(cur_seg_idx, s_);
      return true;
    }
    return false;
  };


  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    auto delt = s2.org - s1.org;
    if ((s1.shift.x < delt.x)||(s2.shift.x < -delt.x))
      return false;
    REAL prod = s1.shift % s2.shift, mul;
    if (((mul = s1.shift % delt) > 0) ^ (mul + prod > 0))
      if (((mul = delt % s2.shift) > 0) ^ (mul - prod > 0))
      {
        if constexpr (_RegistrationType::point & IntersectionRegistrator::reg_type)
          registrator->register_pair_and_point(cur_seg_idx, s_, s1.org + (mul / prod) * s1.shift);
        else
          registrator->register_pair(cur_seg_idx, s_);
        return true;
      };
    return false;
  };

  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    auto delt = s2.org - s1.org;
    REAL prod = s1.shift % s2.shift, mul;
    if (((mul = s1.shift % delt) > 0) ^ (mul + prod > 0))
      if (((mul = delt % s2.shift) > 0) ^ (mul - prod > 0))
      {
        if constexpr (_RegistrationType::point & IntersectionRegistrator::reg_type)
          registrator->register_pair_and_point(cur_seg_idx, s_, s1.org + (mul / prod) * s1.shift);
        else
          registrator->register_pair(cur_seg_idx, s_);
        return true;
      };
    return false;
  };

  bool FindCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    if constexpr (!(_RegistrationType::point & IntersectionRegistrator::reg_type))
      return (search_dir_down ^ UnderActiveEnd(s_)) ?
      registrator->register_pair(cur_seg_idx, s_), true : false;
  
    return FindIntWith(curB, curE, s_);
  };

  bool IsIntersectsCurSegDown(int4 s_) const//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return !UnderActiveEnd(s_);
  };

  bool IsIntersectsCurSegUp(int4 s_) const//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return UnderActiveEnd(s_);
  };

  bool UnderCurPoint(int4 s_) const { return collection[s_].under(cur_point); };//returns true if s is under current point 
  bool UnderActiveEnd(int4 s_) const { return collection[s_].under(active_end); };//returns true if s is under current point 

  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N * 2;
    for (uint4 i = 0; i < NN; ++i)  epoints[i] = i*2;


    for (uint4 i = 0; i < N; ++i) {
      collection[i].shift.x += collection[i].org.x;
    };

    static_assert(sizeof(TLineSegment1) == 4 * sizeof(REAL), "one segment should occupy 4 reals space");
    static_assert(offsetof(TLineSegment1, org.x) == 0, "TLineSegment1.org.x should have 0 offset");
    static_assert(offsetof(TLineSegment1, shift.x) == 2 * sizeof(REAL), "TLineSegment1.shift.x should have 2 reals offset");

    std::sort(epoints, epoints + NN,
      [x = reinterpret_cast<REAL*>(collection)](uint4 pt1, uint4 pt2) {
        return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
      }
    );
    for (uint4 i = 0; i < N; ++i) {
      collection[i].shift.x -= collection[i].org.x;
    };
    for (uint4 i = 0; i < NN; ++i)  epoints[i] /= 2;

    // more memory consume variant
    /*for (uint4 i = 0; i < NN; ++i)  epoints[i] = i;

    struct X1X2 { REAL x1, x2; };
    std::unique_ptr<X1X2[]> buffer(new X1X2[N]);
    X1X2 *x1x2 = buffer.get();

    for (uint4 i = 0; i < N; ++i) {
      x1x2[i].x1 = collection[i].org.x;
      x1x2[i].x2 = collection[i].org.x + collection[i].shift.x;
    };

    std::sort(epoints, epoints + NN,
      [x = reinterpret_cast<REAL*>(x1x2)](uint4 pt1, uint4 pt2) {
        return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
      }
    );
    */
  };

  void clone(CLine1SegmentCollection * c, IntersectionRegistrator *r)
  {
      clone_of = c;
      Init(c->N, c->collection,r);
  };
  //IntersectionRegistrator *GetRegistrator() { return registrator; };

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
  void Init(uint4 n, void * c, IntersectionRegistrator *r)
  {
    N = n;
    collection = reinterpret_cast<TLineSegment1*>(c);
    SetRegistrator(r);

  };

  CLine1SegmentCollection(uint4 n, void* c, IntersectionRegistrator *r)
  {
    Init(n, c, r);
  }
  CLine1SegmentCollection() {};
  
  void SetSearchDirDown(bool dir) { search_dir_down = dir; };

private:
  inline auto GetX(uint4 pt) const
  {
    return is_last(pt) ? collection[get_segm(pt)].org.x + collection[get_segm(pt)].shift.x : collection[get_segm(pt)].org.x;
  };


  IntersectionRegistrator *registrator = nullptr;
  CLine1SegmentCollection *clone_of = nullptr;
  uint4 N=0;
  REAL B, E, curB, curE;
  TPlaneVect cur_point, active_end;

  uint4 cur_seg_idx = 0xFFFFFFFF;
  TLineSegment1 cur_seg;
  TLineSegment1 *collection = nullptr;
  bool search_dir_down = true;
};

