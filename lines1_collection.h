#pragma once
#include "segments.h"
template<class IntersectionRegistrator>
class CLine1SegmentCollection
{
public:
  static const bool is_line_segments = true;

  static inline bool is_last(uint4 pt)
  {
    return pt & 2;
  };
  static inline uint4 get_segm(uint4 pt)
  {
    return pt >> 2;
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
    curB = max(B, curB);
    curE = min(E, curE);
  };
  void SetCurSegCutBeg(uint4 s)
  {
    SetCurSeg(s);
    curB = max(B, curB);
  };
  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    curE = min(E, curE);
  };
  inline void SetCurSeg(uint4 s)
  {
    cur_seg_idx = s;
    cur_seg = collection[s];
    curB = cur_seg.org.x;
    curE = cur_seg.org.x + cur_seg.shift.x;
  };
  bool LBelow(int4 s_1, int4 s_2) //retuns if s1 below s2 at current vertical line
  {
    auto s1 = collection + s_1;
    auto s2 = collection + s_2;
    auto dx1 = s1->shift.x, dx2 = s2->shift.x;
    auto dy1 = s1->shift.y, dy2 = s2->shift.y;
    auto y1 = (s1->org.y*dx1 + (B - s1->org.x)*dy1)*dx2;
    auto y2 = (s2->org.y*dx2 + (B - s2->org.x)*dy2)*dx1;
    return y1<y2;
  };
  bool RBelow(int4 s_1, int4 s_2) //retuns if s1 below s2 at current vertical line
  {
    auto s1 = collection + s_1;
    auto s2 = collection + s_2;
    auto dx1 = s1->shift.x, dx2 = s2->shift.x;
    auto dy1 = s1->shift.y, dy2 = s2->shift.y;
    auto y1 = (s1->org.y*dx1 + (E - s1->org.x)*dy1)*dx2;
    auto y2 = (s2->org.y*dx2 + (E - s2->org.x)*dy2)*dx1;
    return y1<y2;
  };

  template <bool register_int>
  bool FindIntWith(REAL x1, REAL x2, uint4 s_)
  {
    auto s2 = &cur_seg;
    auto s1 = collection + s_;
    TPlaneVect delt = s2->org - s1->org;
    REAL prod = s1->shift%s2->shift, mul, xc;
    mul = s1->shift%delt;
    if ((mul>0) ^ (mul + prod>0))
    {
      if (prod == 0) return false;
      mul = mul / prod;
      xc = s2->org.x - mul*s2->shift.x;
      if (((xc <= x1) || (xc>x2))) return false;
      if (register_int)
      {
        if (_RegistrationType::point&IntersectionRegistrator::reg_type)
        {
          TPlaneVect p(xc, s2->org.y - mul*s2->shift.y);
          registrator->begin_registration(1);
          registrator->register_segments(cur_seg_idx, s_);
          registrator->register_points(&p);
          registrator->end_registration();
        }
        else
        {
          registrator->begin_registration(1);
          registrator->register_segments(cur_seg_idx, s_);
          registrator->end_registration();
        }
      }
      return true;
    }
    return false;
  };


  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = max(curB, s->org.x);
    auto x2 = min(curE, s->org.x + s->shift.x);
    if (x1 >= x2)return false;
    return FindIntWith<true>(x1, x2, s_);
  };

  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    return FindIntWith<true>(curB, curE, s_);
  };

  bool FindCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    return FindIntWith<true>(curB, curE, s_);
  };

  bool IsIntersectsCurSeg(int4 s_)//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return FindIntWith<false>(curB, curE, s_);
  };

  bool UnderCurPoint(int4 s_) { return collection[s_].under(cur_point); };//returns true if s is under current point 
  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N << 1;
    for (uint4 i = 0; i < NN; ++i)     epoints[i] = i << 1;
    std::sort(epoints, epoints + NN,
      [this](uint4 pt1, uint4 pt2) {
      auto x1 = GetX(pt1);
      auto x2 = GetX(pt2);
      return ((x1 < x2) || ((x1 == x2) && (pt1 < pt2)));
    }
    );
  };
  void clone(CLine1SegmentCollection * c, IntersectionRegistrator *r)
  {
    clone_of = c;
    N = c->N;
    registrator = r;
    collection = c->collection;
  };
  void SetRegistrator(IntersectionRegistrator *r) { registrator = r; };
  //IntersectionRegistrator *GetRegistrator() { return registrator; };

  void unclone() { if (clone_of == nullptr)return; collection = nullptr; clone_of = nullptr; };
  void SortAt(uint4 pt, uint4 n, int4 *L)
  {
    SetCurStripeLeft(pt);
    std::sort(L, L + n, [this](int4 s1, int4 s2) {return LBelow(s1, s2); });
  };

  void Init(uint4 n, void * c)
  {
    N = n;
    collection = reinterpret_cast<TLineSegment1*>(c);
  };

private:
  inline auto GetX(uint4 pt)
  {
    return is_last(pt) ? collection[get_segm(pt)].org.x + collection[get_segm(pt)].shift.x : collection[get_segm(pt)].org.x;
  };


  IntersectionRegistrator *registrator = nullptr;
  CLine1SegmentCollection *clone_of = nullptr;
  uint4 N;
  REAL B, E, curB, curE;
  TPlaneVect cur_point;

  uint4 cur_seg_idx = 0xFFFFFFFF;
  TLineSegment1 cur_seg;
  TLineSegment1 *collection;
};

