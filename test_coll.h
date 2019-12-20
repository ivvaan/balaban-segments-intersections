#ifndef TEST_COLL_H
#define TEST_COLL_H
/*
*
*      Copyright (c)  2011  Ivan Balaban 
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

#include "utils.h"
#include "intersection_finder.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include "math.h"
#include <cstddef>

enum _Algorithm
  {
  triv=1,
  simple_sweep=2,
  fast=4,
  optimal=8,
  fast_parallel=16,
  bentley_ottmann=32,
  fast_no_ip=64,
  mem_save=128
  };

enum _Segment
  {
  line1=0,line2,arc
  };
enum _Distribution
  {
  random=0,parallel,mixed,small,param_defined,circul
  };
const int4 n_threads = 6;

PSeg* create_test_collection(int4 seg_type,int4 n,int4 distr,REAL par);
void  delete_test_collection(int4 seg_type,PSeg* );
double find_intersections(int4 seg_type, int4 SN, PSeg* colls, int4 alg, double* counters, bool dont_need_ip = false);

double new_find_int(int4 seg_type, int4 n, PSeg segs,int4 alg);
//#define PRINT_SEG_AND_INT
class TLineSegment1
{
public:
  static const int4 is_line = 1;
  TLineSegment1() {};
  // int4 parnum,intpar;
  TPlaneVect org;
  TPlaneVect shift;
  //TPlaneVect cash_pt;
  void Refine()
  {
    if ((shift.x<0) || ((shift.x == 0) && (shift.y<0)))
    {
      org = org + shift; shift = -shift;
    }
  };
  void Init(TLineSegment1 &s) { org = s.org; shift = s.shift; };
  void InitRandom(CRandomValueGen &rv, int4 seg_n, int4 type, REAL par)
  {
    if (circul == type)
    {
      //Segments ends are on the opposite sides of unit circul, each segment intesect each
      double angle = 1.0 / 16.0+ 3.*M_PI*rv.GetRandomDouble() / 8.0;
      org.x = 1. - cos(angle);
      org.y = 1. - sin(angle);
      shift.x = 1. + cos(angle) - org.x;
      shift.y = 1. + sin(angle) - org.y;
      Refine();

#ifdef PRINT_SEG_AND_INT
      printf("[%6.3f,%6.3f,%6.3f,%6.3f],\n", org.x, org.y, org.x + shift.x, org.y + shift.y);
#endif
      return;
    }
    org.x = rv.GetRandomDouble();
    org.y = rv.GetRandomDouble();
    shift.x = rv.GetRandomDouble() - org.x;
    shift.y = rv.GetRandomDouble() - org.y;
    Refine();
    if (parallel == type)
      shift.y *= par;
    if (mixed == type)
    {
      org.y = org.y + 0.5*(1.0 - par)*shift.y;
      shift.y *= par;
      if (seg_n % 3 != 0)
      {
        org.x = org.x + 0.5*(1.0 - par)*shift.x;
        shift.x *= par;
      }
        
    }
    if (type>mixed)
    {
      org = org + 0.5*(1.0 - par)*shift;
      shift = par*shift;
    }
#ifdef PRINT_SEG_AND_INT
    printf("[%6.3f,%6.3f,%6.3f,%6.3f],\n", org.x, org.y, org.x + shift.x, org.y + shift.y);
#endif
  };
  TPlaneVect BegPoint()
  {
    return org;
  };
  TPlaneVect EndPoint()
  {
    return org + shift;
  };
  void BegPoint(REAL &x, REAL &y)
  {
    x = org.x; y = org.y;
  };
  void EndPoint(REAL &x, REAL &y)
  {
    x = org.x + shift.x; y = org.y + shift.y;
  };
  bool under(const TPlaneVect &v) //segment placed under point v
  {
    REAL res = (v - org) % shift;
    //   if (fabs(res)==0.0) throw 1;
    return (res <= 0);
  };
  bool upper(TPlaneVect &v)
  {
    REAL res = (v - org) % shift;
    //   if (fabs(res)==0.0) throw 1;
    return (res>0);
  };
  REAL YAtX(REAL X)
  {
    return org.y + shift.y*(X - org.x) / shift.x;
  };
  friend bool below(REAL X, TLineSegment1 *s1, TLineSegment1 *s2);
  template <bool _ret_ip> friend int4 IntPoint(TLineSegment1 *s1, TLineSegment1 *s2, TPlaneVect *p);
  template <bool _ret_ip> friend int4 StripePoint(REAL b, REAL e, TLineSegment1 *s1, TLineSegment1 *s2, TPlaneVect *p);
  //for optimal alg only
  friend int4 IntInside(REAL b, REAL e, TLineSegment1 *s1, TLineSegment1 *s2);
};

class TLineSegment2
{
public:
  //int4 seg_n;
  static const int4 is_line = 1;
  REAL x1, a, x2, b;
  TLineSegment2() { x1 = 0; x2 = 0; a = 0; b = 0; };
  void Init(TLineSegment1 &s) { x1 = s.org.x; x2 = s.org.x + s.shift.x; a = s.shift.y / s.shift.x; b = s.org.y - a*s.org.x; };
  void InitRandom(CRandomValueGen &rv, int4 s_n, int4 type, REAL par)
  {
    TLineSegment1 seg;
    //seg_n=s_n;
    seg.InitRandom(rv, s_n, type, par);
    Init(seg);
  };


  TPlaneVect BegPoint()
  {
    return TPlaneVect(x1, a*x1 + b);
  };
  TPlaneVect EndPoint()
  {
    return TPlaneVect(x2, a*x2 + b);
  };
  void BegPoint(REAL &x, REAL &y)
  {
    x = x1; y = a*x1 + b;
  };
  void EndPoint(REAL &x, REAL &y)
  {
    x = x2; y = a*x2 + b;
  };
  bool under(const TPlaneVect &v)
  {
    return v.y >= (v.x*a + b);
  };
  bool upper(TPlaneVect &v)
  {
    return v.y<(v.x*a + b);
  };
  REAL YAtX(REAL X)
  {
    return a*X + b;
  };
  friend bool below(REAL X, TLineSegment2 *s1, TLineSegment2 *s2);
  template <bool _ret_ip> friend int4 IntPoint(TLineSegment2 *s1, TLineSegment2 *s2, TPlaneVect *p);
  template <bool _ret_ip> friend int4 StripePoint(REAL b, REAL e, TLineSegment2 *s1, TLineSegment2 *s2, TPlaneVect *p);
  //for optimal alg only
  friend int4 IntInside(REAL b, REAL e, TLineSegment2 *s1, TLineSegment2 *s2);

};

class TArcSegment// arc formed by intersection of a vertical strip [x1,x2] with a circle having center org and square of radius r2 
{
public:
  static const int4 is_line = 0;
  REAL x1=0;
  TPlaneVect org; 
  REAL x2=0, r2=0;
  
  uint4 is_upper;// upper or lower part of the intersection
  uint4 reserved;// just to have record size 48 bytes
  //TArcSegment() {};
  //void Init(TLineSegment1 &s){x1=s.org.x;x2=s.org.x+s.shift.x;a=s.shift.y/s.shift.x;b=s.org.y-a*s.org.x;};
  void InitRandom(CRandomValueGen &rv, int4 seg_n, int4 type, REAL par)
  {
    REAL tmp, r;
    x1 = rv.GetRandomDouble();
    tmp = rv.GetRandomDouble();
    if (x1<tmp)x2 = tmp;
    else { x2 = x1; x1 = tmp; }
    if (type<2) par = 1.0;
    REAL delta = x2 - x1;
    x1 = 0.5*(x1*(1.0 + par) + x2*(1.0 - par));
    x2 = x1 + par*delta;
    org.x = rv.GetRandomDouble();
    org.y = rv.GetRandomDouble();
    tmp = max(fabs(x1 - org.x), fabs(x2 - org.x));
    par = max(1.0, par);
    while ((r = par*rv.GetRandomDouble()) <= tmp);
    r2 = r*r;
    is_upper = rv.RandomChoose();
  }


  TPlaneVect BegPoint()
  {
    if (is_upper)
      return TPlaneVect(x1, org.y + sqrt(r2 - sq(x1 - org.x)));
    else
      return TPlaneVect(x1, org.y - sqrt(r2 - sq(x1 - org.x)));
  };
  TPlaneVect EndPoint()
  {
    if (is_upper)
      return TPlaneVect(x2, org.y + sqrt(r2 - sq(x2 - org.x)));
    else
      return TPlaneVect(x2, org.y - sqrt(r2 - sq(x2 - org.x)));
  };
  TPlaneVect PointAtX(REAL x)
  {
    if (is_upper)
      return TPlaneVect(x, org.y + sqrt(r2 - sq(x - org.x)));
    else
      return TPlaneVect(x, org.y - sqrt(r2 - sq(x - org.x)));
  };
  REAL YAtX(REAL X)
  {
    if (is_upper)
      return org.y + sqrt(r2 - sq(X - org.x));
    else
      return org.y - sqrt(r2 - sq(X - org.x));
  };
  void BegPoint(REAL& x, REAL& y)
  {
    x = x1; y = is_upper ? org.y + sqrt(r2 - sq(x1 - org.x)) : org.y - sqrt(r2 - sq(x1 - org.x));
  };
  void EndPoint(REAL &x, REAL &y)
  {
    x = x2; y = is_upper ? org.y + sqrt(r2 - sq(x2 - org.x)) : org.y - sqrt(r2 - sq(x2 - org.x));
  };

  bool under(const TPlaneVect &v)//arc placed under point v
  {
    if (is_upper)
      return (v.y>org.y) && ((v - org).get_norm()>r2);
    return (v.y>org.y) || ((v - org).get_norm()<r2);
  };
  bool upper(const TPlaneVect& v)
  {
    return !under(v);
  };
  int4 IsTheSamePart(TPlaneVect &v)
  {
    if (is_upper)return v.y >= org.y;
    return v.y <= org.y;
  };
  friend bool below(REAL X, TArcSegment *s1, TArcSegment *s2);
  friend int4 IntPointsInStripe(REAL x1, REAL x2, TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  template <bool _ret_ip> friend int4 IntPoint(TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  template <bool _ret_ip> friend int4 StripePoint(REAL b, REAL e, TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  //for optimal alg only
  friend int4 IntInside(REAL b, REAL e, TArcSegment *s1, TArcSegment *s2);

};

enum _RegistrationType
{
  count = 1,
  segments = 2,
  point = 4,
  count_and_segments =3,
  full=7
};
template <class ipoint>
class JustCountingRegistrator
{
public:
  static const _RegistrationType reg_type = count;
  double counter=0;
  inline void begin_registration(uint4 inters_numb) { counter += inters_numb; };
  inline void register_segments(uint4 s1, uint4 s2) 
  {
#ifdef PRINT_SEG_AND_INT
      if (s1<s2)
      printf("alt int %i %i\n", s1,s2);
      else
      printf("alt int %i %i\n", s2, s1);

#endif 
  };
  inline void register_points(ipoint*) {};
  inline void end_registration() {};


};

typedef JustCountingRegistrator<TPlaneVect> SimpleCounter;

template<class IntersectionRegistrator>
class CLine2SegmentCollection
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
    E =  GetX(right);
  };
  inline void SetCurStripeRight(uint4 right) {  E = GetX(right); };
  inline void SetCurStripeLeft(uint4 left) { B = GetX(left); };
  void SetCurPoint(uint4 pt) 
  { 
    if(is_last(pt))
      collection[get_segm(pt)].EndPoint(cur_point.x,cur_point.y);
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
  inline void SetCurVLine(uint4 pt) { X = GetX(pt); };
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
  bool Below(int4 s_1, int4 s_2)//retuns if s1 below s2 at current vertical line
  {
    auto s1 = collection + s_1;
    auto s2 = collection + s_2;
    return ((X*(s2->a-s1->a) + s2->b-s1->b)>0);
  };


  int4 TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection+s_;
    auto x1 = max(cur_seg.x1, s->x1);
    auto x2 = min(cur_seg.x2, s->x2);
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
        registrator->begin_registration(1);
        registrator->register_segments(cur_seg_idx, s_);
        registrator->register_points(&p);
        registrator->end_registration();
        return true;
      }
    }
    else
    {
      auto x = (s->b - cur_seg.b) / da;
      if ((x >= x1) && (x <= x2))
      {
        registrator->begin_registration(1);
        registrator->register_segments(cur_seg_idx, s_);
        registrator->end_registration();
        return true;
      }
    }
    return false;

  };

  int4 SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto da = cur_seg.a - s->a;
    if (da == 0)return false;
    auto x1 = max(cur_seg.x1, s->x1);
    auto x2 = min(cur_seg.x2, s->x2);
    if (_RegistrationType::point&IntersectionRegistrator::reg_type)
    {
      TPlaneVect p;
      p.x = (s->b - cur_seg.b) / da;
      if ((p.x >= x1) && (p.x <= x2))
      {
        p.y = p.x*cur_seg.a + cur_seg.b;
        registrator->begin_registration(1);
        registrator->register_segments(cur_seg_idx, s_);
        registrator->register_points(&p);
        registrator->end_registration();
        return true;
      }
    }
    else
    {
      auto x = (s->b - cur_seg.b) / da;
      if ((x >= x1) && (x <= x2))
      {
        registrator->begin_registration(1);
        registrator->register_segments(cur_seg_idx, s_);
        registrator->end_registration();
        return true;
      }
    }
    return false;

  };

  int4 FindCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
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
        registrator->begin_registration(1);
        registrator->register_segments(cur_seg_idx, s_);
        registrator->register_points(&p);
        registrator->end_registration();
        return true;
      }
    }
    else
    {
      auto x = (s->b - cur_seg.b) / da;
      if ((x >= cur_seg.x1) && (x <= cur_seg.x2))
      {
        registrator->begin_registration(1);
        registrator->register_segments(cur_seg_idx, s_);
        registrator->end_registration();
        return true;
      }
    }
    return false;

  };

  int4 IsIntersectsCurSeg(int4 s_)//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    auto s = collection+s_;
    auto da = cur_seg.a - s->a;
    if (da == 0)return false;
    auto x = (s->b - cur_seg.b) / da;
    return ((x >= cur_seg.x1) && (x <= cur_seg.x2));
  };

  bool UnderCurPoint(int4 s_) { auto s = collection + s_; return s->a*cur_point.x+s->b < cur_point.y; };//returns true if s is under current point 
  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N << 1;
      for (uint4 i = 0; i < NN; ++i)     epoints[i] = i<<1;
      std::sort (epoints, epoints + NN, 
            [x = ends](uint4 pt1, uint4 pt2) {
            return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
            }
      );
  };
  void clone(CLine2SegmentCollection * c, IntersectionRegistrator *r)
  {
    clone_of = c;
    N=c->N;
    registrator = r;
    collection = c->collection;
    ends = reinterpret_cast<REAL*>(collection);
  };
  void SetRegistrator(IntersectionRegistrator *r) { registrator = r; };
  //IntersectionRegistrator *GetRegistrator() { return registrator; };

  void unclone() { if (clone_of == nullptr)return; collection = nullptr; clone_of = nullptr; };
  void SortAt(uint4 pt, uint4 n, int4 *L) 
  { 
    SetCurVLine(pt); 
    std::sort(L, L+n, [this](int4 s1, int4 s2) {return Below(s1, s2); });
  };
  
  void Init(uint4 n, void * c)
  {
      N = n;
      collection = reinterpret_cast<TLineSegment2*>(c);
      static_assert(sizeof(TLineSegment2) == 4 * sizeof(REAL),"one segment should occupy 4 reals space");
      static_assert(offsetof(TLineSegment2, x1) == 0,"TLineSegment2.x1 should have 0 offset");
      static_assert(offsetof(TLineSegment2, x2) == 2 * sizeof(REAL),"TLineSegment2.x2 should have 2 reals offset");

      ends = reinterpret_cast<REAL*>(collection);
  };

  private:
  inline auto GetX(uint4 pt){return ends[pt];}; 
  //{ return is_last(pt) ? collection[get_segm(pt)].x2 :collection[get_segm(pt)].x1; };


  IntersectionRegistrator *registrator=nullptr;
  CLine2SegmentCollection *clone_of = nullptr;
  uint4 N;
  REAL B, E, X;
  TPlaneVect cur_point;

  uint4 cur_seg_idx=0xFFFFFFFF;
  TLineSegment2 cur_seg;
  TLineSegment2 *collection;
  REAL *ends;
};
//******************************************************************************************************

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
  inline void SetCurVLine(uint4 pt) { X = GetX(pt); };
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

  bool Below(int4 s_1, int4 s_2)//retuns if s1 below s2 at current vertical line
  {
    return (s_1 != s_2) && (collection[s_1].under(collection[s_2].PointAtX(X)));
  };

#define REGISTER_INTERSECTION(s1,s2,p)\
  if (_RegistrationType::point&IntersectionRegistrator::reg_type)\
  {\
    registrator->begin_registration(1);\
    registrator->register_segments(s1, s2-collection);\
    registrator->register_points(&p);\
    registrator->end_registration();\
  }\
  else\
  {\
    registrator->begin_registration(1);\
    registrator->register_segments(s1, s2-collection);\
    registrator->end_registration();\
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



  int4 TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = max(cur_seg.x1, s->x1);
    auto x2 = min(cur_seg.x2, s->x2);
    if (x1 >= x2)return false;
    return IntPointsInStripe<false>(x1, x2, s);

  };

  int4 SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto s = collection + s_;
    auto x1 = max(cur_seg.x1, s->x1);
    auto x2 = min(cur_seg.x2, s->x2);
    return IntPointsInStripe<false>(x1, x2, s);
  };

  int4 FindCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    return IntPointsInStripe<false>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  int4 IsIntersectsCurSeg(int4 s_)//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return IntPointsInStripe<true>(cur_seg.x1, cur_seg.x2, collection + s_);
  };

  bool UnderCurPoint(int4 s_) { return collection[s_].under(cur_point); };//returns true if s is under current point 
  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = N << 1;
    for (uint4 i = 0; i < NN; ++i)     epoints[i] = i * 3 ;
    std::sort(epoints, epoints + NN,
      [x = ends](uint4 pt1, uint4 pt2) {
      return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
    }
    );
  };
  void clone(CArcSegmentCollection * c, IntersectionRegistrator *r)
  {
    clone_of = c;
    N = c->N;
    registrator = r;
    collection = c->collection;
    ends = reinterpret_cast<REAL*>(collection);
  };
  void SetRegistrator(IntersectionRegistrator *r) { registrator = r; };
  //IntersectionRegistrator *GetRegistrator() { return registrator; };

  void unclone() { if (clone_of == nullptr)return; collection = nullptr; clone_of = nullptr; };
  void SortAt(uint4 pt, uint4 n, int4 *L)
  {
    SetCurVLine(pt);
    std::sort(L, L + n, [this](int4 s1, int4 s2) {return Below(s1, s2); });
  };

  void Init(uint4 n, void* c)
  {
    N = n;
    collection = reinterpret_cast<TArcSegment*>(c);

    static_assert(sizeof(TArcSegment) == 6 * sizeof(REAL), "one segment should occupy 6 reals space");
    static_assert(offsetof(TArcSegment, x1) == 0, "TArcSegment.x1 should have 0 offset");
    static_assert(offsetof(TArcSegment, x2) == 3 * sizeof(REAL), "TArcSegment.x2 should have 3 reals offset");
    ends = reinterpret_cast<REAL*>(collection);
  };

private:
  inline auto GetX(uint4 pt) { return ends[pt]; };
  //{ return is_last(pt) ? collection[get_segm(pt)].x2 :collection[get_segm(pt)].x1; };


  IntersectionRegistrator *registrator = nullptr;
  CArcSegmentCollection *clone_of = nullptr;
  uint4 N;
  REAL B, E, X;
  TPlaneVect cur_point;

  uint4 cur_seg_idx = 0xFFFFFFFF;
  TArcSegment cur_seg;
  TArcSegment *collection;
  REAL *ends;
};


#endif