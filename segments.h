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
#include "utils.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include "math.h"
#include <cstddef>



class TLineSegment1
{
public:
  static const int4 is_line = 1;
  TLineSegment1():org(),shift() {};
  TLineSegment1(TPlaneVect &p1, TPlaneVect &p2) :org(p1), shift(p2-p1) {};
  // int4 parnum,intpar;
  TPlaneVect org;
  TPlaneVect shift;
  //TPlaneVect cash_pt;
  void Refine()
  {
    if ((shift.x<0) || ((shift.x == 0) && (shift.y<0)))
    {
      org = org + shift; 
      shift = -shift;
    }
  };
  void Init(TLineSegment1 &s) { org = s.org; shift = s.shift; };
  void InitRandom(CRandomValueGen &rv, int4 seg_n, int4 type, REAL par);
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
  void write_SVG(int4 id,chostream *SVG_text) {
    if (SVG_text) {
      auto bp = BegPoint();
      auto ep = EndPoint();
      *SVG_text << "<line id='seg" << id;
      *SVG_text << "' x1='" << bp.x;
      *SVG_text << "' y1='" << bp.y;
      *SVG_text << "' x2='" << ep.x;
      *SVG_text << "' y2='" << ep.y;
      *SVG_text << "' class='line1'/>\n";
    }
  
  };
};
  bool below(REAL X, TLineSegment1 *s1, TLineSegment1 *s2);
  template <bool _ret_ip>  int4 IntPoint(TLineSegment1 *s1, TLineSegment1 *s2, TPlaneVect *p);
  template <bool _ret_ip>  int4 StripePoint(REAL b, REAL e, TLineSegment1 *s1, TLineSegment1 *s2, TPlaneVect *p);
  //for optimal alg only
  int4 IntInside(REAL b, REAL e, TLineSegment1 *s1, TLineSegment1 *s2);


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
  bool upper(const TPlaneVect &v)
  {
    return v.y<(v.x*a + b);
  };
  REAL YAtX(REAL X)
  {
    return a*X + b;
  };
  void write_SVG(int4 id,chostream* SVG_text) {
    if (SVG_text) {
      auto bp = BegPoint();
      auto ep = EndPoint();
      *SVG_text << "<line id='seg" << id;
      *SVG_text << "' x1='" << bp.x;
      *SVG_text << "' y1='" << bp.y;
      *SVG_text << "' x2='" << ep.x;
      *SVG_text << "' y2='" << ep.y;
      *SVG_text << "' class='line2'/>\n";
    }
  };
  };
  bool below(REAL X, TLineSegment2 *s1, TLineSegment2 *s2);
  template <bool _ret_ip>  int4 IntPoint(TLineSegment2 *s1, TLineSegment2 *s2, TPlaneVect *p);
  template <bool _ret_ip>  int4 StripePoint(REAL b, REAL e, TLineSegment2 *s1, TLineSegment2 *s2, TPlaneVect *p);
  //for optimal alg only
  int4 IntInside(REAL b, REAL e, TLineSegment2 *s1, TLineSegment2 *s2);



class TArcSegment// arc formed by intersection of a vertical strip [x1,x2] with a circle having center org and square of radius r2 
{
public:
  static const int4 is_line = 0;
  REAL x1 = 0;
  TPlaneVect org;
  REAL x2 = 0, r2 = 0;

  uint4 is_upper;// upper or lower part of the intersection
  uint4 reserved;// just to have record size 48 bytes
                 //TArcSegment() {};
                 //void Init(TLineSegment1 &s){x1=s.org.x;x2=s.org.x+s.shift.x;a=s.shift.y/s.shift.x;b=s.org.y-a*s.org.x;};
  void InitRandom(CRandomValueGen &rv, int4 seg_n, int4 type, REAL par);


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

  bool under(const TPlaneVect &v) const//arc placed under point v
  {
    bool is_lower = !is_upper;
    if (is_lower ^ (v.y < org.y))return is_lower;
    return is_lower ^ ((v - org).get_norm() > r2);
    /*if (is_upper)
      return (v.y>org.y) && ((v - org).get_norm()>r2);
    return (v.y>org.y) || ((v - org).get_norm()<r2);*/ 
  };
  bool upper(const TPlaneVect& v)
  {
    return !under(v);
  };
  bool IsTheSamePart(const TPlaneVect &v)
  {
    return (!is_upper) ^ (v.y > org.y);
  };
  void write_SVG(int4 id, chostream* SVG_text) {
    if (SVG_text) {
      auto bp = BegPoint();
      auto ep = EndPoint();
      auto r = sqrt(r2);
      *SVG_text << "<path id='arc" << id;
      *SVG_text << "' d='M " << bp.x<<" "<<bp.y<<" ";
      *SVG_text << "A " << r << " " << r << " 0 0 ";
      *SVG_text << (is_upper?"0 ":"1 ");
      *SVG_text << ep.x << " " << ep.y << "'";
      *SVG_text << "' class='arc'/>\n";
    }

  }

};
  bool below(REAL X, TArcSegment *s1, TArcSegment *s2);
  int4 IntPointsInStripe(REAL x1, REAL x2, TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  template <bool _ret_ip>  int4 IntPoint(TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  template <bool _ret_ip>  int4 StripePoint(REAL b, REAL e, TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  //for optimal alg only
  int4 IntInside(REAL b, REAL e, TArcSegment *s1, TArcSegment *s2);

