#pragma once
/*
*
*      Copyright (c)  2011-2026  Ivan Balaban
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

constexpr const double crd_scaler = 1.e-3;


class TLineSegment1
{
public:
  constexpr static const uint4 is_line = 1;
  TLineSegment1():org(),shift() {};
  TLineSegment1(const TPlaneVect &p1, const TPlaneVect &p2) :org(p1), shift(p2-p1) {
    Refine();
  };
  TLineSegment1(REAL xb, REAL yb, REAL xs, REAL ys) :org(xb, yb), shift(xs, ys) {
    Refine();
  };

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
  void Rotate(double s,double c) {// s = sin(angle), c = cos(angle);
    assert(fabs(s * s + c * c - 1.) < 1.e-12);
    org = { c * org.x - s * org.y,  
            s * org.x + c * org.y };
    shift = { c * shift.x - s * shift.y,  
              s * shift.x + c * shift.y };
    Refine();
  };

  TLineSegment1 get_rotated(double s, double c) const {// s = sin(angle), c = cos(angle);
    assert(fabs(s * s + c * c - 1.) < 1.e-12);
    TLineSegment1 res;
    res.org = { c * org.x - s * org.y,
            s * org.x + c * org.y };
    res.shift = { c * shift.x - s * shift.y,
              s * shift.x + c * shift.y };
    res.Refine();
    return res;
  };
  void Init(const TLineSegment1 &s) { org = s.org; shift = s.shift; };
  void InitRandom(CRandomValueGen &rv, uint4 seg_n, uint4 type, REAL par);
  TPlaneVect BegPoint() const {
    return org;
  };
  auto bx() const {
    return org.x;
  }
  auto by() const {
    return org.y;
  }
  auto ex() const {
    return org.x + shift.x;
  }
  auto ey() const {
    return org.y + shift.y;
  }
  auto sx() const {
    return shift.x;
  }
  auto sy() const {
    return shift.y;
  }
  TPlaneVect EndPoint() const {
    return org + shift;
  };
  bool under(const TPlaneVect &v) const { //segment placed under point v
    return (v - org) % shift <= 0;
  };
  bool upper(TPlaneVect &v) const {
    return (v - org) % shift > 0;
  };
  REAL YAtX(REAL X) const {
    return org.y + shift.y*(X - org.x) / shift.x;
  };   
  REAL YAtX_Numerator(REAL X) const {
    return shift.x * org.y + shift.y * (X - org.x);
  };
  void write_SVG(int4 id,chostream *SVG_text) const {
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
  template <bool _ret_ip>  uint4 IntPoint(TLineSegment1 *s1, TLineSegment1 *s2, TPlaneVect *p);
  template <bool _ret_ip>  uint4 StripePoint(REAL b, REAL e, TLineSegment1 *s1, TLineSegment1 *s2, TPlaneVect *p);
  //for optimal alg only
  uint4 IntInside(REAL b, REAL e, TLineSegment1 *s1, TLineSegment1 *s2);


class TLineSegment2
{
public:
  //int4 seg_n;
  constexpr static const uint4 is_line = 1;
  REAL x1, a, x2, b;
  TLineSegment2() { x1 = 0; x2 = 0; a = 0; b = 0; };
  void Init(const TLineSegment1 &s) { x1 = s.org.x; x2 = s.org.x + s.shift.x; a = s.shift.y / s.shift.x; b = s.org.y - a*s.org.x; };
  void InitRandom(CRandomValueGen &rv, uint4 s_n, uint4 type, REAL par)
  {
    TLineSegment1 seg;
    //seg_n=s_n;
    seg.InitRandom(rv, s_n, type, par);
    Init(seg);
  };


  TPlaneVect BegPoint() const {
    return TPlaneVect(x1, a*x1 + b);
  };
  TPlaneVect EndPoint() const {
    return TPlaneVect(x2, a*x2 + b);
  };
  void BegPoint(REAL &x, REAL &y) const {
    x = x1; y = a*x1 + b;
  };
  void EndPoint(REAL &x, REAL &y) const {
    x = x2; y = a*x2 + b;
  };
  bool under(const TPlaneVect &v) const {
    return v.y >= (v.x*a + b);
  };
  bool upper(const TPlaneVect &v) const {
    return v.y<(v.x*a + b);
  };
  REAL YAtX(REAL X) const {
    return a*X + b;
  };
  void write_SVG(int4 id,chostream* SVG_text) const {
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
  template <bool _ret_ip>  uint4 IntPoint(TLineSegment2 *s1, TLineSegment2 *s2, TPlaneVect *p);
  template <bool _ret_ip>  uint4 StripePoint(REAL b, REAL e, TLineSegment2 *s1, TLineSegment2 *s2, TPlaneVect *p);
  //for optimal alg only
  uint4 IntInside(REAL b, REAL e, TLineSegment2 *s1, TLineSegment2 *s2);



class TArcSegment// arc formed by intersection of a vertical strip [x1,x2] with a circle having center org and square of radius r2 
{
public:
  constexpr static const uint4 is_line = 0;
  REAL x1 = 0;
  TPlaneVect org;
  REAL x2 = 0, r2 = 0;

  uint4 is_upper;// upper or lower part of the intersection
  uint4 reserved;// just to have record size 48 bytes
                 //TArcSegment() {};
                 //void Init(TLineSegment1 &s) const {x1=s.org.x;x2=s.org.x+s.shift.x;a=s.shift.y/s.shift.x;b=s.org.y-a*s.org.x;};
  void InitRandom(CRandomValueGen &rv, uint4 seg_n, uint4 type, REAL par);


  TPlaneVect BegPoint() const {
    if (is_upper)
      return TPlaneVect(x1, org.y + sqrt(r2 - sq(x1 - org.x)));
    else
      return TPlaneVect(x1, org.y - sqrt(r2 - sq(x1 - org.x)));
  };
  TPlaneVect EndPoint() const {
    if (is_upper)
      return TPlaneVect(x2, org.y + sqrt(r2 - sq(x2 - org.x)));
    else
      return TPlaneVect(x2, org.y - sqrt(r2 - sq(x2 - org.x)));
  };
  TPlaneVect PointAtX(REAL x) const {
    if (is_upper)
      return TPlaneVect(x, org.y + sqrt(r2 - sq(x - org.x)));
    else
      return TPlaneVect(x, org.y - sqrt(r2 - sq(x - org.x)));
  };
  REAL YAtX(REAL X) const {
    if (is_upper)
      return org.y + sqrt(r2 - sq(X - org.x));
    else
      return org.y - sqrt(r2 - sq(X - org.x));
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
  bool upper(const TPlaneVect& v) const {
    return !under(v);
  };
  bool IsTheSamePart(const TPlaneVect &v) const {
    return (!is_upper) ^ (v.y > org.y);
  };
  void write_SVG(uint4 id, chostream* SVG_text) const {
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
  uint4 IntPointsInStripe(REAL x1, REAL x2, TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  template <bool _ret_ip>  uint4 IntPoint(TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  template <bool _ret_ip>  uint4 StripePoint(REAL b, REAL e, TArcSegment *s1, TArcSegment *s2, TPlaneVect *p);
  //for optimal alg only
  uint4 IntInside(REAL b, REAL e, TArcSegment *s1, TArcSegment *s2);

  struct TIntegerSegment
  {
    constexpr static uint4 is_line = 1;
#ifdef PRINT_SEG_AND_INT
    static TIntegerSegment* coll_begin;
#endif

    TIntegerSegment() :org(), shift() {};
    TIntegerSegment(TIntegerVect p1, TIntegerVect p2) :org(p1), shift(p2 - p1) {
      Refine();
    };
    TIntegerSegment(int4 xb, int4 yb, int4 xs, int4 ys) :org( xb,yb ), shift(xs,ys) {
      Refine();
    };

    TIntegerVect org;
    TIntegerVect shift;

    void Refine() noexcept
    {
      if ((shift.x < 0) || ((shift.x == 0) && (shift.y < 0)))
      {
        org = org + shift;
        shift = -shift;
      }
    };
    void Init(TIntegerSegment& s) { org = s.org; shift = s.shift; };
    void Init(const TLineSegment1& s, const transf1D& x_transf, const transf1D& y_transf) {
      org.from_real(s.org, x_transf, y_transf);
      shift= TIntegerVect(s.org + s.shift, x_transf, y_transf) - org;
      Refine();
      //assert((shift.x > 0)||((shift.x==0)&&(shift.y>0)));
#ifdef PRINT_SEG_AND_INT
      if (print_at_lineseg1_init)
        printf("s%i[%i,%i,%i,%i],\n", org.x, org.y, org.x + shift.x, org.y + shift.y);
#endif
    };

    void Init(const TLineSegment1& seg, double si, double co,const transf1D& x_transf, const transf1D& y_transf) {
      auto s = seg.get_rotated(si, co);
      org.from_real(s.org, x_transf, y_transf);
      shift = TIntegerVect(s.org + s.shift, x_transf, y_transf) - org;
      Refine();
      //assert((shift.x > 0)||((shift.x==0)&&(shift.y>0)));
#ifdef PRINT_SEG_AND_INT
      if (print_at_lineseg1_init) {
        uint4 sn = this - TIntegerSegment::coll_begin;
        //printf("is%i[%i,%i,%i,%i],\n", sn, org.x, org.y, org.x + shift.x, org.y + shift.y);
        printf("is%i[%i,%i,%i,%i],\n", sn, org.x, org.y, shift.x,shift.y);
      }
#endif
    };

    TIntegerVect BegPoint() const {
      return org;
    };
    TIntegerVect EndPoint() const {
      return org + shift;
    };
    bool operator==(const TIntegerSegment& s) const = default;

    auto PointX(bool is_last) const {
      return is_last ? org.x + shift.x : org.x;
    };
    auto bx() const {
      return org.x;
    }
    auto by() const {
      return org.y;
    }
    auto ex() const {
      return org.x+shift.x;
    }
    auto ey() const {
      return org.y+shift.y;
    }
    auto sx() const {
      return shift.x;
    }
    auto sy() const {
      return shift.y;
    }

    bool is_vertical() const {
      return sx() == 0;
    }

    auto point_pos(TIntegerVect v) const { //return negative if segment placed under point v
      return (v - org) % shift;
    };
    auto under(TIntegerVect v) const { //segment placed under point v, true == -1, false == 1
      return point_pos(v) <=> 0;
    };
    auto upper(TIntegerVect v) const {//segment placed above point v, true == -1, false == 1
      return 0 <=> point_pos(v);
    };
    static bool exact_true(std::strong_ordering o) {
      return o == std::strong_ordering::less;
    }
    static bool exact_false(std::strong_ordering o) {
      return o == std::strong_ordering::greater;
    }
    static bool no_true_no_false(std::strong_ordering o) {
      return o == std::strong_ordering::equal;
    }
    bool exact_under(TIntegerVect v) const { //segment placed under point v
      return point_pos(v) < 0;
    };
    bool exact_on(TIntegerVect v) const { //segment placed under point v
      return point_pos(v) == 0;
    };
    bool exact_upper(TIntegerVect v) const {
      return 0 < point_pos(v);
    };

    //function to check if segment's coordinates lies inside full_int_range
    bool is_inside_int_range() const {
      return (bx() >= -full_int_range) && (abs(by()) <= full_int_range) &&
        (ex() <= full_int_range) && (abs(ey()) <= full_int_range);
    }

    bool no_common_y(const TIntegerSegment& s2) const {
      auto [min1, max1] = std::pair<int4,int4>(std::minmax(by(), ey()));
      auto [min2, max2] = std::pair<int4, int4>(std::minmax(s2.by(), s2.ey()));
      return (min2 > max1) || (min1 > max2);
    }

    void write_SVG(int4 id, chostream* SVG_text,int4 type=0) const {
      if (SVG_text) {
        *SVG_text << "<line id='seg" << id;
        *SVG_text << "' x1='" << (double)bx() * crd_scaler;
        *SVG_text << "' y1='" << (double)by() * crd_scaler;
        *SVG_text << "' x2='" << (double)ex() * crd_scaler;
        *SVG_text << "' y2='" << (double)ey() * crd_scaler;
        if (type == 0)
          *SVG_text << "' class='int_line'/>\n";
        else {
          *SVG_text << "' class='int_line";
          *SVG_text << type;
          *SVG_text << "'/>\n";
        }
      }
    };

    auto get_int_type(const TIntegerSegment& s) const {
      assert(shift.is_non_zero() && s.shift.is_non_zero());
      uint4 res= _IntType::common_int;
      auto d = s.BegPoint() - BegPoint();
      if (d % shift == 0) res += _IntType::s2_beg_int;
      if (d % s.shift == 0) res += _IntType::s1_beg_int;
      d = s.EndPoint() - EndPoint();
      if (d % shift == 0) res += _IntType::s2_end_int;
      if (d % s.shift == 0) res += _IntType::s1_end_int;
      return res;
    }

    auto get_int_type_beg(const TIntegerSegment& s) const {
      assert(shift.is_non_zero() || s.shift.is_non_zero());

      auto d = s.BegPoint() - BegPoint();

      // `get_int_type_beg()` is a cheap "intersection-kind fingerprint" used to filter duplicates
      // when expanding remapped (normalized) segment pairs back into original segment pairs.
      //
      // For non-degenerate segments it relies on orientation tests `d % shift == 0` that detect
      // whether `s.BegPoint()` lies on the supporting line of `*this` and vice versa.
      //
      // Degenerate (zero-length) segments have `shift == {0,0}` and do not define a supporting line.
      // For them we treat the segment as a single point P = BegPoint() and only report a "beg-touch"
      // when the other segment's begin point is exactly the same point (d == 0). This makes the
      // duplicate-filter logic consistent between original and remapped segments for zero segments.
      if (shift.is_zero()) {
        return _IntType::common_int + _IntType::s1_beg_int +
          d.is_zero() * _IntType::s2_beg_int;
      }

      // Symmetric case: `s` is a point segment.
      if (s.shift.is_zero()) {
        return _IntType::common_int + _IntType::s2_beg_int +
          d.is_zero() * _IntType::s1_beg_int;
      }
      return _IntType::common_int + (d % shift == 0) * _IntType::s2_beg_int + (d % s.shift == 0) * _IntType::s1_beg_int;
    }

    friend std::ostream& operator<<(std::ostream& os, const TIntegerSegment& s) {
      os <<"{"<<s.bx() << "," << s.by() << "," << s.sx() << "," << s.sy() << "}"; //" " <<double(s.shift.y)/ s.shift.x<< std::endl;
      return os;
    };

    auto YAtX_Numerator(int4 X) const {
      assert(X >= org.x);
      return (int8)org.y * (int8)shift.x + ((int8)X - (int8)org.x) * (int8)shift.y;
    };

    frac64 YAtX_frac(int4 X) const {
      return { YAtX_Numerator(X),shift.x };
    };

  };

