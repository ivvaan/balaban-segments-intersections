#ifndef UTILS_FOR_SEGMENT_INTERSECTION
#define UTILS_FOR_SEGMENT_INTERSECTION
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
#include <ostream>

#include <math.h>
typedef int int4;
typedef unsigned int uint4;
typedef double REAL;
//typedef float SINGLE;
typedef int  BOOL;
#define TRUE 1
#define FALSE 0


#ifndef MAX
#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)            (((a) < (b)) ? (a) : (b))
#endif


//#define register
#ifndef NDEBUG
//#define PRINT_SEG_AND_INT
extern bool print_at_lineseg1_init;
#endif

const uint4 max_SVG_items = 20000;
const uint4 max_SVG_points = 150000;


inline REAL sq(REAL x) { return x * x; }

using chostream = std::basic_ostream<char>;

typedef void(*PSeg);
typedef void(*PRegObj);


enum _Algorithm
{
  triv = 1,
  simple_sweep = 2,
  fast = 4,
  optimal = 8,
  fast_parallel = 16,
  bentley_ottmann = 32,
  fast_no_ip = 64,
  mem_save = 128
};

enum _Segment
{
  line1 = 0, line2, arc, graph
};

enum _Distribution
{
  random = 0, parallel, mixed, small, param_defined, circul
};

enum _Registrator
{
  per_segm_reg_just_count_stat = 0,
  per_segm_reg_max_per_segm_stat = 1,
  just_count = 2,
  store_pairs_and_ints_just_count_stat = 3
};


class CRandomValueGen
{
public:
  CRandomValueGen();
  CRandomValueGen(unsigned seed);
  void SetSeeed(unsigned seed);
  double GetRandomDouble();
  //BOOL RandomChoose(double prop);
  bool RandomChoose();
private:
  double (*stduniform)();
};

template<class real>
class couple {
public:
  real x = 0;
  real y = 0;
  couple() {};
  couple(real xc, real yc) :x(xc), y(yc) {};
  couple(const couple& c) = default;
  real getX() const { return x; };
  real getY() const { return y; };
  bool operator<(const couple<real>& v2) const
  {
    return ((x < v2.x) || (x == v2.x) && (y < v2.y));
  };
  bool operator<=(const couple<real>& v2)  const
  {
    return ((x < v2.x) || (x == v2.x) && (y <= v2.y));
  };
  bool operator>(const couple<real>& v2)  const
  {
    return ((x > v2.x) || (x == v2.x) && (y > v2.y));
  };
  couple<real> operator-() const
  {
    return { -x, -y };
  };
  real get_norm() const
  {
    return x * x + y * y;
  }
  real get_length() const
  {
    return sqrt(x * x + y * y);
  };
  couple<real>& normalize()
  {
    real revl = 1.0 / get_length();
    x *= revl; y *= revl;
    return *this;
  };
  couple<real> get_normalized() const
  {
    real revl = 1.0 / get_length();
    return { x * revl, y * revl };
  };
  void InitRandom(CRandomValueGen& rv)
  {
    x = rv.GetRandomDouble();
    y = rv.GetRandomDouble();
  };

};



template <class real> inline real operator%(const couple<real>& v1, const couple<real>& v2)
{
  return v1.x * v2.y - v2.x * v1.y;
}

template <class real> inline couple<real> operator%(real   r, const couple<real>& t)
{
  return { -t.y * r, t.x * r };
}

template <class real> inline couple<real> operator%(const couple<real>& t, real r)
{
  return { t.y * r, -t.x * r };
}

template <class real> inline couple<real> operator*(real   r, const couple<real>& t)
{
  return { t.x * r, t.y * r };
}

template <class real> inline couple<real> operator*(const couple<real>& t, real r)
{
  return { t.x * r, t.y * r };
}

template <class real> inline real operator*(const couple<real>& t, const couple<real>& v)
{
  return t.x * v.x + t.y * v.y;
}

template <class real> inline couple<real> operator+(const couple<real>& v1, const couple<real>& v2)
{
  return { v1.x + v2.x,v1.y + v2.y };
}

template <class real> inline couple<real> operator-(const couple<real>& v1, const couple<real>& v2)
{
  return { v1.x - v2.x,v1.y - v2.y };
}

template<class out, class real>
out& operator <<(out& o, const couple<real>& v)
{
  //return o << "x=" << v.x << ",y=" << v.y;
  return o << "[" << v.x << "," << v.y << "]";
}



typedef couple<REAL> TPlaneVect;

struct minmaxrect { TPlaneVect ld, rt; };
template<class SegArr>
minmaxrect get_mmrect01(SegArr c[], int4 N) {
  REAL xmin = 0, ymin = 0, xmax = 1, ymax = 1;
  for (int4 i = 0; i < N; ++i) {
    TPlaneVect bp = c[i].BegPoint();
    TPlaneVect ep = c[i].EndPoint();
    xmin = MIN(xmin, bp.x);
    xmax = MAX(xmax, ep.x);
    ymin = MIN(ymin, MIN(bp.y, ep.y));
    ymax = MAX(ymax, MAX(bp.y, ep.y));
  }
  return { {xmin,ymin},{xmax,ymax} };
};


#endif
