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
#include <cassert>
#include <type_traits>
#include <math.h>
typedef int int4;
typedef unsigned int uint4;
typedef long long int8;
typedef unsigned long long uint8;
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
#define PRINT_SEG_AND_INT
extern bool print_at_lineseg1_init;
#endif

#define DECL_RAII_ARR(L,len)\
auto L##_ = std::make_unique<std::remove_pointer_t<decltype(L)>[]>(len);\
L = L##_.get();

#define MY_FREE_ARR_MACRO(a) if (a!=nullptr) {delete[] a; a = nullptr;}

constexpr uint4 max_SVG_items = 5000;
constexpr uint4 max_SVG_points = 150000;
constexpr uint4 max_truereg_items = 20000;


template< class REAL> REAL sq(REAL x) { return x * x; }


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

constexpr int4 alg_list[] = { triv, simple_sweep, fast, optimal, fast_parallel, bentley_ottmann,fast_no_ip,mem_save };
constexpr const char* alg_names[] = { "trivial","simple_sweep","fast","optimal","fast_parallel","bentley_ottmann","fast no inters points","fast 'no R'" };

enum _Segment
{
  line1 = 0, line2, arc, graph, intline
};

constexpr uint4 full_int_range = 0x40000000;
constexpr uint4 degenerate_range = 0x4;

enum _Distribution
{
  random = 0, parallel, mixed, small, param_defined, circul, degenerate
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
  auto getX() const { return x; };
  auto getY() const { return y; };
  bool operator<(const couple<real>& v2) const {
    return ((x < v2.x) || (x == v2.x) && (y < v2.y));
  };
  bool operator<=(const couple<real>& v2)  const {
    return ((x < v2.x) || (x == v2.x) && (y <= v2.y));
  };
  bool operator>(const couple<real>& v2)  const {
    return ((x > v2.x) || (x == v2.x) && (y > v2.y));
  };
  couple<real> operator-() const {
    return { -x, -y };
  };
  auto get_norm() const {
    return x * x + y * y;
  }
  auto get_length() const {
    return sqrt(x * x + y * y);
  };
  couple<real>& normalize() {
    real revl = 1.0 / get_length();
    x *= revl; y *= revl;
    return *this;
  };
  couple<real> get_normalized() const {
    real revl = 1.0 / get_length();
    return { x * revl, y * revl };
  };

  void InitRandom(CRandomValueGen& rv)
  {
    x = rv.GetRandomDouble();
    y = rv.GetRandomDouble();
  };

};

template <class real>  auto operator%(const couple<real>& v1, const couple<real>& v2) {
  return v1.x * v2.y - v2.x * v1.y;
}

template <class real>  couple<real> operator%(real   r, const couple<real>& t) {
  return { -t.y * r, t.x * r };
}

template <class real>  couple<real> operator%(const couple<real>& t, real r) {
  return { t.y * r, -t.x * r };
}

template <class real>  couple<real> operator*(real   r, const couple<real>& t) {
  return { t.x * r, t.y * r };
}

template <class real>  couple<real> operator*(const couple<real>& t, real r) {
  return { t.x * r, t.y * r };
}

template <class real>  auto operator*(const couple<real>& t, const couple<real>& v) {
  return t.x * v.x + t.y * v.y;
}

template <class real>  couple<real> operator+(const couple<real>& v1, const couple<real>& v2) {
  return { v1.x + v2.x,v1.y + v2.y };
}

template <class real>  couple<real> operator-(const couple<real>& v1, const couple<real>& v2) {
  return { v1.x - v2.x,v1.y - v2.y };
}

template<class out, class real>
out& operator <<(out& o, const couple<real>& v) {
  //return o << "x=" << v.x << ",y=" << v.y;
  return o << "[" << v.x << "," << v.y << "]";
}


struct TIntegerVect {
  int4 x = 0;
  int4 y = 0;
  TIntegerVect() {};
  TIntegerVect(int4 xc, int4 yc) :x(xc), y(yc) {};
  template <class real>
  TIntegerVect(const couple<real>& v, const couple<real>& x_transf, const couple<real>& y_transf) :
    x((v.x + x_transf.y) * x_transf.x + 0.5),  y((v.y + y_transf.y) * y_transf.x + 0.5){};
  TIntegerVect(const TIntegerVect& c) = default;
  auto getX() const { return x; };
  auto getY() const { return y; };
  bool operator<(const TIntegerVect& v2) const {
    return ((x < v2.x) || (x == v2.x) && (y < v2.y));
  };
  bool operator<=(const TIntegerVect& v2)  const {
    return ((x < v2.x) || (x == v2.x) && (y <= v2.y));
  };
  bool operator>(const TIntegerVect& v2)  const {
    return ((x > v2.x) || (x == v2.x) && (y > v2.y));
  };
  TIntegerVect operator-() const {
    return { -x, -y };
  };
  auto get_norm() const {
    return ((int8)x) * ((int8)x) + ((int8)y) * ((int8)y);
  }
  void InitRandom(CRandomValueGen& rv)
  {
    x = rv.GetRandomDouble() * INT_MAX;
    y = rv.GetRandomDouble() * INT_MAX;
  };

  template <class real>
  void from_real(const couple<real>& v, const couple<real>& x_transf, const couple<real>& y_transf) {
    x = (v.x + x_transf.y) * x_transf.x + 0.5;
    y = (v.y + y_transf.y) * y_transf.x + 0.5;
  };

  template <class real>
  void from_real(const couple<real>& v, real x_scale, real y_scale) {
    x = v.x * x_scale + 0.5;
    y = v.y * y_scale + 0.5;
  };

};

int8 operator%(const TIntegerVect& v1, const TIntegerVect& v2);

int8 operator*(const TIntegerVect& v1, const TIntegerVect& v2);

TIntegerVect operator+(const TIntegerVect& v1, const TIntegerVect& v2);

TIntegerVect operator-(const TIntegerVect& v1, const TIntegerVect& v2);

TIntegerVect operator*(int4 r, const TIntegerVect& t);

TIntegerVect operator*(const TIntegerVect& t, int4 r);


typedef couple<REAL> TPlaneVect;

struct minmaxrect { 
  TPlaneVect ld, rt; 
  minmaxrect get_scaled(REAL factor) {
    auto fp = 0.5 * (1. + factor);
    auto fm = 0.5 * (1. - factor);
    return { fp * ld + fm * rt , fp * rt + fm * ld };
  };
};
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

namespace {
  template<typename SegmentsColl>
  struct has_get_sentinel
  {
    template<typename U, int4(U::*)(bool)> struct SFINAE {};
    template<typename U> static char test(SFINAE<U, &U::get_sentinel>*);
    template<typename U> static int test(...);
    constexpr static bool value = sizeof(test<SegmentsColl>(nullptr)) == sizeof(char);
  };
}

template< class T >
constexpr bool has_sentinels = has_get_sentinel<T>::value;

#define THIS_HAS_SENTINELS (has_sentinels<std::remove_pointer_t<decltype(this)> >)

#endif
