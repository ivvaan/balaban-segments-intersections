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
#include <intrin.h>
#include <algorithm>
#include <vector>
//#define DEBUG_INTERSECTION_SET

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
//#define PRINT_SEG_AND_INT
extern bool print_at_lineseg1_init;
#endif

#define DECL_RAII_ARR(L,len)\
auto __##L##__ = std::make_unique<std::remove_pointer_t<decltype(L)>[]>(len);\
L = __##L##__.get();

#define RAII_ARR(type,arr,len)\
auto __##arr##__ = std::make_unique<type[]>(len);\
type *arr = __##arr##__.get();

#define MY_FREE_ARR_MACRO(a) if (a!=nullptr) {delete[] a; a = nullptr;}

constexpr uint4 max_SVG_items = 5000;
constexpr uint4 max_SVG_points = 150000;
constexpr uint4 max_truereg_items = 20000;


template< class REAL> REAL sq(REAL x) { return x * x; }

using chostream = std::basic_ostream<char>;

typedef void(*PSeg);
typedef void(*PRegObj);


enum _Algorithm{
  triv = 1,
  simple_sweep = 2,
  fast = 4,
  optimal = 8,
  fast_parallel = 16
};

enum _IntType {
  common_int = 0,
  s1_beg_int = 2,
  s2_beg_int = 4,
  s1_end_int = 8,
  s2_end_int = 16,
};


constexpr int4 alg_list[] = { triv, simple_sweep, fast, optimal, fast_parallel };
constexpr const char* alg_names[] = { "trivial","simple_sweep","fast","optimal","fast_parallel" };

enum _Segment
{
  line1 = 0, line2, arc, graph, intline
};

constexpr int4 range_margin = 16;
// sum of 2 coords should not exceed int4 max
constexpr int4 full_int_range = std::numeric_limits<int4>::max() / 2 - range_margin;

enum _Distribution
{
  random = 0, parallel, mixed, small, param_defined, circle
};

enum _Registrator
{
  per_segm_reg_just_count_stat = 0,
  per_segm_reg_max_per_segm_stat = 1,
  just_count = 2,
  store_pairs_and_ints_just_count_stat = 3,
  just_count_require_intersections = 4
};

enum class _Coll_flag_state
{
  state_false = 0,
  state_true = 1,
  state_unimplemented = 2
};

enum class _Coll_flags
{
  line_segments,
  needs_SetCurSegCutBE_at_start
};


class CRandomValueGen
{
public:
  CRandomValueGen();
  CRandomValueGen(unsigned seed);
  void SetSeed(unsigned seed);
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
  couple<real>& operator/=(REAL r) {
    x /= r;
    y /= r;
    return *this;
  };

  void rotate_prod(REAL r) {
    auto t = -y * r;
    y = x * r;
    x = t;
  }

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
out& operator <<(out& o, const couple<real>& v)
{
  //return o << "x=" << v.x << ",y=" << v.y;
  return o << "[" << v.x << "," << v.y << "]";
}

typedef couple<REAL> TPlaneVect;

struct transf1D {
  REAL shift = 0.0;
  REAL scale = 1.0;
  REAL operator()(REAL v) const {
    return (v + shift) * scale;
  };
  REAL doscale(REAL v) const { return v * scale; };
};


struct TIntegerVect {
  int4 x = 0;
  int4 y = 0;
  TIntegerVect() {};
  constexpr TIntegerVect(int4 xc, int4 yc) :x(xc), y(yc) {};
  template <class real>
  TIntegerVect(const couple<real>& v, const transf1D& x_transf, const transf1D& y_transf) :
    x(x_transf(v.x) + 0.5),  
    y(y_transf(v.y) + 0.5)
  {};
  TIntegerVect(const TIntegerVect & c) = default;
  auto getX() const { return x; };
  auto getY() const { return y; };
  bool operator<(TIntegerVect v2) const {
    return ((x < v2.x) || (x == v2.x) && (y < v2.y));
  };
  bool operator!=(TIntegerVect v2) const {
   // *reinterpret_cast<int8 *>((void*)this)!= *reinterpret_cast<int8*>((void*)&v2)
    return (x != v2.x) || (y != v2.y);
  };

  bool operator==(TIntegerVect v2) const {
    // *reinterpret_cast<int8 *>((void*)this)== *reinterpret_cast<int8*>((void*)&v2)
    return (x == v2.x) && (y == v2.y);
  };

  bool is_zero() const {
    return (x == 0) && (y == 0);
  }

  bool is_non_zero() const {
    return !is_zero();
  }

  bool operator<=(TIntegerVect v2)  const {
    return ((x < v2.x) || (x == v2.x) && (y <= v2.y));
  };
  bool operator>(TIntegerVect v2)  const {
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
  void from_real(const couple<real>& v, const transf1D& x_transf, const transf1D& y_transf) {
    x = x_transf(v.x) + 0.5;
    y = y_transf(v.y) + 0.5;
  };

  template <class real>
  void from_real(const couple<real>& v, real x_scale, real y_scale) {
    x = v.x * x_scale + 0.5;
    y = v.y * y_scale + 0.5;
  };
  friend std::ostream& operator<<(std::ostream& os, TIntegerVect v) {
    os << "(" << v.x << "," << v.y<<")"; //" " <<double(s.shift.y)/ s.shift.x<< std::endl;
    return os;
  };
};

inline constexpr TIntegerVect int_max_vect{ std::numeric_limits<decltype(TIntegerVect::x)>::max(),std::numeric_limits<decltype(TIntegerVect::y)>::max() };
inline constexpr TIntegerVect int_zero_vect{ 0,0 };

int8 operator%(TIntegerVect v1, TIntegerVect v2);

int8 operator*(TIntegerVect v1, TIntegerVect v2);

TIntegerVect operator+(TIntegerVect v1, TIntegerVect v2);

TIntegerVect operator-(TIntegerVect v1, TIntegerVect v2);

TIntegerVect operator*(int4 r, TIntegerVect t);

TIntegerVect operator*(TIntegerVect t, int4 r);




struct minmaxrect { 
  TPlaneVect ld, rt; // left down, right top coners
  auto get_diag() const {
    return rt - ld;
  }
  auto get_width() const {
    return rt.x - ld.x;
  }
  auto get_height() const {
    return rt.x - ld.x;
  }

  minmaxrect get_scaled(REAL factor) {// enlage rectangle into factor times 
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

template <class T>
minmaxrect get_minmax(int4 n, T c[])
{
  TPlaneVect bp = c[0].BegPoint();
  TPlaneVect ep = c[0].EndPoint();
  REAL xmin = bp.x;
  REAL xmax = ep.x;
  REAL ymin = MIN(bp.y, ep.y);
  REAL ymax = MAX(bp.y, ep.y);
  for (int4 i = 1; i < n; ++i) {
    bp = c[i].BegPoint();
    ep = c[i].EndPoint();
    xmin = MIN(xmin, bp.x);
    xmax = MAX(xmax, ep.x);
    ymin = std::min({ ymin,bp.y,ep.y });
    ymax = std::max({ ymax,bp.y,ep.y });
  }

  return { {xmin,ymin},{xmax,ymax} };
};

template <class T>
minmaxrect get_rot_minmax(int4 n, T coll[], double s, double c)
{
  //auto s = sin(angle), c = cos(angle);
  auto r = coll[0].get_rotated(s, c);
  TPlaneVect bp = r.BegPoint();
  TPlaneVect ep = r.EndPoint();
  REAL xmin = bp.x;
  REAL xmax = ep.x;
  REAL ymin = MIN(bp.y, ep.y);
  REAL ymax = MAX(bp.y, ep.y);
  for (int4 i = 1; i < n; ++i) {
    r = coll[i].get_rotated(s, c);
    bp = r.BegPoint();
    ep = r.EndPoint();
    xmin = MIN(xmin, bp.x);
    xmax = MAX(xmax, ep.x);
    ymin = std::min({ ymin,bp.y,ep.y });
    ymax = std::max({ ymax,bp.y,ep.y });
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

std::strong_ordering comp_prod(int64_t a, int64_t b, int64_t c, int64_t d);

struct prod64 {// represents product a*b
  int64_t a, b;

  auto operator<=>(const prod64& other) const {
    return comp_prod(a, b, other.a, other.b);
  };
};

struct frac64 {// represents proper fraction num/denum (denum>0)
  int64_t num, denum;

  auto operator<=>(const frac64& other) const {
    assert((denum > 0) == (other.denum > 0));
    return comp_prod(num, other.denum, other.num, denum);
  };
};

inline frac64 div64(int64_t n, int64_t d) {// represents division n/d convert it
      //to proper fraction n'/d', where d'>0
  if (d < 0)
    return { -n, -d };
  return { n, d };
}

template <typename Itterator, typename Compare>
void sort_by_spaceship(Itterator f, Itterator l, Compare cmp) {
  std::sort(f, l, [=](auto&& a, auto&& b) -> bool {
    return cmp(std::forward<decltype(a)>(a), std::forward<decltype(b)>(b)) < 0;
    }
  );
};

template <typename T>
void removeElement(std::vector<T>& vec, const T& value_to_remove) {
  vec.erase(std::remove(vec.begin(), vec.end(), value_to_remove), vec.end());
}

template <typename T>
void remove_by_val(std::vector<T>& vec, T value_to_remove) {
  vec.erase(std::remove(vec.begin(), vec.end(), value_to_remove), vec.end());
};

template <typename T>
T* set_size(std::vector<T>& v, size_t new_size) {
  v.reserve(new_size);
  v.resize(new_size);
  return v.data();
};

struct PrepareResult {
  uint4 ends_len = 0;
  uint4 max_segm_on_vline = 0;
  double avr_segm_on_vline = 0.0;
};


#endif
