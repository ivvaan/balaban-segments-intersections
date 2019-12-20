#ifndef UTILS_FOR_SEGMENT_INTERSECTION
#define UTILS_FOR_SEGMENT_INTERSECTION
/*
*
*      Copyright (c)  2011-2017  Ivan Balaban 
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

#include <math.h>
typedef int int4;
typedef unsigned int uint4;
typedef double REAL;
//typedef float SINGLE;
typedef int  BOOL;
#define TRUE 1
#define FALSE 0

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif
inline REAL sq(REAL x) { return x*x; }

template<class real> 
class couple{
public:
    real x=0;
    real y=0;
    couple() :x(0), y(0) {};
    couple(const real &xc,const real &yc):x(xc),y(yc){};
    int operator<(couple<real>   & v2)
      { return ((x<v2.x)||(x==v2.x) && (y<v2.y));};
    int operator>(couple<real>   & v2)
      { return ((x>v2.x)||(x==v2.x) && (y>v2.y));};
    couple<real> operator-()
      {return couple<real>(-x,-y);};
    real get_norm()
      {return x*x+y*y;}
    real get_length()
      {return sqrt(x*x+y*y);};
    couple<real>& normalize()
      {
      real l=get_length();
      x/=l;y/=l;
      return *this;
      };
    couple<real> get_normalized()
      {
      real l=get_length();
      return couple<real>(x/l,y/l);
      };
  };



template <class real> inline real operator%(const couple<real>    &v1,const couple<real>  &v2)
  { return v1.x*v2.y-v2.x*v1.y;}

template <class real> inline couple<real> operator*(const real   & r,const couple<real>   & t)
  {return couple<real>(t.x*r,t.y*r);}

template <class real> inline couple<real> operator*(const couple<real>   & t,const real   & r)
  {return couple<real>(t.x*r,t.y*r);}

template <class real> inline real operator*(const couple<real>   &t, const couple<real>   &   v)
  {return t.x*v.x+t.y*v.y;}

template <class real> inline couple<real> operator+(const couple<real>   &v1, const couple<real>   &v2)
  {return couple<real>(v1.x+v2.x,v1.y+v2.y);}

template <class real> inline couple<real> operator-(const couple<real>   &v1, const couple<real>   &v2)
  {return couple<real>(v1.x-v2.x,v1.y-v2.y);}


class CRandomValueGen
  {
  public:
    double GetRandomDouble();
    //BOOL RandomChoose(double prop);
    bool RandomChoose();
  };


typedef couple<REAL> TPlaneVect;

#endif
