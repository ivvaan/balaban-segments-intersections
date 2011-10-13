#ifndef UTILS_FOR_SEGMENT_INTERSECTION
#define UTILS_FOR_SEGMENT_INTERSECTION
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

#include <math.h>
typedef int int4;
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

#define SIZECOND(a,b) (a+8<b)

const int4 maxStack=31;
const int4 Har [10]={1,4,13,41,121,365,1097,3293,9889,29667};

template<class real> void ssort(real *List,int4 n)
  {
  register int4 i,j,h;
  int k;
  real temp;
  if (n==1) return;
  for(k=0;Har[k]<n;k++);
  while (k>0)
    {
    h=Har[--k];
    for (j=h;j<n;j++)
      {
      temp=List[i=j];
      while ((i>=h)&&(temp<List[i-h]))
        {
        List[i]=List[i-h];
        i-=h;
        }
      List[i]=temp;
      }
    }
  };

template<class order,class real> void ssort(order &o,real *List,int4 n)
  {
  register int4 i,j,h;
  int k;
  real temp;
  if (n==1) return;
  for(k=0;Har[k]<n;k++);
  while (k>0)
    {
    h=Har[--k];
    for (j=h;j<n;j++)
      {
      temp=List[i=j];
      while ((i>=h)&&(o.lt(temp,List[i-h])))
        {
        List[i]=List[i-h];
        i-=h;
        }
      List[i]=temp;
      }
    }
  };

template<class real> void fastsort(real *List,int4 n)
  {
  int4 LStack[maxStack],RStack[maxStack];
  int4 StackTop;

  int4 ll,lr,lm,nl,nr;
  real current,temp;

  if (n<2) return;
  ll = 0;
  lr = n-1;
  StackTop = 0 ;

  while(1)
    {
    while(SIZECOND(ll,lr)) 
      {
      nl = ll;
      nr = lr;
      lm = (ll+lr)>>1 ;
      current =List[lm];
      // find keys for exchange
      while(1)
        {
        while(List[nl]<current) nl++;
        while(current<List[nr]) nr--;
        if((nl+2)>nr) break ;
        temp=List[nl];
        List[nl]=List[nr];
        List[nr]=temp;
        nl++;
        nr--;
        };
      if( nl <= nr )
        {
        if( nl <  nr )
          {
          temp = List[nl] ;
          List[nl] = List[nr] ;
          List[nr] = temp ;
          }
        nl++;
        nr--;
        }
      // select sub-list to be processed next
      if( nr <  lm )
        {
        LStack[StackTop] = nl ;
        RStack[StackTop] = lr ;
        lr=nr;
        }
      else
        {
        LStack[StackTop] = ll ;
        RStack[StackTop] = nr ;
        ll = nl ;
        }
      StackTop++;
      }; //  while(ll<lr) 
    ssort(List+ll,lr-ll+1);

    if (StackTop == 0) return;
    StackTop-- ;
    ll = LStack[StackTop] ;
    lr = RStack[StackTop] ;
    }
  }

template<class order,class real> void fastsort(order &o,real *List,int4 n)
  {
  int4 LStack[maxStack],RStack[maxStack];
  int4 StackTop;

  int4 ll,lr,lm,nl,nr;
  real current,temp;

  if (n<2) return;
  ll = 0;
  lr = n-1;
  StackTop = 0 ;

  while(1)
    {
    while(SIZECOND(ll,lr)) 
      {
      nl = ll;
      nr = lr;
      lm = (ll+lr)>>1 ;
      current =List[lm];
      while(1)
        {
        while(o.lt(List[nl],current)) nl++;
        while(o.lt(current,List[nr])) nr--;
        if((nl+2)>nr) break ;
        temp=List[nl];
        List[nl]=List[nr];
        List[nr]=temp;
        nl++;
        nr--;
        };
      if( nl <= nr )
        {
        if( nl <  nr )
          {
          temp = List[nl] ;
          List[nl] = List[nr] ;
          List[nr] = temp ;
          }
        nl++;
        nr--;
        }
      if( nr <  lm )
        {
        LStack[StackTop] = nl ;
        RStack[StackTop] = lr ;
        lr=nr;
        }
      else
        {
        LStack[StackTop] = ll ;
        RStack[StackTop] = nr ;
        ll = nl ;
        }
      StackTop++;
      };  
    ssort(o,List+ll,lr-ll+1);

    if (StackTop == 0) return;
    StackTop-- ;
    ll = LStack[StackTop] ;
    lr = RStack[StackTop] ;
    }
  }


template<class real> 
class tuple{
public:
  union
    {
    real arr[1];
    real x;};
    real y;
    tuple()
      { arr[0]=0;
    arr[1]=0;}
    tuple(const real &xc,const real &yc):x(xc),y(yc){};
    int operator<(tuple<real>   & v2)
      { return ((x<v2.x)||(x==v2.x) && (y<v2.y));};
    int operator>(tuple<real>   & v2)
      { return ((x>v2.x)||(x==v2.x) && (y>v2.y));};
    tuple<real> operator-()
      {return tuple<real>(-x,-y);};
    real get_norm()
      {return x*x+y*y;}
    real get_length()
      {return sqrt(x*x+y*y);};
    tuple<real>& normalize()
      {
      real l=get_length();
      x/=l;y/=l;
      return *this;
      };
    tuple<real> get_normalized()
      {
      real l=get_length();
      return tuple<real>(x/l,y/l);;
      };
  };



template <class real> inline real operator%(tuple<real>    v1,tuple<real>   v2)
  { return v1.x*v2.y-v2.x*v1.y;}

template <class real> inline tuple<real> operator*(const real   & r,tuple<real>   & t)
  {return tuple<real>(t.x*r,t.y*r);}

template <class real> inline tuple<real> operator*(tuple<real>   & t,real   & r)
  {return tuple<real>(t.x*r,t.y*r);}

template <class real> inline real operator*(tuple<real>  t,tuple<real>   v)
  {return t.x*v.x+t.y*v.y;}

template <class real> inline tuple<real> operator+(tuple<real>   &v1,tuple<real>   &v2)
  {return tuple<real>(v1.x+v2.x,v1.y+v2.y);}

template <class real> inline tuple<real> operator-(tuple<real>   &v1,tuple<real>   &v2)
  {return tuple<real>(v1.x-v2.x,v1.y-v2.y);}


class CRandomValueGen
  {
  public:
    double GetRandomDouble();
    //BOOL RandomChoose(double prop);
    BOOL RandomChoose();
  };


typedef tuple<REAL> TPlaneVect;

#endif
