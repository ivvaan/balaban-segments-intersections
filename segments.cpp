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
#include "segments.h"

void TLineSegment1::InitRandom(CRandomValueGen &rv, int4 seg_n, int4 type, REAL par)
{
  if (circul == type)
  {
    //Segments ends are on the opposite sides of unit circul, each segment intesect each
    double angle = 1.0 / 16.0 + 3.*M_PI*rv.GetRandomDouble() / 8.0;
    org.x = 1. - cos(angle);
    org.y = 1. - sin(angle);
    shift.x = 1. + cos(angle) - org.x;
    shift.y = 1. + sin(angle) - org.y;
    Refine();

#ifdef PRINT_SEG_AND_INT
    if(print_at_lineseg1_init)
      printf("[%6.3f,%6.3f,%6.3f,%6.3f],\n", org.x, org.y, org.x + shift.x, org.y + shift.y);
#endif
    return;
  }
  org.x = rv.GetRandomDouble();
  org.y = rv.GetRandomDouble();
  shift.x = rv.GetRandomDouble() - org.x;
  shift.y = rv.GetRandomDouble() - org.y;
  Refine();
  shift.x += 1e-6; //to be sure the segment non vertical
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
  if(print_at_lineseg1_init)
    printf("[%6.3f,%6.3f,%6.3f,%6.3f],\n", org.x, org.y, org.x + shift.x, org.y + shift.y);
#endif
};

bool below(REAL X, TLineSegment1* s1, TLineSegment1* s2)
{
  if (s1 == s2) return false;
  REAL dx1 = s1->shift.x, dx2 = s2->shift.x;
  REAL dy1 = s1->shift.y, dy2 = s2->shift.y;
  REAL y1, y2;
  if ((dx1 == 0.0) || (s1->org.x + dx1 == X))
  {
    REAL res = (s1->org + s1->shift - s2->org) % s2->shift;
    //if (fabs(res) == 0.0) throw 1;
    return (res>0);
  }
  if (s2->org.x + dx2 == X)
  {
    REAL res = (s2->org + s2->shift - s1->org) % s1->shift;
    //if (fabs(res) == 0.0) throw 1;
    return (res<0);
  }
  y1 = (s1->org.y*dx1 + (X - s1->org.x)*dy1)*dx2;
  y2 = (s2->org.y*dx2 + (X - s2->org.x)*dy2)*dx1;
  return y1<y2;
}

template <bool _ret_ip>
int4 IntPoint(TLineSegment1* s1, TLineSegment1* s2, TPlaneVect *p)
{
  REAL b1 = max(s1->org.x, s2->org.x);
  REAL b2 = min(s1->org.x + s1->shift.x, s2->org.x + s2->shift.x);
  if (b1 >= b2)return 0;
  TPlaneVect delt = s2->org - s1->org;
  REAL prod = s1->shift%s2->shift, mul;
  if (((mul = s1->shift%delt)>0) ^ (mul + prod>0))
    if (((mul = delt%s2->shift)>0) ^ (mul - prod>0))
    {
      if (prod == 0) return 0;
      if (_ret_ip)*p = s1->org + ((REAL)fabs(mul / prod))*s1->shift;
      return 1;
    }
  return 0;
}

template int4 IntPoint<true>(TLineSegment1* s1, TLineSegment1* s2, TPlaneVect *p);
template int4 IntPoint<false>(TLineSegment1* s1, TLineSegment1* s2, TPlaneVect *p);

template <bool _ret_ip>
int4 StripePoint(REAL b, REAL e, TLineSegment1* s1, TLineSegment1* s2, TPlaneVect *p)
{
  TPlaneVect delt = s2->org - s1->org;
  REAL prod = s1->shift%s2->shift, mul, xc;
  if (((mul = s1->shift%delt)>0) ^ (mul + prod>0))
  {
    if (prod == 0) return 0;
    mul = mul / prod;
    xc = s2->org.x - mul*s2->shift.x;
    if ((xc <= b) || (xc>e)) return 0;
    if (_ret_ip)
    {
      p->x = xc;
      p->y = s2->org.y - mul*s2->shift.y;
    }
    return 1;
  }
  return 0;
}

template int4 StripePoint<true>(REAL b, REAL e, TLineSegment1* s1, TLineSegment1* s2, TPlaneVect *p);
template int4 StripePoint<false>(REAL b, REAL e, TLineSegment1* s1, TLineSegment1* s2, TPlaneVect *p);

int4 IntInside(REAL b, REAL e, TLineSegment1* s1, TLineSegment1* s2)
{
  TPlaneVect delt = s2->org - s1->org;
  REAL prod = s1->shift%s2->shift, mul = s1->shift%delt, xc;
  if ((mul>0) ^ (mul + prod>0))
  {
    if (prod == 0) return 0;
    mul = mul / prod;
    xc = s2->org.x - mul*s2->shift.x;
    if ((xc <= b) || (xc>e)) return 0;
    return 1;
  }
  return 0;
}



bool below(REAL X, TLineSegment2* s1, TLineSegment2* s2)
{

  return (s1 != s2) && ((X*s1->a + s1->b)<(X*s2->a + s2->b));
}

template <bool _ret_ip>
int4 IntPoint(TLineSegment2* s1, TLineSegment2* s2, TPlaneVect *p)
{
  REAL x1 = max(s1->x1, s2->x1);
  REAL x2 = min(s1->x2, s2->x2);
  if (x1 >= x2)return 0;
  REAL da = s1->a - s2->a;
  if (da == 0)return 0;
  if (_ret_ip)
  {
    p->x = (s2->b - s1->b) / da;
    if ((p->x >= x1) && (p->x <= x2)) { p->y = p->x*s1->a + s1->b; return 1; }
  }
  else
  {
    REAL x = (s2->b - s1->b) / da;
    return ((x >= x1) && (x <= x2));
  }
  return 0;
}

template int4 IntPoint<true>(TLineSegment2* s1, TLineSegment2* s2, TPlaneVect *p);
template int4 IntPoint<false>(TLineSegment2* s1, TLineSegment2* s2, TPlaneVect *p);

template <bool _ret_ip>
int4 StripePoint(REAL b, REAL e, TLineSegment2* s1, TLineSegment2* s2, TPlaneVect *p)
{
  REAL x1 = max(s1->x1, s2->x1);
  REAL x2 = min(s1->x2, s2->x2);
  x1 = max(x1, b); x2 = min(x2, e);
  REAL da = s1->a - s2->a;
  if (da == 0)return 0;
  if (_ret_ip)
  {
    p->x = (s2->b - s1->b) / da;
    if ((p->x >= x1) && (p->x <= x2)) { p->y = p->x*s1->a + s1->b; return 1; }
  }
  else
  {
    REAL x = (s2->b - s1->b) / da;
    return ((x >= x1) && (x <= x2));
  }
  return 0;
}


template int4 StripePoint<true>(REAL b, REAL e, TLineSegment2* s1, TLineSegment2* s2, TPlaneVect *p);
template int4 StripePoint<false>(REAL b, REAL e, TLineSegment2* s1, TLineSegment2* s2, TPlaneVect *p);

int4 IntInside(REAL b, REAL e, TLineSegment2 *s1, TLineSegment2 *s2)
{
  REAL x1 = max(s1->x1, s2->x1);
  REAL x2 = min(s1->x2, s2->x2);
  x1 = max(x1, b); x2 = min(x2, e);
  REAL da = s1->a - s2->a;
  if (da == 0)return 0;
  REAL x = (s2->b - s1->b) / da;
  if ((x >= x1) && (x <= x2)) return 1;
  return 0;
}




bool below(REAL X, TArcSegment* s1, TArcSegment* s2)
{
  return (s1 != s2) && (s1->under(s2->PointAtX(X)));
}


int4 IntPointsInStripe(REAL x1, REAL x2, TArcSegment *s1, TArcSegment *s2, TPlaneVect *p)
{
  TPlaneVect oo = s2->org - s1->org;
  REAL org_dist2 = oo.get_norm();
  REAL delta = 0.5*(s1->r2 - s2->r2 + org_dist2);
  REAL t = delta / org_dist2;
  REAL l2 = s1->r2 - delta*t;
  if (l2<0)return 0;
  TPlaneVect m = s1->org + t*oo;
  if (l2 == 0)
  {
    if ((m.x >= x1) && (m.x <= x2) && (s1->IsTheSamePart(m)) && (s2->IsTheSamePart(m)))
    {
      *p = m; return 1;
    }
    return 0;
  }
  t = sqrt(l2 / org_dist2);
  oo = TPlaneVect(oo.y*t, -oo.x*t);
  int4 npoints = 0;
  TPlaneVect res = m + oo;
  if ((res.x >= x1) && (res.x <= x2) && (s1->IsTheSamePart(res)) && (s2->IsTheSamePart(res)))
  {
    if (p)p[npoints] = res; npoints++;
  }
  res = m - oo;
  if ((res.x >= x1) && (res.x <= x2) && (s1->IsTheSamePart(res)) && (s2->IsTheSamePart(res)))
  {
    if (p)p[npoints] = res; npoints++;
  }
  return npoints;
};

template <bool _ret_ip>
int4 IntPoint(TArcSegment* s1, TArcSegment* s2, TPlaneVect *p)
{
  REAL x1 = max(s1->x1, s2->x1);
  REAL x2 = min(s1->x2, s2->x2);
  if (x1 >= x2)return 0;
  return IntPointsInStripe(x1, x2, s1, s2, p);
}

template int4 IntPoint<true>(TArcSegment* s1, TArcSegment* s2, TPlaneVect *p);
template int4 IntPoint<false>(TArcSegment* s1, TArcSegment* s2, TPlaneVect *p);

template <bool _ret_ip>
int4 StripePoint(REAL b, REAL e, TArcSegment* s1, TArcSegment* s2, TPlaneVect *p)
{
  REAL x1 = max(s1->x1, s2->x1);
  REAL x2 = min(s1->x2, s2->x2);
  x1 = max(x1, b); x2 = min(x2, e);
  return IntPointsInStripe(x1, x2, s1, s2, p);
}
template int4 StripePoint<true>(REAL b, REAL e, TArcSegment* s1, TArcSegment* s2, TPlaneVect *p);
template int4 StripePoint<false>(REAL b, REAL e, TArcSegment* s1, TArcSegment* s2, TPlaneVect *p);

int4 IntInside(REAL b, REAL e, TArcSegment* s1, TArcSegment* s2)
{
  REAL x1 = max(s1->x1, s2->x1);
  REAL x2 = min(s1->x2, s2->x2);
  x1 = max(x1, b); x2 = min(x2, e);
  TPlaneVect p[2];
  return IntPointsInStripe(x1, x2, s1, s2, p);
}

void TArcSegment::InitRandom(CRandomValueGen &rv, int4 seg_n, int4 type, REAL par)
{
  if (circul == type)
  {
    //Segments ends are on the opposite sides of unit circul, each segment intesect each
    double angle = 1.0 / 16.0 + 3. * M_PI * rv.GetRandomDouble() / 8.0;
    org.x = cos(angle);
    org.y = sin(angle);
    is_upper = false;
    r2 = 1.0;
    x1 = org.x - 1.+ 1. / 1024.;
    x2 = org.x + 1.- 1. / 1024.;

#ifdef PRINT_SEG_AND_INT
    printf("[%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,0,0,%i],\n", org.x, org.y, r2, x1, x2, (int)is_upper);
#endif
    return;
  }
#ifdef PRINT_SEG_AND_INT
  print_at_lineseg1_init = false;
#endif

  
  TLineSegment1 seg;
  do
  seg.InitRandom(rv, 0, type, par);
  while (seg.shift.x<0.001); //excluding too vertical case they have big radius
  auto slope = seg.shift.y / seg.shift.x;
  TPlaneVect norm(0.5*slope, -0.5);
  auto coef = fabs(seg.shift.x)*(1.+rv.GetRandomDouble())+fabs(seg.shift.y)+1./1024.;
  org = seg.org + 0.5 * seg.shift + (rv.RandomChoose() ? coef * norm : -coef * norm);
  r2 = (org - seg.org).get_norm();
  is_upper = seg.org.y > org.y;
  x1 = seg.org.x;
  x2 = x1 + seg.shift.x;
#ifdef PRINT_SEG_AND_INT
  printf("[%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%i],\n", org.x, org.y, r2, x1, x2, seg.org.y, seg.org.y + seg.shift.y,(int)is_upper);
#endif
}

