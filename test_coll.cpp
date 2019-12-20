/*																 fou
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
#include "test_coll.h"
#include <stdio.h>
#include "new_int_finder.h"
#include "fast_finder.h"
#include "mem_save_finder.h"
#include "optimal_finder.h"




//PSeg first_segment_ptr=NULL;



bool below(REAL X,TLineSegment1* s1,TLineSegment1* s2)
  { 
  if(s1==s2) return 0;
  REAL dx1=s1->shift.x,dx2=s2->shift.x;
  REAL dy1=s1->shift.y,dy2=s2->shift.y;
  REAL y1,y2;
  if ((dx1==0.0)||(s1->org.x+dx1==X))
    {
    REAL res=(s1->org+s1->shift-s2->org)%s2->shift;
    if (fabs(res)==0.0) throw 1;
    return (res>0);
    }
  if (s2->org.x+dx2==X)
    {
    REAL res=(s2->org+s2->shift-s1->org)%s1->shift;
    if (fabs(res)==0.0) throw 1;
    return (res<0);
    }
  y1=(s1->org.y*dx1+(X-s1->org.x)*dy1)*dx2;
  y2=(s2->org.y*dx2+(X-s2->org.x)*dy2)*dx1;
  return y1<y2;
  }

template <bool _ret_ip>
int4 IntPoint(TLineSegment1* s1,TLineSegment1* s2,TPlaneVect *p)
  { 
  REAL b1=max(s1->org.x,s2->org.x);
  REAL b2=min(s1->org.x+s1->shift.x,s2->org.x+s2->shift.x);
  if(b1>=b2)return 0;
  TPlaneVect delt=s2->org-s1->org;
  REAL prod=s1->shift%s2->shift,mul;
  if (((mul=s1->shift%delt)>0)^(mul+prod>0))
    if (((mul=delt%s2->shift)>0)^(mul-prod>0))
      {
      if (prod==0) return 0;
      if(_ret_ip)*p=s1->org+((REAL)fabs(mul/prod))*s1->shift;
      return 1;
      }
    return 0;
  } 

template <bool _ret_ip>
int4 StripePoint(REAL b,REAL e,TLineSegment1* s1,TLineSegment1* s2,TPlaneVect *p)
  { 
  TPlaneVect delt=s2->org-s1->org;
  REAL prod=s1->shift%s2->shift,mul,xc;
  if (((mul=s1->shift%delt)>0)^(mul+prod>0))
    {
    if (prod==0) return 0;
    mul=mul/prod;
    xc=s2->org.x-mul*s2->shift.x;
    if ((xc<=b)||(xc>e)) return 0;
    if(_ret_ip)
      {
        p->x=xc;
        p->y=s2->org.y-mul*s2->shift.y;
      }
    return 1;
    }
  return 0;
  }

int4 IntInside(REAL b,REAL e,TLineSegment1* s1,TLineSegment1* s2)
  { 
  TPlaneVect delt=s2->org-s1->org;
  REAL prod=s1->shift%s2->shift,mul=s1->shift%delt,xc;
  if ((mul>0)^(mul+prod>0))
    {
    if (prod==0) return 0;
    mul=mul/prod;
    xc=s2->org.x-mul*s2->shift.x;
    if ((xc<=b)||(xc>e)) return 0;
    return 1;
    }
  return 0;
  }



bool below(REAL X,TLineSegment2* s1,TLineSegment2* s2)
  { 

  return (s1!=s2)&&((X*s1->a+s1->b)<(X*s2->a+s2->b));
  }

template <bool _ret_ip>
int4 IntPoint(TLineSegment2* s1,TLineSegment2* s2,TPlaneVect *p)
  { 
  REAL x1=max(s1->x1,s2->x1);
  REAL x2=min(s1->x2,s2->x2);
  if(x1>=x2)return 0;
  REAL da=s1->a-s2->a;
  if(da==0)return 0;
  if(_ret_ip)
    {
       p->x=(s2->b-s1->b)/da;
       if ((p->x>=x1)&&(p->x<=x2)){p->y=p->x*s1->a+s1->b; return 1;}
    }
  else
    {
      REAL x=(s2->b-s1->b)/da;
      return ((x>=x1)&&(x<=x2));
    }  
  return 0;
  } 

template <bool _ret_ip>
int4 StripePoint(REAL b,REAL e,TLineSegment2* s1,TLineSegment2* s2,TPlaneVect *p)
  { 
  REAL x1=max(s1->x1,s2->x1);
  REAL x2=min(s1->x2,s2->x2);
  x1=max(x1,b);x2=min(x2,e);
  REAL da=s1->a-s2->a;
  if(da==0)return 0;
  if(_ret_ip)
    {
       p->x=(s2->b-s1->b)/da;
       if ((p->x>=x1)&&(p->x<=x2)){p->y=p->x*s1->a+s1->b; return 1;}
    }
  else
    {
      REAL x=(s2->b-s1->b)/da;
      return ((x>=x1)&&(x<=x2));
    }  
  return 0;
  }

int4 IntInside(REAL b,REAL e,TLineSegment2 *s1,TLineSegment2 *s2)
  { 
  REAL x1=max(s1->x1,s2->x1);
  REAL x2=min(s1->x2,s2->x2);
  x1=max(x1,b);x2=min(x2,e);
  REAL da=s1->a-s2->a;
  if(da==0)return 0;
  REAL x=(s2->b-s1->b)/da;
  if ((x>=x1)&&(x<=x2)) return 1;
  return 0;
  }




bool below(REAL X,TArcSegment* s1,TArcSegment* s2)
  { 
  return (s1!=s2)&&(s1->under(s2->PointAtX(X)));
  }


int4 IntPointsInStripe(REAL x1,REAL x2,TArcSegment *s1,TArcSegment *s2,TPlaneVect *p)
  {
  TPlaneVect oo=s2->org-s1->org;
  REAL org_dist2=oo.get_norm();
  REAL delta=0.5*(s1->r2-s2->r2+org_dist2);
  REAL t=delta/org_dist2;
  REAL l2=s1->r2-delta*t;
  if(l2<0)return 0;
  TPlaneVect m=s1->org+t*oo;
  if(l2==0)
    {
    if((m.x>=x1)&&(m.x<=x2)&&(s1->IsTheSamePart(m))&&(s2->IsTheSamePart(m)))
      {*p=m; return 1;}
    return 0;
    }
  t=sqrt(l2/org_dist2);
  oo=TPlaneVect(oo.y*t,-oo.x*t);
  int4 npoints=0;
  TPlaneVect res=m+oo;
  if((res.x>=x1)&&(res.x<=x2)&&(s1->IsTheSamePart(res))&&(s2->IsTheSamePart(res)))
    { if(p)p[npoints]=res;npoints++;}
  res=m-oo;
  if((res.x>=x1)&&(res.x<=x2)&&(s1->IsTheSamePart(res))&&(s2->IsTheSamePart(res)))
     { if(p)p[npoints]=res;npoints++;}
  return npoints;
  };

template <bool _ret_ip>
int4 IntPoint(TArcSegment* s1,TArcSegment* s2,TPlaneVect *p)
  { 
  REAL x1=max(s1->x1,s2->x1);
  REAL x2=min(s1->x2,s2->x2);
  if(x1>=x2)return 0;
  return IntPointsInStripe(x1,x2,s1,s2,p);
  } 

template <bool _ret_ip>
int4 StripePoint(REAL b,REAL e,TArcSegment* s1,TArcSegment* s2,TPlaneVect *p)
  { 
  REAL x1=max(s1->x1,s2->x1);
  REAL x2=min(s1->x2,s2->x2);
  x1=max(x1,b);x2=min(x2,e);
  return IntPointsInStripe(x1,x2,s1,s2,p);
  }
int4 IntInside(REAL b,REAL e,TArcSegment* s1,TArcSegment* s2)
  { 
  REAL x1=max(s1->x1,s2->x1);
  REAL x2=min(s1->x2,s2->x2);
  x1=max(x1,b);x2=min(x2,e);
  TPlaneVect p[2];
  return IntPointsInStripe(x1,x2,s1,s2,p);
  }

PSeg first_segment_ptr=NULL;
template <class SEGMENT>
class SegmentFunctions
  {
  public:
    static void register_intersection(PRegObj r,PSeg s1,PSeg s2,int4 nInt,TPlaneVect *points)
      {
#ifdef PRINT_SEG_AND_INT
      int4 n1=(SEGMENT*)s1-(SEGMENT*)first_segment_ptr,n2=(SEGMENT*)s2-(SEGMENT*)first_segment_ptr;
      if (n1<n2)
		  printf("found int %i, %i\n",n1,n2);
      else
      printf("found int %i, %i\n", n2, n1);
#endif
      (*(double*)r) += nInt;
      };

    static REAL __YAtX(PSeg s,REAL X)
      {
        return ((SEGMENT *)s)->YAtX(X);
      };

    static bool __Below(REAL x,PSeg s1,PSeg s2)
      {
      return below(x,(SEGMENT*)s1,(SEGMENT*)s2);
      };
    template<bool  _reg_ip>
    static  int4 __FindAndRegIPoints(PSeg s1,PSeg s2,PRegObj intersection_registrator)
      {
      TPlaneVect P[2];
      int4 n=IntPoint<_reg_ip>((SEGMENT*)s1,(SEGMENT*)s2,P);
      if(n)
		  register_intersection(intersection_registrator,s1,s2,n, _reg_ip?P:NULL);
      return n; 
      };

    static  int4 __FindAndRegIPointOnRightOfSWL(REAL swl_x,PSeg s1,PSeg s2,TPlaneVect *p1,TPlaneVect *p2)
      {
          TPlaneVect p[2];
          int4 n=IntPoint<true>((SEGMENT*)s1,(SEGMENT*)s2,p);
//              if(n&&(p->x>swl_x)) *p1=*p;
              
          if(n==0)return 0;
          TPlaneVect *pt=p1;
          int4 res=0;
          int4 f=0,l=1;
          if((n==2)&&(p[0].x>p[1].x)){f=1;l=0;};
          if(p[f].x>swl_x)
            {
                *pt=p[f];
                pt=p2;
                res++;
            }
          if((n==2)&&(p[l].x>swl_x))
            {
                *pt=p[l];
                res++;
            } 
            return res;   
      };
    
    template<bool  _reg_ip>
    static int4 __FindAndRegIPointsInStripe(REAL b,REAL e, PSeg s1,PSeg s2,PRegObj intersection_registrator)
      {
      TPlaneVect P[2];
      int4 n=StripePoint<_reg_ip>(b,e,(SEGMENT*)s1,(SEGMENT*)s2,P);
      if(n)register_intersection(intersection_registrator,s1,s2,n, _reg_ip ? P : NULL);
      return n; 
      };
    static void __BegPoint(PSeg s,REAL &x,REAL &y)
      {
      ((SEGMENT*)s)->BegPoint(x,y);
      };
    static  void __EndPoint(PSeg s,REAL &x,REAL &y)
      {
      ((SEGMENT*)s)->EndPoint(x,y);
      };
    static  int4 __Under(PSeg s1,PSeg s2)
      {
      return ((SEGMENT*)s1)->under(((SEGMENT*)s2)->BegPoint());
      };

    static int4 __IntInside(REAL b,REAL e,PSeg s1,PSeg s2)
      {
      return IntInside(b,e,(SEGMENT*)s1,(SEGMENT*)s2);
      };  

  };

PSeg* create_test_collection(int4 seg_type,int4 n,int4 distr,REAL par)
  {
  PSeg *colls=NULL;
  CRandomValueGen random_gen;
  int4 i;
  colls=new PSeg[n];
  if(distr!=param_defined)par/=n/33.0;
  if(distr==parallel)par/=2.4;
  switch (seg_type)
    {
    case line1:
      {
      TLineSegment1 *Colls=new TLineSegment1[n];
      for(i=0;i<n;i++)colls[i]=Colls+i;
      for (i=0;i<n;i++)
        Colls[i].InitRandom(random_gen,i,distr,par);
      };break;    
    case line2:
      {
      TLineSegment2 *Colls=new TLineSegment2[n];
      for(i=0;i<n;i++)colls[i]=Colls+i;
      for (i=0;i<n;i++)
        Colls[i].InitRandom(random_gen,i,distr,par);
      };break;    
    case arc:
      {
      TArcSegment *Colls=new TArcSegment[n];
      for(i=0;i<n;i++)colls[i]=Colls+i;
      for (i=0;i<n;i++)
        Colls[i].InitRandom(random_gen,i,distr,par);
      };break;    
    }
  first_segment_ptr=colls[0];
  return colls;  
  };

void  delete_test_collection(int4 seg_type,PSeg *colls)
  {
  switch (seg_type)
    {
    case line1:delete[]  (TLineSegment1 *)colls[0];break;
    case line2:delete[]  (TLineSegment2 *)colls[0];break;
    case arc:delete[]  (TArcSegment *)colls[0];break;
    }
  delete[]  colls; 
  };

const int4 reg_obj_margin = 32;// for reg objects to be in different CPU cash blocks
double reg_objects[reg_obj_margin*n_threads];




double find_intersections(int4 seg_type, int4 SN, PSeg* colls, int4 alg, double* counters, bool dont_need_ip )
  {
  if(colls==NULL) return 0;
  double int_numb=0;
  
//  bool dont_need_ip = (alg == fast_no_ip);
 
  
  BOOL is_line;
  if ((seg_type == line1) || (seg_type == line2))
	  {
		  CLineFinder intersection_finder;
		  switch (seg_type)
			  {
			  case line1:
				  intersection_finder.set_segm_fuctions(
					  SegmentFunctions<TLineSegment1>::__Below,
					  dont_need_ip ? SegmentFunctions<TLineSegment1>::__FindAndRegIPoints<false> : SegmentFunctions<TLineSegment1>::__FindAndRegIPoints<true>,
					  dont_need_ip ? SegmentFunctions<TLineSegment1>::__FindAndRegIPointsInStripe<false> : SegmentFunctions<TLineSegment1>::__FindAndRegIPointsInStripe<true>,
					  SegmentFunctions<TLineSegment1>::__IntInside,
					  SegmentFunctions<TLineSegment1>::__BegPoint,
					  SegmentFunctions<TLineSegment1>::__EndPoint,
					  SegmentFunctions<TLineSegment1>::__Under,
					  SegmentFunctions<TLineSegment1>::__FindAndRegIPointOnRightOfSWL,
					  SegmentFunctions<TLineSegment1>::__YAtX,
					  SegmentFunctions<TLineSegment1>::register_intersection,
					  &int_numb,
					  is_line = TLineSegment1::is_line
				  );
				  break;
			  case line2:
				  intersection_finder.set_segm_fuctions(
					  SegmentFunctions<TLineSegment2>::__Below,
					  dont_need_ip ? SegmentFunctions<TLineSegment2>::__FindAndRegIPoints<false> : SegmentFunctions<TLineSegment2>::__FindAndRegIPoints<true>,
					  dont_need_ip ? SegmentFunctions<TLineSegment2>::__FindAndRegIPointsInStripe<false> : SegmentFunctions<TLineSegment2>::__FindAndRegIPointsInStripe<true>,
					  SegmentFunctions<TLineSegment2>::__IntInside,
					  SegmentFunctions<TLineSegment2>::__BegPoint,
					  SegmentFunctions<TLineSegment2>::__EndPoint,
					  SegmentFunctions<TLineSegment2>::__Under,
					  SegmentFunctions<TLineSegment2>::__FindAndRegIPointOnRightOfSWL,
					  SegmentFunctions<TLineSegment2>::__YAtX,
					  SegmentFunctions<TLineSegment2>::register_intersection,
					  &int_numb,
					  is_line = TLineSegment2::is_line
				  );
				  break;
			  }
		  switch (alg)
			  {
				  case triv:intersection_finder.trivial(SN, colls); break;
				  case simple_sweep:intersection_finder.simple_sweep(SN, colls); break;
				  case fast: intersection_finder.balaban_fast(SN, colls); break;
				  case mem_save: intersection_finder.balaban_memory_save(SN, colls); break;
				  case fast_no_ip: intersection_finder.balaban_no_ip(SN, colls); break;
				  case optimal:intersection_finder.balaban_optimal(SN, colls); break;
				  case fast_parallel:
					  {
						  PRegObj additional_reg_obj[n_threads];
						  for (int i = 0; i < n_threads - 1; i++)
						  {
							  reg_objects[i*reg_obj_margin] = 0;
							  additional_reg_obj[i] = (PRegObj *)reg_objects + i*reg_obj_margin;
						  }
						  intersection_finder.fast_parallel(SN, colls, n_threads, additional_reg_obj);
						  for (int i = 0; i < n_threads - 1; i++)int_numb += *(double*)additional_reg_obj[i];
					  }
				  break;
				  case bentley_ottmann:intersection_finder.bentley_ottmann(SN, colls); break;
			  };
		  memcpy(counters, intersection_finder.my_counter, sizeof(intersection_finder.my_counter));
      }

  if (seg_type == arc)
  {
	  CCommonFinder intersection_finder;
	  intersection_finder.set_segm_fuctions(
		  SegmentFunctions<TArcSegment>::__Below,
		  SegmentFunctions<TArcSegment>::__FindAndRegIPoints<true>,
		  SegmentFunctions<TArcSegment>::__FindAndRegIPointsInStripe<true>,
		  SegmentFunctions<TArcSegment>::__IntInside,
		  SegmentFunctions<TArcSegment>::__BegPoint,
		  SegmentFunctions<TArcSegment>::__EndPoint,
		  SegmentFunctions<TArcSegment>::__Under,
		  SegmentFunctions<TArcSegment>::__FindAndRegIPointOnRightOfSWL,
		  SegmentFunctions<TArcSegment>::__YAtX,
		  SegmentFunctions<TArcSegment>::register_intersection,
		  &int_numb,
		  is_line = TArcSegment::is_line
	  );
	  switch (alg)
	  {
		  case triv:intersection_finder.trivial(SN, colls); break;
		  case simple_sweep:intersection_finder.simple_sweep(SN, colls); break;
		  case fast: intersection_finder.balaban_fast(SN, colls); break;
      case mem_save:intersection_finder.balaban_memory_save(SN, colls);break;
		  case fast_no_ip:  break;//incompatible
		  case optimal:intersection_finder.balaban_optimal(SN, colls); break;
		  case fast_parallel:
		  {
			  PRegObj additional_reg_obj[n_threads];
			  for (int i = 0; i < n_threads - 1; i++)
			  {
				  reg_objects[i*reg_obj_margin] = 0;
				  additional_reg_obj[i] = (PRegObj *)reg_objects + i*reg_obj_margin;
			  }
			  intersection_finder.fast_parallel(SN, colls, n_threads, additional_reg_obj);
			  for (int i = 0; i < n_threads - 1; i++)int_numb += *(double*)additional_reg_obj[i];
		  }
		  break;
		  case bentley_ottmann:intersection_finder.bentley_ottmann(SN, colls); break;
	  }

	  memcpy(counters, intersection_finder.my_counter, sizeof(intersection_finder.my_counter));
  }

  
  return int_numb;
  };

  using CLines2CollSimpleReg = CLine2SegmentCollection<SimpleCounter>;
  using CArcCollSimpleReg = CArcSegmentCollection<SimpleCounter>;


 //const uint4 n_threads=16;
  template<class RegistratingSegmentsCollection>
  double find_int( int4 n, PSeg segs,int4 alg)
  {
    auto seg_type = line2;
    SimpleCounter reg;
    reg.counter = 0;
    RegistratingSegmentsCollection coll;
    coll.Init(n, segs);
    coll.SetRegistrator(&reg);
    switch (alg) {
        case triv: {
            CTrivialIntFinder<RegistratingSegmentsCollection> fi;
            fi.find_intersections(&coll);
        } break;
        case simple_sweep:{
            CSimpleSweepIntFinder<RegistratingSegmentsCollection> fi;
            fi.find_intersections(&coll);
        };
            break;
        case fast: {
            CFastIntFinder<RegistratingSegmentsCollection> fi;
            fi.find_intersections(&coll);
        } break;
        case mem_save:{
            CMemSaveIntFinder<RegistratingSegmentsCollection> fi;
            fi.find_intersections(&coll);
        }
            break;
        case optimal: {
          COptimalIntFinder<RegistratingSegmentsCollection> fi;
          fi.find_intersections(&coll);
        }
            break;
        case fast_parallel: {
            SimpleCounter reg_objects[reg_obj_margin * n_threads];
            SimpleCounter* additional_reg_obj[n_threads];
            for (int i = 0; i < n_threads - 1; i++) {
                //reg_objects[i * reg_obj_margin] = 0;
                additional_reg_obj[i] = reg_objects + i * reg_obj_margin;
            }
            CFastIntFinder<RegistratingSegmentsCollection> fi;
            fi.find_intersections(&coll, n_threads, additional_reg_obj);
            for (int i = 0; i < n_threads - 1; i++)
                reg.counter +=additional_reg_obj[i]->counter;

        } break;

    }
    return reg.counter;
  };

  double new_find_int(int4 seg_type, int4 n, PSeg segs, int4 alg)
  {
    switch (seg_type) {
//    case line1:return find_int(n, segs, alg);
    case line2:return find_int<CLines2CollSimpleReg>(n, segs, alg);
    case arc:return find_int<CArcCollSimpleReg>(n, segs, alg);
    }
  };
