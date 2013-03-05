/*
*
*      Copyright (c)  2011-2013  Ivan Balaban 
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
#include "seg_coll.h"
#include <stddef.h>
#include <stdio.h>

/*
Contains test implementation of CRobustSegmentCollection interface. 
CTestRobustSegmentCollection is a set of vertises connected by the edges exposing this graph
as a collection of segments. 
*/


REAL sq(REAL x) {return x*x;}
//typedef void (*PSeg);
//int4 first_segment_ptr=NULL;

#define FIRST_END_POINT(e) (endpoint_idx[FIRST_POINT_END(e)]) 
#define LAST_END_POINT(e) (endpoint_idx[LAST_POINT_END(e)])
#define OTHER_END_POINT(e) (endpoint_idx[OTHER_POINT_END(e)])

 struct SSegmEndIndexes
   {
   int4 B,E;
   bool operator<(SSegmEndIndexes si){return (B<si.B)||((B==si.B)&&(E<si.E));};
   bool operator!=(SSegmEndIndexes si){return (B!=si.B)||(E!=si.E);};
   }; 

class TLineSegment
  {
  public:
    static int4 is_line;
    TLineSegment(){};
    TLineSegment(TPlaneVect *ends,int4 B,int4 E){org=ends[B];shift=ends[E]-ends[B];};
    TLineSegment(TPlaneVect *ends,SSegmEndIndexes *s){org=ends[s->B];shift=ends[s->E]-ends[s->B];};
    // int4 parnum,intpar;
    TPlaneVect org;
    TPlaneVect shift;
/*    void Refine()
      {
      if((shift.x<0)||((shift.x==0)&&(shift.y<0)))
        {org=org+shift;shift=-shift;}
      };*/
    void Init(TLineSegment &s){org=s.org;shift=s.shift;};
    friend int4 below(REAL X, TLineSegment &s1, TLineSegment &s2);
    friend int4 IntPoint( TLineSegment &s1, TLineSegment &s2,TPlaneVect *p);
    friend int4 StripePoint(REAL b,REAL e, TLineSegment &s1, TLineSegment &s2,TPlaneVect *p);
    //for optimal alg only
    friend int4 IntInside(REAL b,REAL e, TLineSegment &s1, TLineSegment &s2);
  };

int4 TLineSegment::is_line=1;


int4 below(REAL X, TLineSegment & s1, TLineSegment & s2,int4 r)
  { 
  REAL dx1=s1.shift.x,dx2=s2.shift.x;
  REAL dy1=s1.shift.y,dy2=s2.shift.y;
  REAL y1,y2;
  y1=(s1.org.y*dx1+(X-s1.org.x)*dy1)*dx2;
  y2=(s2.org.y*dx2+(X-s2.org.x)*dy2)*dx1;
  if(y1==y2)
  {
    REAL z=s1.shift%s2.shift;
    if(z==0) return r;
    return z>0;
  }
  return y1<y2;
  }


int4 IntPoint(TLineSegment & s1,TLineSegment & s2,TPlaneVect *p)
  { 
  REAL b1=max(s1.org.x,s2.org.x);
  REAL b2=min(s1.org.x+s1.shift.x,s2.org.x+s2.shift.x);
  if(b1>=b2)return 0;
  REAL prod=s1.shift%s2.shift,mul;
  if(prod==0) return 0;
  TPlaneVect delt=s2.org-s1.org;
  if (((mul=s1.shift%delt)>0)^(mul+prod>0))
    if (((mul=delt%s2.shift)>0)^(mul-prod>0))
      {
      if (prod==0) return 0;
      *p=s1.org+((REAL)fabs(mul/prod))*s1.shift;
      return 1;
      }
    return 0;
  } 

int4 StripePoint(REAL b,REAL e,TLineSegment & s1,TLineSegment & s2,TPlaneVect *p)
  { 
  REAL prod=s1.shift%s2.shift,mul,xc;
  if(prod==0) return 0;
  TPlaneVect delt=s2.org-s1.org;
  if (((mul=s1.shift%delt)>0)^(mul+prod>0))
    {
    if (prod==0) return 0;
    mul=mul/prod;
    xc=s2.org.x-mul*s2.shift.x;
    if ((xc<b)||(xc>e)) return 0;
    p->x=xc;
    p->y=s2.org.y-mul*s2.shift.y;
    return 1;
    }
  return 0;
  }

int4 IntInside(REAL b,REAL e,TLineSegment & s1,TLineSegment & s2)
  { 
  REAL prod=s1.shift%s2.shift;
  if(prod==0) return 0;
  TPlaneVect delt=s2.org-s1.org;
  REAL mul=s1.shift%delt,xc;
  if ((mul>0)^(mul+prod>0))
    {
    if (prod==0) return 0;
    mul=mul/prod;
    xc=s2.org.x-mul*s2.shift.x;
    if ((xc<=b)||(xc>e)) return 0;
    return 1;
    }
  return 0;
  }
 


class CTestRobustSegmentCollection:public CRobustSegmentCollection
  {
  public:
    CTestRobustSegmentCollection()
    {
    nPoints=nSegments=0;
    end_points=NULL;
    endpoint_idx=NULL;
    se_index=NULL;
    int_numb=0;
    };

    ~CTestRobustSegmentCollection()
    {
      if(end_points){delete[] end_points; end_points=NULL;};
      if(se_index){delete[] se_index; se_index=NULL;endpoint_idx=NULL;};
    };




 virtual int4 GetSegmNumb(){return nSegments;}; 
 virtual int4 Below(Pt x,int4 s1,int4 s2)
   {
     int4 point_idx=endpoint_idx[x];
     int4 B1,E1,B2,E2;
     int4 flag=0;
     if((B1=se_index[s1].B)==point_idx) flag+=1;
     if((B2=se_index[s2].B)==point_idx) flag+=2;
     if((E1=se_index[s1].E)==point_idx) flag+=4;
     if((E2=se_index[s2].E)==point_idx) flag+=8;
     REAL z;
     TPlaneVect v1,v2;
     switch(flag)
       {
         case 3://begin points coincide with x
           {
             v1= end_points[E1]-end_points[B1];
             v2= end_points[E2]-end_points[B2];
             z=v1%v2;
             if(z==0)//segments (patially) coincide
               return s1<s2;
             return z>0;
           }; 
         case 12://end points coincide with x
           {
             v1= end_points[E1]-end_points[B1];
             v2= end_points[E2]-end_points[B2];
             z=v1%v2;
             if(z==0) //segments (patially) coincide
               return s1<s2;
             return z<0;
           }; break;
         case 4://end point of s1 coincides with x
           {
             TPlaneVect v2= end_points[E2]-end_points[B2];
             z=(end_points[E1]-end_points[B2])%v2;
             if(z==0)//end point of s1 lies on s2
             {
               z=v2%(end_points[E1]-end_points[B1]);
               if(z==0)//segments (patially) coincide
                 return s1<s2;
             }
             return z>0;
           };break;
         case 8://end point of s2 coincides with x
           {
             TPlaneVect v1= end_points[E1]-end_points[B1];
             z=(end_points[E2]-end_points[B1])%v1;
             if(z==0)//end point of s2 lies on s1
             {
               z=v1%(end_points[E2]-end_points[B2]);
               if(z==0)//segments (patially) coincide
                 return s1<s2;
             }
             return z<0;
           };break;
         case 0:
         case 1:
         case 2:break;
         default: throw "unsupported case";
       };
       //return res^^exchanged;
       return below(end_points[point_idx].x,TLineSegment(end_points,B1,E1),TLineSegment(end_points,B2,E2),s1<s2);
   };//retuns if s1 below s2 at x
 virtual int4 FindAndRegIPoints(int4 s1,int4 s2)//finds all intersections of s1 and s2 and register it
   {
     int4 B1,E1,B2,E2;
     B1=se_index[s1].B;
     B2=se_index[s2].B;
     E1=se_index[s1].E;
     E2=se_index[s2].E;
     if((B1==B2)||(E1==B2))return 0;
     if((B1==E2)||(E1==E2))return 0; 
      TPlaneVect p;
      int4 res=IntPoint(TLineSegment(end_points,B1,E1),TLineSegment(end_points,B2,E2), &p);
      if(res) RegisterIntersection(s1,s2,&p);
      return res;
   };
 virtual int4 FindAndRegIPointsInStripe(Pt b,Pt e, int4 s1,int4 s2)//finds all intersection points of s1 and s2 in the stripe b,e and register them
   {
     int4 B1,E1,B2,E2;
     B1=se_index[s1].B;
     B2=se_index[s2].B;
     E1=se_index[s1].E;
     E2=se_index[s2].E;
     if((B1==B2)||(E1==E2))
       {  
      //if((end_points[E1]-end_points[B1])%(end_points[E2]-end_points[B2])==0) {    }
      return 0; 
       };
      TPlaneVect p;
      int4 res=StripePoint(end_points[endpoint_idx[b]].x,end_points[endpoint_idx[e]].x,TLineSegment(end_points,B1,E1),TLineSegment(end_points,B2,E2), &p);
      if(res) RegisterIntersection(s1,s2,&p);
      return res;
   }
 virtual int4 Under(int4 s1,int4 s2)//returns TRUE if s2 begin point is above s1
 {
   int4 B1,E1,B2,E2;
     B1=se_index[s1].B;
     B2=se_index[s2].B;
     E1=se_index[s1].E;
     E2=se_index[s2].E;
   TPlaneVect v1= end_points[E1]-end_points[B1];
   REAL z=v1%(end_points[B2]-end_points[B1]);
   if(z==0)//s2 begin point lies on s1
   {
     z=v1%(end_points[E2]-end_points[B2]);
     if(z==0)//segments (patially) coincide
       return s1<s2;
   }
   return z>0;
 };
 virtual int4 IsLineSegments(){return true;};
 virtual bool lt(Pt p1, Pt p2)
   {
     if(END_SEGMENT(p1)==END_SEGMENT(p2))
       return p1<p2;// if points belong to the same segment, then having lower number first (segment must be directed from lef to right)
     int4 e1=endpoint_idx[p1];
     int4 e2=endpoint_idx[p2];
     TPlaneVect *pe1=end_points+e1;
     TPlaneVect *pe2=end_points+e2;
     if((e1==e2)||((pe1->x==pe2->x)&&(pe1->y==pe2->y)))//the same phisical points
       { // if points are same we find other ends indexes
         int4 o1=OTHER_END_POINT(p1); 
         int4 o2=OTHER_END_POINT(p2);
         //IS_LAST(p2)?FIRST_END_POINT(p2):LAST_END_POINT(p2);
         TPlaneVect *po1=end_points+o1;// other ends coords
         TPlaneVect *po2=end_points+o2;
         if((o1==o2)||((po1->x==po2->x)&&(po1->y==po2->y)))//the same phisical points
           return END_SEGMENT(p1)<END_SEGMENT(p2);// if other ends are same, then order according segment number (for both 0 length segments doesn't work correcly!)
         if((o1==e1)||((po1->x==pe1->x)&&(po1->y==pe1->y)))//first segm is 0 length
           return !IS_LAST(p2);
         if((o2==e2)||((po2->x==pe2->x)&&(po2->y==pe2->y)))//second segm is 0 length
           return IS_LAST(p1);
           
         if(IS_LAST(p1))
           {
             if(!IS_LAST(p2))return true;// p2- last, ð1 - first then p1<p2
             // both last points of segments
             double r=(*po1-*pe1)%(*po2-*pe1);
             if(r==0) //segments are parallel
               return END_SEGMENT(p1)<END_SEGMENT(p2);
             return r<0;
             
           }
         else
           {
             if(IS_LAST(p2))return false;// p1- first, ð2 - last then p2<p1
             // both first points of segments
             double r=(*po1-*pe1)%(*po2-*pe1);
             if(r==0) //segments are parallel
               return END_SEGMENT(p1)<END_SEGMENT(p2);
             return r>0;
           }  
       } 
     if((pe1->x==pe2->x))return(pe1->y<pe2->y);
     return (pe1->x<pe2->x);
   };   

  void CreateRandom (int4 nPt,int4 nSeg,CRandomValueGen *rg)
  {
    nPoints=nPt;
    nSegments=nSeg;
    end_points=new TPlaneVect[nPoints];
    se_index=new SSegmEndIndexes[nSegments];
    endpoint_idx=(int4 *)se_index;
    int4 i,n;
    for(i=0;i<nPoints;i++)
      {
        end_points[i].x=rg->GetRandomDouble();
        end_points[i].y=rg->GetRandomDouble();
      }
    
    //int4 ns=2*nSegments
    //SSegmEndIndexes *tmp=new SSegmEndIndexes[ns];
    int4 ns=nSegments;
    SSegmEndIndexes *tmp=se_index;
    for(i=0;i<ns;i++)
      {
        int4 l,f=rg->GetRandomInt(nPoints);
        do
         l=rg->GetRandomInt(nPoints);
        while(end_points[f].x==end_points[l].x);
        if(end_points[f].x<end_points[l].x)//segments MUST be directed from left to right (first point x<last point x)
        {tmp[i].B=f;tmp[i].E=l;}
        else
        {tmp[i].B=l;tmp[i].E=f;}

      }
     /*fastsort(tmp,ns);
     n=0;se_index[n]=tmp[0];
     for(i=1;i<ns;i++)
     {
       if(se_index[n]!=tmp[i])
         se_index[++n]=tmp[i];
       if(n==(nSegments-1))break;
     }
     delete [] tmp;
     nSegments=n+1;*/
  };

  double int_numb;
  void RegisterIntersection(int4 s1,int4 s2,TPlaneVect *p)
  {
/*    int4 tmp;
    if(s1>s2){tmp=s1;s1=s2;s2=tmp;}
    printf("%i,%i\n",s1,s2);*/
    int_numb++;
  };

  private:
    int4 nPoints,nSegments;
    TPlaneVect *end_points;
    int4 *endpoint_idx;
    SSegmEndIndexes *se_index;

  };
  
CRobustSegmentCollection *create_segment_collection (int4 n_ends, int4 n_seg)
  {
  CRandomValueGen random_gen;
  CTestRobustSegmentCollection *tc=new CTestRobustSegmentCollection();
  tc->CreateRandom(n_ends,n_seg,&random_gen);
  return tc;
  };
double get_int_numb_and_reset(CRobustSegmentCollection *rsc)
{
  CTestRobustSegmentCollection *tc=(CTestRobustSegmentCollection *)rsc;
  double res=tc->int_numb;
  tc->int_numb=0;
  return res;
};
void delete_segment_collection (CRobustSegmentCollection *rsc)
  {
  delete (CTestRobustSegmentCollection *)rsc;
  };


/*double  find_intersections(int4 seg_type,int4 SN,int4 *colls,int4 alg, double *counters)
  {
  if(colls==NULL) return 0;
  double int_numb=0;
  CIntersectionFinder intersection_finder;
  double int_numb1=0;
    }

  return int_numb+int_numb1;
  };*/

