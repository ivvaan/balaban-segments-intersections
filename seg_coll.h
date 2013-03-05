#ifndef TEST_COLL_H
#define TEST_COLL_H
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

#include "robust_utils.h"

#define END_SEGMENT(e) (e>>1)
#define IS_LAST(e) (e&1) 
#define FIRST_SEGM_END(s) (s<<1) 
#define LAST_SEGM_END(s) ((s<<1)|1) 

//finds FIRST_SEGM_END(END_SEGMENT(e))
#define FIRST_POINT_END(e) (e&0xFFFFFFFE) 

//funds LAST_SEGM_END(END_SEGMENT(e))
#define LAST_POINT_END(e) (e|1)

//finds another endpoint
#define OTHER_POINT_END(e) (e^1)


enum _Algorithm
  {
  triv=1,
  simple_sweep=2,
  fast=4,
  optimal=8,
  fast_parallel=16
  };

enum _Segment
  {
  line1=0,line2,arc
  };
enum _Distribution
  {
  random=0,mixed,small,param_defined
  };
 
// typedef int4 int4;
 typedef int4 Pt;

class CRobustSegmentCollection
  {
  public:
 virtual int4 GetSegmNumb()=0; 
 virtual int4 Below(Pt x,int4 s1,int4 s2)=0;//retuns if s1 below s2 at x
 virtual int4 FindAndRegIPoints(int4 s1,int4 s2)=0;//finds all intersections of s1 and s2 and register it
 virtual int4 FindAndRegIPointsInStripe(Pt b,Pt e, int4 s1,int4 s2)=0;//finds all intersection points of s1 and s2 in the stripe b,e and register them
 virtual int4 Under(int4 s1,int4 s2)=0;//returns TRUE if s2 begin point is above s1
 virtual int4 IsLineSegments()=0;
 virtual bool lt(Pt e1, Pt e2)=0;   
  };
  
CRobustSegmentCollection *create_segment_collection (int4 n_ends, int4 n_seg);
double get_int_numb_and_reset(CRobustSegmentCollection *rsc);
void delete_segment_collection (CRobustSegmentCollection *rsc);


#endif