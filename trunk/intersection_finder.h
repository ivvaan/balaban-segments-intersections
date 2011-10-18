#ifndef INTSEGM_H
#define INTSEGM_H
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

#include "utils.h"
#include <math.h>
#include <memory.h>
#include <time.h>


typedef void (*PSeg);
typedef void (*PRegObj);


typedef int4 (*FBelow)(REAL x,PSeg s1,PSeg s2);//retuns if s1 below s2 at x
typedef int4 (*FFindAndRegIPoints)(PSeg s1,PSeg s2,PRegObj intersection_registrator);//finds all intersections of s1 and s2 and register it
typedef int4 (*FFindAndRegIPointsInStripe)(REAL b,REAL e, PSeg s1,PSeg s2,PRegObj intersection_registrator);//finds all intersection points of s1 and s2 in the stripe b,e and register them
typedef int4 (*FIsIntersetInStripe)(REAL b,REAL e, PSeg s1,PSeg s2);//finds if s1 and s2 in the stripe b,e have an intersection point (no registration)
typedef void (*FGetPoint)(PSeg s,REAL &x,REAL &y);//puts a segment s end point coordinates into x and y
typedef int4 (*FUnder)(PSeg s1,PSeg s2);//returns TRUE if s2 begin point is above s1

const int4 undef_loc=-2147483648;
const int4 inherit_each=32; // defines how often in optimal algorithm stairs are inherited,  one inherited per inherit_each


struct EndIndexes
  {
  int4 B,E;
  };

struct IStackRec
  {
  int4 Q_pos, right_bound;
  inline void Set(int4 Q,int4 R){Q_pos=Q;right_bound=R;}
  };

/*struct SegmentInfo
  {
  PSeg segment;
  EndIndexes * end_indexes;
  SegmentInfo(PSeg s=NULL,EndIndexes *ei=NULL):segment(s),end_indexes(ei){};
  };
 
#define GETPSEG(s) (s.segment)
#define GETEIP(s) (s.end_indexes)*/

typedef int4 SegmentInfo;

#define GETPSEG(s) (Scoll[s])
#define GETEIP(s) (Pairs+s)

//typedef int4 SegmentInfo;

struct TEnd
  {
  public:
    REAL x;
    // REAL y;
    int4 islast;
    bool operator<(const TEnd &p)
      {
      return ((x<p.x)||
        ((x==p.x)&&(islast<p.islast)));
      };
  };

class TSEnd:public TEnd
  {
  public:
    PSeg s;
  };

class TSegmEnd:public TEnd
  {
  public:
    SegmentInfo s;
  };


class CSegmCompare
  {
  public:
    PSeg *Scoll;
    FBelow _Below;
    REAL x;
    CSegmCompare(){};
    CSegmCompare(PSeg *sc,FBelow _b,REAL _x){Scoll=sc;_Below=_b;x=_x;};
    bool lt(SegmentInfo s1,SegmentInfo s2){return _Below(x,GETPSEG(s1),GETPSEG(s2));}
  };


class CIntersectionFinder
  {
  int4 is_line_segments;
  int4 nTotSegm;
  PRegObj _reg_obj;
  FBelow _Below;
  FFindAndRegIPoints _FindAndRegIPoints;
  FFindAndRegIPointsInStripe _FindAndRegIPointsInStripe;
  FGetPoint _BegPoint;
  FGetPoint _EndPoint;
  FUnder _Under;
  // for optimal algorithm
  FIsIntersetInStripe _IsIntersetInStripe;


  SegmentInfo * L,* R,* Q;
  EndIndexes *Pairs;
  TSegmEnd *ENDS;
  int4 *Loc;
  PSeg *Scoll;

  IStackRec *I_Stack;
  int4 I_Stack_size;

  long int_numb;
  int4 LBoundIdx,RBoundIdx;
  REAL B,E;

  //for optimal algorithm
  int4 *father_loc;
  //for parallel version  
  int4 run_to;

  //functions

  //checks if IStak full
  void CheckIStack(int4 stack_pos)
    {
    while(stack_pos>I_Stack_size-2)
      {
      IStackRec *tmp=new IStackRec[2*I_Stack_size];
      memcpy(tmp,I_Stack,sizeof(I_Stack[0])*I_Stack_size);
      delete[] I_Stack;
      I_Stack=tmp;
      I_Stack_size*=2;
      };
    };

  //helpers for segment functions
  inline int4 IntersectionsInCurStripe(SegmentInfo  s1,PSeg s2)
    {
    return _FindAndRegIPointsInStripe(B,E,GETPSEG(s1),s2,_reg_obj);
    };

  inline int4 IsIntersectInCurStripe(SegmentInfo  s1,PSeg s2)
    {
    return _IsIntersetInStripe(B,E,GETPSEG(s1),s2);
    };

  inline int4 Intersections(SegmentInfo  s1,PSeg s2)
    {
    return _FindAndRegIPoints(GETPSEG(s1),s2,_reg_obj);
    };

  inline void ExchangeLR()
    {
    register SegmentInfo* tmp=L;
    L=R; R=tmp;
    };
  //common  functions
  void SearchInStrip(int4 QP,int4 LSize);
  void FreeMem();
  void prepare_ends(int4 n);
  void prepare_ends_and_pairs(int4 n);

  //functions for fast algorithm
  void FindInt(int4 QB,int4 QE,int4 l,PSeg s);
  void FindIntI(int4 r_index,int4 stack_pos,PSeg seg);
  void FindIntL(int4 QB,int4 QE,int4 segm_numb);
  int4 InsDel(int4 n,int4 stack_pos,int4 LSize);
  int4 Merge(int4 QB,int4 QE,int4 LSize);
  int4 Split(int4 &step_index,int4 LSize);
  void AllocMem(int4 n);

  //same for optimal algorithm
  void optFindInt(int4 QB,int4 QE,int4 l,PSeg s);
  void optFindIntI(int4 r_index,int4 stack_pos,PSeg seg);
  int4 optInsDel(int4 n,int4 stack_pos,int4 LSize);
  int4 optMerge(int4 QB,int4 QE,int4 LSize);
  int4 optSplit(int4 father_first_step, int4 &step_index,int4 LSize);
  void optAllocMem(int4 n);

  //additional functions for fast parallel algorithm
  CIntersectionFinder *clone(PRegObj robj);
  void unclone();
  int4 CalcLAt(int4 end_index);


  public:
    // some counters to explore algorithm (not nessesary for functioning)
    double my_counter[16];


    CIntersectionFinder();
    ~CIntersectionFinder(){FreeMem();};

    void set_segm_fuctions(FBelow _below, 
      FFindAndRegIPoints findAndRegIPoints,
      FFindAndRegIPointsInStripe findAndRegIPointsInStripe,
      FIsIntersetInStripe isIntersetInStripe,
      FGetPoint begPoint,FGetPoint endPoint,
      FUnder under,PRegObj reg_obj,int4 is_line);

    void balaban_fast(int4 n,PSeg _Scoll[]);
    int4 FindR(int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,int4 stack_pos,int4 LSize);

    void balaban_optimal(int4 n,PSeg _Scoll[]);
    int4 optFindR(int4 father_first_step,int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,
      int4 stack_pos,int4 LSize);

    void fast_parallel(int4 n,PSeg _Scoll[],PRegObj add_reg);
    void Run(){ FindR(-1,0,   run_to  ,-1,1); }

    // trivial algorithm
    void trivial(int4 n,PSeg sgm[]);

    // simple sweepline algorithm (not Bently-Ottmann!)
    void simple_sweep(int4 n,PSeg Scoll[]);

  };


#endif

