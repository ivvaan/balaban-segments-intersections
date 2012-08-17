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

typedef int4 (*FFindAndRegIPointOnRightOfSWL)(REAL swl_x,PSeg s1,PSeg s2,PRegObj intersection_registrator,TPlaneVect *p1,TPlaneVect *p2); 
//finds first intersection on right of sweep line (swl_x), register it and place x coord to int_point_x


const int4 undef_loc=-2147483648;
const int4 inherit_each=32; // defines how often in optimal algorithm stairs are inherited,  one inherited per inherit_each
const int4 max_call=5; //max number of sequential recursive call (opt)FindR before dividing current strip


struct EndIndexes
  {
  int4 B,E;
  };

struct ProgramStackRec
  {
  ProgramStackRec *prev;
  int4 Q_pos, right_bound;
  inline ProgramStackRec(int4 qp):Q_pos(qp){};
  inline ProgramStackRec(int4 qp,int4 rb):prev(NULL),Q_pos(qp),right_bound(rb){};
  inline ProgramStackRec *Set(ProgramStackRec *p,int4 rb){ prev=p,right_bound=rb;return this;};
  };

typedef int4 SegmentInfo;


//uncomment line bellow to switch on counters
//#define COUNTERS_ON

#ifdef COUNTERS_ON
#define INC_COUNTER(i) (my_counter[i]++)
#define SET_COUNTER(i,value) my_counter[i]=value
#define ADD_COUNTER(i,value) my_counter[i]+=value
#define SET_COUNTER_IFGT(i,value) if(value>my_counter[i])my_counter[i]=value
#else
#define INC_COUNTER(i) 
#define SET_COUNTER(i,value)
#define ADD_COUNTER(i,value)
#define SET_COUNTER_IFGT(i,value)
#endif

struct TEnd
  {
    REAL x;
    // REAL y;
    int4 islast;
    bool operator<(const TEnd &p)
      {
      return ((x<p.x)||
        ((x==p.x)&&(islast<p.islast)));
      };
  };

struct TSEnd:public TEnd
  {
    PSeg s;
  };

struct TSegmEnd:public TEnd
  {
    SegmentInfo s;
  };


struct CSegmCompare
  {
    PSeg *Scoll;
    FBelow _Below;
    REAL x;
    CSegmCompare(){};
    CSegmCompare(PSeg *sc,FBelow _b,REAL _x){Scoll=sc;_Below=_b;x=_x;};
    bool lt(SegmentInfo s1,SegmentInfo s2){return _Below(x,Scoll[s1],Scoll[s2]);}
  };
  
  
//for Bentley & Ottmann algorithm
  
struct TBinTreeNode// for sweep line support 
  {
    int4 left,right,next,prev;
    int4 segment;
    int4 seg_node_idx;
            // always should be 
            //SweepLine[SweepLine[s].seg_node_idx].segment==s
            //SweepLine[SweepLine[node].segment].seg_node_idx==node
    TBinTreeNode()
      {
         segment=seg_node_idx=left=right=next=prev=-1;
      };
  }; 

struct TEvent
  {
    TPlaneVect pt;
    int4 s1,s2;
  } ;  


class CIntersectionFinder
  {
  int4 is_line_segments;//shows if segments is straitline and more fast calculations are possible
  
  //begin functions to deal with segments
  PRegObj _reg_obj;
  FBelow _Below;
  FFindAndRegIPoints _FindAndRegIPoints;
  FFindAndRegIPointsInStripe _FindAndRegIPointsInStripe;
  FGetPoint _BegPoint;
  FGetPoint _EndPoint;
  FUnder _Under;
  // for optimal algorithm
  FIsIntersetInStripe _IsIntersetInStripe;
  
  //for Bentley & Ottmann alg
  FFindAndRegIPointOnRightOfSWL _FindAndRegIPointOnRightOfSWL;
  // end functions to deal with segments
  
  //begin strutures for Balaban algorithms
  int4 nTotSegm;
  SegmentInfo * L,* R,* Q;
  EndIndexes *SegmentBE;
  TSegmEnd *ENDS;
  int4 *Loc;
  PSeg *Scoll;

  long int_numb;
  int4 LBoundIdx,RBoundIdx;
  REAL B,E;

  //for optimal algorithm
  int4 *father_loc;
  //for parallel version  
  int4 run_to;
  //end strutures for Balaban algorithms

  //begin strutures for Bentley & Ottmann algorithm
  REAL sweep_line_x;
  int4 swl_root;
  TBinTreeNode *SweepLine;
  int4 events_n,events_max;
  TEvent *Events;
//  int4 *SWL_stack;
  //end strutures for Bentley & Ottmann algorithm
  
  //private functions for Balaban algorithms

  //helpers for segment functions
  inline int4 IntersectionsInCurStripe(SegmentInfo  s1,PSeg s2)
    {
    return _FindAndRegIPointsInStripe(B,E,Scoll[s1],s2,_reg_obj);
    };

  inline int4 IsIntersectInCurStripe(SegmentInfo  s1,PSeg s2)
    {
    return _IsIntersetInStripe(B,E,Scoll[s1],s2);
    };

  inline int4 Intersections(SegmentInfo  s1,PSeg s2)
    {
    return _FindAndRegIPoints(Scoll[s1],s2,_reg_obj);
    };

  inline void ExchangeLR()
    {
    register SegmentInfo* tmp=L;
    L=R; R=tmp;
    };
  //common  functions
  int4 SearchInStrip(int4 QP,int4 Size);
  void AllocMem(int4 n);
  void FreeMem();
  int4 prepare_ends(int4 n);

  //functions for fast algorithm
  void FindInt(int4 QB,int4 QE,int4 l,PSeg s);
  void FindIntI(int4 r_index,ProgramStackRec * stack_pos,PSeg seg);
  void FindIntL(int4 QB,int4 QE,int4 segm_numb);
  int4 InsDel(int4 n,ProgramStackRec * stack_pos,int4 Size);
  int4 Merge(int4 QB,int4 QE,int4 Size);
  int4 Split(int4 &step_index,int4 Size);

  //same for optimal algorithm
  void optFindInt(int4 QB,int4 QE,int4 l,PSeg s);
  void optFindIntI(int4 r_index,ProgramStackRec *stack_pos,PSeg seg);
  int4 optInsDel(int4 n,ProgramStackRec *stack_pos,int4 Size);
  int4 optMerge(int4 QB,int4 QE,int4 Size);
  int4 optSplit(int4 father_first_step, int4 &step_index,int4 Size);

  //additional functions for fast parallel algorithm
  CIntersectionFinder *clone(PRegObj robj);
  void unclone();
  int4 CalcLAt(int4 end_index);
  
  //private functions for Bentley & Ottmann algorithm
 void SweepLineInsert(int4 s);
 void SweepLineDelete(int4 s);
 void SweepLineExchange(int4 s1,int4 s2);
 void IntOnRightOfSWL(int4 s1,int4 s2);
 void PrepareEvents();
 void EventsDelMin();
 void EventsAddNew();
 void CheckEventsMem();
 void AllocBOMem();
 void FreeBOMem();
// void SweepLinePrint();


  public:
    // some counters to explore algorithm (not nessesary for functioning)
    double my_counter[8];


    CIntersectionFinder();
    ~CIntersectionFinder(){FreeMem();};

    void set_segm_fuctions(FBelow _below, 
      FFindAndRegIPoints findAndRegIPoints,
      FFindAndRegIPointsInStripe findAndRegIPointsInStripe,
      FIsIntersetInStripe isIntersetInStripe,
      FGetPoint begPoint,FGetPoint endPoint,
      FUnder under,PRegObj reg_obj,int4 is_line,
     FFindAndRegIPointOnRightOfSWL findAndRegIPointOnRightOfSWL);

    void balaban_fast(int4 n,PSeg _Scoll[]);
    int4 FindR(int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,ProgramStackRec *stack_pos,int4 Size,int4 call_numb);

    void balaban_optimal(int4 n,PSeg _Scoll[]);
    int4 optFindR(int4 father_first_step,int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,
      ProgramStackRec *stack_pos,int4 Size,int4 call_numb);

    void fast_parallel(int4 n,PSeg _Scoll[],PRegObj add_reg);
    void Run(){ ProgramStackRec stack_rec(-1,2*nTotSegm);  FindR(-1,0,   run_to  ,&stack_rec,1,0); }

    // trivial algorithm
    void trivial(int4 n,PSeg sgm[]);

    // simple sweepline algorithm (not Bentley & Ottmann!)
    void simple_sweep(int4 n,PSeg Scoll[]);

    //  Bentley & Ottmanns
    void bentley_ottmann(int4 n,PSeg _Scoll[]);

  };


#endif

