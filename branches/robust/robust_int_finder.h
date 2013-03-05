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

#include <math.h>
#include <memory.h>
#include <time.h>
#include "robust_utils.h"
#include "seg_coll.h"






struct ProgramStackRec
  {
  ProgramStackRec *prev;
  int4 Q_pos, right_bound;
  inline ProgramStackRec(int4 qp):Q_pos(qp){};
  inline ProgramStackRec(int4 qp,int4 rb):prev(NULL),Q_pos(qp),right_bound(rb){};
  inline ProgramStackRec *Set(ProgramStackRec *p,int4 rb){ prev=p,right_bound=rb;return this;};
  };


typedef int4 int4;

//#define GETEIP(s) (Pairs+s)

#define INC_COUNTER(i) 
//uncoment line bellow to switch on counters
//#define INC_COUNTER(i) (my_counter[i]++)

//typedef int4 int4;

class CSegmCompare
  {
  public:
    CRobustSegmentCollection *collection;
    Pt x;
    CSegmCompare(){};
    CSegmCompare(CRobustSegmentCollection *sc,Pt _x){collection=sc;x=_x;};
    bool lt(int4 s1,int4 s2){return collection->Below(x,(s1),(s2));}
  };


class CIntersectionFinder
  {
  int4 is_line_segments;
  int4 nTotSegm;
  CRobustSegmentCollection *collection;

  int4 * L,* R,* Q;
  int4 *EndRank;
  Pt *ENDS;

  long int_numb;
  int4 LBoundRank,RBoundRank;
  Pt B,E;

  //for parallel version  
  int4 run_to;

  //functions

  //helpers for segment functions
  inline int4 IntersectionsInCurStripe(int4  s1,int4 s2)
    {
    return collection->FindAndRegIPointsInStripe(B,E,(s1),s2);
    };

  inline int4 Intersections(int4  s1,int4 s2)
    {
    return collection->FindAndRegIPoints((s1),s2);
    };

  inline void ExchangeLR()
    {
    register int4* tmp=L;
    L=R; R=tmp;
    };
  //common  functions
  int4 SearchInStrip(int4 QP,int4 Size);
  void FreeMem();
  void prepare_ends(int4 n);
  void prepare_ends_and_pairs(int4 n);

  //functions for fast algorithm
  void FindInt(int4 QB,int4 QE,int4 l,int4 s);
  void FindIntI(ProgramStackRec * stack_pos,int4 seg);
  void FindIntL(int4 QB,int4 QE,int4 segm_numb);
  int4 InsDel(int4 n,ProgramStackRec * stack_pos,int4 Size);
  int4 Merge(int4 QB,int4 QE,int4 Size);
  int4 Split(int4 &step_index,int4 Size);
  void AllocMem(int4 n);

  //additional functions for fast parallel algorithm
/*  CIntersectionFinder *clone(PRegObj robj);
  void unclone();
  int4 CalcLAt(int4 end_rank);*/


  public:
    // some counters to explore algorithm (not nessesary for functioning)
    double my_counter[8];


    CIntersectionFinder();
    ~CIntersectionFinder(){FreeMem();};

    void balaban_fast(CRobustSegmentCollection * _Scoll);
    int4 FindR(int4 ladder_start_index,int4 interval_left_rank,int4 interval_right_rank,ProgramStackRec *stack_pos,int4 Size,int4 call_numb);

/*    void fast_parallel(int4 n,int4 _Scoll[],PRegObj add_reg);
    void Run(){ ProgramStackRec stack_rec(-1,2*nTotSegm);  FindR(-1,0,   run_to  ,&stack_rec,1); }*/

    // simple sweepline algorithm (not Bently-Ottmann!)
    void simple_sweep(CRobustSegmentCollection * _Scoll);
    void trivial(CRobustSegmentCollection * _Scoll);

  };


#endif

