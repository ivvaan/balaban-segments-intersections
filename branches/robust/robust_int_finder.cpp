/*
*
*      Copyright (c)  2011-2013 Ivan Balaban 
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

/*
Contains robust implementation of suboptimal algorithm. 
The implementation uses only integers and performs "logical" calculations. 
All float point operations moved out from the implementation to CRobustSegmentCollection. 
So if you implement the CRobustSegmentCollection functions robust way the whole procedure would be robust. 
Function balaban_fast works with logical N segments represented by integers 0,1..N-1. 
The ends for the segment s represented by integers 2s (begin point) and 2s+1 (end point). 
The function assumes that all ends logically different and begin point lies "on the left" from end point. 
Ends also represents logical x-coordinates from 0 to 2N-1.
*/

#include "robust_utils.h"
#include "robust_int_finder.h"

const int4 max_call=5; //max number of sequential recursive call FindR before dividing current strip

#define FIRST_END_RANK(s) (EndRank[FIRST_SEGM_END(s)]) 
#define LAST_END_RANK(s) (EndRank[LAST_SEGM_END(s)])
//#define IS_LINE_SEGMENTS is_line_segments
//#define IS_LINE_SEGMENTS true


// FindInt finds and reports intersection of segment whith staircase QB<stair_index<=QE, given location of the segment begin point - l
void CIntersectionFinder::FindInt(int4 QB,int4 QE,int4 l,int4 s)
  {
  int4 c=l;
  while ((c>QB) && (IntersectionsInCurStripe(Q[c],s))) //first get intersections below
    c--;
  if(is_line_segments&&(c!=l))return; //if found and segment is line it can't be any more
  c=l+1;
  while ((c<=QE) && (IntersectionsInCurStripe(Q[c],s))) // get intersections above
    c++;
  };

// FindIntL finds and reports intersection of L with current staircase QB<stair_index<=QE for fast algorithm
// locations of the segments should be placed to Loc array before calling the function
void CIntersectionFinder::FindIntL(int4 QB,int4 QE,int4 segm_numb)
  {
  for (int4 i=0;i<segm_numb;i++)
    {
    int4 s=L[i];
    int4 loc=R[i];
    if (LAST_END_RANK(s)<RBoundRank) //if segment is not covering the stripe
      {
      FindInt(QB,QE,loc,s);// find intersections
      }
    else // segment is covering the stripe
      if(!is_line_segments)// for line segments all intersections already found 
        {// if segment is not line then find intersections above (below found in Split)
        int4 c=++loc;
        while ((c<=QE) && (IntersectionsInCurStripe(Q[c],s)))  c++;
        int_numb+=c-loc;
        };
    }

  }


// FindIntI - finds and reports intersections of the segment seg with ancestor staircases (listed in program stack)
// from top of the hierarchy to bottom stopping on last staircase for which seg is internal
void CIntersectionFinder::FindIntI(ProgramStackRec * stack_pos,int4 seg)
  {
  int4 r_rank=LAST_END_RANK(seg);
  while(stack_pos->right_bound<=r_rank)stack_pos=stack_pos->prev;// go from bottom to top and find staircase to start
  int4 l,r,m,QB,QE=stack_pos->Q_pos; 
  for(stack_pos=stack_pos->prev;stack_pos;stack_pos=stack_pos->prev)
    {
      l=QB=stack_pos->Q_pos; 
      r=QE+1;
      while((r-l)>1)// binary search
        {
        m=(r+l)>>1;//        m=(r+l)/2;
        INC_COUNTER(0);
        if(collection->Under((Q[m]),seg))
          l=m;
        else
          r=m;
        }
      FindInt(QB,QE,l,seg);
      QE=QB; // move staircase bound to parent staircase
    };
  }//*/
  

// InsDel - inserts or deletes segment associated with end_rank to/from L
//retuns new L size
int4 CIntersectionFinder::InsDel(int4 end_rank,ProgramStackRec * stack_pos,int4 Size)
  {
  int4 j,i;
  Pt current_end=ENDS[end_rank];
  ExchangeLR();  
  int4 seg=END_SEGMENT(current_end);
  if (IS_LAST(current_end)) // if last point - delete
    {
    L[Size]=seg;
    i=0; 
    while(L[i]!=seg)i++;
    Size--;
    for (j=i;j<Size;j++)L[j]=L[j+1];
    }
  else// if first point - insert
    {
    for (i=Size-1;(i>-1)&& (!collection->Under(L[i],seg));i--)
      L[i+1]=L[i];
    L[i+1]=seg;
    Size++;
    E=LAST_POINT_END(current_end);// adjust right bound for segment functions helpers
    FindIntI(stack_pos,seg);// get internal intersections
    }
    return Size;
  }


// Merge merges R and current staircase QB<stair_index<=QE on right bound of the stripe and place it in L, 
// than exchange it back to R and returns new R size.  While merging finds locations of R segments
// in current staircase and appropriate intersections.
int4 CIntersectionFinder::Merge(int4 QB,int4 QE,int4 Size)
  {
  int4 LB=0,new_size=0;
  int4 cur_stair=QB;
  while ((cur_stair<QE)&&(LB<Size))
    {
    int4 s=R[LB];
    if (collection->Below(E,s,Q[cur_stair+1]))
      {
      if (FIRST_END_RANK(R[LB])>LBoundRank) FindInt(QB,QE,cur_stair,s);      
      L[new_size++]=R[LB++];
      }
    else
      L[new_size++]=Q[++cur_stair];
    }
  while (LB<Size)
    {
    if (FIRST_END_RANK(R[LB])>LBoundRank) FindInt(QB,QE,QE,(R[LB]));      
    L[new_size++]=R[LB++];
    }
  while (cur_stair<QE)
    L[new_size++]=Q[++cur_stair];
  ExchangeLR();
  return new_size;
  };

// Split - splits L into new staircase and new L, new L temporary placed to R, at the end L and R exchanged.
// While executing finds location of new L segments in new staicase and finds intersection of the new L segments
// covering current strip with the staircase stairs below.   
int4 CIntersectionFinder::Split(int4 &step_index,int4 Size)
  {
  int4 father_last_step,new_L_size,j;

  father_last_step=step_index; new_L_size=0;
  for (j=0;j<Size;j++)
    {
    if(LAST_END_RANK(L[j])>=RBoundRank)//segment is covering current stripe
      {
      int4 step=step_index;
      int4 s=L[j];
      while((father_last_step<step)&&(IntersectionsInCurStripe(Q[step],s)))
        step--;
      if(step_index!=step)
        {
        int_numb+=step_index-step;
        L[new_L_size]=L[j];
        R[new_L_size++]=step_index;
        }
      else
        {
        step_index++;
        Q[step_index]=L[j];
        }  
      }
    else
      {
      L[new_L_size]=L[j];
      R[new_L_size++]=step_index;
      }
    }
  return new_L_size;
  };


// SearchInStrip finds and reports all intersection in a stripe containing no segment ends
int4 CIntersectionFinder::SearchInStrip(int4 QP,int4 Size)
  {
  if(Size<2){ExchangeLR(); return Size;}
  int4 NQP=QP;
  int4 size=Size;int4 *q=Q+QP+1;
  if(size=Split(NQP,size))
    {
      do
      {
        FindIntL(QP,NQP,size);
        QP=NQP;
      }
      while (size=Split(NQP,size));
      // at this point we can just place Q starting from QP+1 to R and sort it
      for(QP=0;QP<Size;QP++)R[QP]=q[QP];
      fastsort(CSegmCompare(collection ,E),R,Size);
    }
  else
    for(QP=0;QP<Size;QP++)R[QP]=q[QP];
    return Size;
  }

// main function to find intersections for fast algorithm
int4 CIntersectionFinder::FindR( int4 ladder_start_index,int4 interval_left_rank,int4 interval_right_rank,
                                ProgramStackRec *stack_pos, int4 Size,int4 call_numb)
  {
//  int4 RSize;
  B=ENDS[interval_left_rank];E=ENDS[RBoundRank=interval_right_rank]; 
  if (interval_right_rank-interval_left_rank==1)
     return SearchInStrip(ladder_start_index,Size);
  ProgramStackRec stack_rec(ladder_start_index);
  int_numb=0;// variable to count intersections on Split stage
  if (Size>0)
    {
    Size=Split(stack_rec.Q_pos,Size);
    if((ladder_start_index<stack_rec.Q_pos))
      {
      FindIntL(ladder_start_index,stack_rec.Q_pos,Size);
      stack_pos=stack_rec.Set(stack_pos,interval_right_rank);
      }
    };
  if ((int_numb>Size)&&(call_numb<max_call)) //if found a lot of intersections repeat FindR
    Size=FindR(stack_rec.Q_pos,interval_left_rank,interval_right_rank,stack_pos,Size,call_numb+1);
  else //cut stripe on the middle
    {
    int4 m=(interval_left_rank+interval_right_rank)/2; 
    Size=FindR(stack_rec.Q_pos,interval_left_rank,m,stack_pos,Size,0);
    Size=InsDel(m,stack_pos,Size);
    Size=FindR(stack_rec.Q_pos,m,interval_right_rank,stack_pos,Size,0);
    }
  if (ladder_start_index>=stack_rec.Q_pos) return Size;
  B=ENDS[LBoundRank=interval_left_rank];   E=ENDS[interval_right_rank];
  Size=Merge(ladder_start_index,stack_rec.Q_pos,Size);
  return Size;
  };


void CIntersectionFinder::FreeMem()
  {
#define MY_FREE_ARR_MACRO(a) if(a){delete[] a;a=NULL;}
  MY_FREE_ARR_MACRO(Q);
  MY_FREE_ARR_MACRO(EndRank);
  MY_FREE_ARR_MACRO(L);
  MY_FREE_ARR_MACRO(R);
  MY_FREE_ARR_MACRO(ENDS);

  };

void CIntersectionFinder::AllocMem(int4 n)
  {
  EndRank=new int4 [2*n];
  L=new int4  [n+1];
  R=new int4  [n+1];
  ENDS=new int4[2*n];
  Q=new int4  [n];

  };



void CIntersectionFinder::prepare_ends_and_pairs(int4 n)
  {
  int4 i,j;
  for (i=0;i<2*n;i++)
    ENDS[i]=i;
/*  for (i=0,j=0;i<n;i++)
    {
    ENDS[j++]=FIRST_SEGM_END(i);
    ENDS[j++]=LAST_SEGM_END(i);
    }*/
  fastsort(*collection,ENDS,2*n);
  for (int4 i=0;i<2*n;i++) EndRank[ENDS[i]]=i;
  };

CIntersectionFinder::CIntersectionFinder()
  {
  collection=NULL;
  EndRank=NULL;
  Q=NULL;
  L=NULL;
  R=NULL;
  ENDS=NULL;
  memset(my_counter,0,sizeof(my_counter));
  };

void CIntersectionFinder::balaban_fast(CRobustSegmentCollection * _Scoll)
  {
  int4 Size;
  int4 n=_Scoll->GetSegmNumb();
  collection=_Scoll;
  nTotSegm=n;
  AllocMem(n);
  prepare_ends_and_pairs(n);

  if (END_SEGMENT(ENDS[0])==END_SEGMENT(ENDS[2*n-1]))
    {
    int4 s=END_SEGMENT(ENDS[0]);
    for(int4 i=0;i<n;i++)  if (s!=i) _Scoll->FindAndRegIPoints(s,i);
    Size=0;
    }
  else
    {
    L[0]=END_SEGMENT(ENDS[0]); 
    Size=1;
    }
  ProgramStackRec stack_rec(-1,2*n);  //need to be initialized this way
  Size=FindR(-1,0,  n  ,&stack_rec,Size,0);
  Size=InsDel(n,&stack_rec,Size);
  Size=FindR(-1,n,2*n-1,&stack_rec, Size,0);

  FreeMem();
  }
  
void CIntersectionFinder::simple_sweep(CRobustSegmentCollection * _Scoll)
  {
  int4 i,j;
  int4 s;
  int4 pos,k;
  int4 n=_Scoll->GetSegmNumb();
  int4 *Ends=new int4[2*n];
  for (i=0;i<2*n;i++)
    Ends[i]=i;
/*  for (i=0,j=0;i<n;i++)
    {
    Ends[j++]=FIRST_SEGM_END(i);
    Ends[j++]=LAST_SEGM_END(i);
    }*/
  fastsort(*_Scoll,Ends,2*n);
  int4 *sgm=new int4[n+1];
  pos=0;
  for (i=0;i<2*n;i++)
    {
    s=END_SEGMENT(Ends[i]);
    if (IS_LAST(Ends[i]))
      {
      j=0; sgm[pos]=s;
      while (sgm[j]!=s) j++;
      pos--;
      for (;j<pos;j++)sgm[j]=sgm[j+1];
      }
    else
      {
      for (j=0;j<pos;j++)
        _Scoll->FindAndRegIPoints(sgm[j],s);
      sgm[pos]=s;
      pos++;
      }
    }
  delete [] Ends;
  delete []sgm;
  }  
  
void CIntersectionFinder::trivial(CRobustSegmentCollection * _Scoll)
{
  int4 i,j,n=_Scoll->GetSegmNumb();
  int4 m=n-1;
  for(i=0;i<m;i++)
    for(j=i+1;j<n;j++)
      _Scoll->FindAndRegIPoints(i,j);
};
