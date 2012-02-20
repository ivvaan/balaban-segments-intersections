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
#include "intersection_finder.h"


//#define IS_LINE_SEGMENTS is_line_segments
//#define IS_LINE_SEGMENTS true

#define IS_ORIGINAL(v) (v<1)

// FindInt finds and reports intersection of segment whith staircase QB<stair_index<=QE, given location of the segment begin point - l
void CIntersectionFinder::FindInt(int4 QB,int4 QE,int4 l,PSeg s)
  {
  int4 c=l;
  while ((c>QB) && (IntersectionsInCurStripe(Q[c],s))) //first get intersections below
    c--;
  if(is_line_segments&&(c!=l))return; //if found and segment is line it can't be any more
  c=l+1;
  while ((c<=QE) && (IntersectionsInCurStripe(Q[c],s))) // get intersections above
    c++;
  };

// optFindInt finds and reports intersection of segment whith staircase QB<stair_index<=QE, given location of the segment begin point - l
// for optimal algorithm. Intersections whith inherited stairs are not reported.
void CIntersectionFinder::optFindInt(int4 QB,int4 QE,int4 l,PSeg s)
  {
  int4 c=l;
  while (c>QB)
    {
    if(IS_ORIGINAL(father_loc[c])) 
      {
      if(!IntersectionsInCurStripe(Q[c],s)) break;
      }
    else //if inherited stair just check intersection
      if(!IsIntersectInCurStripe(Q[c],s)) break;
    c--;  
    }
  int_numb+=l-c;
  if(is_line_segments&&(c!=l))return;
  c=++l;
  while (c<=QE)
    {
    if(IS_ORIGINAL(father_loc[c])) 
      {
      if(!IntersectionsInCurStripe(Q[c],s)) break;
      }
    else
      if(!IsIntersectInCurStripe(Q[c],s)) break;
    c++;  
    }
  int_numb+=c-l;
  };

// FindIntL finds and reports intersection of L with current staircase QB<stair_index<=QE for fast algorithm
// locations of the segments should be placed to Loc array before calling the function
void CIntersectionFinder::FindIntL(int4 QB,int4 QE,int4 segm_numb)
  {
  for (int4 i=0;i<segm_numb;i++)
    {
    PSeg s=GETPSEG(L[i]);
    int4 loc=Loc[i];
    if (GETEIP(L[i])->E<RBoundIdx) //if segment is not covering the stripe
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

// optFindIntI - finds and reports intersections of the segment seg with ancestor staircases (listed in I_Stack)
// from bottom of the hierarchy to top starting from first staircase for which seg is internal
void CIntersectionFinder::optFindIntI(int4 r_index,ProgramStackRec *stack_pos,PSeg seg)
  {
  while(stack_pos->right_bound<=r_index)stack_pos=stack_pos->prev;// go from bottom to top and find staircase to start
    E=ENDS[r_index].x;// adjust right bound for segment functions helpers
    int4  l=undef_loc,r,m,QB,QE=stack_pos->Q_pos; 
    for(stack_pos=stack_pos->prev;stack_pos;stack_pos=stack_pos->prev)// for all staircases
      {
      QB=stack_pos->Q_pos; 
      if(l==undef_loc)// search all range to find seg location for first or small staircase 
        {l=QB;r=QE+1;}
      else
        { // else search from l to l+inherit_each
        l=abs(l);
        if(l<QB)l=QB;
        r=l+inherit_each;           
        if(r>QE+1)r=QE+1;
        }
      while((r-l)>1)// binary search
        {
        m=(r+l)>>1;//        m=(r+l)/2;
        //my_counter[0]++;
        if(_Under(GETPSEG(Q[m]),seg))
          l=m;
        else
          r=m;
        }
      optFindInt(QB,QE,l,seg);
      QE=QB; // move staircase bound to parent staircase
      l=father_loc[l];  //using father_loc get approximate location in parent staircase 
      };
  }


// FindIntI - finds and reports intersections of the segment seg with ancestor staircases (listed in I_Stack)
// from top of the hierarchy to bottom stopping on last staircase for which seg is internal
void CIntersectionFinder::FindIntI(int4 r_index,ProgramStackRec * stack_pos,PSeg seg)
  {
  while(stack_pos->right_bound<=r_index)stack_pos=stack_pos->prev;// go from bottom to top and find staircase to start
  E=ENDS[r_index].x;// adjust right bound for segment functions helpers
  int4 l,r,m,QB,QE=stack_pos->Q_pos; 
  for(stack_pos=stack_pos->prev;stack_pos;stack_pos=stack_pos->prev)
    {
      l=QB=stack_pos->Q_pos; 
      r=QE+1;
      while((r-l)>1)// binary search
        {
        m=(r+l)>>1;//        m=(r+l)/2;
        //my_counter[0]++;
        if(_Under(GETPSEG(Q[m]),seg))
          l=m;
        else
          r=m;
        }
      FindInt(QB,QE,l,seg);
      QE=QB; // move staircase bound to parent staircase
    };
  }//*/
  

// InsDel - inserts or deletes segment associated with end_index to/from L
//retuns new L size
int4 CIntersectionFinder::InsDel(int4 end_index,ProgramStackRec * stack_pos,int4 Size)
  {
  int4 j,i;
  TSegmEnd *current_end=ENDS+end_index;
  ExchangeLR();  
  if (current_end->islast) // if endpoint - delete
    {
    int4 sn=L[Size]=current_end->s;
    i=0; 
    while(L[i]!=sn)i++;
    Size--;
    for (j=i;j<Size;j++)L[j]=L[j+1];
    }
  else// if startpoint - insert
    {
    PSeg seg=GETPSEG(current_end->s);
    for (i=Size-1;(i>-1)&& (!_Under(GETPSEG(L[i]),seg));i--)
      L[i+1]=L[i];
    L[i+1]=current_end->s;
    Size++;

    FindIntI(GETEIP(current_end->s)->E,stack_pos,seg);// get internal intersections
    }
    return Size;
  }

// see InsDel comments
int4 CIntersectionFinder::optInsDel(int4 end_index,ProgramStackRec *stack_pos,int4 Size)
  {
  int4 j,i;
  TSegmEnd *current_end=ENDS+end_index;
  ExchangeLR();
  if (current_end->islast)
    {
    int4 sn=L[Size]=current_end->s;
    i=0; 
    while(L[i]!=sn)i++;
    Size--;
    for (j=i;j<Size;j++)L[j]=L[j+1];
    }
  else
    {
    PSeg seg=GETPSEG(current_end->s);
    for (i=Size-1;(i>-1)&& (!_Under(GETPSEG(L[i]),seg));i--)
      L[i+1]=L[i];
    L[i+1]=current_end->s;
    Size++;

    optFindIntI(GETEIP(current_end->s)->E,stack_pos,seg);
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
    PSeg s=GETPSEG(R[LB]);
    if (_Below(E,s,GETPSEG(Q[cur_stair+1])))
      {
      if (GETEIP(R[LB])->B>LBoundIdx) FindInt(QB,QE,cur_stair,s);      
      L[new_size++]=R[LB++];
      }
    else
      L[new_size++]=Q[++cur_stair];
    }
  while (LB<Size)
    {
    if (GETEIP(R[LB])->B>LBoundIdx) FindInt(QB,QE,QE,GETPSEG(R[LB]));      
    L[new_size++]=R[LB++];
    }
  while (cur_stair<QE)
    L[new_size++]=Q[++cur_stair];
  ExchangeLR();
  return new_size;
  };

// Same as Merge but merges only original stairs
int4 CIntersectionFinder::optMerge(int4 QB,int4 QE,int4 Size)
  {
  int4 LB=0,new_size=0;
  int4 cur_stair=QB;
  while ((cur_stair<QE)&&(LB<Size))
    {
    PSeg s=GETPSEG(R[LB]);
    if (_Below(E,s,GETPSEG(Q[cur_stair+1])))
      {
      if (GETEIP(R[LB])->B>LBoundIdx) optFindInt(QB,QE,cur_stair,s);      
      L[new_size++]=R[LB++];
      }
    else
      {
      cur_stair++;
      if(IS_ORIGINAL(father_loc[cur_stair]))
        L[new_size++]=Q[cur_stair];
      }
    }  
    while (LB<Size)
      {
      if (GETEIP(R[LB])->B>LBoundIdx) optFindInt(QB,QE,QE,GETPSEG(R[LB]));      
      L[new_size++]=R[LB++];
      }
    while (cur_stair<QE)
      {
      cur_stair++;
      if(IS_ORIGINAL(father_loc[cur_stair]))
        L[new_size++]=Q[cur_stair];
      }
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
    if(GETEIP(L[j])->E>=RBoundIdx)//segment is covering current stripe
      {
      int4 step=step_index;
      PSeg s=GETPSEG(L[j]);
      while((father_last_step<step)&&(IntersectionsInCurStripe(Q[step],s)))
        step--;
      if(step_index!=step)
        {
        int_numb+=step_index-step;
        R[new_L_size]=L[j];
        Loc[new_L_size++]=step_index;
        }
      else
        {
        step_index++;
        Q[step_index]=L[j];
        }  
      }
    else
      {
      R[new_L_size]=L[j];
      Loc[new_L_size++]=step_index;
      }
    }
  ExchangeLR();
  return new_L_size;
  };

// same as Split but  add a part of the parent staircase stairs (inherited stairs) to the new one.
int4 CIntersectionFinder::optSplit(int4 cur_father_pos,int4 &step_index,int4 Size)
  {
  int4 father_last_step,new_L_size,cur_L_pos;

  father_last_step=step_index; 
  new_L_size=cur_L_pos=0;

  while ((cur_L_pos<Size)&&(cur_father_pos<=father_last_step))
    {
    PSeg s=GETPSEG(L[cur_L_pos]);
    if(_Below(B,s,GETPSEG(Q[cur_father_pos]))) // current segment below candidate to inherit
      {
      if(GETEIP(L[cur_L_pos])->E>=RBoundIdx) // current segment covering strip
        {
        if(
          ( (step_index>father_last_step)&&(IsIntersectInCurStripe(Q[step_index],s) ))
          || (IsIntersectInCurStripe(Q[cur_father_pos],s))
          ) // current segment intersect last stair or candidate to inherit 
          { //add it to new L (R)
          R[new_L_size]=L[cur_L_pos];
          Loc[new_L_size++]=step_index;
          }
        else
          { //add it to Q as original
          step_index++;
          // if segment is over this stair it is over over previous inherited segment in father staircase
          // negative value of  father_loc[step_index] means it is original segment ( always cur_father_pos  >  inherit_each )              
          father_loc[step_index]=inherit_each-cur_father_pos; 
          Q[step_index]=L[cur_L_pos];
          }  
        }
      else
        {
        R[new_L_size]=L[cur_L_pos];
        Loc[new_L_size++]=step_index;
        }
      cur_L_pos++; // move to next in L
      }
    else
      {  // add candidate to inherit to Q
      step_index++;
      // if segment is over this stair it is over the same inherited segment in father staircase  
      // positive value of  father_loc[step_index] means it is inherited segment                  
      father_loc[step_index]=cur_father_pos;
      Q[step_index]=Q[cur_father_pos];
      cur_father_pos+=inherit_each;
      };
    }    
  for(;cur_L_pos<Size;cur_L_pos++)  // father stair case is over. we are adding rest segments the same manner as in suboptimal case
    if(GETEIP(L[cur_L_pos])->E>=RBoundIdx)//segment is covering current stripe
      {
      if((step_index>father_last_step)&&(IsIntersectInCurStripe(Q[step_index],GETPSEG(L[cur_L_pos]))))// if it intesects last stair
        {  // add it to new L
        R[new_L_size]=L[cur_L_pos];
        Loc[new_L_size++]=step_index;
        }
      else
        {  // add it to Q
        step_index++;
        // negative value of  father_loc[step_index] means it is original segment ( always cur_father_pos  >  inherit_each )              
        father_loc[step_index]=inherit_each-cur_father_pos;
        Q[step_index]=L[cur_L_pos];
        }  
      }
    else
      { // add it to new L
      R[new_L_size]=L[cur_L_pos];
      Loc[new_L_size++]=step_index;
      }
    for(;cur_father_pos<=father_last_step;cur_father_pos+=inherit_each) // L segments is over. We are adding rest inherited father segments to Q
      {
      step_index++;
      // positive value of  father_loc[step_index] means it is inherited segment                  
      father_loc[step_index]=cur_father_pos;
      Q[step_index]=Q[cur_father_pos];
      };
    ExchangeLR();
    return new_L_size;
  };

// SearchInStrip finds and reports all intersection in a stripe containing no segment ends
int4 CIntersectionFinder::SearchInStrip(int4 QP,int4 Size)
  {
  if(Size<2){ExchangeLR(); return Size;}
  int4 NQP=QP;
  int4 size=Size;
  if(size=Split(NQP,size))
    {
    SegmentInfo *q=Q+QP+1;
      do
      {
        FindIntL(QP,NQP,size);
        QP=NQP;
      }
      while (size=Split(NQP,size));
      // at this point we can just place Q starting from QP+1 to L and sort it
      for(QP=0;QP<Size;QP++)R[QP]=q[QP];
      fastsort(CSegmCompare(Scoll,_Below ,E),R,Size);
    }
    return Size;
  }

// main function to find intersections for fast algorithm
int4 CIntersectionFinder::FindR( int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,
                                ProgramStackRec *stack_pos, int4 Size)
  {
//  int4 RSize;
  B=ENDS[interval_left_index].x;E=ENDS[RBoundIdx=interval_right_index].x; 
  if (interval_right_index-interval_left_index==1)
     return SearchInStrip(ladder_start_index,Size);
  ProgramStackRec stack_rec(ladder_start_index);
  int_numb=0;
  if (Size>0)
    {
    Size=Split(stack_rec.Q_pos,Size);
    if((ladder_start_index<stack_rec.Q_pos))
      {
      FindIntL(ladder_start_index,stack_rec.Q_pos,Size);
      stack_pos=stack_rec.Set(stack_pos,interval_right_index);
      }
    };
  if (int_numb>Size)
    Size=FindR(stack_rec.Q_pos,interval_left_index,interval_right_index,stack_pos,Size);
  else
    {
    int4 m=(interval_left_index+interval_right_index)/2; 
    Size=FindR(stack_rec.Q_pos,interval_left_index,m,stack_pos,Size);
    Size=InsDel(m,stack_pos,Size);
    Size=FindR(stack_rec.Q_pos,m,interval_right_index,stack_pos,Size);
    }
  if (ladder_start_index>=stack_rec.Q_pos) return Size;
  B=ENDS[LBoundIdx=interval_left_index].x;   E=ENDS[interval_right_index].x;
  Size=Merge(ladder_start_index,stack_rec.Q_pos,Size);
  return Size;
  };


// main function to find intersections for optimal algorithm
int4 CIntersectionFinder::optFindR(int4 father_first_step, int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,
                                   ProgramStackRec *stack_pos,int4 Size)
  {
  B=ENDS[interval_left_index].x;E=ENDS[RBoundIdx=interval_right_index].x; 
  if (interval_right_index-interval_left_index==1)
     return SearchInStrip(ladder_start_index,Size);
  ProgramStackRec stack_rec(ladder_start_index);
  int_numb=0;
  if(father_first_step+inherit_each<stack_rec.Q_pos)
    {// use optimal variant if new staircase is big
    Size=optSplit(father_first_step,stack_rec.Q_pos,Size);
    for (int4 i=0;i<Size;i++)  
      optFindInt(ladder_start_index,stack_rec.Q_pos,Loc[i],GETPSEG(L[i]));
    }
  else
    {// use fast variant if new staircase is small 
    Size=Split(stack_rec.Q_pos,Size);
    if((ladder_start_index<stack_rec.Q_pos))
      FindIntL(ladder_start_index,stack_rec.Q_pos,Size);
    for (int4 i=ladder_start_index+1;i<=stack_rec.Q_pos;i++)
      father_loc[i]=undef_loc;
    }           
  if((ladder_start_index<stack_rec.Q_pos))
    stack_pos=stack_rec.Set(stack_pos,interval_right_index);
  if (int_numb>Size+1)
    Size=optFindR(ladder_start_index+1, stack_rec.Q_pos,interval_left_index,interval_right_index,stack_pos,Size);
  else
    {
    int4 m=(interval_left_index+interval_right_index)/2; 
    Size=optFindR(ladder_start_index+1,stack_rec.Q_pos,interval_left_index,m,stack_pos,Size);
    Size=optInsDel(m,stack_pos,Size);
    Size=optFindR(ladder_start_index+1,stack_rec.Q_pos,m,interval_right_index,stack_pos,Size);
    }
  if (ladder_start_index>=stack_rec.Q_pos) return Size;
  LBoundIdx=interval_left_index; 
  B=ENDS[interval_left_index].x;   E=ENDS[interval_right_index].x;
  if(father_first_step+inherit_each<stack_rec.Q_pos)
    Size=optMerge(ladder_start_index,stack_rec.Q_pos,Size);
  else // use fast variant if new staircase is small 
    Size=Merge(ladder_start_index,stack_rec.Q_pos,Size);
    return Size;
  };





void CIntersectionFinder::FreeMem()
  {
#define MY_FREE_ARR_MACRO(a) if(a){delete[] a;a=NULL;}
  MY_FREE_ARR_MACRO(Q);
  MY_FREE_ARR_MACRO(Pairs);
  MY_FREE_ARR_MACRO(L);
  MY_FREE_ARR_MACRO(R);
  MY_FREE_ARR_MACRO(ENDS);
  MY_FREE_ARR_MACRO(Loc);
  MY_FREE_ARR_MACRO(father_loc);

  };

void CIntersectionFinder::AllocMem(int4 n)
  {
  Pairs=new EndIndexes [n];
  L=new SegmentInfo  [n+1];
  R=new SegmentInfo  [n+1];
  ENDS=new TSegmEnd[2*n];
  Loc=new int4[n];
  Q=new SegmentInfo  [n];

  };

void CIntersectionFinder::optAllocMem(int4 n)
  {
  Pairs=new EndIndexes [n];
  L=new SegmentInfo  [n+1];
  R=new SegmentInfo  [n+1];
  ENDS=new TSegmEnd[2*n];
  // for optimal alg it needs additional space reserve n/(inherit_each-1)+inherit_each+1 to store inherited stairs
  Loc=new int4[n+n/(inherit_each-1)+inherit_each+1];
  Q=new SegmentInfo  [n+n/(inherit_each-1)+inherit_each+1];
  father_loc=new int4[n+n/(inherit_each-1)+inherit_each+1];
  };


CIntersectionFinder *CIntersectionFinder::clone(PRegObj robj)
  {
  CIntersectionFinder *res=new CIntersectionFinder();
  res->nTotSegm=nTotSegm;
  res->Scoll=Scoll;
  res->Pairs=Pairs;
  res->ENDS=ENDS;
  res->L=new SegmentInfo  [nTotSegm+1];
  res->L[0]=L[0];


  res->Q=new SegmentInfo  [nTotSegm];
  res->R=new SegmentInfo  [nTotSegm+1];
  res->Loc=new int4[nTotSegm];

  res->set_segm_fuctions(_Below,
    _FindAndRegIPoints,
    _FindAndRegIPointsInStripe,
    _IsIntersetInStripe,
    _BegPoint,
    _EndPoint,
    _Under,
    robj,
    is_line_segments);
  return res;       
  };

void CIntersectionFinder::unclone()
  {

  Pairs=NULL;
  ENDS=NULL;

  MY_FREE_ARR_MACRO(Q);
  MY_FREE_ARR_MACRO(L)
  MY_FREE_ARR_MACRO(R);
  MY_FREE_ARR_MACRO(Loc);

  Scoll=NULL; 
  nTotSegm=0;
  };

void CIntersectionFinder::prepare_ends_and_pairs(int4 n)
  {
  int4 i,j,k;
  REAL x,y;
  for (i=0,j=0;i<n;i++)
    {
    ENDS[j].s=i;
//    GETPSEG(ENDS[j].s)=Scoll[i];
//    ENDS[j].s.end_indexes=Pairs+i;
    _BegPoint(Scoll[i],ENDS[j].x,y);
    //      _BegPoint(Scoll[i],ENDS[j].x,ENDS[j].y);
    ENDS[j].islast=0;
    j++;
    _EndPoint(Scoll[i],ENDS[j].x,y);
    //      _EndPoint(Scoll[i],ENDS[j].x,ENDS[j].y);
    ENDS[j].s=i;
//    GETPSEG(ENDS[j].s)=Scoll[i];
//    ENDS[j].s.end_indexes=Pairs+i;
    ENDS[j].islast=1;
    j++;
    }
  fastsort(ENDS,j);
  for (int4 i=0;i<2*n;i++)
    if (ENDS[i].islast)
      GETEIP(ENDS[i].s) -> E=i;
    else
      GETEIP(ENDS[i].s) -> B=i;
  };

CIntersectionFinder::CIntersectionFinder()
  {
  Pairs=NULL;
  _reg_obj=NULL;
  Q=NULL;
  L=NULL;
  R=NULL;
  ENDS=NULL;
  Scoll=NULL;
  Loc=NULL;
  father_loc=NULL;
  memset(my_counter,0,sizeof(my_counter));
  };

void CIntersectionFinder::set_segm_fuctions(
  FBelow below, 
  FFindAndRegIPoints findAndRegIPoints,
  FFindAndRegIPointsInStripe findAndRegIPointsInStripe,
  FIsIntersetInStripe isIntersetInStripe,
  FGetPoint begPoint,
  FGetPoint endPoint,
  FUnder under,
  PRegObj reg_obj,
  int4 is_line)
  {
  _Below=below;
  _FindAndRegIPoints=findAndRegIPoints;
  _FindAndRegIPointsInStripe=findAndRegIPointsInStripe;
  _IsIntersetInStripe=isIntersetInStripe;
  _BegPoint=begPoint;
  _EndPoint=endPoint;
  _Under=under;
  _reg_obj=reg_obj;
  is_line_segments=is_line;
  }

void CIntersectionFinder::balaban_fast(int4 n,PSeg _Scoll[])
  {
  int4 Size;

  Scoll=_Scoll;
  nTotSegm=n;
  AllocMem(n);
  prepare_ends_and_pairs(n);

  if (ENDS[0].s==ENDS[2*n-1].s)
    {
    PSeg s=GETPSEG(ENDS[0].s);
    for(int4 i=0;i<n;i++)  if (s!=Scoll[i]) _FindAndRegIPoints(s,Scoll[i],_reg_obj);
    Size=0;
    }
  else
    {
    L[0]=ENDS[0].s; 
    Size=1;
    }
  ProgramStackRec stack_rec(-1,2*n);  //need to be initialized this way
  Size=FindR(-1,0,  n  ,&stack_rec,Size);
  Size=InsDel(n,&stack_rec,Size);
  Size=FindR(-1,n,2*n-1,&stack_rec, Size);

  FreeMem();
  }

void CIntersectionFinder::balaban_optimal(int4 n,PSeg _Scoll[])
  {
  int4 Size;

  Scoll=_Scoll;
  nTotSegm=n;
  optAllocMem(n);
  prepare_ends_and_pairs(n);

  if (ENDS[0].s==ENDS[2*n-1].s)
    {
    PSeg s=GETPSEG(ENDS[0].s);
    for(int4 i=0;i<n;i++)  if (s!=Scoll[i]) _FindAndRegIPoints(s,Scoll[i],_reg_obj);
    Size=0;
    }
  else
    {
    L[0]=ENDS[0].s; 
    Size=1;
    }

  ProgramStackRec stack_rec(inherit_each,2*n);  //need to be initialized this way  
  Size=optFindR(inherit_each+1,inherit_each,0,  n ,&stack_rec,Size);
  Size=optInsDel(n,&stack_rec,Size);
  Size=optFindR(inherit_each+1,inherit_each,n,2*n-1,&stack_rec,Size);

  FreeMem();
  }

int4 CIntersectionFinder::CalcLAt(int4 end_index)
  {
  int4 i,j,k;
  int4 Size=0;
  EndIndexes * cur_ends;
  for(i=0;i<nTotSegm;i++)
    {

    cur_ends=Pairs+i;
    if((cur_ends->B<=end_index)&&(cur_ends->E>end_index))
      L[Size++]=i;
    }

  fastsort(CSegmCompare(Scoll,_Below ,ENDS[end_index].x),L,Size);
  return Size;

  };   

void CIntersectionFinder::trivial(int4 n,PSeg sgm[])
  { 
  int4 i,j,m=n-1;
  for(i=0;i<m;i++)
    for(j=i+1;j<n;j++)
      _FindAndRegIPoints(sgm[i],sgm[j],_reg_obj);
  };


void CIntersectionFinder::simple_sweep(int4 n,PSeg Scoll[])
  {
  int4 i,j;
  PSeg s;
  int4 pos,k;

  TSEnd *Ends=new TSEnd[2*n];
  PSeg *sgm=new PSeg[n+1];
  j=0;
  for (i=0;i<n;i++)
    {
    REAL y;
    _BegPoint(Scoll[i],Ends[j].x,y);
    Ends[j].s=Scoll[i];
    Ends[j].islast=0;
    j++;
    _EndPoint(Scoll[i],Ends[j].x,y);
    Ends[j].s=Scoll[i];
    Ends[j].islast=1;
    j++;
    }
  fastsort(Ends,2*n);
  pos=0;
  for (i=0;i<2*n;i++)
    {
    s=Ends[i].s;
    if (Ends[i].islast)
      {
      j=0; sgm[pos]=s;
      while (sgm[j]!=s) j++;
      pos--;
      for (;j<pos;j++)sgm[j]=sgm[j+1];
      }
    else
      {
      for (j=0;j<pos;j++)
        _FindAndRegIPoints(s,sgm[j],_reg_obj);
      sgm[pos]=s;
      pos++;
      }
    }
  delete [] Ends;
  delete []sgm;
  }

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
DWORD WINAPI thread_routine(void *param)
  {
  ((CIntersectionFinder *)param)->Run();
  return 0;
  }    
#endif

void CIntersectionFinder::fast_parallel(int4 n,PSeg _Scoll[],PRegObj add_reg)
  {
#ifdef WIN32 //compiling for Windows CreateThread supported
  int4 Size,m;
  HANDLE hThread;
  DWORD tID;    

  Scoll=_Scoll;
  nTotSegm=n;
  AllocMem(n);
  prepare_ends_and_pairs(n);


  CIntersectionFinder *intersection_finder=clone(add_reg);
  intersection_finder->L[0]=ENDS[0].s;
  m=n*1.05;
  intersection_finder->run_to=m;
  hThread=CreateThread(NULL,0,thread_routine,(void *)intersection_finder,0,&tID);
  if(hThread) 
    {
    Size=CalcLAt(m);
    ProgramStackRec stack_rec(-1,2*nTotSegm); 
    FindR(-1,m,2*n-1,&stack_rec, Size);    
    WaitForSingleObject(hThread,INFINITE);
    CloseHandle(hThread);
    }
  intersection_finder->unclone();
  delete intersection_finder;
  FreeMem();

#else //compiling for some another enviroment just call balaban_fast. Rewrite thread creation to use parallel version 
  // somethere else (not Windows)
  balaban_fast( n, _Scoll);
#endif    
  }

