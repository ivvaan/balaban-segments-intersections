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
  if(is_line_segments&&(c!=l))return; //if found and segment is line it cant'nt be any more
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
    PSeg s=L[i].segment;
    int4 loc=Loc[i];
    if (L[i].end_indexes->E<RBoundIdx) //if segment is not covering the stripe
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
void CIntersectionFinder::optFindIntI(int4 r_index,int4 stack_pos,PSeg seg)
  {
  while(I_Stack[stack_pos].right_bound<=r_index)stack_pos--;// go from bottom to top and find staircase to start
  if(stack_pos)// if found
    {
    E=ENDS[r_index].x;// adjust right bound for segment functions helpers
    int4 l=undef_loc,r,m,QB,QE=I_Stack[stack_pos].Q_pos; 
    while(stack_pos)// for all staircases
      {
      QB=I_Stack[--stack_pos].Q_pos; 
      if(l==undef_loc)// for first or small staircase to find seg location search all range
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
        if(_Under(Q[m].segment,seg))
          l=m;
        else
          r=m;
        }
      optFindInt(QB,QE,l,seg);
      QE=QB; // move staircase bound to parent staircase
      l=father_loc[l];  //using father_loc get approximate location in parent staircase 
      };
    }
  }


// FindIntI - finds and reports intersections of the segment seg with ancestor staircases (listed in I_Stack)
// from top of the hierarchy to bottom stopping on last staircase for which seg is internal
void CIntersectionFinder::FindIntI(int4 r_index,int4 stack_pos,PSeg seg)
  {
  if(I_Stack->right_bound<=r_index)return;// if not internal for topmost staircase exit
  E=ENDS[r_index].x;// adjust right bound for segment functions helpers
  IStackRec *stack_ptr;
  int4 l,r,m,QB,QE=-1; 
  I_Stack[stack_pos+1].right_bound=r_index;//set this way circle below to stop somewhere
  for(stack_ptr=I_Stack;stack_ptr->right_bound>r_index;stack_ptr++)
    {
    l=QB=QE;
    QE=stack_ptr->Q_pos;
    r=QE+1;// always search all staircase QB,QE
    while((r-l)>1)// binary search
      {
      m=(r+l)>>1;//        m=(r+l)/2;
      //my_counter[0]++;
      if(_Under(Q[m].segment,seg))
        l=m;
      else
        r=m;
      }
    FindInt(QB,QE,l,seg);
    };
  }

// InsDel - inserts or deletes segment associated with end_index to/from L
//retuns new L size
int4 CIntersectionFinder::InsDel(int4 end_index,int4 stack_pos,int4 LSize)
  {
  int4 j,i;
  TSegmEnd *current_end=ENDS+end_index;
  PSeg seg=current_end->s.segment;
  if (current_end->islast) // if endpoint - delete
    {
    i=0; L[LSize].segment=seg;
    while(L[i].segment!=seg)i++;
    LSize--;
    for (j=i;j<LSize;j++)L[j]=L[j+1];
    }
  else// if startpoint - insert
    {
    for (i=LSize-1;(i>-1)&& (!_Under(L[i].segment,seg));i--)
      L[i+1]=L[i];
    L[i+1]=current_end->s;
    LSize++;

    FindIntI(current_end->s.end_indexes->E,stack_pos,seg);// get internal intersections
    }
    return LSize;
  }

// see InsDel comments
int4 CIntersectionFinder::optInsDel(int4 end_index,int4 stack_pos,int4 LSize)
  {
  int4 j,i;
  TSegmEnd *current_end=ENDS+end_index;
  PSeg seg=current_end->s.segment;
  if (current_end->islast)
    {
    i=0; L[LSize].segment=seg;
    while(L[i].segment!=seg)i++;
    LSize--;
    for (j=i;j<LSize;j++)L[j]=L[j+1];
    }
  else
    {
    for (i=LSize-1;(i>-1)&& (!_Under(L[i].segment,seg));i--)
      L[i+1]=L[i];
    L[i+1]=current_end->s;
    LSize++;

    optFindIntI(current_end->s.end_indexes->E,stack_pos,seg);
    }
    return LSize;
  }

// Merge merges L and current staircase QB<=stair_index<=QE on right bound of the stripe and place it in L
// returns new L size. Previous L is placed to R. While merging finds locations of R segments
// in current staircase and place them in Loc.
int4 CIntersectionFinder::Merge(int4 QB,int4 QE,int4 LSize)
  {
  int4 LB=0,new_L_size=0;
  while ((QB<=QE)&&(LB<LSize))
    if (_Below(E,L[LB].segment,Q[QB].segment))
      {
      Loc[LB]=QB-1;
      R[new_L_size++]=L[LB++];
      }
    else
      R[new_L_size++]=Q[QB++];
  while (LB<LSize)
    {
    Loc[LB]=QE; 
    R[new_L_size++]=L[LB++];
    }
  while (QB<=QE)
    R[new_L_size++]=Q[QB++];
  ExchangeLR();
  return new_L_size;
  };

// Same as Merge but merges only original stairs
int4 CIntersectionFinder::optMerge(int4 QB,int4 QE,int4 LSize)
  {
  int4 LB=0,new_L_size=0;
  while ((QB<=QE)&&(LB<LSize))
    if (_Below(E,L[LB].segment,Q[QB].segment))
      {
      Loc[LB]=QB-1;
      R[new_L_size++]=L[LB++];
      }
    else
      {
      if(IS_ORIGINAL(father_loc[QB]))
        R[new_L_size++]=Q[QB];
      QB++;
      }
    while (LB<LSize)
      {
      Loc[LB]=QE;
      R[new_L_size++]=L[LB++];
      }
    while (QB<=QE)
      {
      if(IS_ORIGINAL(father_loc[QB]))
        R[new_L_size++]=Q[QB];
      QB++;
      }
    ExchangeLR();
    return new_L_size;
  };

// Split - splits L into new staircase and new L, new L temporary placed to R, at the end L and R exchanged.
// While executing finds location of new L segments in new staicase and finds intersection of the new L segments
// covering current strip with the staircase stairs below.   
int4 CIntersectionFinder::Split(int4 &step_index,int4 LSize)
  {
  int4 father_last_step,new_L_size,j;

  father_last_step=step_index; new_L_size=0;
  for (j=0;j<LSize;j++)
    {
    if(L[j].end_indexes->E>=RBoundIdx)//segment is covering current stripe
      {
      int4 step=step_index;
      PSeg s=L[j] . segment;
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
int4 CIntersectionFinder::optSplit(int4 cur_father_pos,int4 &step_index,int4 LSize)
  {
  int4 father_last_step,new_L_size,cur_L_pos;

  father_last_step=step_index; 
  new_L_size=cur_L_pos=0;

  while ((cur_L_pos<LSize)&&(cur_father_pos<=father_last_step))
    {
    PSeg s=L[cur_L_pos] . segment;
    if(_Below(B,s,Q[cur_father_pos].segment)) // current segment below candidate to inherit
      {
      if(L[cur_L_pos].end_indexes->E>=RBoundIdx) // current segment covering strip
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
  for(;cur_L_pos<LSize;cur_L_pos++)  // father stair case is over. we are adding rest segments the same manner as in suboptimal case
    if(L[cur_L_pos].end_indexes->E>=RBoundIdx)//segment is covering current stripe
      {
      if((step_index>father_last_step)&&(IsIntersectInCurStripe(Q[step_index],L[cur_L_pos] . segment)))// if it intesects last stair
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
void CIntersectionFinder::SearchInStrip(int4 QP,int4 lsize)
  {
  int4 NQP,*stack_ptr=SIS_Stack;
  int4 LB,NL,QB;
  *stack_ptr=NQP=QP;
  while (lsize=Split(NQP,lsize))
    {
    *(++stack_ptr)=NQP;
    FindIntL(QP,NQP,lsize);
    QP=NQP;
    }
  // at this point we can just sort Q starting from *stack_ptr and place result to L, but we use sequential merging;  
  ExchangeLR();
  lsize=NQP-QP;
  NQP=QP;
  while (stack_ptr!=SIS_Stack)
    {
    QP=*(--stack_ptr);
    //merging Q and L
    QB=QP+1;
    LB=NL=0;
    while ((QB<=NQP)&&(LB<lsize))
      R[NL++]=_Below(E,L[LB].segment,Q[QB].segment)?L[LB++]:Q[QB++];
    while (LB<lsize)
      R[NL++]=L[LB++];
    while (QB<=NQP)
      R[NL++]=Q[QB++];
    ExchangeLR();
    lsize=NL;
    NQP=QP;
    }
  }

// main function to find intersections for fast algorithm
int4 CIntersectionFinder::FindR( int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,
                                int4 stack_pos, int4 LSize)
  {
  int4 m,ladder_cur_index,RSize;
  LBoundIdx=interval_left_index; RBoundIdx=interval_right_index;
  B=ENDS[interval_left_index].x;E=ENDS[interval_right_index].x; 
  if (interval_right_index-interval_left_index==1)
    {
    if (LSize>1)
      SearchInStrip(ladder_start_index,LSize);
    return LSize;
    };
  ladder_cur_index=ladder_start_index;
  int_numb=0;
  if (LSize!=0)
    {
    LSize=Split(ladder_cur_index,LSize);
    if((ladder_start_index<ladder_cur_index))
      {
      FindIntL(ladder_start_index,ladder_cur_index,LSize);
      CheckIStack(++stack_pos); 
      I_Stack[stack_pos].Set(ladder_cur_index,interval_right_index);
      };
    };
  ladder_start_index++;
  if (int_numb>LSize)
    LSize=FindR(ladder_cur_index,interval_left_index,interval_right_index,stack_pos,LSize);
  else
    {
    m=(interval_left_index+interval_right_index)/2; 
    LSize=FindR(ladder_cur_index,interval_left_index,m,stack_pos,LSize);
    LSize=InsDel(m,stack_pos,LSize);
    LSize=FindR(ladder_cur_index,m,interval_right_index,stack_pos,LSize);
    }
  if (ladder_start_index>ladder_cur_index) return LSize;
  LBoundIdx=interval_left_index; RBoundIdx=interval_right_index;
  B=ENDS[interval_left_index].x;   E=ENDS[interval_right_index].x;
  LSize=Merge(ladder_start_index--,ladder_cur_index,RSize=LSize);
  for (int4 i=0;i<RSize;i++)   if (R[i] .end_indexes->B>interval_left_index) 
    FindInt(ladder_start_index,ladder_cur_index,Loc[i],R[i].segment);
  return LSize;
  };


// main function to find intersections for optimal algorithm
int4 CIntersectionFinder::optFindR(int4 father_first_step, int4 ladder_start_index,int4 interval_left_index,int4 interval_right_index,
                                   int4 stack_pos,int4 LSize)
  {
  int4 m,ladder_cur_index,RSize;

  LBoundIdx=interval_left_index; RBoundIdx=interval_right_index;
  B=ENDS[interval_left_index].x;E=ENDS[interval_right_index].x; 
  if (interval_right_index-interval_left_index==1)
    {
    if (LSize>1)
      SearchInStrip(ladder_start_index,LSize);
    return LSize;
    };
  ladder_cur_index=ladder_start_index;
  int_numb=0;
  if(father_first_step+inherit_each<ladder_cur_index)
    {// if new staircase is big use optimal variant
    LSize=optSplit(father_first_step,ladder_cur_index,LSize);
    for (int4 i=0;i<LSize;i++)  
      optFindInt(ladder_start_index,ladder_cur_index,Loc[i],L[i].segment);
    }
  else
    {// if new staircase is small use fast variant
    LSize=Split(ladder_cur_index,LSize);
    if((ladder_start_index<ladder_cur_index))
      FindIntL(ladder_start_index,ladder_cur_index,LSize);
    for (int4 i=ladder_start_index+1;i<=ladder_cur_index;i++)
      father_loc[i]=undef_loc;
    }           
  if((ladder_start_index<ladder_cur_index))
    {
    CheckIStack(++stack_pos); 
    I_Stack[stack_pos].Set(ladder_cur_index,interval_right_index);
    }
  ladder_start_index++;     
  if (int_numb>LSize+1)
    LSize=optFindR(ladder_start_index, ladder_cur_index,interval_left_index,interval_right_index,stack_pos,LSize);
  else
    {
    m=(interval_left_index+interval_right_index)/2; 
    LSize=optFindR(ladder_start_index,ladder_cur_index,interval_left_index,m,stack_pos,LSize);
    LSize=optInsDel(m,stack_pos,LSize);
    LSize=optFindR(ladder_start_index,ladder_cur_index,m,interval_right_index,stack_pos,LSize);
    }
  if (ladder_start_index>ladder_cur_index) return LSize;
  LBoundIdx=interval_left_index; RBoundIdx=interval_right_index;
  B=ENDS[interval_left_index].x;   E=ENDS[interval_right_index].x;
  if(father_first_step+inherit_each<ladder_cur_index)
    {// if new staircase is big use optimal variant
    LSize=optMerge(ladder_start_index--,ladder_cur_index,RSize=LSize);
    for (int4 i=0;i<RSize;i++)   if (R[i] .end_indexes->B>interval_left_index) 
      optFindInt(ladder_start_index,ladder_cur_index,Loc[i],R[i].segment);
    }
  else
    { // if new staircase is small use fast variant
    LSize=Merge(ladder_start_index--,ladder_cur_index,RSize=LSize);
    for (int4 i=0;i<RSize;i++)   if (R[i] .end_indexes->B>interval_left_index) 
      FindInt(ladder_start_index,ladder_cur_index,Loc[i],R[i].segment);
    }
    return LSize;
  };





void CIntersectionFinder::FreeMem()
  {
#define MY_FREE_ARR_MACRO(a) if(a){delete[] a;a=NULL;}
  MY_FREE_ARR_MACRO(Q);
  MY_FREE_ARR_MACRO(Pairs);
  MY_FREE_ARR_MACRO(L);
  MY_FREE_ARR_MACRO(R);
  MY_FREE_ARR_MACRO(ENDS);
  MY_FREE_ARR_MACRO(SIS_Stack);
  MY_FREE_ARR_MACRO(I_Stack);
  MY_FREE_ARR_MACRO(Loc);
  MY_FREE_ARR_MACRO(father_loc);

  };

void CIntersectionFinder::AllocMem(int4 n)
  {
  Pairs=new EndIndexes [n];
  L=new SegmentInfo  [n+1];
  R=new SegmentInfo  [n+1];
  ENDS=new TSegmEnd[2*n];
  I_Stack=new IStackRec[I_Stack_size=128+n/100];
  SIS_Stack=new int4[n];
  Loc=new int4[n];
  Q=new SegmentInfo  [n];

  };

void CIntersectionFinder::optAllocMem(int4 n)
  {
  Pairs=new EndIndexes [n];
  L=new SegmentInfo  [n+1];
  R=new SegmentInfo  [n+1];
  ENDS=new TSegmEnd[2*n];
  SIS_Stack=new int4[n];
  I_Stack=new IStackRec[I_Stack_size=128+n/100];
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
  res->SIS_Stack=new int4[nTotSegm]; 
  res->I_Stack=new IStackRec[res->I_Stack_size=128+nTotSegm/100];

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
  MY_FREE_ARR_MACRO(SIS_Stack);
  MY_FREE_ARR_MACRO(I_Stack);

  Scoll=NULL; 
  nTotSegm=0;
  };

void CIntersectionFinder::prepare_ends_and_pairs(int4 n)
  {
  int4 i,j,k;
  REAL x,y;
  for (i=0,j=0;i<n;i++)
    {
    ENDS[j].s.segment=Scoll[i];
    ENDS[j].s.end_indexes=Pairs+i;
    _BegPoint(Scoll[i],ENDS[j].x,y);
    //      _BegPoint(Scoll[i],ENDS[j].x,ENDS[j].y);
    ENDS[j].islast=0;
    j++;
    _EndPoint(Scoll[i],ENDS[j].x,y);
    //      _EndPoint(Scoll[i],ENDS[j].x,ENDS[j].y);
    ENDS[j].s.segment=Scoll[i];
    ENDS[j].s.end_indexes=Pairs+i;
    ENDS[j].islast=1;
    j++;
    }
  fastsort(ENDS,j);
  for (int4 i=0;i<2*n;i++)
    if (ENDS[i].islast)
      ENDS[i].s.end_indexes -> E=i;
    else
      ENDS[i].s.end_indexes -> B=i;
  };

CIntersectionFinder::CIntersectionFinder()
  {
  Pairs=NULL;
  SIS_Stack=NULL;
  _reg_obj=NULL;
  Q=NULL;
  L=NULL;
  R=NULL;
  ENDS=NULL;
  Scoll=NULL;
  Loc=NULL;
  father_loc=NULL;
  I_Stack=NULL;
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
  int4 LSize;

  Scoll=_Scoll;
  nTotSegm=n;
  AllocMem(n);
  prepare_ends_and_pairs(n);

  if (ENDS[0].s.segment==ENDS[2*n-1].s.segment)
    {
    PSeg s=ENDS[0].s.segment;
    for(int4 i=0;i<n;i++)  if (s!=Scoll[i]) _FindAndRegIPoints(s,Scoll[i],_reg_obj);
    LSize=0;
    }
  else
    {
    L[0]=ENDS[0].s; 
    LSize=1;
    }

  LSize=FindR(-1,0,  n  ,-1,LSize);
  LSize=InsDel(n,0,LSize);
  LSize=FindR(-1,n,2*n-1,-1, LSize);

  FreeMem();
  }

void CIntersectionFinder::balaban_optimal(int4 n,PSeg _Scoll[])
  {
  int4 LSize;

  Scoll=_Scoll;
  nTotSegm=n;
  optAllocMem(n);
  prepare_ends_and_pairs(n);

  if (ENDS[0].s.segment==ENDS[2*n-1].s.segment)
    {
    PSeg s=ENDS[0].s.segment;
    for(int4 i=0;i<n;i++)  if (s!=Scoll[i]) _FindAndRegIPoints(s,Scoll[i],_reg_obj);
    LSize=0;
    }
  else
    {
    L[0]=ENDS[0].s; 
    LSize=1;
    }

  I_Stack[0].Set(inherit_each,2*n ); //need to be initialized this way
  LSize=optFindR(inherit_each+1,inherit_each,0,  n ,0,LSize);
  LSize=optInsDel(n,0,LSize);
  LSize=optFindR(inherit_each+1,inherit_each,n,2*n-1,0,LSize);

  FreeMem();
  }

int4 CIntersectionFinder::CalcLAt(int4 end_index)
  {
  int4 i,j,k;
  int4 LSize=0;
  EndIndexes * cur_ends;
  for(i=0;i<nTotSegm;i++)
    {

    cur_ends=Pairs+i;
    if((cur_ends->B<=end_index)&&(cur_ends->E>end_index))
      {
      L[LSize].segment=Scoll[i];
      L[LSize].end_indexes=cur_ends;
      LSize++;
      }
    }

  fastsort(CSegmCompare(_Below ,ENDS[end_index].x),L,LSize);
  return LSize;

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
  int4 LSize,m;
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
    LSize=CalcLAt(m);
    FindR(-1,m,2*n-1,-1, LSize);    
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

