/*
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

#include "utils.h"
#include "intersection_finder.h"
//#include <stdio.h>
#include <crtdbg.h>
#include <algorithm>
//#include <limits>


#define IS_ORIGINAL(v) (v<1)
#define INC_POS(p) p++;if(p>=N)p-=N;
#define DEC_POS(p) p--;if(p<0)p+=N;
#define MOVE_POS(p,m) p+=m;if(p>=N)p-=N;
#define POS_DIFF(p1,p2) (p1>=p2?p1-p2:p1+N-p2);


// FindInt finds and reports intersection of segment whith staircase QB<stair_index<=QE, given location of the segment begin point - l
template<bool is_line_seg> 
void CIntersectionFinder<is_line_seg>::FindInt(int4 QB, int4 QE, int4 l, PSeg s)
{
	int4 c = l;
	while ((c > QB) && (IntersectionsInCurStripe(Q[c], s))) //first get intersections below
		--c;
	if (is_line_seg && (c != l))return; //if found and segment is line it can't be any more
	c = l + 1;
	while ((c <= QE) && (IntersectionsInCurStripe(Q[c], s))) // get intersections above
		++c;
};

// no_ip case FindInt finds and reports intersection of segment whith staircase QB<stair_index<=QE, given location of the segment begin point - l
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::no_ipFindInt(int4 QB, int4 QE, int4 l, PSeg s)
{
	int4 c = l;
	while ((c>QB) && (IntersectionsInCurStripe(no_ipQ[c].s, s))) //first get intersections below
		c--;
	if (c != l)return; //if found  it can't be any more
	c = l + 1;
	while ((c <= QE) && (IntersectionsInCurStripe(no_ipQ[c].s, s))) // get intersections above
		c++;
};

// optFindInt finds and reports intersection of segment whith staircase QB<stair_index<=QE, given location of the segment begin point - l
// for optimal algorithm. Intersections whith inherited stairs are not reported.
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::optFindInt(int4 QB, int4 QE, int4 l, PSeg s)
{
	int4 c = l;
	while (c>QB)
	{
		if (IS_ORIGINAL(father_loc[c]))
		{
			if (!IntersectionsInCurStripe(Q[c], s)) break;
		}
		else //if inherited stair just check intersection
			if (!IsIntersectInCurStripe(Q[c], s)) break;
		c--;
	}
	int_numb += l - c;
	if (is_line_seg && (c != l))return;
	c = ++l;
	while (c <= QE)
	{
		if (IS_ORIGINAL(father_loc[c]))
		{
			if (!IntersectionsInCurStripe(Q[c], s)) break;
		}
		else
			if (!IsIntersectInCurStripe(Q[c], s)) break;
		c++;
	}
	int_numb += c - l;
};

// FindIntL finds and reports intersection of L with current staircase QB<stair_index<=QE for fast algorithm
// locations of the segments should be placed to R array before calling the function
/*
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::FindIntL(int4 qe, int4 segm_numb)
{
	auto Q_tail = Q + len_of_Q;
	for (int4 i = 0; i<segm_numb; ++i)
	{
		PSeg s = Scoll[L[i]]; 
		int4 loc = *(--Q_tail) + 1; 
		int4 c=loc;
		while ((c <= qe) && (IntersectionsInCurStripe(Q[c++], s))) ;
		int_numb += c - loc;
	}
}
*/


template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::ring_bufFindIntL(int4 QB, int4 QE, int4 segm_numb)
{
	auto LQ = Q + nTotSegm;
	register int4 _ringb_end = ringb_end;
	register int4 ringb_pos = ringb_beg;
	register int4 N = nTotSegm;
	auto _L = L;

	while (ringb_pos!= _ringb_end)
	{
		int4 s_idx = _L[ringb_pos];
		PSeg s = Scoll[s_idx];
		int4 loc = *(--LQ);
		if (SegmentBE[s_idx].E<RBoundIdx) //if segment is not covering the stripe
			FindInt(QB, QE, loc, s);// find intersections
		else // segment is covering the stripe
			if (!is_line_seg)// for line segments all intersections already found 
			{// if segment is not line then find intersections above (below found in Split)
				int4 c = ++loc;
				while ((c <= QE) && (IntersectionsInCurStripe(Q[c], s)))  c++;
				int_numb += c - loc;
			};
		INC_POS(ringb_pos);
	}

}


// FindIntL finds and reports intersection of L with current staircase QB<stair_index<=QE for fast algorithm
// locations of the segments should be placed to Loc array before calling the function
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::no_ipFindIntL(int4 QB, int4 QE, int4 segm_numb)
{
	for (int4 i = 0; i<segm_numb; i++)
	{
		int4 s_idx = no_ipL[i].s;
		if (SegmentBE[s_idx].E<RBoundIdx) //if segment is not covering the stripe
			no_ipFindInt(QB, QE, Loc[i], Scoll[s_idx]);// find intersections
	}
}

// optFindIntI - finds and reports intersections of the segment seg with ancestor staircases (listed in program stack)
// from bottom of the hierarchy to top starting from first staircase for which seg is internal
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::optFindIntI(uint4 r_index, ProgramStackRec *stack_pos, PSeg seg)
{
    while (stack_pos->right_bound <= r_index)stack_pos = stack_pos->prev;// go from bottom to top and find staircase to start
    if (!stack_pos->prev)return;
    E = ENDS[r_index].x;// adjust right bound for segment functions helpers
    int4  l = undef_loc, r, m,  QE = stack_pos->Q_pos + 1;
    int4 QB = undef_loc+1;//just to differ from l at first loop
    for (stack_pos = stack_pos->prev; stack_pos; stack_pos = stack_pos->prev)// for all staircases above
    {
        if ((big_staircase_threshold>inherit_each) && (l == QB)) //if location is below first stair we don't have data in father_loc for that location
            // because the number of locations is greater (by one) than number of stairs 
            // (for example one stair has two locations (below and above), but we keep in father_loc data for only one (above))
            // so we use data from next location QE==l+1
        {
            l = QB = stack_pos->Q_pos;// set lower location bound to the location below first step of current staircase
            r = father_loc[QE] ? QB + inherit_each : QE ;// father_loc[QE] !=undef_loc (i.e. 0) means staircase was created by optSplit 
            //and we use range [QB;QB+inherit_each] to search, otherwise it was created by Split and full range search [QB;QE] is used
        }
        else
        {
            QB = stack_pos->Q_pos;
            m=abs(father_loc[l]);//using father_loc get approximate location in parent staircase; 
            // line above sometimes executed when l == undef_loc, so we must initialize father_loc[undef_loc]=undef_loc in AllocMem()
            //to keep l unchanged in this case; undef_loc was chosen to be zero so abs(undef_loc)==undef_loc
            //otherwise we need additional checks
            l=max(QB,m);//here we use a fact that always undef_loc < QB and if m==undef_loc or m<QB l should be QB
            if (m) { r = l + inherit_each; if (r > QE)r = QE; }
            else r = QE;//if m==0 i.e. undef_loc we use [QB;QE] as a range
        }

        while ((r - l)>1)// binary search
        {
            m = (r + l) >> 1;//        m=(r+l)/2;
            if (_Under(Scoll[Q[m]], seg))
                l = m;
            else
                r = m;
        }

        optFindInt(QB, QE - 1, l, seg);
        QE = QB + 1; //set upper location bound to the location above last step of prev (which will be current) staircase
    };
 
}


// FindIntI - finds and reports intersections of the segment seg with ancestor staircases (listed in program stack)
// from top of the hierarchy to bottom stopping on last staircase for which seg is internal
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::FindIntI(uint4 r_index, ProgramStackRec * stack_pos, PSeg seg)
{
	while (stack_pos->right_bound <= r_index)stack_pos = stack_pos->prev;// go from bottom to top and find staircase to start
	E = ENDS[r_index].x;// adjust right bound for segment functions helpers
	int4 l, r, m, QB, QE = stack_pos->Q_pos;
	for (stack_pos = stack_pos->prev; stack_pos; stack_pos = stack_pos->prev)
	{
		l = QB = stack_pos->Q_pos;
		r = QE + 1;
		while ((r - l)>1)// binary search
		{
			m = (r + l) >> 1;//        m=(r+l)/2;
			if (_Under(Scoll[Q[m]], seg))
				l = m;
			else
				r = m;
		}
		FindInt(QB, QE, l, seg);
		QE = QB; // move staircase bound to parent staircase
	};
}//*/

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::no_ipFindIntI(uint4 r_index, ProgramStackRec * stack_pos, PSeg seg)
{
	while (stack_pos->right_bound <= r_index)stack_pos = stack_pos->prev;// go from bottom to top and find staircase to start
	E = ENDS[r_index].x;// adjust right bound for segment functions helpers
	int4 l, r, m, QB, QE = stack_pos->Q_pos;
	for (stack_pos = stack_pos->prev; stack_pos; stack_pos = stack_pos->prev)
	{
		l = QB = stack_pos->Q_pos;
		r = QE + 1;
		while ((r - l)>1)// binary search
		{
			m = (r + l) >> 1;//        m=(r+l)/2;
			if (_Under(Scoll[no_ipQ[m].s], seg))
				l = m;
			else
				r = m;
		}
		no_ipFindInt(QB, QE, l, seg);
		QE = QB; // move staircase bound to parent staircase
	};
}//*/



// InsDel - inserts or deletes segment associated with end_index to/from L
//retuns new L size

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::InsDel(uint4 end_index, ProgramStackRec * stack_pos, int4 Size)
{
	int4 i;
    auto sn = ENDS[end_index].s();
	if (ENDS[end_index].islast()) // if endpoint - delete
	{
		i = --Size;
		auto cur = L[i];
		while (cur != sn)
		{
			--i;
			auto buf = cur;
			cur = L[i];
			L[i] = buf;
		}
	}
	else// if startpoint - insert
	{
 		PSeg seg = Scoll[sn];
        for (i = Size - 1; (i > -1) && (!_Under(Scoll[L[i]], seg)); --i)
			L[i + 1] = L[i];
		L[i + 1] = sn;
		Size++;

		FindIntI(SegmentBE[sn].E, stack_pos, seg);// get internal intersections
	}
	return Size;
}



template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::msInsDel(uint4 end_index, ProgramStackRec * stack_pos, int4 Size)
{
	int4 i;
	auto sn = ENDS[end_index].s();
	if (from_begin)
	{
            if (ENDS[end_index].islast()) // if endpoint - delete
		{
			i = --Size;
			auto cur = L[i];
			while (cur != sn)
			{
				--i;
				auto buf = cur;
				cur = L[i];
				L[i] = buf;
			}
		}
		else// if startpoint - insert
		{
			PSeg seg = Scoll[sn];
			for (i = Size - 1; (i>-1) && (!_Under(Scoll[L[i]], seg)); --i)
				L[i + 1] = L[i];
			L[i + 1] = sn;
			Size++;

			FindIntI(SegmentBE[sn].E, stack_pos, seg);// get internal intersections
		}

	}
	else
	{
        if (ENDS[end_index].islast()) // if endpoint - delete
		{
			i = nTotSegm - Size;
			--Size;
			auto cur = L[i];
			while (cur != sn)
			{
				++i;
				auto buf = cur;
				cur = L[i];
				L[i] = buf;
			}
		}
		else// if startpoint - insert
		{
			PSeg seg = Scoll[sn];
            
			for (i = nTotSegm - Size; (i<nTotSegm) && (_Under(Scoll[L[i]], seg)); ++i)
				L[i - 1] = L[i];
			L[i - 1] = sn;
            ++Size;

			FindIntI(SegmentBE[sn].E, stack_pos, seg);// get internal intersections
		}

	}
	return Size;

}

/*template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::InsDelFT(uint4 end_index, ProgramStackRec * stack_pos, int4 Size)
{
	int4 i;
	TSegmEnd *current_end = ENDS + end_index;
	if (current_end->islast()) // if endpoint - delete
	{
		auto sn = current_end->s();
		i = nTotSegm - Size;
		--Size;
		auto cur = L[i];
		while (cur != sn)
		{
			++i;
			auto buf = cur;
			cur = L[i];
			L[i] = buf;
		}
	}
	else// if startpoint - insert
	{
		PSeg seg = Scoll[current_end->s()];
		for (i = nTotSegm - Size; (i<nTotSegm) && (_Under(Scoll[L[i]], seg)); ++i)
			L[i - 1] = L[i];
		L[i - 1] = current_end->s();
		++Size;

		FindIntI(SegmentBE[current_end->s()].E, stack_pos, seg);// get internal intersections
	}
	return Size;
}
*/

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::ring_bufInsDel(uint4 end_index, ProgramStackRec * stack_pos, int4 Size)
{
	register int4 new_ringb_pos;
	register int4 N = nTotSegm;
	register int4 _ringb_end ;
	register int4 ringb_pos;	
	auto _L = L;
	int4 sn;

	TSegmEnd *current_end = ENDS + end_index;
	if (current_end->islast()) // if endpoint - delete
	{
		_ringb_end = ringb_end;
		ringb_pos = ringb_beg;		
		sn =current_end->s();
		//L[_ringb_end] = sn;
		while (_L[ringb_pos] != sn) { INC_POS(ringb_pos); };
		DEC_POS(ringb_end);
		new_ringb_pos = ringb_pos;
		INC_POS(new_ringb_pos);
		
		while (new_ringb_pos != _ringb_end) 
		{ 
			_L[ringb_pos] = _L[new_ringb_pos]; 
			ringb_pos = new_ringb_pos;
			INC_POS(new_ringb_pos);
		};
		return Size - 1;
	}
	_ringb_end = ringb_beg;
	DEC_POS(_ringb_end);
	ringb_pos = ringb_end;
	new_ringb_pos = ringb_pos;
	DEC_POS(new_ringb_pos);
	PSeg seg = Scoll[current_end->s()];
	sn = _L[new_ringb_pos];
	while ((new_ringb_pos != _ringb_end) && (!_Under(Scoll[sn], seg)))
	{
		_L[ringb_pos] = sn;
		ringb_pos = new_ringb_pos;
		DEC_POS(new_ringb_pos);
		sn = _L[new_ringb_pos];
	}
			
	_L[ringb_pos] = current_end->s();
	INC_POS(ringb_end);

	FindIntI(SegmentBE[current_end->s()].E, stack_pos, seg);// get internal intersections
	return Size + 1;
}


// InsDel - inserts or deletes segment associated with end_index to/from L
//retuns new L size
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::no_ipInsDel(uint4 end_index, ProgramStackRec * stack_pos, int4 Size)
{
	int4 i;
	TSegmEnd *current_end = ENDS + end_index;
	no_ipExchangeLR();
	if (current_end->islast()) // if endpoint - delete
	{

		int4 sn = no_ipL[Size].s = current_end->s();
		for (i = 0; no_ipL[i].s != sn; i++);
		//_ASSERTE(i!=Size);
		Size--;
		for (; i<Size; i++)no_ipL[i] = no_ipL[i + 1];
	}
	else// if startpoint - insert
	{
		PSeg seg = Scoll[current_end->s()];
		REAL X, Y;
		_BegPoint(seg, X, Y);
		for (i = Size - 1; (i>-1) && (no_ipL[i].y >= Y); i--)
			no_ipL[i + 1] = no_ipL[i];
		no_ipL[i + 1].s = current_end->s();
		Size++;

		no_ipFindIntI(SegmentBE[current_end->s()].E, stack_pos, seg);// get internal intersections
	}
	return Size;
}

// see InsDel comments
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::optInsDel(uint4 end_index, ProgramStackRec *stack_pos, int4 Size)
{
	int4 i;
	TSegmEnd *current_end = ENDS + end_index;
//	ExchangeLR();
	if (current_end->islast())
	{
		int4 sn = L[Size] = current_end->s();
		for (i = 0; L[i] != sn; i++);
		//_ASSERTE(i!=Size);
		Size--;
		for (; i<Size; i++)L[i] = L[i + 1];
	}
	else
	{
		PSeg seg = Scoll[current_end->s()];
		for (i = Size - 1; (i>-1) && (!_Under(Scoll[L[i]], seg)); i--)
			L[i + 1] = L[i];
		L[i + 1] = current_end->s();
		Size++;

		optFindIntI(SegmentBE[current_end->s()].E, stack_pos, seg);
	}
	return Size;
}

// Merge merges R and current staircase QB<stair_index<=QE on right bound of the stripe and place it in L, 
// than exchange it back to R and returns new R size.  While merging finds locations of R segments
// in current staircase and appropriate intersections.

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::Merge(uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
{
	int4 cur_R_pos = 0, new_size = 0, cur_stair = QB;
    PSeg s;
    ExchangeLR();
    int4 cur_seg = R[cur_R_pos];
	while ((cur_stair<QE) && (cur_R_pos<Size))
	{
        
        s = Scoll[cur_seg];
		if (_Below(E, s, Scoll[Q[cur_stair + 1]]))
		{
			if (SegmentBE[cur_seg].B>LBoundIdx) 
				FindInt(QB, QE, cur_stair, s);
			L[new_size++] = cur_seg;
            cur_seg = R[++cur_R_pos];
		}
		else
			L[new_size++] = Q[++cur_stair];
	}
	while (cur_R_pos<Size)
	{
		if (SegmentBE[cur_seg].B>LBoundIdx) FindInt(QB, QE, QE, Scoll[cur_seg]);
		L[new_size++] = cur_seg;
        cur_seg = R[++cur_R_pos];
	}
	while (cur_stair<QE)
		L[new_size++] = Q[++cur_stair];
	return new_size;
};


template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::msMerge(uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
{
	int4 cur_R_pos = 0, new_size = 0;
	int4 cur_seg;
	PSeg s;
	if (from_begin)
	{
		auto _R = L + Size - 1;
		auto _L = L + (nTotSegm - 1);
        auto cur_stair = QE, _Size = -Size;
        cur_seg = _R[cur_R_pos];
		while ((cur_stair>QB) && (cur_R_pos>_Size))
		{
			s = Scoll[cur_seg];
			if (_Below(E, Scoll[Q[cur_stair]], s))
			{
				if (SegmentBE[cur_seg].B>LBoundIdx)
					FindInt(QB, QE, cur_stair, s);
				_L[new_size--] = cur_seg;
                cur_seg = _R[--cur_R_pos];
			}
			else
				_L[new_size--] = Q[cur_stair--];
		}
		while (cur_R_pos>_Size)
		{
			if (SegmentBE[cur_seg].B>LBoundIdx) FindInt(QB, QE, QB, Scoll[cur_seg]);
			_L[new_size--] = cur_seg;
            cur_seg = _R[--cur_R_pos];
		}
		while (cur_stair>QB)
			_L[new_size--] = Q[cur_stair--];
		from_begin = false;
		return -new_size;
	}


	auto _R = L + (nTotSegm - Size);
	int4 cur_stair = QB;
    cur_seg = _R[cur_R_pos];
	while ((cur_stair<QE) && (cur_R_pos<Size))
	{
		s = Scoll[cur_seg];
		if (_Below(E, s, Scoll[Q[cur_stair + 1]]))
		{
			if (SegmentBE[cur_seg].B>LBoundIdx)
				FindInt(QB, QE, cur_stair, s);
			L[new_size++] = cur_seg;
            cur_seg = _R[++cur_R_pos];
		}
		else
			L[new_size++] = Q[++cur_stair];
	}
	while (cur_R_pos<Size)
	{
		if (SegmentBE[cur_seg].B>LBoundIdx) FindInt(QB, QE, QE, Scoll[cur_seg]);
		L[new_size++] = cur_seg;
        cur_seg = _R[++cur_R_pos];
	}
	while (cur_stair<QE)
		L[new_size++] = Q[++cur_stair];
	from_begin = true;
	return new_size;
};


template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::ring_bufMerge(uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
{
	int4 cur_R_pos = 0, new_size = 0;
	int4 cur_stair = QB;
	int4 cur_seg;
	PSeg s;
	register int4 _ringb_end = ringb_end;
	register int4 ringb_pos = ringb_beg;
	register int4 new_ringb_pos = _ringb_end;
	register int4 N = nTotSegm;
	auto _L = L;
	while ((cur_stair<QE) && (ringb_pos!=_ringb_end))
	{
		cur_seg = _L[ringb_pos];
		s = Scoll[cur_seg];
		if (_Below(E, s, Scoll[Q[cur_stair + 1]]))
		{
			if (SegmentBE[cur_seg].B>LBoundIdx) FindInt(QB, QE, cur_stair, s);
			_L[new_ringb_pos] = cur_seg;
			INC_POS(ringb_pos);
		}
		else
			_L[new_ringb_pos] = Q[++cur_stair];
		INC_POS(new_ringb_pos)
	}
	while (ringb_pos != _ringb_end)
	{
		cur_seg = _L[ringb_pos];
		if (SegmentBE[cur_seg].B>LBoundIdx) FindInt(QB, QE, QE, Scoll[cur_seg]);
		_L[new_ringb_pos] = cur_seg;
		INC_POS(ringb_pos);
		INC_POS(new_ringb_pos);
	}
	while (cur_stair < QE)
	{
		_L[new_ringb_pos] = Q[++cur_stair];
		INC_POS(new_ringb_pos);
	}
	new_size = POS_DIFF(new_ringb_pos, _ringb_end);
	ringb_beg = ringb_end;
	ringb_end = new_ringb_pos;
	return new_size;
};



// Merge merges R and current staircase QB<stair_index<=QE on right bound of the stripe and place it in L, 
// than exchange it back to R and returns new R size.  While merging finds locations of R segments
// in current staircase and appropriate intersections.
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::no_ipMerge(uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
{
	int4 new_size = 0;
	no_ipSegmentInfo *cur_R_pos = no_ipR, *last_R_pos = no_ipR + Size;
	int4 cur_stair = QB;
	while ((cur_stair<QE) && (cur_R_pos<last_R_pos))
	{
		if (cur_R_pos->y<no_ipQ[cur_stair + 1].y)
		{
			if (SegmentBE[cur_R_pos->s].B>LBoundIdx) no_ipFindInt(QB, QE, cur_stair, Scoll[cur_R_pos->s]);
			no_ipL[new_size++] = *cur_R_pos;
			cur_R_pos++;
		}
		else
			no_ipL[new_size++] = no_ipQ[++cur_stair];
	}
	while (cur_R_pos<last_R_pos)
	{
		if (SegmentBE[cur_R_pos->s].B>LBoundIdx) no_ipFindInt(QB, QE, QE, Scoll[cur_R_pos->s]);
		no_ipL[new_size++] = *cur_R_pos;
		cur_R_pos++;
	}
	while (cur_stair<QE)
		no_ipL[new_size++] = no_ipQ[++cur_stair];
	no_ipExchangeLR();
	return new_size;
};

// Same as Merge but merges only original stairs
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::optMerge(uint4 LBoundIdx, int4 QB, int4 QE, int4 Size)
{
	int4 cur_R_pos = 0, new_size = 0;
	int4 cur_stair = QB;
	int4 cur_seg;
	PSeg s;
    ExchangeLR();
	while ((cur_stair<QE) && (cur_R_pos<Size))
	{
		s = Scoll[cur_seg = R[cur_R_pos]];
		if (_Below(E, s, Scoll[Q[cur_stair + 1]]))
		{
			if (SegmentBE[cur_seg].B>LBoundIdx) optFindInt(QB, QE, cur_stair, s);
			L[new_size++] = cur_seg;
			cur_R_pos++;
		}
		else
		{
			cur_stair++;
			if (IS_ORIGINAL(father_loc[cur_stair]))
				L[new_size++] = Q[cur_stair];
		}
	}
	while (cur_R_pos<Size)
	{
		cur_seg = R[cur_R_pos++];
		if (SegmentBE[cur_seg].B>LBoundIdx) optFindInt(QB, QE, QE, Scoll[cur_seg]);
		L[new_size++] = cur_seg;
	}
	while (cur_stair<QE)
	{
		cur_stair++;
		if (IS_ORIGINAL(father_loc[cur_stair]))
			L[new_size++] = Q[cur_stair];
	}

	return new_size;
};

// Split - splits L into new staircase and new L, new L temporary placed to R, at the end L and R exchanged.
// While executing finds location of new L segments in new staicase and finds intersection of the new L segments
// covering current strip with the staircase stairs below.   


template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::Split(int4 &step_index, int4 Size)
{
	int4 father_last_step = step_index, new_L_size = 0;
	auto Q_tail = Q + len_of_Q;
  int4 *R_pos;
  if (is_line_seg)R_pos=R;
	for (int4 cur_L_pos = 0; cur_L_pos<Size; cur_L_pos++)
	{
		int4 cur_seg = L[cur_L_pos];
		int4 step = step_index;
		PSeg ps = Scoll[cur_seg];
		while ((father_last_step<step) && (IntersectionsInCurStripe(Q[step], ps)))
			step--;
		int4 n_int= step_index - step;
		
		if (n_int == 0)
		{
			if(SegmentBE[cur_seg].E >= RBoundIdx)//segment is covering current stripe and doesn't intersect ladder stairs
				Q[++step_index] = cur_seg;
			else
			{ 
				//place segment in L
				L[new_L_size++] = cur_seg;
				//storing segment position in Q_tail
				*(--Q_tail) = step_index+1;
        //storing segment number in R
        if (is_line_seg)*(R_pos++) = cur_seg;
			}
		}
		else
		{
			int_numb += n_int;
			//place segment in L
			L[new_L_size++] = cur_seg;
			//storing segment position in Q_tail, only for nonline segments
			//*(--Q_tail) = is_line_seg? INT_MAX :step_index+1;
      if (!is_line_seg) *(--Q_tail) = step_index + 1;
		}
	}
  Q_tail = Q + len_of_Q;
  if (!is_line_seg)
      R_pos = L + new_L_size;
  for (auto _R = is_line_seg ? R : L; _R < R_pos; ++_R) 
	{
		PSeg s = Scoll[*_R];             //getting segment which needs additional checks from R
		int4 loc = *(--Q_tail);          // getting ts position from tail of Q
		int4 c = loc;
		while ((c <= step_index) && (IntersectionsInCurStripe(Q[c++], s)));
		int_numb += c - loc;
	}

	return new_L_size;
};



template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::msSplit(int4 &step_index, int4 Size)
{
	auto _L = from_begin ? L : L + (nTotSegm - Size);
	int4 father_last_step = step_index, new_L_size = 0;
	auto Q_tail = Q + len_of_Q;
	for (int4 cur_L_pos = 0; cur_L_pos<Size; cur_L_pos++)
	{
		int4 cur_seg = _L[cur_L_pos];
		int4 step = step_index;
		PSeg ps = Scoll[cur_seg];
		while ((father_last_step<step) && (IntersectionsInCurStripe(Q[step], ps)))
			step--;
		int4 n_int = step_index - step;

		if (n_int == 0)
		{
			if (SegmentBE[cur_seg].E >= RBoundIdx)//segment is covering current stripe and doesn't intersect ladder stairs
				Q[++step_index] = cur_seg;
			else
			{
				//place segment in L
				L[new_L_size++] = cur_seg;
				//storing segment position in Q_tail
				*(--Q_tail) = step_index+1;
			}
		}
		else
		{
			int_numb += n_int;
			//place segment in L
			L[new_L_size++] = cur_seg;
			//storing segment position in Q_tail
			*(--Q_tail) = is_line_seg ? INT_MAX : step_index+1;
		}
	}
    Q_tail = Q + len_of_Q;
    for (int4 i = 0; i < new_L_size; ++i) 
	{
        PSeg s = Scoll[L[i]];
        int4 loc = *(--Q_tail);
        int4 c = loc;
        while ((c <= step_index) && (IntersectionsInCurStripe(Q[c++], s)));
        int_numb += c - loc;
    }

	from_begin = true;
	return new_L_size;
};


template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::ring_bufSplit(int4 &step_index, int4 Size)
{
	int4 father_last_step = step_index, new_L_size = 0;
	register int4 _ringb_end = ringb_end;
	register int4 ringb_pos = ringb_beg;
	register int4 new_ringb_pos = ringb_beg;
	register int4 N = nTotSegm;
	auto _L = L;

	auto LQ = Q + nTotSegm;
	while (ringb_pos!= _ringb_end)
	{
		int4 cur_seg = _L[ringb_pos];
		if (SegmentBE[cur_seg].E >= RBoundIdx)//segment is covering current stripe
		{
			int4 step = step_index;
			PSeg ps = Scoll[cur_seg];
			while ((father_last_step<step) && (IntersectionsInCurStripe(Q[step], ps)))
				step--;
			if (step_index != step)
			{
				int_numb += step_index - step;
				//place segment in L
				_L[new_ringb_pos] = cur_seg;
				INC_POS(new_ringb_pos);
				//storing segment position in LQ
				*(--LQ) = step_index;
			}
			else
			{
				++step_index;
				Q[step_index] = cur_seg;
			}
		}
		else
		{
			//place segment in L
			_L[new_ringb_pos] = cur_seg;
			INC_POS(new_ringb_pos);
			//storing segment position in LQ
			*(--LQ) = step_index;
		}
		INC_POS(ringb_pos);
	}
	int4 new_size = POS_DIFF(new_ringb_pos, ringb_beg);
	ringb_end = new_ringb_pos;
	return new_size;
};



// same as Split but more fast intersection checking without finding intersection point
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::no_ipSplit(int4 &step_index, int4 Size, int4 stripe_divided)
{
	int4 father_last_step = step_index, new_L_size = 0;
	no_ipSegmentInfo *cur_L_pos = no_ipL, *last_L_pos = no_ipL + Size;

	for (; cur_L_pos<last_L_pos; cur_L_pos++)
	{
		if (SegmentBE[cur_L_pos->s].E >= RBoundIdx)//segment is covering current stripe
		{
			int4 step = step_index;
			PSeg ps = Scoll[cur_L_pos->s];
			//  at the left edge of the stripe current segment must have Y coord greater then all segments of the current staircase 
			if (stripe_divided)cur_L_pos->y = _YAtX(ps, E); // Y coordinate of current segment at right edge of the stripe
			while ((father_last_step<step) && (no_ipQ[step].y>cur_L_pos->y))// current segment intersects all segments of the staircase having bigger Y coords at the stripe right edge  
			{
				_RegIntersection(_reg_obj, Scoll[no_ipQ[step].s], ps, 1, NULL); step--;
			}// intersection point remains unknown, so we register only pair passing NULL as last param to _RegIntersection
			if (step_index != step)
			{
				int_numb += step_index - step;
				no_ipR[new_L_size] = *cur_L_pos;
				Loc[new_L_size++] = step_index;
			}
			else
			{
				step_index++;
				no_ipQ[step_index] = *cur_L_pos;
			}
		}
		else
		{
			no_ipR[new_L_size].s = cur_L_pos->s;
			Loc[new_L_size++] = step_index;
		}
	}
	no_ipExchangeLR();
	return new_L_size;
};


// same as Split but  add a part of the parent staircase stairs (inherited stairs) to the new one.
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::optSplit(int4 cur_father_pos, int4 &step_index, int4 Size)
{
	int4 father_last_step = step_index, new_L_size = 0, cur_L_pos = 0;
	int4 cur_seg;
	while ((cur_L_pos<Size) && (cur_father_pos <= father_last_step))
	{
		cur_seg = L[cur_L_pos];
		PSeg s = Scoll[cur_seg];
		if (_Below(B, s, Scoll[Q[cur_father_pos]])) // current segment below candidate to inherit
		{
			if ((SegmentBE[cur_seg].E < RBoundIdx)// current segment not covering strip
				    || (step_index>father_last_step) && IsIntersectInCurStripe(Q[step_index], s)// or intersects last stair
					|| IsIntersectInCurStripe(Q[cur_father_pos], s)) //or intersects candidate to inherit 
			{
					//place segment in L
					L[new_L_size] = cur_seg;
					//storing segment position in R
					R[new_L_size++] = step_index;
			}
			else
			{ //add it to Q as original
				step_index++;
				// if segment is over this stair it is over previous inherited segment in father staircase
				// negative value of  father_loc[step_index] means it is original segment ( always cur_father_pos  >  inherit_each )              
				father_loc[step_index] = inherit_each - cur_father_pos;
				Q[step_index] = cur_seg;
			}
			cur_L_pos++; // move to next in L
		}
		else
		{  // add candidate to inherit to Q
			step_index++;
			// if segment is over this stair it is over the same inherited segment in father staircase  
			// positive value of  father_loc[step_index] means it is inherited segment                  
			father_loc[step_index] = cur_father_pos;
			Q[step_index] = Q[cur_father_pos];
			cur_father_pos += inherit_each;
		};
	}
	for (; cur_L_pos<Size; cur_L_pos++)  // father stair case is over. we are adding rest segments the same manner as in suboptimal case
	{
		cur_seg = L[cur_L_pos];
		if ((SegmentBE[cur_seg].E < RBoundIdx)||//segment is not covering current stripe
			 (step_index>father_last_step) && IsIntersectInCurStripe(Q[step_index], Scoll[cur_seg]))// if it intesects last stair
		{  
				//place segment in L
				L[new_L_size] = cur_seg;
				//storing segment position in R
				R[new_L_size++] = step_index;
		}
		else
		{  // add it to Q
			step_index++;
			// negative value of father_loc[step_index] means it is original segment ( always cur_father_pos  >  inherit_each )              
			father_loc[step_index] = inherit_each - cur_father_pos;
			Q[step_index] = cur_seg;
		}
	};
	for (; cur_father_pos <= father_last_step; cur_father_pos += inherit_each) // L segments is over. We are adding rest inherited father segments to Q
	{
		step_index++;
		// positive value of  father_loc[step_index] means it is inherited segment                  
		father_loc[step_index] = cur_father_pos;
		Q[step_index] = Q[cur_father_pos];
	};
	//  ExchangeLR();
	return new_L_size;
};

// SearchInStrip finds and reports all intersection in a stripe containing no segment ends

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::SearchInStrip(int4 QP, int4 Size)
{
    auto _L = L, q = Q + QP + 1, last = L + Size;
    int4 n_split = 0, size = Size;
    while (size = Split(QP, size))++n_split; //it change QP must be after q = Q + QP + 1 
    // at this point we can just place Q starting from QP+1 to L and sort it
    for (; _L < last; ++_L, ++q)  *_L = *q;
    if(n_split)std::sort(L, last, CSegmCompare(Scoll, _Below, E));
    return Size;
}

template <bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::msSearchInStrip(int4 QP, int4 Size)
{
    auto _L = L, q = Q + QP + 1, last = L + Size;
    int4 n_split = 0, size = Size;
    while (size = msSplit(QP, size))
        ++n_split; //it change QP must be after q = Q + QP + 1
    // at this point we can just place Q starting from QP+1 to L and sort it
    for (; _L < last; ++_L, ++q)
        *_L = *q;
    if (n_split)
        std::sort(L, last, CSegmCompare(Scoll, _Below, E));
    return Size;
}

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::ring_bufSearchInStrip(int4 QP, int4 Size)
{
	if (Size<2)  return Size; 
	int4 NQP = QP;
	int4 size = Size; int4 *q = Q + QP + 1;
	if (size = ring_bufSplit(NQP, size))
	{

		do
		{
			ring_bufFindIntL(QP, NQP, size);
			QP = NQP;
		} while (size = ring_bufSplit(NQP, size));
		// at this point we can just place Q starting from QP+1 to R and sort it
		for (QP = 0; QP<Size; QP++)L[QP] = q[QP];
		ringb_beg = 0;
		ringb_end = Size;
		std::sort(L, L + Size, CSegmCompare(Scoll, _Below, E));

	}
	else
	{
		// at this point we can just place Q starting from QP+1 to R and sort it
		for (QP = 0; QP<Size; QP++)L[QP] = q[QP];
		ringb_beg = 0;
		ringb_end = Size;
	}
	return Size;
}

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::no_ipSearchInStrip(int4 QP, int4 Size)
{
	if (!Size)return 0;
	if (Size<2)
	{
		no_ipR->y = _YAtX(Scoll[no_ipR->s = no_ipL->s], E);
		return 1;
	}
	int4 NQP = QP;
	int4 size = Size;
	if (size = no_ipSplit(NQP, size, TRUE))
	{
		no_ipSegmentInfo *q = no_ipQ + QP + 1;
		do
		{
			no_ipFindIntL(QP, NQP, size);
			QP = NQP;
		} while (size = no_ipSplit(NQP, size, FALSE));
		// at this point we can just place Q starting from QP+1 to R and sort it
		for (QP = 0; QP<Size; QP++)no_ipR[QP] = q[QP];
		std::sort(no_ipR, no_ipR+Size);
	}
	return Size;
}

// main function to find intersections for fast algorithm
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::FindR(int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index,
	ProgramStackRec *stack_pos, int4 Size, int4 call_numb)
{
	B = ENDS[interval_left_index].x; E = ENDS[RBoundIdx = interval_right_index].x;
    if (interval_right_index - interval_left_index == 1)
        return Size>1?SearchInStrip(ladder_start_index, Size):Size;
	ProgramStackRec stack_rec(ladder_start_index);
	int_numb = 0;// variable to count intersections on Split stage
	if (Size>0)
	{
		Size = Split(stack_rec.Q_pos, Size);
		if ((ladder_start_index<stack_rec.Q_pos))
			stack_pos = stack_rec.Set(stack_pos, interval_right_index);
	};
	if ((int_numb>Size) && (call_numb<max_call)) //if found a lot of intersections repeat FindR
		Size = FindR(stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
	else //cut stripe 
	{

      uint4 m = (interval_left_index + interval_right_index) / 2;
      if(call_numb>1){// if L contains a lot of segments then cut on two parts
		       Size = FindR(stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
		       Size = InsDel(m, stack_pos, Size);
		       Size = FindR(stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
      }
      else{// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
          uint4 q = (interval_left_index + m) / 2;
          if (interval_left_index != q) {
              Size = FindR(stack_rec.Q_pos, interval_left_index, q, stack_pos, Size, 0);
              Size = InsDel(q, stack_pos, Size);
          }
          if (q != m) {
              Size = FindR(stack_rec.Q_pos, q, m, stack_pos, Size, 0);
              Size = InsDel(m, stack_pos, Size);
          }
          q = (interval_right_index + m) / 2;
          if (q != m) {
              Size = FindR(stack_rec.Q_pos, m, q, stack_pos, Size, 0);
              Size = InsDel(q, stack_pos, Size);
          }
          Size = FindR(stack_rec.Q_pos, q, interval_right_index, stack_pos, Size, 0);
      }
	}
	if (ladder_start_index >= stack_rec.Q_pos) return Size;
	B = ENDS[interval_left_index].x;   E = ENDS[interval_right_index].x;
	return Merge(interval_left_index,ladder_start_index, stack_rec.Q_pos, Size);
};

template <bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::msFindR(int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index,
    ProgramStackRec* stack_pos, int4 Size, int4 call_numb)
{
    B = ENDS[interval_left_index].x;
    E = ENDS[RBoundIdx = interval_right_index].x;
    if (interval_right_index - interval_left_index == 1)
        return Size > 1 ? msSearchInStrip(ladder_start_index, Size) : Size;
    ProgramStackRec stack_rec(ladder_start_index);
    int_numb = 0; // variable to count intersections on Split stage
    if (Size > 0) {
        Size = msSplit(stack_rec.Q_pos, Size);
        if ((ladder_start_index < stack_rec.Q_pos))
            stack_pos = stack_rec.Set(stack_pos, interval_right_index);
    };
    if ((int_numb > Size) && (call_numb < max_call)) //if found a lot of intersections repeat FindR
        Size = msFindR(stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
    else //cut stripe
    {

        uint4 m = (interval_left_index + interval_right_index) / 2;
        if (Size > nTotSegm / 100) { // if L contains a lot of segments then cut on two parts
            Size = msFindR(stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
            Size = msInsDel(m, stack_pos, Size);
            Size = msFindR(stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
        } else { // if L contains not so many segments than cut on four parts (works faster for some segment distributions)
            uint4 q = (interval_left_index + m) / 2;
            if (interval_left_index != q) {
                Size = msFindR(stack_rec.Q_pos, interval_left_index, q, stack_pos, Size, 0);
                Size = msInsDel(q, stack_pos, Size);
            }
            if (q != m) {
                Size = msFindR(stack_rec.Q_pos, q, m, stack_pos, Size, 0);
                Size = msInsDel(m, stack_pos, Size);
            }
            q = (interval_right_index + m) / 2;
            if (q != m) {
                Size = msFindR(stack_rec.Q_pos, m, q, stack_pos, Size, 0);
                Size = msInsDel(q, stack_pos, Size);
            }
            Size = msFindR(stack_rec.Q_pos, q, interval_right_index, stack_pos, Size, 0);
        }
    }
    if (ladder_start_index >= stack_rec.Q_pos)
        return Size;
    B = ENDS[interval_left_index].x;
    E = ENDS[interval_right_index].x;
    return msMerge(interval_left_index, ladder_start_index, stack_rec.Q_pos, Size);
};

template <bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::ring_bufFindR(int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index,
	ProgramStackRec *stack_pos, int4 Size, int4 call_numb)
{
	//  int4 RSize;
	B = ENDS[interval_left_index].x; E = ENDS[RBoundIdx = interval_right_index].x;
	if (interval_right_index - interval_left_index == 1)
		return ring_bufSearchInStrip(ladder_start_index, Size);
	ProgramStackRec stack_rec(ladder_start_index);
	int_numb = 0;// variable to count intersections on Split stage
	if (Size>0)
	{
		Size = ring_bufSplit(stack_rec.Q_pos, Size);
		if ((ladder_start_index<stack_rec.Q_pos))
		{
			if (Size)ring_bufFindIntL(ladder_start_index, stack_rec.Q_pos, Size);
			stack_pos = stack_rec.Set(stack_pos, interval_right_index);
		}
	};
	if ((int_numb>Size) && (call_numb<max_call)) //if found a lot of intersections repeat FindR
		Size = ring_bufFindR(stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
	else //cut stripe 
	{

		uint4 m = (interval_left_index + interval_right_index) / 2;
		if (Size>nTotSegm / 100) {// if L contains a lot of segments then cut on two parts
			Size = ring_bufFindR(stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
			Size = ring_bufInsDel(m, stack_pos, Size);
			Size = ring_bufFindR(stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
		}
		else {// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
			uint4 q = (interval_left_index + m) / 2;
			if (interval_left_index != q) {
				Size = ring_bufFindR(stack_rec.Q_pos, interval_left_index, q, stack_pos, Size, 0);
				Size = ring_bufInsDel(q, stack_pos, Size);
			}
			if (q != m) {
				Size = ring_bufFindR(stack_rec.Q_pos, q, m, stack_pos, Size, 0);
				Size = ring_bufInsDel(m, stack_pos, Size);
			}
			q = (interval_right_index + m) / 2;
			if (q != m) {
				Size = ring_bufFindR(stack_rec.Q_pos, m, q, stack_pos, Size, 0);
				Size = ring_bufInsDel(q, stack_pos, Size);
			}
			Size = ring_bufFindR(stack_rec.Q_pos, q, interval_right_index, stack_pos, Size, 0);
		}
	}
	if (ladder_start_index >= stack_rec.Q_pos) return Size;
	B = ENDS[interval_left_index].x;   E = ENDS[interval_right_index].x;
	return ring_bufMerge(interval_left_index, ladder_start_index, stack_rec.Q_pos, Size);
};


// main function to find intersections for optimal algorithm
template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::optFindR(int4 father_first_step, int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index,
	ProgramStackRec *stack_pos, int4 Size, int4 call_numb)
{
	B = ENDS[interval_left_index].x; E = ENDS[RBoundIdx = interval_right_index].x;
    if (interval_right_index - interval_left_index == 1)
        return Size > 1 ? SearchInStrip(ladder_start_index, Size):Size;
	ProgramStackRec stack_rec(ladder_start_index);
	int_numb = 0;// variable to count intersections on Split stage
  bool use_opt = (ladder_start_index- father_first_step>big_staircase_threshold);
  if (use_opt)
  {// use optimal variant if father staircase is big
      Size = optSplit(father_first_step, stack_rec.Q_pos, Size);
      for (int4 i = 0; i<Size; i++)
          optFindInt(ladder_start_index, stack_rec.Q_pos, R[i], Scoll[L[i]]);//R[i] should contain the location of the segment Scoll[L[i]]
      stack_pos = stack_rec.Set(stack_pos, interval_right_index);
      father_first_step = ladder_start_index + inherit_offset;
  }
  else
  {// use fast variant if new staircase is small 
      Size = Split(stack_rec.Q_pos, Size);
      if ((ladder_start_index < stack_rec.Q_pos))
      {
          //if(Size)FindIntL(stack_rec.Q_pos, Size);
          for (int4 i = ladder_start_index + 1; i <= stack_rec.Q_pos; i++)
              father_loc[i] = undef_loc;
          stack_pos = stack_rec.Set(stack_pos, interval_right_index);
          father_first_step = ladder_start_index + inherit_offset;
      }
  }
	if ((int_numb>Size) && (call_numb<max_call)) //if found a lot of intersections repeat optFindR
		Size = optFindR(father_first_step, stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
	else //cut stripe 
	{
      uint4 m = (interval_left_index + interval_right_index) / 2;
      if (Size>nTotSegm / 100) {// if L contains a lot of segments then cut on two parts
          Size = optFindR(father_first_step, stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
          Size = optInsDel(m, stack_pos, Size);
          Size = optFindR(father_first_step, stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
      }
      else {// if L contains not so many segments than cut on four parts (works faster for some segment distributions)
          uint4 q = (interval_left_index + m) / 2;
          if (interval_left_index != q) {
              Size = optFindR(father_first_step, stack_rec.Q_pos, interval_left_index, q, stack_pos, Size, 0);
              Size = optInsDel(q, stack_pos, Size);
          }
          if (q != m) {
              Size = optFindR(father_first_step, stack_rec.Q_pos, q, m, stack_pos, Size, 0);
              Size = optInsDel(m, stack_pos, Size);
          }
          q = (interval_right_index + m) / 2;
          if (q != m) {
              Size = optFindR(father_first_step, stack_rec.Q_pos, m, q, stack_pos, Size, 0);
              Size = optInsDel(q, stack_pos, Size);
          }
          Size = optFindR(father_first_step, stack_rec.Q_pos, q, interval_right_index, stack_pos, Size, 0);
      }

	}
	if (ladder_start_index >= stack_rec.Q_pos) return Size;
	B = ENDS[interval_left_index].x;   E = ENDS[interval_right_index].x;
  if(use_opt)
      Size= optMerge(interval_left_index, ladder_start_index, stack_rec.Q_pos, Size);
  else
      Size=Merge(interval_left_index,ladder_start_index, stack_rec.Q_pos, Size);
  return Size;
};

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::no_ipFindR(int4 ladder_start_index, uint4 interval_left_index, uint4 interval_right_index,
	ProgramStackRec *stack_pos, int4 Size, int4 call_numb)
{
	//  int4 RSize;
	B = ENDS[interval_left_index].x; E = ENDS[RBoundIdx = interval_right_index].x;
	if (interval_right_index - interval_left_index == 1)
		return no_ipSearchInStrip(ladder_start_index, Size);
	ProgramStackRec stack_rec(ladder_start_index);
	int_numb = 0;// variable to count intersections on Split stage
	if (Size>0)
	{
		Size = no_ipSplit(stack_rec.Q_pos, Size, call_numb ? FALSE : TRUE);
		if ((ladder_start_index<stack_rec.Q_pos))
		{
			if(Size)no_ipFindIntL(ladder_start_index, stack_rec.Q_pos, Size);
			stack_pos = stack_rec.Set(stack_pos, interval_right_index);
		}
	};
	if ((int_numb>Size) && (call_numb<max_call)) //if found a lot of intersections repeat FindR
		Size = no_ipFindR(stack_rec.Q_pos, interval_left_index, interval_right_index, stack_pos, Size, call_numb + 1);
	else //cut stripe on the middle
	{
		uint4 m = (interval_left_index + interval_right_index) / 2;
		Size = no_ipFindR(stack_rec.Q_pos, interval_left_index, m, stack_pos, Size, 0);
		Size = no_ipInsDel(m, stack_pos, Size);
		Size = no_ipFindR(stack_rec.Q_pos, m, interval_right_index, stack_pos, Size, 0);
	}
	if (ladder_start_index >= stack_rec.Q_pos) return Size;
	B = ENDS[interval_left_index].x;   E = ENDS[interval_right_index].x;
	Size = no_ipMerge(interval_left_index,ladder_start_index, stack_rec.Q_pos, Size);
	return Size;
};




template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::FreeMem()
{
#define MY_FREE_ARR_MACRO(a) if(a){delete[] a;a=NULL;}
	MY_FREE_ARR_MACRO(SegmentBE);
  MY_FREE_ARR_MACRO(ENDS);
  MY_FREE_ARR_MACRO(Q);
	MY_FREE_ARR_MACRO(L);
	MY_FREE_ARR_MACRO(R);
	MY_FREE_ARR_MACRO(Loc);
	MY_FREE_ARR_MACRO(father_loc);
	MY_FREE_ARR_MACRO(no_ipL);
	MY_FREE_ARR_MACRO(no_ipR);
	MY_FREE_ARR_MACRO(no_ipQ);
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::AllocMem(BOOL for_optimal)
{
	SegmentBE = new EndIndexes[nTotSegm];
	L = new int4[nTotSegm];
	R = new int4[nTotSegm];
	ENDS = new TSegmEnd[2 * nTotSegm];
	
	len_of_Q = nTotSegm;
	if (for_optimal)
	{
		// for optimal alg it needs additional space reserve nTotSegm/(inherit_each-1)+inherit_each+1 to store inherited stairs
		len_of_Q = nTotSegm + nTotSegm / (inherit_each - 1) + inherit_each + 1;
		father_loc = new int4[len_of_Q];
		father_loc[undef_loc] = undef_loc;//nesessary precondition for optFindIntI
	}
	Q = new int4[len_of_Q];
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::no_rAllocMem()
{
	SegmentBE = new EndIndexes[nTotSegm];
	L = new int4[nTotSegm];
	Q = new int4[nTotSegm];
	len_of_Q = nTotSegm;
	R = NULL;
	ENDS = new TSegmEnd[2 * nTotSegm];
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::no_ipAllocMem()
{
	SegmentBE = new EndIndexes[nTotSegm];
	no_ipL = new no_ipSegmentInfo[nTotSegm + 1];
	no_ipR = new no_ipSegmentInfo[nTotSegm + 1];
	ENDS = new TSegmEnd[2 * nTotSegm];
	Loc = new int4[nTotSegm];
	no_ipQ = new no_ipSegmentInfo[nTotSegm];
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::clone(CIntersectionFinder *master, PRegObj robj)
{
	
	nTotSegm = master->nTotSegm;
	len_of_Q = nTotSegm;
	from_begin = true;
	Scoll = master->Scoll;
	SegmentBE = master->SegmentBE;
	ENDS = master->ENDS;
	L = new int4[nTotSegm];


	Q = new int4[nTotSegm];
	R = new int4[nTotSegm];
	//  res->Loc=new int4[nTotSegm+1];

	set_segm_fuctions(master->_Below,
      master->_FindAndRegIPoints,
      master->_FindAndRegIPointsInStripe,
      master->_IsIntersectInStripe,
      master->_BegPoint,
      master->_EndPoint,
      master->_Under,
		NULL,
		NULL,
		NULL,
		robj,
      master->is_line_segments);
  clone_of= master;
 };

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::unclone()
{
    // to be done: my_counter summation should be here
    // lock(clone_of->my_counter) for(i...)clone_of->my_counter[i]+=my_counter[i];unlock(clone_of->my_counter);
	SegmentBE = NULL;
	ENDS = NULL;
	MY_FREE_ARR_MACRO(Q);
	MY_FREE_ARR_MACRO(L)
	MY_FREE_ARR_MACRO(R);
  MY_FREE_ARR_MACRO(Loc);
	Scoll = NULL;
	nTotSegm = 0;
  clone_of = NULL;
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::prepare_ends(uint4 n)
{
	uint4 i, j;
	REAL x, y;
	for (i = 0, j = 0; i<n; i++)
	{
		ENDS[j].segm = 2*i;// even value means first (leftmost) end
		_BegPoint(Scoll[i], ENDS[j].x, y);
		//ENDS[j].islast() = false;
		j++;
		_EndPoint(Scoll[i], ENDS[j].x, y);
		ENDS[j].segm =  2*i+1;// odd value means last (rightmost) end
		//ENDS[j].islast() = true;
		j++;
	}
  std::sort(ENDS, ENDS + j);
	for (uint4 i = 0; i<2 * n; i++)
		if (ENDS[i].islast())
			SegmentBE[ENDS[i].s()].E = i;
		else
			SegmentBE[ENDS[i].s()].B = i;
};

template<bool is_line_seg>
CIntersectionFinder<is_line_seg>::CIntersectionFinder()
{
	SegmentBE = NULL;
	_reg_obj = NULL;
	Q = NULL;
	L = NULL;
	R = NULL;
	ENDS = NULL;
	Scoll = NULL;
	Loc = NULL;
	father_loc = NULL;
	no_ipL = NULL;
	no_ipR = NULL;
	no_ipQ = NULL;
  clone_of = NULL;
	memset(my_counter, 0, sizeof(my_counter));
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::set_segm_fuctions(
	FBelow below,
	FFindAndRegIPoints findAndRegIPoints,
	FFindAndRegIPointsInStripe findAndRegIPointsInStripe,
	FIsIntersectInStripe isIntersectInStripe,
	FGetPoint begPoint,
	FGetPoint endPoint,
	FUnder under,
	FFindAndRegIPointOnRightOfSWL findAndRegIPointOnRightOfSWL,
	FYAtX yAtX,
	FRegIntersection regIntersection,
	PRegObj reg_obj,
	int4 is_line
	)
{
	_Below = below;
	_FindAndRegIPoints = findAndRegIPoints;
	_FindAndRegIPointsInStripe = findAndRegIPointsInStripe;
	_IsIntersectInStripe = isIntersectInStripe;
	_BegPoint = begPoint;
	_EndPoint = endPoint;
	_Under = under;
	_FindAndRegIPointOnRightOfSWL = findAndRegIPointOnRightOfSWL;
	_YAtX = yAtX;
	_RegIntersection = regIntersection;
	_reg_obj = reg_obj;
	is_line_segments = is_line;
}


template <bool is_line_seg>
void CIntersectionFinder<is_line_seg>::balaban_fast(uint4 n, PSeg _Scoll[])
{
    Scoll = _Scoll;
    nTotSegm = n;
    AllocMem(FALSE);
    prepare_ends(n);
    from_begin = true;
    ProgramStackRec stack_rec(-1, 2 * n); //need to be initialized this way
    L[0] = ENDS[0].s();
    FindR(-1, 0, 2 * n - 1, &stack_rec, 1, 0);
    FreeMem();
}

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::balaban_memory_save(uint4 n, PSeg _Scoll[])
{
	Scoll = _Scoll;  nTotSegm = n;
	no_rAllocMem();
	prepare_ends(n);
	
	ProgramStackRec stack_rec(-1, 2 * n);  //need to be initialized this way
    L[0] = ENDS[0].s();
	from_begin = true;
    msFindR(-1, 0, 2*n-1, &stack_rec, 1, 0);
	FreeMem();
}

template <bool is_line_seg>
void CIntersectionFinder<is_line_seg>::balaban_ring_buf(uint4 n, PSeg _Scoll[])
{
	Scoll = _Scoll;  nTotSegm = n;
	no_rAllocMem();
	prepare_ends(n);
	ProgramStackRec stack_rec(-1, 2 * n);  //need to be initialized this way
	L[0] = ENDS[0].s();
	ringb_beg = 0;
	ringb_end = 1;
	ring_bufFindR(-1, 0, 2 * n - 1, &stack_rec, 1, 0);
	FreeMem();
}

//same as balaban_fast but without recursion (just for the case, it doesn't work faster)
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::balaban_no_recursion(uint4 n, PSeg _Scoll[])
{
	Scoll = _Scoll;  nTotSegm = n;
	AllocMem(FALSE);
  prepare_ends(n);
  uint4 interval_left_index = 0, interval_right_index = 2 * n - 1,  m;
  int4 ladder_start_index = -1, call_numb = 0;


	struct SNoRecursionRec :public ProgramStackRec
	{
      uint4 left_bound, m;
      int4   Q_pos_prev;
		SNoRecursionRec() :ProgramStackRec(-1,0){ left_bound = m = 0; };
	} *stack = new SNoRecursionRec[32 * (max_call + 1)];

	SNoRecursionRec *stack_pos = stack, *stack_bottom = stack;
	stack->right_bound = 2 * n;
  int4 Size = 1;
  L[0] = ENDS[0].s();
	while (1)
	{
		B = ENDS[interval_left_index].x;
		while ((interval_right_index - interval_left_index)>1)
		{
			stack++;
			stack->left_bound = interval_left_index; stack->right_bound = interval_right_index;
			E = ENDS[RBoundIdx = interval_right_index].x;
			stack->Q_pos_prev = ladder_start_index;
			int_numb = 0;// variable to count intersections on Split stage
			stack->prev = stack_pos;
			if (Size>0)
			{
				Size = Split(ladder_start_index, Size);
				if ((stack->Q_pos_prev<ladder_start_index))
				{
					//if(Size)FindIntL(stack->Q_pos_prev, ladder_start_index, Size);
					//if(Size)FindIntL(ladder_start_index, Size);
                    stack_pos = stack;
				}
			}
			stack->m = m = (interval_right_index + interval_left_index) >> 1;
			stack->Q_pos = ladder_start_index;
			if ((int_numb>Size) && (call_numb<max_call))
				call_numb++;
			else
			{
				call_numb = 0;
				interval_right_index = m;
			}
		}
		E = ENDS[RBoundIdx = interval_right_index].x;
		if(Size>1)SearchInStrip(ladder_start_index, Size);
		while (interval_right_index == stack->right_bound)
		{
			if (stack->Q_pos_prev<stack->Q_pos)
			{
				B = ENDS[stack->left_bound].x;
				Size = Merge(stack->left_bound,stack->Q_pos_prev, stack->Q_pos, Size);
			};
			stack--;
		}
		if (stack == stack_bottom)break;
		ladder_start_index = stack->Q_pos;
		stack_pos = stack;
		Size = InsDel(interval_left_index = stack->m, stack, Size);
		interval_right_index = stack->right_bound;
	}
	delete[] stack_bottom;
	FreeMem();
}


template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::balaban_no_ip(uint4 n, PSeg _Scoll[])
{
	Scoll = _Scoll;  nTotSegm = n;
	no_ipAllocMem();
	prepare_ends(n);
	ProgramStackRec stack_rec(-1, 2 * n);  //need to be initialized this way
  no_ipL[0].s = ENDS[0].s();
	no_ipFindR(-1, 0, 2*n-1, &stack_rec, 1, 0);
	FreeMem();
}



template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::balaban_optimal(uint4 n, PSeg _Scoll[])
{
	Scoll = _Scoll;  nTotSegm = n;
	AllocMem(TRUE);
	prepare_ends(n);
	ProgramStackRec stack_rec(inherit_each, 2 * n);  //need to be initialized this way 
    L[0] = ENDS[0].s();
    optFindR(inherit_each + 1, inherit_each, 0, 2*n-1, &stack_rec, 1, 0);
	FreeMem();
}

template<bool is_line_seg>
int4 CIntersectionFinder<is_line_seg>::CalcLAt(uint4 end_index)
{
	int4 i, Size = 0;
	EndIndexes * cur_ends;
	for (i = 0; i<nTotSegm; i++)
	{

		cur_ends = SegmentBE + i;
		if ((cur_ends->B <= end_index) && (cur_ends->E>end_index))
			L[Size++] = i;
	}

  std::sort(L, L + Size, CSegmCompare(Scoll, _Below, ENDS[end_index].x));
	return Size;

};



#include <thread>
#include <vector>
template <bool is_line_seg>
void CIntersectionFinder<is_line_seg>::fast_parallel(uint4 n, PSeg _Scoll[], int4 n_threads,PRegObj add_reg[])
{
	Scoll = _Scoll;
	nTotSegm = n;
	AllocMem(FALSE);
	prepare_ends(n);
  using  namespace std;
  
  vector<thread> wrk_threads;
  auto thread_func = [](CIntersectionFinder *master, uint4 from, uint4 to, PRegObj add_reg)
        {
            CIntersectionFinder<is_line_seg> i_f;
            i_f.clone(master, add_reg); 
            i_f.FindR(-1, from, to, &ProgramStackRec(-1, 2 *  i_f.nTotSegm), i_f.CalcLAt(from), 0);
            i_f.unclone();
        };
  double part = 2 * n /(double) n_threads;
  int4 i = 1;
  uint4 start_from = part;
  uint4 from,to = start_from;
  while ( i < n_threads)
  {
      from = to;
      i++;
      to = (i == n_threads) ? 2 * n - 1 : part*i;
      wrk_threads.emplace_back(thread_func, this, from,to, add_reg[i-2]);// starts intersection finding in a stripe <from,to>
  }
  L[0] = ENDS[0].s();
  FindR(-1, 0, start_from, &ProgramStackRec(-1, 2 * n), 1, 0);
  auto cur_thread = wrk_threads.begin();
  for (; cur_thread != wrk_threads.end(); cur_thread++) cur_thread->join();//waiting for calculation of all threads are finished
  FreeMem();
}




template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::trivial(uint4 n, PSeg sgm[])
{
	uint4 i, j, m = n - 1;
	for (i = 0; i<m; i++)
		for (j = i + 1; j<n; j++)
			_FindAndRegIPoints(sgm[i], sgm[j], _reg_obj);
};


template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::simple_sweep(uint4 n, PSeg Scoll[])
{
	uint4 i, j;
	PSeg s;
	uint4 pos, k;
	PSeg *sgm = new PSeg[n + 1];
	struct TSEnd :public TEnd
	{
		PSeg s;
	} *Ends = new TSEnd[2 * n];

	j = 0;
	for (i = 0; i<n; i++)
	{
		REAL y;
		_BegPoint(Scoll[i], Ends[j].x, y);
		Ends[j].s = Scoll[i];
		Ends[j].islast = false;
		j++;
		_EndPoint(Scoll[i], Ends[j].x, y);
		Ends[j].s = Scoll[i];
		Ends[j].islast = true;
		j++;
	}
	std::sort(Ends, Ends+2 * n);
	pos = 0;
	for (i = 0; i<2 * n; i++)
	{
		s = Ends[i].s;
		if (Ends[i].islast)
		{
      j = --pos;
      auto cur = sgm[j];
      while (cur != s) {
          --j;
          auto buf = cur;
          cur = sgm[j];
          sgm[j] = buf;
      }
		}
		else
		{
			for (j = 0; j<pos; j++)
				_FindAndRegIPoints(s, sgm[j], _reg_obj);
			sgm[pos] = s;
			pos++;
		}
	}
	delete[] Ends;
	delete[]sgm;
}


// Bentley & Ottmann algorithm
template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::SweepLineInsert(int4 s)
{
	TBinTreeNode *inserted = SweepLine + s;
	inserted->segment = inserted->seg_node_idx = s;
	if (SweepLine[swl_root].segment == -1){ swl_root = s; return; }
	int4 cur_node = swl_root;
	int4 onright;
	int4 cur_seg = SweepLine[cur_node].segment;
	PSeg ps = Scoll[s];
	while (cur_seg != -1)
	{
		if (onright = _Under(Scoll[cur_seg], ps))
		{
			inserted->prev = cur_node;
			cur_node = SweepLine[cur_node].right;
		}
		else
		{
			inserted->next = cur_node;
			cur_node = SweepLine[cur_node].left;
		};
		cur_seg = SweepLine[cur_node].segment;
	};
	if (onright)
		SweepLine[inserted->prev].right = s;
	else
		SweepLine[inserted->next].left = s;

	SweepLine[inserted->next].prev = SweepLine[inserted->prev].next = s;

	IntOnRightOfSWL(SweepLine[inserted->prev].segment, s);
	IntOnRightOfSWL(s, SweepLine[inserted->next].segment);
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::SweepLineDelete(int4 s)
{
	static bool del_prev = true;
	int4 cur_node = SweepLine[s].seg_node_idx;
	TBinTreeNode *deleted = SweepLine + cur_node;
	IntOnRightOfSWL(SweepLine[deleted->prev].segment, SweepLine[deleted->next].segment);
	bool left_null, right_null;
	while (deleted->segment != -1)
	{
		left_null = SweepLine[deleted->left].segment == -1;
		right_null = SweepLine[deleted->right].segment == -1;
		if (left_null&&right_null)
		{
			SweepLine[deleted->prev].next = deleted->next;
			SweepLine[deleted->next].prev = deleted->prev;
			deleted->segment = -1; break;
		}
		else
		{
			if ((!left_null&&del_prev) || (!del_prev&&right_null))
			{
				deleted->segment = SweepLine[deleted->prev].segment;
				SweepLine[deleted->segment].seg_node_idx = cur_node;
				cur_node = deleted->prev;
			}
			else
			{
				deleted->segment = SweepLine[deleted->next].segment;
				SweepLine[deleted->segment].seg_node_idx = cur_node;
				cur_node = deleted->next;
			}
			deleted = SweepLine + cur_node;
		}
	}
	del_prev = !del_prev;
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::SweepLineExchange(int4 s1, int4 s2)
{
	int4 node1 = SweepLine[s1].seg_node_idx;
	int4 node2 = SweepLine[s2].seg_node_idx;
	IntOnRightOfSWL(s1, SweepLine[SweepLine[node2].next].segment);
	IntOnRightOfSWL(SweepLine[SweepLine[node1].prev].segment, s2);
	SweepLine[node1].segment = s2;
	SweepLine[s2].seg_node_idx = node1;
	SweepLine[node2].segment = s1;
	SweepLine[s1].seg_node_idx = node2;
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::EventsDelMin()
{
	Events[0] = Events[--events_n];
	if (events_n<2) return;
	uint4 cur = 0;
	uint4 son = 2;
	TEvent tmp;
	do
	{
		if (Events[son - 1].pt.x<Events[son].pt.x)
			son--;
		if (Events[cur].pt.x>Events[son].pt.x)
		{
			tmp = Events[cur]; Events[cur] = Events[son]; Events[son] = tmp;
			cur = son;
			son = 2 + (cur << 1);
		}
		else break;
	} while (son<events_n);
	if (son == events_n)
	{
		son--;
		if (Events[cur].pt.x>Events[son].pt.x)
		{
			tmp = Events[cur]; Events[cur] = Events[son]; Events[son] = tmp;
		}
	}
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::EventsAddNew()
{
	int4 cur = events_n++;
	if (events_n == 1)return;
	int4 father;
	TEvent tmp;
	while (cur)
	{
		father = (cur - 1) >> 1;
		if (Events[cur].pt.x<Events[father].pt.x)
		{
			tmp = Events[cur]; Events[cur] = Events[father]; Events[father] = tmp;
		}
		else
			break;
		cur = father;
	}
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::IntOnRightOfSWL(int4 s1, int4 s2)// important s1 below s2 at sweep_line_x 
{
	if ((s1 == -1) || (s2 == -1))return;
	if (events_n + 3>events_max)
	{
		events_max *= 2;
		TEvent *newEvents = new TEvent[events_max];
		for (uint4 i = 0; i<events_n; i++)
			newEvents[i] = Events[i];
		delete[] Events;
		Events = newEvents;
	}
	TEvent *ev = Events + events_n;
	int4 n;
	if (n = _FindAndRegIPointOnRightOfSWL(sweep_line_x, Scoll[s1], Scoll[s2], &(ev->pt), &(ev[1].pt)))// find and place to ev->pt
	{
		ev->s1 = s1;
		ev->s2 = s2;
		EventsAddNew();
		if (n == 2) //for arcs intersection case
		{
			ev++;
			ev->s1 = s2;
			ev->s2 = s1;
			EventsAddNew();
		}
	}
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::PrepareEvents()
{
	TEvent *ev = Events;
	for (uint4 i = 0; i<nTotSegm; i++)
	{
		_BegPoint(Scoll[i], ev->pt.x, ev->pt.y);
		ev->s1 = i;
		ev->s2 = -1;
		EventsAddNew();
		ev++;
		_EndPoint(Scoll[i], ev->pt.x, ev->pt.y);
		ev->s1 = -1;
		ev->s2 = i;
		EventsAddNew();
		ev++;
	}
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::AllocBOMem()
{
	swl_root = -1;
	SweepLine = new TBinTreeNode[nTotSegm + 1];
	SweepLine++;
	events_n = 0;
	Events = new TEvent[events_max = nTotSegm * 3];
};

template<bool is_line_seg>
void CIntersectionFinder<is_line_seg>::FreeBOMem()
{
	if (SweepLine){ SweepLine--; delete[] SweepLine; SweepLine = NULL; }
	MY_FREE_ARR_MACRO(Events);
	events_max = 0;
};

template<bool is_line_seg> 
void CIntersectionFinder<is_line_seg>::bentley_ottmann(uint4 n, PSeg _Scoll[])
{
	Scoll = _Scoll;  nTotSegm = n;
	AllocBOMem();
	PrepareEvents();
	int4 s1, s2;
	while (events_n>0)
	{
		sweep_line_x = Events->pt.x;
		s1 = Events->s1;
		s2 = Events->s2;
		if ((s1 != -1) && (s2 != -1))
			_RegIntersection(_reg_obj, Scoll[s1], Scoll[s2], 1, &(Events->pt));// register found intersection
		EventsDelMin();
		while ((s1 == Events->s1) && (s2 == Events->s2))
		{
			s1 = Events->s1;
			s2 = Events->s2;
			EventsDelMin();
		}
		if (s1 != -1)
		{
			if (s2 != -1)
				SweepLineExchange(s1, s2);
			else
				SweepLineInsert(s1);
		}
		else
			SweepLineDelete(s2);
	}
	FreeBOMem();
};

template class CIntersectionFinder<true>;
template class CIntersectionFinder<false>;
