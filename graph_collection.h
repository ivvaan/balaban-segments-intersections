#pragma once
#pragma once
/*
*
*      Copyright (c)  2011-2020  Ivan Balaban
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
// NEW IMPLEMENTATION

/*struct SegmEndIndexes
{
  int4 B, E;
  bool operator<(SegmEndIndexes si) { return (B<si.B) || ((B == si.B) && (E<si.E)); };
  bool operator!=(SegmEndIndexes si) { return (B != si.B) || (E != si.E); };
};*/

#include "segments.h"
#include <iostream>
template<class IntersectionRegistrator>
class CGraphSegmentCollection
{
public:
  static constexpr bool is_line_segments = true;

  static inline bool is_last(uint4 pt)
  {
    return pt & 1;
  };
  static inline uint4 get_segm(uint4 pt)
  {
    return pt >> 1;
  };
  static inline uint4 get_first_pt(uint4 s)
  {
    return s << 1;
  };
  static inline uint4 get_last_pt(uint4 s)
  {
    return (s << 1) + 1;
  };
  static inline uint4 get_other_pt(uint4 pt)
  {
    return pt ^ 1;
  };
 inline uint4 get_first_idx(uint4 s) const
  {
    return vertex_idx[get_first_pt(s)];
  };
 inline uint4 get_last_idx(uint4 s) const
  {
    return vertex_idx[get_last_pt(s)];
  };
 inline uint4 get_other_idx(uint4 pt) const
  {
    return vertex_idx[get_other_pt(pt)];
  };

  //TPlaneVect
  uint4  GetSegmNumb() const { return nEdges; };
  inline void SetCurStripe(uint4 left, uint4 right)
  {

    SetCurStripeLeft(left);
    SetCurStripeRight(right);
  };
  inline void SetCurStripeLeft(uint4 left) { ptB = left, idxB = vertex_idx[left]; B = vertices[idxB].x; };
  inline void SetCurStripeRight(uint4 right) { ptE = right, idxE = vertex_idx[right]; E = vertices[idxE].x; };
  inline void SetCurPoint(uint4 pt)
  {
    cur_pt = pt;
    cur_pt_idx = vertex_idx[pt];
    cur_point = vertices[cur_pt_idx];
  };
  inline void SetCurPointAtBeg(uint4 s)
  {
    SetCurPoint(get_first_pt(s));
//    collection[s].BegPoint(cur_point.x, cur_point.y);
  };
  void SetCurPointAtEnd(uint4 s)
  {
    SetCurPoint(get_last_pt(s));
  };
  void SetCurSegCutBE(uint4 s)
  {
    SetCurSeg(s);
    curB = MAX(curB, B);
    curE = MIN(curE, E);
  };
  void SetCurSegCutBeg(uint4 s)
  {
    SetCurSeg(s);
    curB = MAX(curB, B);
  };
  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    curE = MIN(curE, E);
  };

  void SetCurSeg(uint4 s)
  {
    cur_seg_idx = s;
    cur_seg_beg_idx = get_first_idx(s);
    cur_seg_end_idx = get_last_idx(s);

    cur_seg = TLineSegment1(vertices[cur_seg_beg_idx],vertices[cur_seg_end_idx]);
    curB = cur_seg.org.x;
    curE = cur_seg.org.x + cur_seg.shift.x;
  };

  void SetCurSegAE(uint4 s)
  {
    SetCurSeg(s);
  };

  bool LBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    auto beg1 = get_first_idx(s_1);
    auto beg2 = get_first_idx(s_2);
    auto end1 = get_last_idx(s_1);
    auto end2 = get_last_idx(s_2);
    auto org1 = vertices[beg1];
    if ((beg1 == beg2)&&(beg1==idxB)) 
      return (vertices[end1] - org1) % (vertices[end2] - org1) > 0;
    auto org2 = vertices[beg2];
    auto shift1 = vertices[end1] - org1;
    auto shift2 = vertices[end2] - org2;
    auto y1 = (org1.y*shift1.x + (B - org1.x)*shift1.y)*shift2.x;
    auto y2 = (org2.y*shift2.x + (B - org2.x)*shift2.y)*shift1.x;
    return y1<y2;
  };
  bool RBelow(int4 s_1, int4 s_2) const //retuns if s1 below s2 at current vertical line
  {
    auto beg1 = get_first_idx(s_1);
    auto beg2 = get_first_idx(s_2);
    auto end1 = get_last_idx(s_1);
    auto end2 = get_last_idx(s_2);
    auto org1 = vertices[end1];
    if ((end1 == end2) && (end1 == idxE))
      return (vertices[beg2] - org1) % (vertices[beg1] - org1) > 0;
    auto org2 = vertices[end2];
    auto shift1 = vertices[beg1] - org1;
    auto shift2 = vertices[beg2] - org2;
    auto y1 = (org1.y*shift1.x + (E - org1.x)*shift1.y)*shift2.x;
    auto y2 = (org2.y*shift2.x + (E - org2.x)*shift2.y)*shift1.x;
    return y1<y2;
  };

  template <bool register_int>
  bool FindIntWith(REAL x1, REAL x2, uint4 s_)
  {
    auto beg = get_first_idx(s_);
    auto end = get_last_idx(s_);
    if ((cur_seg_beg_idx == beg) || (cur_seg_end_idx == end))return false;
    auto s2 = TLineSegment1(vertices[beg], vertices[end]);
    TPlaneVect delt = s2.org - cur_seg.org;
    REAL prod = cur_seg.shift%s2.shift, mul, xc;
    mul = cur_seg.shift%delt;
    if ((mul>0) ^ (mul + prod>0))
    {
      if (prod == 0) return false;
      mul = mul / prod;
      xc = s2.org.x - mul*s2.shift.x;
      if (((xc <= x1) || (xc>x2))) return false;
      if constexpr(register_int)
      {
        if constexpr((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0)
        {
          TPlaneVect p(xc, s2.org.y - mul*s2.shift.y);
          registrator->register_pair_and_point(cur_seg_idx, s_,p);
        }
        else
          registrator->register_pair(cur_seg_idx, s_);
      }
      return true;
    }
    return false;
  };


  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto x1 = std::max(curB, vertices[get_first_idx(s_)].x);
    auto x2 = std::min(curE, vertices[get_last_idx(s_)].x);
    if (x1 >= x2)return false;
    return FindIntWith<true>(x1, x2, s_);
  };

  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    return FindIntWith<true>(curB, curE, s_);
  };

  bool FindCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    return FindIntWith<true>(curB, curE, s_);
  };

  bool IsIntersectsCurSegDown(int4 s_)//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return FindIntWith<false>(curB, curE, s_);
  };

  bool IsIntersectsCurSegUp(int4 s_)//check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return FindIntWith<false>(curB, curE, s_);
  };

  bool UnderCurPoint(int4 s_) const
  {
    auto beg= get_first_idx(s_);
    auto org = vertices[beg];
    if (beg != cur_pt_idx)
      return (vertices[get_last_idx(s_)] - org) % (cur_point - org) >0;
    return (vertices[get_last_idx(s_)] - org) % (vertices[get_other_idx(cur_pt)] - org) >0;
  };//returns true if s is under current point 

  void PrepareEndpointsSortedList(uint4 *epoints)// endpoints allocated by caller and must contain space for at least 2*GetSegmNumb() points 
  {
    auto NN = nEdges * 2;
    for (uint4 i = 0; i < NN; ++i)     epoints[i] = i;
    std::sort(epoints, epoints + NN, 
      [this](uint4 pt1, uint4 pt2) {
        auto idx1 = vertex_idx[pt1];
        auto idx2 = vertex_idx[pt2];
        if (idx1 != idx2)
          return vertices[idx1] < vertices[idx2];
        bool l1 = is_last(pt1);
        bool l2 = is_last(pt2);
        if (l1 < l2)return false;
        if (l1 > l2)return true;
    //    return (bool)(l1 ^ ((vertices[get_other_idx(pt1)] - vertices[idx1]) % (vertices[get_other_idx(pt2)] - vertices[idx2])>0));
        return (vertices[get_other_idx(pt1)] - vertices[idx1]) % (vertices[get_other_idx(pt2)] - vertices[idx2])>0;
    }
    );
  };
  void clone(CGraphSegmentCollection &c, IntersectionRegistrator *r)
  {
    clone_of = &c;
    nVertices = c.nVertices;
    nEdges= c.nEdges;
    vertices = c.vertices;
    vertex_idx = c.vertex_idx;
    SetRegistrator(r);
  };
  
  void unclone() 
  { 
    if (clone_of == nullptr)return; 
    nVertices = 0;
    nEdges = 0;
    vertices = nullptr;
    vertex_idx = nullptr;
    clone_of = nullptr;
  };

  void SortAt(uint4 pt, uint4 n, int4 *L)
  {
    SetCurStripeLeft(pt);
    std::sort(L, L + n, [this](int4 s1, int4 s2) {return LBelow(s1, s2); });
  };

  void SetRegistrator(IntersectionRegistrator *r)
  {
    registrator = r;
    r->Alloc(GetSegmNumb());
  };
  IntersectionRegistrator *GetRegistrator() { return registrator; };

  void Init(uint4 n, void * c, IntersectionRegistrator *r)
  {
    nVertices = n; nEdges=n;
    vertices = reinterpret_cast<TPlaneVect*>(c);
    //vertices = new TPlaneVect[n];
    vertex_idx = new uint4[2 * n];
    for (uint4 i = 0; i < n; i++)
    {
      auto ii = i << 1;
      auto j = (i + 1) % n;
      // adding edge connecting current and next points
      if (vertices[i] < vertices[j])//first point of the edge should be "on the left" from second
        {
          vertex_idx[ii] = i;
          vertex_idx[ii + 1] = j;
        }
      else
        {
          vertex_idx[ii] = j;
          vertex_idx[ii + 1] = i;
        }
    }
    /*std::cout << "vertices *********************\n[";
      for (uint4 i = 0; i < nVertices; ++i)
        std::cout <<  vertices[i] << ",\n";
      std::cout << "]\n edges *********************\n[";
      for (uint4 i = 0; i < nEdges; ++i)
        std::cout << "[" << vertex_idx[i*2] <<"," << vertex_idx[i * 2+1] << "]\n";
      std::cout << "]\n";*/
    SetRegistrator(r);

  };

  CGraphSegmentCollection(uint4 n, void* c, IntersectionRegistrator *r)
  {
    Init(n, c, r);
  }
  CGraphSegmentCollection() {};
  ~CGraphSegmentCollection() 
  {
    unclone();
    //MY_FREE_ARR_MACRO(vertices);
    MY_FREE_ARR_MACRO(vertex_idx);

  };

  void SetSearchDirDown(bool dir) { };

  void coll_to_SVG(chostream* SVG_text) {
    if (!SVG_text)return;
    int4 n = MIN(max_SVG_items, GetSegmNumb());
    REAL xmin = 0, ymin = 0, xmax = 1, ymax = 1;
    for (int4 i = 0; i < nVertices; ++i) {
      xmin = MIN(xmin, vertices[i].x);
      xmax = MAX(xmax, vertices[i].x);
      ymin = MIN(ymin, vertices[i].y);
      ymax = MAX(ymax, vertices[i].y);
    }
    *SVG_text << "<svg height='100%' width='100%' viewBox='";
    *SVG_text << xmin << " " << ymin << " "
      << xmax - xmin << " " << ymax - ymin << "'>\n";
    for (uint4 s = 0; s < n; ++s)
    {
      auto& bp = vertices[get_first_idx(s)];
      auto& ep = vertices[get_last_idx(s)];
      *SVG_text << "<line id='seg" << s;
      *SVG_text << "' x1='" << bp.x;
      *SVG_text << "' y1='" << bp.y;
      *SVG_text << "' x2='" << ep.x;
      *SVG_text << "' y2='" << ep.y;
      *SVG_text << "' class='edge'/>\n";

    }

  };


private:
  inline auto GetX(uint4 pt) const
  {
    return vertices[vertex_idx[pt]].x;
  };


  IntersectionRegistrator *registrator = nullptr;
  CGraphSegmentCollection *clone_of = nullptr;
  uint4 nVertices, nEdges;
  TPlaneVect *vertices=nullptr;
  uint4 *vertex_idx=nullptr;
  //SegmEndIndexes *se_index;
  REAL B, E, curB, curE;
  uint4 ptB, ptE, idxB, idxE;
  uint4 cur_pt= 0xFFFFFFFF;
  uint4 cur_pt_idx = 0xFFFFFFFF;
  TPlaneVect cur_point;

  uint4 cur_seg_idx = 0xFFFFFFFF;
  uint4 cur_seg_beg_idx = 0xFFFFFFFF;
  uint4 cur_seg_end_idx = 0xFFFFFFFF;
  TLineSegment1 cur_seg;
};

