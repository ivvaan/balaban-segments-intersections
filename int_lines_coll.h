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

// New template-based implementation.
//
// This collection is the integer-geometry pipeline for the degenerate_cases branch.
// It supports degenerate configurations that break the original Balaban assumptions:
// - collinear overlaps (handled via a preprocessing "remaper" that splits overlaps)
// - multiple endpoints with identical X ("multi-events", i.e., collapsed same-X boundaries)
//
// Zero-length segments (begin == end) are supported (they are treated as degenerate vertical segments).
// NOTE: The implementation assumes there are no two distinct zero segments that coincide as points
// (asserted in `CRemaper::TurnRemapOn()` in debug builds).

#include "segments.h"
#include "registrator.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <numeric>

template<typename Real>
Real get_max_gap_middle(uint4 N, Real arr[])
{
  // Returns the midpoint of the maximum gap between adjacent elements of `arr`.
  // Note: `arr` is sorted in-place.
  assert(N > 1);
  std::sort(arr, arr + N);

  // Find the first maximum gap and return its midpoint.
  uint4 best_i = 1;
  auto _arr = arr - 1;
  Real max_gap = arr[best_i] - _arr[best_i];

  for (uint4 i = 2; i < N; ++i) {
    if (Real g = arr[i] - _arr[i]; g > max_gap) {
      max_gap = g;
      best_i = i;
    }
  }

  return (arr[best_i] + _arr[best_i]) / Real(2);
}


template<class IntersectionRegistrator>
class CIntegerSegmentCollection
{
  using CTHIS = CIntegerSegmentCollection<IntersectionRegistrator>;

  static REAL get_rot_angle(uint4 n, TLineSegment1 sc[])
  {
    REAL* arr;
    DECL_RAII_ARR(arr, n + 2);
    arr[0] = -M_PI / 2.;
    std::transform(sc, sc + n, arr,
      [](const TLineSegment1& s) { return std::atan2(s.shift.y, s.shift.x); });
    arr[n + 1] = M_PI / 2.;
    auto max_gap_mid = get_max_gap_middle(n + 2, arr);
    return max_gap_mid + (max_gap_mid > 0 ? -M_PI / 2. : M_PI / 2.);
  }

public:
  static constexpr _Coll_flag_state get_coll_flag(_Coll_flags flag)
  {
    if (flag == _Coll_flags::line_segments)
      // With degeneracies enabled (exact endpoint-on-stair cases),
      // "line segment monotonicity" optimizations are impossible. Treat as non-linear
      // so algorithms check intersections in both directions.
      return _Coll_flag_state::state_false;

    return _Coll_flag_state::state_unimplemented;
  }

  enum _Stages {
    stage_split = 0, stage_bubble, stage_merge
  };

  // Endpoints are encoded as point indices in `[0, 2*nSegments)`,
  // where `pt = (segment_index << 1) | is_last_bit`.
  static bool is_last(uint4 pt) { return pt & 1; }
  static bool is_first(uint4 pt) { return (pt & 1) == 0; }
  static uint4 get_segm(uint4 pt) { return pt >> 1; }
  static uint4 first_point(uint4 s) { return s << 1; }
  static uint4 last_point(uint4 s) { return (s << 1) + 1; }
  static uint4 other_point(uint4 pt) { return pt ^ 1; }

private:
  struct ListBounds {
    uint4 beg = 0;
    uint4 end = 0;
  };

  //===========================================================================
#include "remaper.h"
  //===========================================================================

public:
  int4* tmp = nullptr;

  std::pair<uint4, uint4> get_dublicate_stat(uint4 len, uint4* arr) const
  {
    // Counts unique X-coordinates and the number of X-values that appear more than once.
    // Example X: {0,1,3,3,3,4,4} -> unique_count=4, duplicated_elements_count=2 (X=3 and X=4).
    uint4 unique_count = 0;
    uint4 duplicated_elements_count = 0;

    for (uint4 i = 0, j; i < len; i = j) {
      j = i + 1;
      while (j < len && pts[arr[j]].x == pts[arr[i]].x)
        ++j;
      ++unique_count;
      if (j > i + 1)
        ++duplicated_elements_count;
    }
    return { unique_count, duplicated_elements_count };
  }
  template<bool for_simple_sweep = true>
  uint4 PrepareEndpointsSortedList(uint4* epoints)
  {
    // Fills `epoints` with all endpoint indices `[0..2*nSegments)` and sorts them by:
    // X, then Y, then (first before last), then angular tie-break for stable ordering.
    // The output is later collapsed into "multi-events" (same-X groups) in `Prepare()`.
    //
    // Unlike the older codepath, we keep zero-length segments too. They participate in multi-events
    // and in vertical-boundary processing (AllIntCurLine).
    uint4 nPoints = 0;
    if (exclude_zero_segs) {
      for (auto NN = GetSegmNumb() << 1, i = 0U; i < NN; ++i)
        if (collection[get_segm(i)].shift.is_non_zero())
          epoints[nPoints++] = i;
    }
    else {
      nPoints = GetSegmNumb() << 1;
      std::iota(epoints, epoints + nPoints, 0);
    }
    auto is_below = [pts = pts, coll = collection](uint4 pt1, uint4 pt2) {
      if (pts[pt1].x != pts[pt2].x)
        return pts[pt1].x < pts[pt2].x;
      if (pts[pt1].y != pts[pt2].y)
        return pts[pt1].y < pts[pt2].y;

      auto le1 = is_last(pt1), le2 = is_last(pt2);
      auto s1i = get_segm(pt1);
      auto s2i = get_segm(pt2);
      auto& s1 = coll[s1i];
      auto& s2 = coll[s2i];

      if constexpr (for_simple_sweep) {
        if (le1 != le2) {
          return le2;
        };
      }
      else{
        auto vert1 = s1.is_vertical(), vert2 = s2.is_vertical();

        if (vert1 != vert2)
          return vert2; // vertical segments are "above" non-vertical ones at shared endpoints.
        // both vertical or both non-vertical
        if (le1 != le2) 
          return vert1 ? (s1i == s2i) == le2 : le2;

        if (vert1) {
          //ends of vertical segments both first or both last
          //first -> zero - first
          //last -> zero - last
          auto is_zero1 = s1.shift.is_zero();
          auto is_zero2 = s2.shift.is_zero();
          if (is_zero1 != is_zero2)
            // it should be: le1 ? (is_zero1 < is_zero2) : (is_zero2 < is_zero1);
            // taking into account that is_zero1 != is_zero2 => is_zero2 == !is_zero1
            // equivalent to: le1 ? is_zero2 : !is_zero2;
            // equivalent to:
            return le1 == is_zero2;
        };
      };
      //ends of non-vertical segments both first or both last
      auto angle_like = s1.shift % s2.shift;
      if (angle_like != 0)
        return le1 != (0 < angle_like);

      return (pt1 < pt2);
      };

    std::sort(epoints, epoints + nPoints, is_below);
    return nPoints;
  }

  // Multi-event encoding:
  // `ENDS[rank]` is either a normal endpoint index (< 2*nSegments),
  // or a marker `mark_multiple(head)` that points into a secondary list stored in ENDS[].
  bool is_multiple(uint4 pt) const { return pt >= (nSegments << 1); }
  auto mark_multiple(uint4 pt) const { return pt + (nSegments << 1); }
  auto unmark_multiple(uint4 pt) const { return pt - (nSegments << 1); }

  constexpr const static uint4 list_stop = std::numeric_limits<uint4>::max();

  PrepareResult Prepare()
  {
    // Preparation does two degenerate-case steps:
    // 1) TurnRemapOn(): optional normalization that splits collinear overlaps via CRemaper.
    // 2) Collapse equal-X endpoints into "multi-events" to avoid zero-width stripes.
    Reset();
    TurnRemapOn();

    uint4 nSeg = GetSegmNumb();
    if (nSeg == 0) return {};
    auto __tmpENDS__ = std::make_unique<uint4[]>(nSeg*2);
    uint4* tmpENDS = __tmpENDS__.get();
    uint4 nEnds = PrepareEndpointsSortedList<false>(tmpENDS);

    auto [nX, nV] = get_dublicate_stat(nEnds, tmpENDS);
    auto nAll = nEnds + nV;

    // tmp[] maps each endpoint (in sorted tmpENDS order) to its final rank in ENDS[].
    // This keeps seg_L_rank / seg_R_rank construction compatible with multi-events.
    tmp = new int4[nEnds];

    if (nV != 0) {
      // ENDS layout when multi-events exist:
      // - ENDS[0..nX) contains either a normal endpoint index or a multi-event marker.
      // - Secondary lists for same-X groups are appended starting from D_last=nX.
      // - Each list is terminated by `list_stop`.
      ENDS = new uint4[nAll + nV];

      uint4 D_last = nX;
      for (uint4 i = 0, ne_pos = 0, j; i < nEnds; i = j) {
        j = i + 1;
        while (j < nEnds && pts[tmpENDS[j]].x == pts[tmpENDS[i]].x)
          ++j;

        auto len = j - i;
        if (len > 1) {
          std::copy(tmpENDS + i, tmpENDS + j, ENDS + D_last);
          ENDS[ne_pos] = mark_multiple(D_last);
          std::fill(tmp + i, tmp + j, ne_pos);

          D_last += len;
          ENDS[D_last++] = list_stop; // Marks the end of the per-X list.
        }
        else {
          ENDS[ne_pos] = tmpENDS[i];
          tmp[i] = ne_pos;
        }
        ++ne_pos;
      }
    }
    else {
      std::iota(tmp, tmp + nEnds, 0);
      ENDS = tmpENDS;
      __tmpENDS__.release();
    }

    seg_L_rank = new uint4[nSeg];
    seg_R_rank = new uint4[nSeg];

    uint4 max_segm_on_vline = 0, nsegm_on_vline = 0;
    double avr = 0;

    for (uint4 i = 0; i < nEnds; ++i) {
      if (is_last(tmpENDS[i])) {
        seg_R_rank[get_segm(tmpENDS[i])] = tmp[i];
        --nsegm_on_vline;
      }
      else {
        seg_L_rank[get_segm(tmpENDS[i])] = tmp[i];
        ++nsegm_on_vline;
        if (nsegm_on_vline > max_segm_on_vline)
          max_segm_on_vline = nsegm_on_vline;
      }
      avr += nsegm_on_vline;
    }

    avr /= (double)nX;
    nTotX = nX;
    nCollideX = nV;
    return { nX, max_segm_on_vline, avr };
  }

  bool is_remapped() const {
    return is_collection_remapped;
  }

  template<class IntersectionFinder, class ProgramStackRec>
  void InsDelOrdinary(IntersectionFinder& i_f, uint4& L_size, int4* L, uint4 pt, ProgramStackRec* stack_pos)
  {
    // Ordinary event: a single endpoint at this X.
    // - If it's a last point: remove the segment from L.
    // - If it's a first point: insert and optionally search internal intersections higher in the stack.
    auto sn = get_segm(pt);
    if (is_last(pt)) { // Endpoint: remove.
      auto last = std::remove(L, L + L_size, sn);
      assert(last + 1 == L + L_size);
      --L_size;
      ReorderStep(E, L_size, L);
    }
    else { // Start point: insert.
      ReorderStep(E, L_size, L);
      if (stack_pos->isnot_top()) {
        SetCurSegAndPoint(sn);
        i_f.FindIntI(*this, sn, stack_pos); // Find internal intersections in ancestor staircases.
      }
      else {
        SetCurPointAtBeg(sn); // Sets current point to the begin of segment sn.
      }

      int4 i = L_size;
      for (auto _L = L - 1; (i != 0) && (!UnderCurPoint(_L[i])); --i)
        L[i] = _L[i];
      L[i] = sn;
      ++L_size;
    }
  }

  template<class IntersectionFinder, class ProgramStackRec>
  void InsDel(IntersectionFinder& i_f, uint4& L_size, int4* L, uint4 end_rank, ProgramStackRec* stack_pos)
  {
    // Degenerate handling: if ENDS[end_rank] is a "multi-event marker", process the whole same-X boundary
    // as a single combined event to avoid zero-width stripes.
    auto pt = PointAtRank(end_rank);
    if (is_multiple(pt)) {
      DelStep(end_rank, L_size, L);

      auto f = unmark_multiple(pt);
      assert(E == XAtRank(end_rank));

      ReorderStep(E, L_size, L);
      AllIntCurLine(f, L_size, L);
      InsStep(f, L_size, L);

      if (stack_pos->isnot_top())
        // Find intersections for all inserted (first endpoints) and vertical segments
        // in staircase levels above (in the stack).
        for (uint4 i = f, cur_pt = ENDS[f]; cur_pt != list_stop; cur_pt = ENDS[++i]) {
          if (is_first(cur_pt)) {
            auto sn = get_segm(cur_pt);
            SetCurSegAndPoint(sn);
            i_f.FindIntI(*this, sn, stack_pos);
          }
        }
    }
    else {
      InsDelOrdinary(i_f, L_size, L, pt, stack_pos);
    }
  }

  void ReorderStep(int4 x, uint4 L_size, int4* L)
  {
#ifdef TURN_ON_ADDITIONAL_ORDER
    // Optional experiment: for segments with equal YAtX at the boundary, reverse order.
    // This is related to the "node on boundary" phenomenon described in Degenerate-cases.md.
    // Left disabled by default because normal split/merge tends to reestablish correct order anyway.
    for (uint4 i = 0, j; i < L_size; i = j) {
      j = i + 1;
      while (j < L_size && 0 == (collection[L[i]].YAtX_frac(x) <=> collection[L[j]].YAtX_frac(x)))
        ++j;
      if (j > i + 1) {
        // Reverse order for sublists of equal-Y segments (passing segments should be in reverse order).
        std::reverse(L + i, L + j);
      }
    }
#endif
  }

  void AllIntCurLine(uint4 f, uint4& L_size, int4* L)
  {
    // Registers intersections on a multi-event vertical boundary X = const.
    // The list ENDS[f..] contains endpoints with that same X, terminated by list_stop.
    //
    // We handle:
    // - vertical-vertical intersections (stacks),
    // - vertical vs segments passing through the boundary (segments already in L),
    // - vertical vs endpoints lying on the boundary,
    // - plus same-Y begin/end combinations via LastFirstInt().
    auto vertical_segments = tmp; // Temporary storage for vertical segment indices.
    uint4 v_size = 0;

    auto non_vertical_pts = tmp + GetSegmNumb(); // Temporary storage for non-vertical endpoints on this X.
    uint4 nv_size = 0;

    for (uint4 i = f, cur_pt = ENDS[f]; cur_pt != list_stop; cur_pt = ENDS[++i]) {
      auto s = get_segm(cur_pt);
      if (collection[s].is_vertical()) {
        if (is_first(cur_pt))
          vertical_segments[v_size++] = s;
      }
      else {
        non_vertical_pts[nv_size++] = cur_pt;
      }
    }

    LastFirstInt(nv_size, non_vertical_pts);

    if (v_size ==0)
      return;

    auto prev_vert = vertical_segments[0];
    for (uint4 i = 1; i < v_size; ++i) { // vertical-vertical intersections
      auto cur_vert = vertical_segments[i];
      // fix me: check for zero-length vertical segments!!
      if (pts[first_point(cur_vert)].y == pts[last_point(prev_vert)].y) {
        register_pair(this, prev_vert, cur_vert);
        if((i>1)&& (pts[first_point(cur_vert)].y == pts[last_point(vertical_segments[i-2])].y))
          register_pair(this, vertical_segments[i - 2], cur_vert);
      }
      prev_vert = cur_vert;
    }

    for (uint4 i = 0, j = 0; j < L_size && i < v_size; ++i) { // vertical vs passing intersections
      auto s = vertical_segments[i];

      // Skip segments below the vertical segment.
      while (j < L_size && collection[L[j]].exact_under(pts[first_point(s)]))
        j++;

      // Register segments intersecting the vertical segment, except those ending exactly at its upper end.
      for (; j < L_size && collection[L[j]].exact_under(pts[last_point(s)]); ++j)
        register_pair(this, L[j], s);

      // Register segments passing exactly through the upper end of the vertical segment.
      // Keep j unchanged for the next vertical segment beginning at the same point.
      for (auto k = j; k < L_size && collection[L[k]].exact_on(pts[last_point(s)]); ++k)
        register_pair(this, L[k], s);
    }

    for (uint4 i = 0, j = 0; j < nv_size && i < v_size; ++i) { // vertical vs endpoint intersections
      auto s = vertical_segments[i];

      // Skip endpoints below the vertical segment.
      while (j < nv_size && pts[non_vertical_pts[j]].y < pts[first_point(s)].y)
        j++;

      // Register endpoints strictly inside the vertical segment.
      for (; j < nv_size && pts[non_vertical_pts[j]].y < pts[last_point(s)].y; ++j)
        register_pair(this, get_segm(non_vertical_pts[j]), s);

      // Register endpoints exactly at the upper end; keep j unchanged for the next vertical segment.
      for (auto k = j; k < nv_size && pts[non_vertical_pts[k]].y == pts[last_point(s)].y; ++k)
        register_pair(this, get_segm(non_vertical_pts[k]), s);
    }

  }

  void LastFirstInt(uint4 len, int4* nv_pts)
  {
    // For endpoints on the same X: if multiple non-vertical endpoints share the same Y,
    // and the block contains first-points followed by last-points, register all (begin,end) pairs.
    for (uint4 i = 0, j; i < len; i = j) {
      j = i + 1;
      while (j < len && pts[nv_pts[j]].y == pts[nv_pts[i]].y)
        ++j;
      if (j > i + 1) {
        auto first_beg = nv_pts + i;
        auto all_ends = nv_pts + j;
        auto first_end = std::find_if(first_beg, all_ends, +[](uint4 pt) { return is_last(pt); });
        assert(all_ends == std::find_if(first_end, all_ends, +[](uint4 pt) { return is_first(pt); }));
        if (first_end != all_ends)
          for (auto beg_pt = first_beg; beg_pt < first_end; ++beg_pt)
            for (auto end_pt = first_end; end_pt < all_ends; ++end_pt)
              register_pair(this, get_segm(*beg_pt), get_segm(*end_pt));
      }
    }
  }


  void DelStep(uint4 end_rank, uint4& L_size, int4* L)
  {
    uint4 new_size = 0;
    for (uint4 i = 0; i < L_size; ++i)
      if (seg_R_rank[L[i]]!= end_rank) {
        L[new_size++] = L[i];
      }
    L_size = new_size;
  }

  template <class Comp, class T>
  static auto merge_two_ranges(
    const T* it1, const T* last1,
    const T* it2, const T* last2,
    T* out,  Comp comp   // comp(x_from_1, y_from_2)
  ) {
    while (it1 != last1 && it2 != last2) 
      *out++ = comp(*it1, *it2) ? *it1++ : *it2++;
    while (it1 != last1)
      *out++ = *it1++;
    while (it2 != last2)
      *out++ = *it2++;
    return out; 
  }

  void InsStep(uint4 f, uint4& L_size, int4* L)
  {
    auto first_to_insert = tmp;
    auto last_to_insert = first_to_insert;
    auto first_L = tmp + GetSegmNumb();
    auto last_L = std::copy(L, L + L_size, first_L);
    for (uint4 i = f, cur_pt = ENDS[f]; cur_pt != list_stop; cur_pt = ENDS[++i]) {
      if (is_first(cur_pt)) {
        auto s = get_segm(cur_pt);
        if (!collection[s].is_vertical())
          *(last_to_insert++)  = s;
      }
    }
    auto is_below = [coll = this->collection, pts=this->pts](int4 sL, int4 sI) {
      auto& s = coll[sL];
      auto cmp = s.under(pts[first_point(sI)]);
      if (cmp != std::strong_ordering::equal)
        return cmp == std::strong_ordering::less;
      return 0 < s.shift % coll[sI].shift;
      };
    
    auto last_res = merge_two_ranges(first_L, last_L, first_to_insert, last_to_insert, L, is_below);
    L_size = last_res - L;
  }

  static void reg_pair(CIntegerSegmentCollection* c, uint4 s1, uint4 s2) {
    c->remaper.registrator->register_pair(s1, s2);
  };
  static void reg_remapped_pair(CIntegerSegmentCollection* c, uint4 s1, uint4 s2) {
    c->remaper.register_pair(s1, s2);
  };
  uint4  GetSegmNumb() const {
    return nSegments;
  };
  uint4 GetSegR(uint4 sn) const {
    return seg_R_rank[sn];
  };
  uint4 GetSegL(uint4 sn) const {
    return seg_L_rank[sn];
  };
  uint4 PointAtRank(uint4 rank) const {
    return ENDS[rank];
  }

  bool is_right_pt_on_bound(uint4 s) const {
    return stripe_right_rank == seg_R_rank[s];
  }
  bool is_left_pt_on_bound(uint4 s) const {
    return stripe_left_rank == seg_L_rank[s];
  }

  void SetCurStripe(uint4 left_rank, uint4 right_rank) {
    stripe_left_rank = left_rank;
    stripe_right_rank = right_rank;
    B = XAtRank(left_rank);
    E = XAtRank(right_rank);
  };

  void SetCurStripeLeft(uint4 left_rank) {
    stripe_left_rank = left_rank;
    B = XAtRank(left_rank);
  };

  void SetCurStripeRight(uint4 right_rank) {
    stripe_right_rank = right_rank;
    E = XAtRank(right_rank);
  };


  void SetCurPointAtBeg(uint4 s) {
    cur_point_seg = s;
    cur_point = collection[s].BegPoint();
  };

  void SetCurSeg(uint4 s){
    cur_seg_idx = s;
    cur_seg = collection[s];
    is_rstump=false;
  };

  void SetCurSegCutBE(uint4 s){
    SetCurSeg(s);
    stage = _Stages::stage_split;
    is_rstump = E < cur_seg.ex();
    active_end = cur_seg.EndPoint();
  };

  void SetCurSegCutBeg(uint4 s) {
    SetCurSeg(s);
    stage = _Stages::stage_split;
    active_end = cur_seg.EndPoint();
  };

  void SetCurSegCutEnd(uint4 s)
  {
    SetCurSeg(s);
    stage = _Stages::stage_merge;
    active_end = cur_seg.BegPoint();
    cur_seg_pt_on_right_bound = is_right_pt_on_bound(s);
  };

  void SetCurSegAndPoint(uint4 s) {
    SetCurPointAtBeg(s);

    SetCurSeg(s);
    stage = _Stages::stage_bubble;
    active_end = pts[last_point(s)];
  };

  constexpr bool static way_up = true;
  constexpr bool static way_down = false;

  template<bool is_left_bound>
  bool BoundIntersect(int4 s_, int4 X) const {
    assert(stage == _Stages::stage_split || stage == _Stages::stage_merge);
    auto cmp = collection[s_].YAtX_frac(X) <=> cur_seg.YAtX_frac(X);
    //  on the left bound always must be registered
    // on the right bound no registration except when the bound is defined by the right point of s_
    return cmp == std::strong_ordering::equal && 
      (is_left_bound || cur_seg_pt_on_right_bound || is_right_pt_on_bound(s_));
  };

  template<bool check_way_up_s_is_under_cur_seg>// s_  is_under or is_upper than cur_seg right stump
  bool RStumpIntersect(int4 s_) const {
    assert(stage == _Stages::stage_split);
    auto& s = collection[s_];

    auto cmpE = s.YAtX_frac(E) <=> cur_seg.YAtX_frac(E);// cmpE==1 is s.YAtX_frac(E) > cur_seg.YAtX_frac(E)
    if (cmpE == std::strong_ordering::equal) {// inresection on the right bound
      return is_right_pt_on_bound(s_);
    }

    if constexpr (check_way_up_s_is_under_cur_seg)
      return (cmpE == std::strong_ordering::less); //s < cur_seg
    else //check_way_down_s_is_upper_cur_seg
      return (cmpE == std::strong_ordering::greater); //s > cur_seg
  };

  constexpr bool static left_bound = true;
  constexpr bool static right_bound = false;

  template<bool is_left>
  bool XBelow(int4 s_1, int4 s_2, int4 X) const {
    auto& s1 = collection[s_1];
    auto& s2 = collection[s_2];//for right bound s2 always a stair
    auto cmp = s1.YAtX_frac(X) <=> s2.YAtX_frac(X);
    if (cmp != std::strong_ordering::equal)
      return cmp == std::strong_ordering::less;
    auto prod = s1.shift % s2.shift;
    assert(prod != 0);
    return is_left == prod > 0;//intersections is on the right of X (prod>0 would be on the left)
  }

  bool LBelow(int4 s_1, int4 s_2)const{  //retuns if s1 below s2 at current vertical line
    return XBelow<left_bound>(s_1, s_2, B);
  };

  bool RBelow(int4 s_1, int4 s_2)const{  //retuns if s1 below s2 at current vertical line
    return XBelow<right_bound>(s_1, s_2, E);
  };

  bool TrivCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
  {
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    if ((s1.ex() < s2.bx()) || (s2.ex() < s1.bx()))
      return false;//segments are not on vertical line and can't intersect
    return SSCurSegIntWith(s_);
  };


  template <bool reg = true>
  bool SSCurSegIntWith(int4 s_)//finds all intersection points of cur_seg and s (both are on some vertical line) and register them
  {
    //precondition segments share some x coords
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    if (s1.no_common_y(s2))
      return false;
    //at this moment, we can be sure that segments' x,y-projections intersect each other (1)
    // let l1 is the line through s1 and l2 is the line through s2
    auto delt = s2.BegPoint() - s1.BegPoint();
    auto prod = s1.shift % s2.shift;
    auto beg = s1.shift % delt;
    if (prod == 0) { // segments parallel
      if(always_exclude_zero_segs() && (s1.shift.is_zero()||s2.shift.is_zero()))
        return false;//zero segments are excluded

      if ((beg == 0) && (s2.shift % delt == 0)) {// they are on the same line
        if constexpr (reg)register_pair(this, cur_seg_idx, s_);
        return true;
      }
      return false;//segments are not on one line and can't intersect
    };
    auto end = beg + prod;
    if (((beg > 0) != (end < 0)) && (beg != 0) && (end != 0))//begin and end of s2 are on the same side of l1 
        return false;

    //auto int_type = _IntType::s2_beg_int * (beg == 0) + _IntType::s2_end_int * (end == 0);

    //at this moment, we can be sure that s2 intesects l1 (2)
    beg = delt % s2.shift;
    end = beg - prod;// if beg==0 && end == 0 -> prod == 0 impossible here  
    if ((beg == 0) ||//begin of s1 lies on l2, taking (1)(2) to account 
        (end == 0) ||//end of s1 lies on l2, taking (2) to account
        ((beg > 0) != (end > 0))) {//begin and end of s1 are on diffenent sides l2

      //int_type += _IntType::s1_beg_int * (beg == 0) + _IntType::s1_end_int * (end == 0);

      if constexpr (reg)register_pair(this, cur_seg_idx, s_);
      return true;
    };
    return false;
  };


  template<bool is_way_up>
  bool IsCurSegIntWith(int4 s_) {//checks if cur_seg and s intersects
    if (is_rstump) // can be in stage_split only
      return RStumpIntersect<is_way_up>(s_) || BoundIntersect<left_bound>(s_, B);
    if (stage == _Stages::stage_split)
      return ActiveEndIntersect<is_way_up>(s_) || BoundIntersect<left_bound>(s_, B);
    if (stage == _Stages::stage_merge)
      return ActiveEndIntersect<is_way_up>(s_) || BoundIntersect<right_bound>(s_, E);
    // no need to check BoundIntersect in stage_bubble but check cur_point intersection
    return ActiveEndIntersect<is_way_up>(s_) || CurPointIntersect(s_);
  };

  template<bool is_way_up>
  auto FindCurSegIntWith(int4 s_) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    auto cs = cur_seg_idx;
    if (is_rstump) {//can be in stage_split only
      if (RStumpIntersect<is_way_up>(s_) || BoundIntersect<left_bound>(s_, B)) {
        register_pair(this, cs, s_);
        return true;
      };
      return false;
    };
    if (stage == _Stages::stage_split) {
      if (ActiveEndIntersect<is_way_up>(s_) || BoundIntersect<left_bound>(s_, B)) {
        register_pair(this, cs, s_);
        return true;
      };
      return false;
    };
    if (stage == _Stages::stage_merge) {
      if (ActiveEndIntersect<is_way_up>(s_) || BoundIntersect<right_bound>(s_, E)) {
        register_pair(this, cs, s_);
        return true;
      };
      return false;
    };
    // no need to check BoundIntersect in stage_bubble but check cur_point intersection
    if (ActiveEndIntersect<is_way_up>(s_) || CurPointIntersect(s_)) {
      register_pair(this, cs, s_);
      return true;
    };
    return false;
  };

  template<bool is_way_up>
  auto FindCurSegIntWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    auto cs = cur_seg_idx;
    int4 increment = is_way_up ? 1 : -1;
    if (is_rstump) {//can be in stage_split only
      auto X = B;
      while ((s_ != last) && (RStumpIntersect<is_way_up>(*s_) || BoundIntersect<left_bound>(*s_, X))) {
        register_pair(this, cs, *s_);
        s_ += increment;
      }
      return s_;
    }
    if (stage == _Stages::stage_split) {
      auto X = B;
      while ((s_ != last) && (ActiveEndIntersect<is_way_up>(*s_) || BoundIntersect<left_bound>(*s_, X))) {
        register_pair(this, cs, *s_);
        s_ += increment;
      }
      return s_;
    }
    if (stage == _Stages::stage_merge) {
      auto X = E;
      while ((s_ != last) && (ActiveEndIntersect<is_way_up>(*s_) || BoundIntersect<right_bound>(*s_, X))) {
        register_pair(this, cs, *s_);
        s_ += increment;
      }
      return s_;
    }
    // no need to check BoundIntersect in stage_bubble
    // but check cur_point intersection
    while ((s_ != last) && (ActiveEndIntersect<is_way_up>(*s_) || CurPointIntersect(*s_))) {
      register_pair(this, cs, *s_);
      s_ += increment;
    }
    return s_;
  };

  bool IsIntersectsCurSegDown(int4 s_) //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return IsCurSegIntWith<way_down>(s_);
  };

  bool IsIntersectsCurSegUp(int4 s_) //check if cur_seg and s intersects (in the stripe b,e if cur_seg set in b,e) 
  {
    return IsCurSegIntWith<way_up>(s_);
  };

  auto FindCurSegIntDownWith(int4 s_) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_down>(s_);
  };

  auto FindCurSegIntUpWith(int4 s_) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_up>(s_);
  };

  auto FindCurSegIntDownWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_down>(s_, last);
  };

  auto FindCurSegIntUpWith(int4* s_, int4* last) {//finds all intersection points of cur_seg and s (in the stripe b,e if cur_seg set in b,e) and register them
    return FindCurSegIntWith<way_up>(s_, last);
  };


  bool UnderCurPoint(int4 s_) const { //returns true if s_ is under current point 
    //write as ActiveEndIntersect for consistency
    auto& s = collection[s_]; 
    auto pp = s.point_pos(cur_point);
    if(pp != 0)
      return pp < 0;//s_ placed under current point
    auto prod = collection[cur_point_seg].shift % s.shift;
    //check me: so far it works, but it can be that for zero cur_point_seg we need more accurate implementation !!!
    //assert((prod != 0));
    return prod < 0;//if prod<0 then s_ is clockwise from cur_seg at cur_point, 
    // taking into account that cur_seg must be on the right of cur_point 
    // (current point is always the begin point of cur_seg),
    //s_ is under cur_seg
  };

  template<bool is_way_up>
  bool ActiveEndIntersect(int4 s_) const { //
    auto pp = collection[s_].point_pos(active_end);
    if constexpr(is_way_up)
      return (pp <= 0);//s_ placed under active_end
    else
      return (pp >= 0);//s_ placed over active_end
  };

  bool CurPointIntersect(int4 s_) const { //
    return (collection[s_].point_pos(cur_point) == 0);//current point placed on  s_
  };


  void ResetRegistration()
  {
    if (factory) {// only main collection has factory, not clones
      if (GetRegistrator())//regstrar of main collection, clones will reset via factory
        GetRegistrator()->Reset();
      factory->Reset();
    }
  }

  void CombineRegData()
  {
    if (factory) {// only main collection has factory, not clones
      factory->combine_reg_data();
    }
  }

  void IntersectionsFindingDone()
  {// to be called after all intersections found
// empty function for the future use
  }

  void Reset()
  {
    MY_FREE_ARR_MACRO(ENDS);
    MY_FREE_ARR_MACRO(seg_L_rank);
    MY_FREE_ARR_MACRO(seg_R_rank);
    MY_FREE_ARR_MACRO(tmp);
  }

  void clone(CIntegerSegmentCollection &c, IntersectionRegistrator *r)
  {
      clone_of = &c;
      nSegments = c.nSegments;
      ENDS = c.ENDS;
      seg_L_rank = c.seg_L_rank;
      seg_R_rank = c.seg_R_rank;
      collection = c.collection;
      pts = c.pts;
      tmp = new int4[nSegments << 1];
      remaper.clone_from(&c.remaper);
      register_pair = c.register_pair;
      SetRegistrator(r);
  };

  void unclone() { 
    if (GetRegistrator())
      GetRegistrator()->Flash();

    if (clone_of == nullptr)return; 
    remaper.unclone();
    collection = nullptr; 
    pts = nullptr;
    clone_of = nullptr; 
    ENDS = nullptr;
    seg_L_rank = nullptr;
    seg_R_rank = nullptr;
    remaper.unclone();
    MY_FREE_ARR_MACRO(tmp);

    //unclone remaper !!!!!!!!

  };

  void SortAt(uint4 rank, uint4 n, int4 *L)
  {
    SetCurStripeLeft(rank);
    std::sort(L, L + n, [this](int4 s1, int4 s2) {return LBelow(s1, s2); });
  };

  void SetRegistrator(IntersectionRegistrator* r)
  {
    assert(("integer collection registers only intersecting pairs:\n\
       intersection points don't have integer coords ",
      (IntersectionRegistrator::reg_type & _RegistrationType::point) == 0));
    remaper.registrator = r;
  };

  IntersectionRegistrator* GetRegistrator() { return  remaper.registrator; };

  uint4 int_seg_from_real(uint4 n, TLineSegment1 sc[],int4 range)
  {
    // Converts real-coordinate segments into integer segments (scaled into [-range, range]).
    REAL si = 0., co = 1.;
    auto mm_rect = get_rot_minmax(n, sc, si, co);
    auto x_transf = transf1D{ .shift = -mm_rect.ld.x - mm_rect.get_width() / 2.,
      .scale = 2. * range / mm_rect.get_width() };
    auto y_transf = transf1D{ .shift = -mm_rect.ld.y - mm_rect.get_height() / 2.,
      .scale = 2. * range / mm_rect.get_height() };
    auto seg = set_size(segments, n);
    set_size(points, 2 * n);
#ifdef PRINT_SEG_AND_INT
    TIntegerSegment::coll_begin = coll_begin;
#endif
    TIntegerSegment s;
    uint4 i = 0;
    uint4 n_zero_seg = 0;
    for (uint4 m = 0; m != n; ++m) {
      s.Init(sc[m], si, co, x_transf, y_transf);
      assert(s.is_inside_int_range());
      if (s.shift.is_zero()) 
          ++n_zero_seg;
      seg[i] = s;
      points[first_point(i)] = s.BegPoint();
      points[last_point(i)] = s.EndPoint();
      ++i;
    }
    return i - n_zero_seg;
  }


  auto NoRemapInit(uint4 n, TLineSegment1* sc, IntersectionRegistrator* r, int4 range) {
    auto initial_SN = n;
    points = std::vector<TIntegerVect>(2 * n);
    auto nonzero_N = int_seg_from_real(n, sc,range);
    SetRegistrator(r);
    remaper.init(initial_SN, nonzero_N);
  };

  void NoRemapInit(uint4 n, TIntegerSegment* sc, IntersectionRegistrator* r, int4 range) {
    auto initial_SN = n;
    points = std::vector<TIntegerVect>(2 * n);
    uint4 i = 0;
    uint4 n_zero_seg = 0;
    auto seg = set_size(segments, n);
    for (uint4 m = 0; m != n; ++m) {
      assert(sc[m].is_inside_int_range());
      if (sc[m].shift.is_zero())
        ++n_zero_seg;
      seg[i] = sc[m];
      points[first_point(i)] = sc[m].BegPoint();
      points[last_point(i)] = sc[m].EndPoint();
      ++i;
    }
    auto nonzero_N = i - n_zero_seg;
    SetRegistrator(r);
    remaper.init(initial_SN, nonzero_N);
  };


  template <class Segment>
  void Init(uint4 n, Segment* c, IntersectionRegistrator* r, int4 range){
    assert(("integer collection registers only intersecting pairs:\n\
       intersection points don't have integer coords ",
      (IntersectionRegistrator::reg_type & _RegistrationType::point) == 0));
    NoRemapInit(n, c, r, range);

    // Precondition: no two distinct zero-length segments coincide as points.
    // This must hold in release too; otherwise event ordering becomes ambiguous
    // and downstream degenerate-case handling can misbehave.
    exclude_zero_segs = false;
    {
      std::vector<TIntegerVect> zero_pts;
      zero_pts.reserve(n);
      for (auto& s : segments) 
        if (s.shift.is_zero())
          zero_pts.push_back(s.BegPoint());

      if (zero_pts.size() > 1) {
        std::sort(zero_pts.begin(), zero_pts.end());

        auto it = std::adjacent_find(zero_pts.begin(), zero_pts.end(),
          [](TIntegerVect a, TIntegerVect b) { return a == b; });

        auto has_coinciding_zero_segs = (it != zero_pts.end());
#ifndef NDEBUG
        //        assert(!has_coinciding_zero_segs);
#endif

        if (has_coinciding_zero_segs)// Fallback: drop zero segments and proceed deterministically.
        {
          if (print_msg_on_coincide_zero())
            printf("CRemaper: coinciding zero-length segments encounter. Dropping all zero-length segments and proceeding.\n");
          exclude_zero_segs = true;
        }
      }
    }
    is_collection_remapped = false;
    register_pair = reg_pair;
    nSegments = remaper.get_N();
    collection = segments.data();
    pts = points.data();

  }

  template <class Segment>
  CIntegerSegmentCollection(const CollectionOptions& co, Segment* c, IntersectionRegistrator *r)
  {
    coincide_zero_segm_flags = co.coincide_zero_segm_flags;
    Init(co.n, c, r, co.range_for_int_seg);
  }

  void InitClone(uint4 n_threads) {
    if (factory)
      factory->InitClone(n_threads);
  }

  template <class Segment>
  CIntegerSegmentCollection(const CollectionOptions& co, Segment* c, CRegistratorFactory<IntersectionRegistrator>* f)
  {
    factory = f;
    factory->PrepareAlloc(co.n);
    coincide_zero_segm_flags = co.coincide_zero_segm_flags;
    Init(co.n, c, factory->GetRegistrator(0), co.range_for_int_seg);
  }

  void FromIntSegVect(std::vector<TIntegerSegment> &v, IntersectionRegistrator* r)
  {
    assert(("integer collection registers only intersecting pairs:\n\
       intersection points don't have integer coords ",
      (IntersectionRegistrator::reg_type & _RegistrationType::point) == 0));
    remaper.FromIntSegVect(v, r, *this);
    is_collection_remapped = false;
    register_pair = reg_pair;
    nSegments = remaper.get_N();
    collection = segments.data();
    pts = points.data();

  }


  bool TurnRemapOn() {
    if(exclude_zero_segs)
      is_collection_remapped = remaper.TurnRemapOn<false>(*this);
    else
      is_collection_remapped = remaper.TurnRemapOn<true>(*this);

#ifdef DEBUG_INTERSECTION_SET
    register_pair = reg_pair;
    is_collection_remapped = false;
#else
    register_pair = is_collection_remapped ? reg_remapped_pair : reg_pair;
#endif

    nSegments = remaper.get_N();
    collection = segments.data();
    pts = points.data();
    return is_collection_remapped;
  }

  TIntegerSegment& get_original_segment(uint4 s) {
    if (is_collection_remapped)
      return remaper.get_original_segment(s);
    return collection[s];
  }

  // Print difference vectors produced by multiset_difference.
  void print_difference_vectors(const std::vector<int_rec>& a_minus_b,
    const std::vector<int_rec>& b_minus_a,
    const char* nameA = "reg1",
    const char* nameB = "reg2")
  {
    printf("Pairs present in %s but not in %s (with multiplicity): %zu\n", nameA, nameB, a_minus_b.size());
    for (const auto& t : a_minus_b) {
      auto n1 = std::get<0>(t);
      auto n2 = std::get<1>(t);
      auto& s1 = get_original_segment(n1);
      auto& s2 = get_original_segment(n2);
      printf("%u %u:[(%i,%i)(%i,%i){%i,%i}] [(%i,%i)(%i,%i){%i,%i}]\n", n1, n2,
        s1.org.x, s1.org.y, s1.org.x + s1.shift.x, s1.org.y + s1.shift.y, s1.shift.x, s1.shift.y,
        s2.org.x, s2.org.y, s2.org.x + s2.shift.x, s2.org.y + s2.shift.y, s2.shift.x, s2.shift.y);
    }

    printf("Pairs present in %s but not in %s (with multiplicity): %zu\n", nameB, nameA, b_minus_a.size());
    for (const auto& t : b_minus_a) {
      auto n1 = std::get<0>(t);
      auto n2 = std::get<1>(t);
      auto& s1 = get_original_segment(n1);
      auto& s2 = get_original_segment(n2);
      printf("%u %u:[(%i,%i)(%i,%i){%i,%i}] [(%i,%i)(%i,%i){%i,%i}]\n", n1, n2,
        s1.org.x, s1.org.y, s1.org.x + s1.shift.x, s1.org.y + s1.shift.y, s1.shift.x, s1.shift.y,
        s2.org.x, s2.org.y, s2.org.x + s2.shift.x, s2.org.y + s2.shift.y, s2.shift.x, s2.shift.y);
    }
  }

  CIntegerSegmentCollection() = default;
 
  CIntegerSegmentCollection(CIntegerSegmentCollection& coll, IntersectionRegistrator* r)
  {
    clone(coll, r);
  }

  CIntegerSegmentCollection(CIntegerSegmentCollection& coll, uint4 thread_index)
  {
    if (coll.factory)
      clone(coll, coll.factory->GetRegistrator(thread_index));
    else
      clone(coll, coll.GetRegistrator());
  }

  ~CIntegerSegmentCollection()
  {
    unclone();
    Reset();
  }


 
  void coll_to_SVG(chostream* SVG_stream) {
    if (!SVG_stream)return;
    int4 n = MIN(max_SVG_items, nSegments);
    auto mmr = get_mmrect01(collection, n);
    *SVG_stream << "<svg height='100%' width='100%' viewBox='";
    *SVG_stream << mmr.ld.x*crd_scaler << " " << mmr.ld.y * crd_scaler << " "
      << (mmr.rt.x - mmr.ld.x) * crd_scaler << " " << (mmr.rt.y - mmr.ld.y) * crd_scaler << "' transform='scale(1, -1)'>\n";
    for (int4 i = 0; i < n; ++i)collection[i].write_SVG(i, SVG_stream);
  };
  auto &get_segments() {
    return segments;
  };

  auto& get_ini_segments() {
    if(is_collection_remapped)
      return remaper.get_seg_v();
    return segments;
  };


  auto XAtRank(uint4 rank) const {
    auto pt = PointAtRank(rank);
    return GetXEx(pt);
  };
  bool always_exclude_zero_segs() const {
    return exclude_zero_segs && 
      ((coincide_zero_segm_flags & coincide_zero_filter_type)== coincide_zero_filter_for_all);
  };

  bool print_msg_on_coincide_zero() const {
    return (coincide_zero_segm_flags & coincide_zero_print_warning) == coincide_zero_print_warning;
  };

private:
  auto GetXEx(uint4 pt) const {
    if (is_multiple(pt)){
      auto f = unmark_multiple(pt);
      return pts[ENDS[f]].x;
    };
    return GetX(pt);
  };

  auto GetX(uint4 pt) const {
    assert(pt < 2 * nSegments);
    return pts[pt].x;
  };

  TIntegerSegment* collection = nullptr;
  TIntegerVect* pts = nullptr;
  uint4* seg_L_rank = nullptr, * seg_R_rank = nullptr, * ENDS = nullptr;
  CRegistratorFactory<IntersectionRegistrator>* factory = nullptr;

  TIntegerSegment cur_seg;
  TIntegerVect cur_point, active_end;
  std::vector<TIntegerSegment> segments;
  std::vector<TIntegerVect> points;
  CRemaper remaper;
  decltype(reg_pair) *register_pair = reg_pair;
  CIntegerSegmentCollection *clone_of = nullptr;
  uint4 nSegments=0,nTotX=0,nCollideX=0;
  int4 B, E,  stripe_left_rank, stripe_right_rank;
  uint4 stage = _Stages::stage_split;

  uint4 cur_seg_idx = 0xFFFFFFFF;// , cur_point_idx = 0xFFFFFFFF;
  uint4 cur_point_seg = 0xFFFFFFFF;// , active_end_idx = 0xFFFFFFFF;
  uint4 coincide_zero_segm_flags = _CoincideZeroSegmFlags::coincide_zero_filter_for_all;
  bool is_rstump = false;
  bool is_collection_remapped = false;
  bool cur_seg_pt_on_right_bound = false;
  bool exclude_zero_segs = false;
//  bool rbelong_to_stripe=true;
};

