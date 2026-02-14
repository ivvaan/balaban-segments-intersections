#pragma once
/*
*
*      Copyright (c)  2011-2026  Ivan Balaban
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

#include "segments.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <vector>
#include <numeric> 
// CRemaper is a degenerate-case preprocessor for integer segments.
//
// It normalizes a set of potentially overlapping / collinear segments by splitting
// shared subsegments into separate pieces. After normalization:
// - no two distinct segments overlap along a 1D interval,
// - only endpoint-level degeneracies remain (shared endpoints, multi-events handled elsewhere).
//
// The sweep algorithms run on the normalized collection, and reported intersecting pairs are
// expanded back into pairs of original segments via `register_pair(...)`.

template<class IntersectionRegistrator>
class CRemaper {

  struct MarkedListBounds {
    uint4 beg = 0;
    uint4 end = 0;
    uint4 flags = 0;

    bool is_mapped_entirely() const
    {
      // Non-zero means: the remapped segment corresponds exactly to a single original segment.
      // (A fast path for one-to-one mapping.)
      return (flags != 0);
    }
  };

  template<typename Collection>
  bool prepare_remap(std::vector<uint4>& indexes, TIntegerVect points[], Collection& collection)
  {
    std::vector<TIntegerSegment>& res_coll = collection.get_segments();
    std::vector<TIntegerVect>& res_pts = collection.get_points();
    // Input:
    // - `indexes`: endpoint ids sorted by the comparator from TurnRemapOn().
    //   Endpoint id is an encoded integer where:
    //     - `Collection::get_segm(pt)` gives the segment index,
    //     - `Collection::is_last(pt)` distinguishes start/end.
    // - `points[pt]`: the integer endpoint coordinate for that encoded endpoint id.
    //
    // Output:
    // - `res_coll` / `res_pts`: normalized (remapped) segments and their endpoints.
    // - `rrec[new_seg]`: [beg,end) slice inside `remap_v` telling which original segments
    //   contributed to the normalized segment `new_seg`.
    //
    // Key invariant relied upon by the sweep code:
    // After normalization, on a given supporting line there are no overlapping 1D intervals.
    // Any maximal collinear overlap is split into disjoint pieces; each piece remembers
    // the list of original segments that cover it.

    uint4 N = indexes.size();
    uint4 remap_size = 0; // total number of back-references that will be stored in `remap_v`
    remapped_SN = 0;      // total count of produced (normalized) segments

    // Pass #1:
    // Sweep over the sorted endpoints and determine:
    // - how many normalized segments will be created (`remapped_SN`),
    // - how many original-segment references we will need in total (`remap_size`).
    //
    // `size` is the current active set size on the current supporting line:
    // it is incremented at "first" endpoints and decremented at "last" endpoints.
    //
    // For a non-degenerate interval between `indexes[i]` and `indexes[i+1]`:
    // - if `size>0` and `points differ`, there is a real 1D piece covered by `size` segments.
    //
    // Degenerate (zero-length) intervals are normally ignored, BUT:
    // - for true zero segments (segment length == 0) we *do* want to keep a piece even though
    //   `beg==end`, because that piece represents an actual input segment.
    //   This is handled by `allow_zero_len_piece` below.
    bool not_remapped = true;// initial_SN == nonzero_N;
    for (uint4 i = 0, size = 0; i < N; ++i) {
      auto pt = indexes[i];
      size += Collection::is_last(pt) ? -1 : 1;
      assert(size >= 0);

      if (size != 0) {
        assert(i + 1 < N);
        auto next_pt = indexes[i + 1];

        // If pt and next_pt belong to the same original segment, allow zero-length piece.
        // For non-zero segments this is redundant (their endpoints differ anyway); for zero segments it's required.
        auto allow_zero_len_piece = (Collection::get_segm(pt) == Collection::get_segm(next_pt));

        // If there is a real-length interval OR we explicitly allow a zero-length one,
        // then it will become a normalized segment piece.
        if (allow_zero_len_piece || (points[pt] != points[next_pt])) {
          // If `size>=2`, then we are splitting/normalizing overlapping parts => remap happens.
          not_remapped = not_remapped && size < 2;

          // This normalized piece will carry `size` original segment ids in `remap_v`.
          remap_size += size;
          ++remapped_SN;
        }
      };
    }

    // If normalization did not change anything (no overlaps and no split points introduced),
    // we can return without producing remap tables.
    not_remapped = not_remapped && (remapped_SN == nonzero_N);
    if (not_remapped) {
      remapped_SN = initial_SN;
      return true;//not remapped
    }

    // Allocate remap tables:
    // - `rrec_v` has one record per produced normalized segment (remapped_SN).
    // - `remap_v` is a flat array of original segment indices; each normalized segment points
    //   to its slice via `rrec[new_seg].beg/end`.
    rrec = set_size(rrec_v, remapped_SN);
    auto remap = set_size(remap_v, remap_size);

    // Allocate normalized geometry buffers.
    set_size(res_coll, remapped_SN + 2);// +2 for sentinels
    set_size(res_pts, 2 * remapped_SN);

    // Pass #2:
    // Re-run the sweep, maintain the active set, and actually emit normalized segments.
    //
    // `stack` holds the currently active original segment indices on the current line.
    // For each emitted normalized piece we:
    // - copy `stack` into `remap_v`,
    // - store the slice bounds into `rrec`,
    // - store the geometry (beg/end points) into `res_coll` / `res_pts`.
    uint4 remaper_pos = 0; // current write position inside `remap_v`
    uint4 new_seg_num = 0; // current normalized segment index

    std::vector<uint4> stack;
    stack.reserve(nonzero_N);

    for (uint4 i = 0; i < N; ++i) {
      auto pt = indexes[i];

      // Update active set.
      if (Collection::is_last(pt)) {
        remove_by_val(stack, Collection::get_segm(pt));
      }
      else {
        auto s = Collection::get_segm(pt);

        // When a segment becomes active on a line, it overlaps the ones already active on that line.
        // That overlap is a "collinear intersection" and is easiest to register right here
        // because we still have the full 1D overlap context.
        for (auto p : stack)
          registrator->register_pair(s, p);

        stack.push_back(s);
      }

      auto size = (uint4)stack.size();
      if (size != 0) {
        assert(i + 1 < N);

        auto next_pt = indexes[i + 1];
        auto beg = points[pt];
        auto end = points[next_pt];

        // Same special-case: keep a zero-length piece only when it is exactly a zero segment's endpoints.
        // All other zero-length "gaps" between unrelated endpoints are ignored.
        auto allow_zero_len_piece = (Collection::get_segm(pt) == Collection::get_segm(next_pt));

        if (allow_zero_len_piece || (beg != end)) {
          // Store back-references to original segments.
          std::copy(stack.begin(), stack.end(), remap + remaper_pos);

          // flags: for fast path in register_pair(s1,s2):
          // If a normalized segment maps to exactly one original segment AND corresponds to
          // that original segment entirely, we can skip list expansion later.
          //
          // For an allowed zero-length piece (which is always "entirely" that segment),
          // this is true whenever only one segment is active.
          uint4 is_mapped_entirely = allow_zero_len_piece && (size == 1);

          rrec[new_seg_num] = { remaper_pos, remaper_pos + size, is_mapped_entirely };
          remaper_pos += size;

          // Store normalized geometry.
          res_coll[new_seg_num] = { beg, end };
          res_pts[Collection::first_point(new_seg_num)] = beg;
          res_pts[Collection::last_point(new_seg_num)] = end;

          ++new_seg_num;
        }
      }
    }

    // Consistency checks.
    assert(remapped_SN == new_seg_num);
    remap_size = remaper_pos;

    // If we got here, remap definitely happened.
    return not_remapped;
  }


public:

  template<typename Collection>
  auto TurnRemapOn(Collection& collection) {
    auto n = initial_SN;

    auto points = std::move(collection.get_points());
    seg_v = std::move(collection.get_segments());
    seg = seg_v.data();
    std::vector<uint4> indexes(2 * n);

    std::iota(indexes.begin(), indexes.end(), 0);

    nonzero_N = n;

    auto comparator = [segments = seg, pts = points.data()](uint4 pt1, uint4 pt2) {
      auto i1 = Collection::get_segm(pt1);
      auto i2 = Collection::get_segm(pt2);
      if (i1 == i2)
        return pt1 < pt2; //same segment

      auto& s1 = segments[i1];
      auto& s2 = segments[i2];

      // Treat "zero segments" (degenerate to a point) as vertical for grouping purposes,
      // and force their endpoints to come first when point coordinates coincide.
      auto is_zero1 = s1.shift.is_zero();
      auto is_zero2 = s2.shift.is_zero();

      TIntegerVect shift1 = is_zero1 ? TIntegerVect(0, 1) : s1.shift;
      TIntegerVect shift2 = is_zero2 ? TIntegerVect(0, 1) : s2.shift;

      auto S = shift1 % shift2;
      if (S != 0)//segments non parallel
        return S < 0;

      //segments parallel
      auto shift = shift1;// +shift2;
      assert(shift.is_non_zero());
      auto oo = pts[pt2] - pts[pt1];
      S = shift % oo;
      if (S != 0)// segments on different lines
        return  S < 0;

      //segments parallel and lies on one line
      auto prod = oo * shift;
      if (prod != 0)// points not coinside
        return  0 < prod;
      // here points coinside, but segments can be different (zero-seg can coincide with non-zero seg)
      bool last1 = Collection::is_last(pt1);
      bool last2 = Collection::is_last(pt2);
      if (last1 != last2) {
        //return (is_zero1 && is_zero2) ? last1 < last2 : last1 > last2; //different segments end first exept when both zero, then start first
        return (is_zero1 && is_zero2) == last2; //equal to previous line, but faster to compute

      }
      // Both endpoints are at the same coordinate and are of the same kind (both first or both last).
      // Among FIRST endpoints: zero first. Among LAST endpoints: non-zero first.
      // here Collection::is_last(pt1) == Collection::is_last(pt2)
      if (is_zero1 != is_zero2)
        //return last2 ? (is_zero1 < is_zero2) : (is_zero1 > is_zero2);
        return last2 == is_zero2;//equal to previous line, but faster to compute
      // Both endpoints are at the same coordinate, are of the same kind, and either both zero or both non-zero.
      return last2 ? pt1 > pt2:pt1 < pt2; //  LAST endpoints in reverse order to FIRST. As a result segments form piramid like structure.
        //      |---------------|
        //   |------------------------|
        // |---------------------------------|
      };

    std::sort(indexes.begin(), indexes.end(), comparator);


    bool not_remapped = prepare_remap(indexes, points.data(), collection);
    if (not_remapped) {
      collection.get_segments() = std::move(seg_v);
      collection.get_points() = std::move(points);
      return false;
    }
    remapped_segs = collection.get_segments().data();
    remapped_segs[remapped_SN] = seg[n];//bottom sentinel
    remapped_segs[remapped_SN + 1] = seg[n + 1];//top sentinel
    remapped_ends = collection.get_points().data();

    return true;
  };


  void register_pair(uint4 s1, uint4 s2, uint4 int_type) {
    if(s1 == s2)
      return;
    auto orig_int_type = seg[s1].get_int_type_beg(seg[s2]);
    if (int_type == orig_int_type)
      registrator->register_pair(s1, s2);
  }

  void register_pair(uint4 s1, uint4 s2) {
    auto rr1 = rrec[s1];
    auto rr2 = rrec[s2];

    if (rr1.is_mapped_entirely() && rr2.is_mapped_entirely()) {
      // If both remapped segments map to exactly one original segment each,
      // we can register the corresponding original pair directly.
      registrator->register_pair(remap_v[rr1.beg], remap_v[rr2.beg]);
      return;
    }

    auto int_type = remapped_segs[s1].get_int_type_beg(remapped_segs[s2]);

    if (int_type == _IntType::common_int) {
      // Common (non-degenerate) intersection of remapped segments.
      // Expand to original pairs via the remap lists.

      if (rr1.end + rr2.end - rr1.beg - rr2.beg < 3) { // one-to-one mapping
        registrator->register_pair(remap_v[rr1.beg], remap_v[rr2.beg]);
        return;
      }

      for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1) { // many-to-many mapping
        for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
          registrator->register_pair(remap_v[seg1], remap_v[seg2]);
      }
      return;
    }

    // Non-common intersection: need to filter by intersection type to avoid double registration.
    if (rr1.end + rr2.end - rr1.beg - rr2.beg < 3) {
      register_pair(remap_v[rr1.beg], remap_v[rr2.beg], int_type);
      return;
    }

    // Most complex case: possible duplicates inside remap lists; mark duplicates to avoid double registration.
    constexpr const uint4 mark = 1 << (sizeof(uint4) * 8 - 1);
    for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
      for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
        if (remap_v[seg1] == remap_v[seg2]) { // same original segment appears in both lists
          remap_v[seg1] |= mark;
          remap_v[seg2] |= mark;
          break;
        }

    //register only intersections with nonmarked segments, otherwise we can have double registration
    for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
      if ((remap_v[seg1] & mark) == 0) {
        for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
          if ((remap_v[seg2] & mark) == 0)
            //also need to check type of intersection to exclude double registration
            register_pair(remap_v[seg1], remap_v[seg2], int_type);
      }

    // erase marks, back to original values
    for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
      remap_v[seg1] &= ~mark;
    for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
      remap_v[seg2] &= ~mark;
  };

  auto init(uint4 _initial_SN, uint4 _nonzero_N) 
  {
    initial_SN = _initial_SN;
    nonzero_N = _nonzero_N;
    remapped_SN = initial_SN;
    return;
  }

  auto& get_seg_v() {
    return seg_v;
  }

  TIntegerSegment& get_original_segment(uint4 s) {
    return seg[s];
  }

  auto get_N() {
    return remapped_SN;
  }

  void clone_from(CRemaper* other) {
    initial_SN = other->initial_SN;
    nonzero_N = other->nonzero_N;
    remapped_SN = other->remapped_SN;
    remap_size = other->remap_size;
    //range = other->range;
    registrator = other->registrator;
    remap_v = other->remap_v;//the data mast be copied when cloning (read\write access)
    //just copy pointers (for read only access)
    rrec = other->rrec;
    seg = other->seg;
    remapped_segs = other->remapped_segs;
    remapped_ends = other->remapped_ends;
  };

  void unclone() {
    rrec = nullptr;
    remapped_segs = nullptr;
    remapped_ends = nullptr;
    registrator = nullptr;
  };

IntersectionRegistrator* registrator = nullptr;

private:
  uint4 initial_SN = 0;
  uint4 nonzero_N = 0;
  uint4 remapped_SN = 0;
  uint4 remap_size = 0;
  //int4 range = 1;
  std::vector<MarkedListBounds> rrec_v;
  std::vector<uint4> remap_v;//mast be copied when cloning
  std::vector<TIntegerSegment> seg_v;//if remapped original int segments are stored here
  MarkedListBounds* rrec = nullptr;
  TIntegerSegment* seg = nullptr;
  TIntegerSegment* remapped_segs = nullptr;
  TIntegerVect* remapped_ends = nullptr;
  //CRemaper* clone_of = nullptr;

};



















