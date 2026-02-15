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
// New implementation: "optimal" variant.
//
// This file implements an optimized version of the recursive "ladder/staircase" approach
// from `CFastIntFinder`.
//
// Key idea (vs. `CFastIntFinder`):
// - When the parent staircase is large enough, the algorithm *inherits* every `inherit_each`-th step
//   from the parent staircase into the child staircase.
// - `father_loc[]` stores an approximate mapping between a location (between/at steps) in a child
//   staircase and a nearby location in the parent staircase, which reduces the search range when
//   bubbling newly inserted segments through higher-level staircases.
//
// See wiki pages: "Optimal-algorithm", "New-implementation", "Sentinels-and-guards", "Split".

#include "fast_finder.h"
#include <algorithm>
#include <vector>
#include <thread>
#include <cassert>


// `COptimalIntFinder` is a drop-in replacement for `CFastIntFinder` for collections that satisfy
// `OptimalSegColl`.
//
// It reuses most of the fast finder infrastructure (buffers `L/R/Q`, stripe recursion, event handling),
// but adds:
// - inherited stairs (sampling from parent staircase),
// - `father_loc` mapping to accelerate internal intersection checks (`FindIntI`).
//
// Typical usage:
// - call `prepare_ends(segments)` inherited from `CFastIntFinder`
// - then call `find_intersections(segments)`.
class COptimalIntFinder : public CFastIntFinder
{
public:
  using CFAST = CFastIntFinder;

  ~COptimalIntFinder() { FreeMem(); };
  void FreeMem() {};

  // Runs the optimal algorithm on the full x-range.
  //
  // Preconditions:
  // - `prepare_ends(segments)` was called (from base class) to initialize `LR_len`, `nTotX`, etc.
  // - `segments` satisfies `OptimalSegColl` and supports the additional operations used by this file
  //   (in particular: `IsIntersectsCurSegUp/Down`, and the same cut/stripe protocol as in the fast finder).
  //
  // Postconditions:
  // - All intersections are reported via `segments` registrars/callbacks.
  // - `segments.IntersectionsFindingDone()` is called exactly once.
  template<OptimalSegColl SegmentsColl>
  void find_intersections(SegmentsColl &segments)
  {
    // Allocate working buffers.
    // `Q` is larger than in `CFastIntFinder` because inherited stairs increase the number of stored steps.
    // The formula below allocates slack for (~1 inherited per `inherit_each-1`) plus some constants.
    len_of_Q = LR_len + LR_len / (inherit_each - 1) + inherit_each + 2;
    DECL_RAII_ARR(L, LR_len+1);
    ++L;//to have one cell before L for sentinel
    DECL_RAII_ARR(R, LR_len+1);
    ++R;//to have one cell before R for sentinel
    DECL_RAII_ARR(Q, len_of_Q);

    // `father_loc` maps a location index in the current staircase to an approximate location in the parent staircase.
    // Sign convention:
    // - `father_loc[i] <= 0` => original step (not inherited).
    // - `father_loc[i]  > 0` => inherited step; value is the parent step index from which it was inherited.
    DECL_RAII_ARR(father_loc, len_of_Q);
    std::fill_n(father_loc, len_of_Q, undef_loc);

    // `bottom_index` reserves a small prefix in `Q` so that the first generated staircase starts
    // at a position greater than zero (important for location math and the `undef_loc` sentinel value).
    constexpr uint4 bottom_index = inherit_each + 1;
    ProgramStackRec stack_rec(bottom_index, nTotX);  // Must be initialized this way (top frame).

    // Initialize stripe to the very first endpoint rank, build initial active set, and start recursion.
    segments.SetCurStripeRight(0);
    InsDel(*this,segments, 0, &stack_rec);

    // Start with `father_first_step = bottom_index + 1` (a valid "first inherited" candidate location).
    FindR(segments, bottom_index + 1, bottom_index, 0, nTotX - 1, &stack_rec, 0);

    InsDel(*this,segments, nTotX - 1, &stack_rec);
    segments.IntersectionsFindingDone();
 }

  // Optimized "bubble up" internal intersection checks, using `father_loc` to shrink search ranges.
  //
  // Compared to `CFastIntFinder::FindIntI`, this method attempts to reuse the *approximate* location
  // found in the previous (child) staircase as a hint for the search range in the parent staircase.
  //
  // Preconditions:
  // - `segments.SetCurSeg4Bubble(cur_seg)` was called by the caller (`InsDel` in the base class).
  // - `stack_pos` points to the current recursion frame corresponding to where the insertion happened.
  //
  // Postconditions:
  // - Intersections between `cur_seg` and staircases on all higher recursion levels are reported.
  template<class SegmentsColl>
  void FindIntI(SegmentsColl& segments, uint4 cur_seg, ProgramStackRec* stack_pos) const
  {
    auto r_index = segments.GetSegR(cur_seg);
    while (stack_pos->right_bound <= r_index)stack_pos = stack_pos->prev;// Go from bottom to top and find staircase to start.
    if (stack_pos->prev == nullptr)return;

    // `l` is a "location index" (between/at steps) carried across levels.
    // `QE` is the upper location bound of the current staircase (one past last step location).
    int4  l = undef_loc, r, m, QE = stack_pos->Q_pos + 1;
    int4 QB = undef_loc + 1;// Just to differ from `l` at the first loop.

    for (stack_pos = stack_pos->prev; stack_pos != nullptr; stack_pos = stack_pos->prev)// For all staircases above.
    {
      if ((big_staircase_threshold > inherit_each) && (l == QB)) // If location is below the first stair, we don't have data in father_loc for that location
        // because the number of locations is greater (by one) than the number of stairs.
        // (For example, one stair has two locations (below and above), but we keep in father_loc data for only one (above).)
        // So we use data from the next location: QE == l + 1.
      {
        l = QB = stack_pos->Q_pos;// Set lower location bound to the location below the first step of current staircase.
        r = father_loc[QE] ? QB + inherit_each : QE;// `father_loc[QE] != undef_loc` means staircase was created by `Split()` (optimal variant with inheritance).
        // If created by optimal `Split`, search only in [QB..QB+inherit_each]; otherwise use full range [QB..QE].
      }
      else
      {
        QB = stack_pos->Q_pos;

        // Use `father_loc` to get an approximate location in the parent staircase.
        // This line can run when `l == undef_loc`, so `father_loc[undef_loc]` must be `undef_loc`.
        m = abs(father_loc[l]);

        // Clamp to staircase bounds. If mapping is missing (`m == undef_loc`) or outside, fall back to QB.
        l = MAX(QB, m);
        r = l + inherit_each;
        if ((m == undef_loc) || (r > QE))r = QE;
        // If `m == undef_loc` we use [QB..QE] as a range.
      }

      // Binary search in [l..r) to locate location near current y-position.
      auto len = r - l;
      auto Ql = Q + l;
      while (len > 1) // Binary search.
      {
        m = len / 2;
        if (segments.UnderCurPoint(Ql[m]))
          Ql += m;
        len -= m;
      }

      // Store the refined location for the next (parent) iteration.
      l = Ql - Q;

      // Scan intersections around the found location in the current staircase slice.
      FindInt(segments, QB, QE - 1, l);

      // Move upper bound to the location above the last step of the previous staircase
      // (which becomes the current staircase on the next iteration).
      QE = QB + 1;
    };
  };



private:
  // `undef_loc` is the "no location hint" marker.
  // It is intentionally 0 to keep `abs(undef_loc) == undef_loc` and minimize special-case checks.
  static constexpr const int4 undef_loc = 0;

  // Limits the number of repeated recursive calls (`dont_cut_stripe` case) before forcing a cut.
  static constexpr const uint4 max_call = 32; // Max number of sequential recursive calls of (opt)FindR before dividing current stripe.

  // Inheritance configuration:
  // - Inherit one parent stair per `inherit_each` steps (sampling rate).
  // - `inherit_offset` picks the first inherited step (must be in [1..inherit_each-1]).
  static constexpr const uint4 inherit_each = 15; // In the optimal algorithm, one inherited step per `inherit_each`.
  static constexpr const uint4 inherit_offset = inherit_each / 2; // First stair to inherit; must be in [1..inherit_each-1].
  // Inherited stair positions are: inherit_offset, inherit_offset + inherit_each, inherit_offset + 2*inherit_each, ...

  // Switch threshold:
  // - Use the optimal (inheritance) variant only if the parent staircase is large enough.
  static constexpr const uint4 big_staircase_threshold = 64; // Used in (opt)FindR.

  // `father_loc` stores mapping/origin metadata per step/location for the current staircase region in `Q`.
  // Must be `int4` to store both negative (original) and positive (inherited) values.
  int4* father_loc = nullptr;

  // Helper: original steps have `father_loc <= 0` (including `undef_loc`).
  static bool is_original(int4 v) { return v < 1; };

  // Intersection scan around a location `l` inside staircase slice `[qb..qe]`.
  //
  // This variant differs from `CFastIntFinder::FindInt`:
  // - For inherited steps, it only *tests* intersection (`IsIntersectsCurSeg*`) without registering,
  //   because inherited steps are duplicates of parent steps and would otherwise cause duplicate reporting.
  // - For original steps, it performs and registers intersections (`FindCurSegInt*With`).
  //
  // Returns:
  // - Number of stairs "touched" (a proxy for number of potential intersections), used by heuristics.
  //
  // Preconditions:
  // - Current segment is set in `segments` (via `SetCurSeg...` protocol).
  template<class SegmentsColl>
  uint4 FindInt(SegmentsColl& segments, uint4 qb, uint4 qe, uint4 l) const
  {

    auto c = l;
    while (c>qb)
    {
      if (is_original(father_loc[c]))
      { // If original step: try to find and register intersection.
        if (!segments.FindCurSegIntDownWith(Q[c])) break;
      }
      else // If inherited step: only test intersection (do not register).
        if (!segments.IsIntersectsCurSegDown(Q[c])) break;
      c--;
    }

    constexpr const bool line_seg = (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) == _Coll_flag_state::state_true);

    if (line_seg && (c != l))
      return l-c;

    ++c; ++l;
    while (l <= qe)
    {
      if (is_original(father_loc[l]))
      {
        if (!segments.FindCurSegIntUpWith(Q[l])) break;
      }
      else
        if (!segments.IsIntersectsCurSegUp(Q[l])) break;
      l++;
    }
    return l-c;
  };

  // Merge step for the optimal variant.
  //
  // Differences from `CFastIntFinder::Merge`:
  // - Only original staircase steps are emitted into the merged active set (`L`),
  //   inherited steps are skipped to avoid duplicating segments in `L`.
  // - Uses `FindInt(segments, QB, QE, ...)` which respects original/inherited semantics.
  //
  // Preconditions:
  // - `Q[QB+1..QE]` contains the staircase for the current recursion level.
  // - `father_loc[QB+1..QE]` matches `Q` and marks steps as original/inherited.
  //
  // Postconditions:
  // - `L` becomes the merged active set at the stripe right boundary, containing no inherited duplicates.
  template<class SegmentsColl>
  void Merge(SegmentsColl &segments, uint4 LBoundIdx, uint4 QB, uint4 QE)
  {
    auto Size = L_size;
    uint4 cur_R_pos = 0, new_size = 0;
    auto cur_stair = QB;
    auto _R = L, _L = R;
    auto cur_seg = _R[cur_R_pos];
    while ((cur_stair<QE) && (cur_R_pos<Size))
    {
      if (segments.RBelow(cur_seg, Q[cur_stair + 1]))
      {
        if (segments.GetSegL(cur_seg) >LBoundIdx)
        {
          segments.SetCurSegCutEnd(cur_seg);
          FindInt(segments, QB, QE, cur_stair);
        }
        _L[new_size++] = cur_seg;
        cur_seg = _R[++cur_R_pos];
      }
      else
        if (is_original(father_loc[++cur_stair]))
          _L[new_size++] = Q[cur_stair];
    }
    while (cur_R_pos<Size)
    {
      if (segments.GetSegL(cur_seg) >LBoundIdx)
      {
        segments.SetCurSegCutEnd(cur_seg);
        FindInt(segments, QB, QE, QE);
      }
      _L[new_size++] = cur_seg;
      cur_seg = _R[++cur_R_pos];
    }
    while (cur_stair<QE)
      if (is_original(father_loc[++cur_stair]))
        _L[new_size++] = Q[cur_stair];
    L = _L; R = _R;
    L_size = new_size;
  };

  // Optimal split: builds a new staircase by mixing:
  // - original steps from current `L` (candidate segments that cover the stripe and do not intersect the ladder),
  // - inherited steps sampled from the parent staircase (`Q[cur_father_pos]` at stride `inherit_each`).
  //
  // It also records, for each segment kept in `L`, the location (step index) in the new staircase
  // into `location[]` (stored in `R` buffer), so that we can quickly call `FindInt` for that segment
  // after staircase construction.
  //
  // Parameters:
  // - `RBoundIdx`: stripe right boundary rank.
  // - `cur_father_pos`: initial step in the parent staircase to consider for inheritance.
  // - `_step_index`: in/out; initial last-step index, updated to new last-step index.
  //
  // Preconditions:
  // - `L[0..L_size)` is sorted at the stripe left boundary (same as in base algorithm).
  // - `Q` currently contains the parent staircase entries used for inheritance.
  //
  // Postconditions:
  // - `Q` is extended to contain the new staircase steps up to `_step_index`.
  // - `father_loc` is populated for each created step.
  // - `L` is compacted to contain only the segments that must remain in the active set.
  // - `dont_cut_stripe` is updated based on observed intersection density.
  template<class SegmentsColl>
  void Split(SegmentsColl& segments, uint4 RBoundIdx, uint4 cur_father_pos, uint4& _step_index)
  {
    auto Size = L_size;
    auto step_index = _step_index;
    auto father_last_step = step_index;
    uint4 new_L_size = 0, cur_L_pos = 0;
    uint4 cur_seg;
    auto location = R;
    while ((cur_L_pos<Size) && (cur_father_pos <= father_last_step))
    {
      cur_seg = L[cur_L_pos];
      segments.SetCurSegCutBE(cur_seg);
      if (segments.LBelow( cur_seg, Q[cur_father_pos])) // Current segment is below candidate step to inherit.
      {
        if ((segments.GetSegR(cur_seg) < RBoundIdx)// Current segment does not cover stripe.
          || (step_index>father_last_step) && segments.IsIntersectsCurSegDown(Q[step_index])// Or intersects the last stair.
          || segments.IsIntersectsCurSegUp(Q[cur_father_pos])) // Or intersects the candidate-to-inherit.
        {
          // Keep segment in `L` (cannot become a stair in the new staircase).
          L[new_L_size] = cur_seg;
          // Store segment's staircase location hint (used later for `FindInt`).
          location[new_L_size++] = step_index;
        }
        else
        { // Add it to `Q` as an original step.
          step_index++;
          // Negative (<= 0) means original step.
          father_loc[step_index] = (int4)inherit_each - (int4)cur_father_pos;
          Q[step_index] = cur_seg;
        }
        cur_L_pos++; // Move to next in `L`.
      }
      else
      {  // Add candidate to inherit into `Q`.
        step_index++;
        // Positive means inherited step; value is parent position.
        father_loc[step_index] = cur_father_pos;
        Q[step_index] = Q[cur_father_pos];
        cur_father_pos += inherit_each;
      };
    }
    for (; cur_L_pos<Size; cur_L_pos++)  // Parent staircase is over. Add remaining segments in the same manner as in the suboptimal case.
    {
      cur_seg = L[cur_L_pos];
      segments.SetCurSegCutBE(cur_seg);
      if ((segments.GetSegR(cur_seg) < RBoundIdx) ||// Segment does not cover stripe.
        (step_index>father_last_step) && segments.IsIntersectsCurSegDown(Q[step_index]))// If it intersects the last stair.
      {
        L[new_L_size] = cur_seg;
        location[new_L_size++] = step_index;
      }
      else
      {  // Add it to `Q` as an original step.
        step_index++;
        father_loc[step_index] = (int4)inherit_each - (int4)cur_father_pos;
        Q[step_index] = cur_seg;
      }
    };
    for (; cur_father_pos <= father_last_step; cur_father_pos += inherit_each) // `L` is over. Add the remaining inherited parent segments to `Q`.
    {
      step_index++;
      father_loc[step_index] = cur_father_pos;
      Q[step_index] = Q[cur_father_pos];
    };

    // After staircase construction, check remaining `L` segments against the new staircase.
    long long n_int = 0;
    for (decltype(new_L_size) i = 0; i != new_L_size; i++)
    {
      segments.SetCurSegCutBE(L[i]);
      // `location[i]` contains the approximate location for the segment `L[i]` in the newly built staircase.
      n_int +=FindInt(segments,father_last_step, step_index, location[i]);
    }

    dont_cut_stripe = n_int > new_L_size + cut_margin;
    _step_index = step_index;
    L_size = new_L_size;

  };

  // Recursive driver for the optimal algorithm.
  //
  // Parameters:
  // - `father_first_step`: first step index in parent staircase eligible for inheritance in this subtree.
  // - `ladder_start_index`: location in `Q` where current staircase begins.
  // - `[interval_left_rank, interval_right_rank]`: stripe interval for this recursion.
  // - `stack_pos`: recursion stack chain for bubbling.
  // - `call_numb`: number of consecutive "repeat without cutting" calls (dense stripe mitigation).
  //
  // Preconditions:
  // - Active set `L` matches `interval_left_rank` before entry (maintained by caller + `InsDel` protocol).
  //
  // Postconditions:
  // - Intersections inside the interval are reported.
  // - On return, if a staircase was built, the merged active set for parent is restored.
  template<class SegmentsColl>
  void FindR(SegmentsColl &segments, uint4 father_first_step, uint4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec *stack_pos, uint4 call_numb,uint4 _max_call=max_call)
  {
    segments.SetCurStripe(interval_left_rank, interval_right_rank);
    if (interval_right_rank - interval_left_rank == 1) {
      // Base case: a unit stripe. Use base implementation (fast path) for in-stripe processing.
      if (L_size > 1) CFAST::SearchInStrip(segments, ladder_start_index);
      return;
    }


    ProgramStackRec stack_rec(ladder_start_index);

    // Decide whether to use the optimal split (inherit) or fallback to the fast split.
    bool use_opt =  (ladder_start_index > big_staircase_threshold + father_first_step);
    if (use_opt)
    {// Use optimal variant if father staircase is big.
      Split(segments, interval_right_rank, father_first_step, stack_rec.Q_pos);
      stack_pos = stack_rec.Set(stack_pos, interval_right_rank);
      father_first_step = ladder_start_index + inherit_offset;
    }
    else
    {// Use fast variant if new staircase is small.
      stack_rec.Q_pos+=CFAST::Split(segments,Q+ stack_rec.Q_pos,interval_right_rank);
      if ((ladder_start_index < stack_rec.Q_pos))
      {
        // Mark all newly created steps as "no inheritance metadata" (treat as original).
        for (auto i = ladder_start_index + 1; i <= stack_rec.Q_pos; ++i)
          father_loc[i] = undef_loc;

        stack_pos = stack_rec.Set(stack_pos, interval_right_rank);
        father_first_step = ladder_start_index + inherit_offset;
      }
    }

    // If this stripe appears intersection-heavy, avoid cutting too early: repeat recursion on same interval.
    if (dont_cut_stripe && (call_numb < _max_call))
      FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_rank, interval_right_rank, stack_pos, call_numb + 1,_max_call);
    else // Cut stripe.
    {
      uint4 m = (interval_left_rank + interval_right_rank) >>1;
      if (call_numb > 1)
      {
        _max_call -=  (_max_call > 4) * 2;
        // If we had repeats, cut into two parts.
        FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_rank, m, stack_pos, 0, _max_call);
        InsDel(*this,segments,m, stack_pos);
        FindR(segments,father_first_step, stack_rec.Q_pos, m, interval_right_rank, stack_pos, 0, _max_call);
      }
      else
      {// Otherwise cut into four parts (often faster for some distributions).
        _max_call -=  (_max_call > 8) * 4;

        uint4 q = (interval_left_rank + m) >>1;
        if (interval_left_rank != q) {
          FindR(segments,father_first_step, stack_rec.Q_pos, interval_left_rank, q, stack_pos, 0, _max_call);
          InsDel(*this,segments, q, stack_pos);
        }
        if (q != m) {
          FindR(segments, father_first_step, stack_rec.Q_pos, q, m, stack_pos, 0, _max_call);
          InsDel(*this,segments, m, stack_pos);
        }
        q = (interval_right_rank + m) / 2;
        if (q != m) {
          FindR(segments, father_first_step, stack_rec.Q_pos, m, q, stack_pos, 0, _max_call);
          InsDel(*this,segments, q, stack_pos);
        }
        FindR(segments, father_first_step, stack_rec.Q_pos, q, interval_right_rank, stack_pos, 0, _max_call);
      }

    }

    // No staircase was created => nothing to merge for parent.
    if (ladder_start_index >= stack_rec.Q_pos) return;

    // Restore stripe-left context before merge (some `SegmentsColl` require this as part of protocol).
    segments.SetCurStripeLeft(interval_left_rank);

    // Merge based on the variant used to build the staircase.
    if(use_opt)
      Merge(segments, interval_left_rank, ladder_start_index, stack_rec.Q_pos);
    else
      CFAST::Merge(segments, interval_left_rank, ladder_start_index, stack_rec.Q_pos);
  };

};

