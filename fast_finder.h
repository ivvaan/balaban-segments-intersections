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
// New implementation (recursive "ladder/staircase" stripe processing).

#include <algorithm>
#include <vector>
#include <thread>
#include <cassert>
#include <type_traits>

// `CFastIntFinder` implements the fast intersection-finding algorithm described in the wiki
// ("New-implementation", "Search-in-strip", "Split", "Sentinels-and-guards").
//
// Glossary (as used by this implementation):
// - "stripe": x-interval between two consecutive endpoint ranks (or a larger interval in recursion).
// - `L`: active set of segments intersecting the current stripe, ordered by y at stripe boundary.
// - `Q`: "ladder"/"staircase" representation (a monotone sequence of segments in y-order)
//        used to reduce intersection checks during recursion.
// - `R`: scratch buffer (used as a temporary `Q` in some procedures and as merge destination).
// - "sentinel": a special segment id inserted into arrays to avoid bounds checks in inner loops.
//
// Main expected workflow:
// 1) `prepare_ends(segments)` to compute ends and statistics needed for sizing internal buffers.
// 2) `find_intersections(...)` to run the algorithm and report intersections via `segments`.
class CFastIntFinder
{
  // RAII helper that temporarily replaces an array cell with a collection-provided sentinel value.
  // For collections without sentinel support, this is a no-op.
  //
  // Usage pattern:
  //   CSentinel s(segments, array[idx]);          // bottom sentinel
  //   CSentinel t(segments, array[idx], true);   // top sentinel
  //
  // Precondition:
  // - `ini` refers to a valid writable cell which may be overwritten for the duration of the guard.
  template<class SegmentsColl, bool HasSentinels = has_sentinels<SegmentsColl> >
  struct CSentinel {
    CSentinel(SegmentsColl& coll, uint4& ini, bool is_top = false) {};
  };

  template<class SegmentsColl>
  struct CSentinel <SegmentsColl, true> {
    CSentinel(SegmentsColl& coll, uint4& ini, bool is_top = false) : to_restore(ini) {
      prev_val = to_restore;
      to_restore = coll.get_sentinel(is_top);
    };
    ~CSentinel() {
      to_restore = prev_val;
    };
  private:
    uint4& to_restore;
    uint4 prev_val;
  };

public:

  CFastIntFinder() = default;
  ~CFastIntFinder() { unclone(); FreeMem(); };
  void FreeMem(){};

  // Find all intersections in a stripe range `[from, to]` (endpoint ranks).
  //
  // If called with default arguments (`to == 0`), processes the full range and finishes
  // by calling `segments.IntersectionsFindingDone()`.
  //
  // Parallel usage:
  // - Worker threads call this overload with explicit `from/to` and MUST NOT call
  //   `segments.IntersectionsFindingDone()` (this is handled by the master).
  //
  // Preconditions:
  // - `prepare_ends(segments)` was called on the same finder (or on the master that was cloned).
  // - `segments` is prepared and provides stable endpoint ranks in `[0, nTotX-1]`.
  // - For `from > 0`, the active set `L` is reconstructed by `CalcLAt()` (O(nTotSegm)).
  //
  // Postconditions:
  // - Intersections inside processed stripes are reported through `segments` callbacks.
  // - `L` is left representing the active set at `to` (after `InsDel` at `to`).
  template<FastSegColl SegmentsColl>
  void find_intersections(SegmentsColl  &segments,uint4 from=0,uint4 to=0)
  {
    // Allocate working buffers sized from statistics computed in `prepare_ends()`.
    // Layout notes:
    // - `L` and `R` are allocated with one extra cell and then incremented by 1 so that
    //   `L[-1]`/`R[-1]` exists for placing a bottom sentinel/guard when needed.
    len_of_Q = LR_len;
    DECL_RAII_ARR(L, LR_len+1);
    ++L;//to have one cell before L for sentinel
    DECL_RAII_ARR(R, LR_len+1);
    ++R;//to have one cell before R for sentinel
    DECL_RAII_ARR(Q, len_of_Q);

    bool not_parallel = false;
    if (to == 0) {
      from = 0;
      to = nTotX - 1;
      not_parallel = true;
    }
    constexpr uint4 bottom_index = 0;
    ProgramStackRec stack_rec(bottom_index, nTotX); // Must be initialized this way (top frame).
    if (from == 0) {
      segments.SetCurStripeRight(0);
      InsDel(*this,segments,0,&stack_rec);
    }
    else
      L_size = CalcLAt(segments, from);
    if (not_parallel) {
      // Heuristic: for sparse vertical density, a simple scan-in-stripe is faster than recursion.
      if (avr_segm_on_vline < 20)
        SISFindR(segments, bottom_index, from, to, &stack_rec);
      else
        MultipleCutting(segments, bottom_index, from, to, &stack_rec, GetDivPow(to - from));
    }
    else // Parallel call: recursively process only `[from,to]` without global completion notification.
      FindR(segments, bottom_index, from, to, &stack_rec);
    InsDel(*this,segments, to, &stack_rec);
    if(not_parallel)
      segments.IntersectionsFindingDone();
  }

  // Parallel dispatcher. Splits the full x-range into `n_threads` stripes and runs
  // independent workers. Each worker uses a cloned `SegmentsColl` instance (per-thread state).
  //
  // Preconditions:
  // - `prepare_ends(segments)` was called.
  // - `segments` supports cloning protocol: `InitClone(n_threads)` and cloning ctor.
  //
  // Postconditions:
  // - All stripes are processed; then `segments.IntersectionsFindingDone()` is called once.
  template<ParallelSegColl SegmentsColl>
  void find_intersections(uint4 n_threads,SegmentsColl& segments)
  {
    using namespace std;
    vector<thread> wrk_threads;
    auto thread_func = [](CFastIntFinder* master, SegmentsColl* segments,uint4 thread_index, uint4 from, uint4 to) {
      SegmentsColl coll(*segments, thread_index);
      CFastIntFinder(master).find_intersections(coll, from, to);
      };
    segments.InitClone(n_threads);
    auto n = segments.GetSegmNumb();
    double part = nTotX / (double)n_threads;
    uint4 start_from = part;
    for (uint4 i = 2, from = start_from, to; i <= n_threads; ++i) {
      to = (i == n_threads) ? nTotX-1 : (uint4)(part * i);
#if defined(DEBUG) || defined(_DEBUG)
      thread_func(this, &segments,i-1, from, to); // Start intersection finding in stripe `[from,to]`.
#else
      wrk_threads.emplace_back(thread_func, this, &segments,i-1, from, to); // Start intersection finding in stripe `[from,to]`.
#endif
      from = to;
    }

    find_intersections(segments, 0, start_from);
#if !(defined(DEBUG) || defined(_DEBUG))
    for (auto& cur_thread : wrk_threads)
      cur_thread.join(); // Wait until all threads finish.
#endif
    segments.IntersectionsFindingDone();
  }


  // Prepares endpoint ranks and computes statistics used for sizing internal arrays.
  //
  // Preconditions:
  // - `segments` contains all segments for which intersections should be found.
  //
  // Postconditions:
  // - `nTotSegm`, `nTotX`, `LR_len`, `avr_segm_on_vline` are initialized.
  // - Must be called before `find_intersections(... )`.
  template <class SegmentsColl>
  void prepare_ends(SegmentsColl& segments)
  {
    // Moved to `SegmentsColl::Prepare()` - it computes ENDS, SegL, SegR and statistics.
    auto pr = segments.Prepare();
    // Set internal common statistics derived from collection `Prepare()` result.
    nTotSegm = segments.GetSegmNumb();
    nTotX = pr.ends_len;
    LR_len = pr.max_segm_on_vline + 3; // +3 for sentinels/guards.
    avr_segm_on_vline = pr.avr_segm_on_vline;
  };


protected:
  // "Clone constructor": used in parallel workers to copy only read-only statistics from master.
  // Working buffers are still allocated per-thread in `find_intersections()`.
  CFastIntFinder(CFastIntFinder* master) {
    clone(master);
  };

  // Represents a chain of recursion frames ("program stack") for the ladder/staircase hierarchy.
  // Each frame stores:
  // - `Q_pos`: the index in `Q` where the current staircase for this recursion level begins.
  // - `right_bound`: the stripe right boundary rank for this level (used for bubbling logic).
  // - `prev`: pointer to parent level frame (nullptr for the top frame).
  //
  // Used by `InsDel()` + `FindIntI()` to "bubble up" and check newly inserted segments
  // against staircases of all higher recursion levels.
  struct ProgramStackRec // Structure to use the program stack to store a list of staircases above the current call.
  {
    uint4 Q_pos; // Starting position in `Q`.
    uint4 right_bound;// Right bound of current staircase.
    ProgramStackRec* prev = nullptr; // Link to parent staircase record.
    ProgramStackRec(uint4 qp) :Q_pos(qp) {};
    ProgramStackRec(uint4 qp, uint4 rb) : Q_pos(qp), right_bound(rb) {};
    ProgramStackRec(uint4 qp, uint4 rb, ProgramStackRec* pr) : Q_pos(qp), right_bound(rb), prev(pr) {};
    ProgramStackRec* Set(ProgramStackRec* p, uint4 rb) { prev = p, right_bound = rb; return this; };
    bool isnot_top() { return prev != nullptr; };
  };

  // Cutting heuristics:
  // - `cut_margin`: extra slack to decide that a stripe is "too intersection-heavy" and should not be cut.
  // - `min_strip_width`: base recursion stop threshold; narrow stripes are processed by `SISFindR()`.
  constexpr static const uint4 cut_margin = 24;
  constexpr static const uint4 min_strip_width = 6;// at least 2

  // Working buffers (allocated per call to `find_intersections()`).
  uint4* Q = nullptr; // Ladder/staircase storage.
  uint4* L = nullptr; // Active set (sorted by y at stripe boundary).
  uint4* R = nullptr; // Scratch / merge destination.

  CFastIntFinder* clone_of = nullptr;

  // Statistics / configuration derived from `prepare_ends()`.
  double avr_segm_on_vline = 0;
  uint4 LR_len = 0;
  uint4 nTotSegm = 0;
  uint4 len_of_Q = 0;
  uint4 L_size = 0;
  uint4 nTotX = 0;

  // Set by `Split*()` when a stripe appears to be "dense" (too many intersections with ladder stairs),
  // in which case `FindRNoChecks()` repeats splitting without further cutting.
  bool dont_cut_stripe;



  // Returns pointer to the tail buffer used by `Split*()` for storing per-segment ladder indices.
  // This reuses `Q` memory by writing indices from the end backwards.
  uint4* GetQTail() {
    return Q + len_of_Q;
  }

  // Computes a division power for `MultipleCutting()` based on stripe width and average density.
  // Larger values => more sub-stripes (2^divide_pow).
  uint4 GetDivPow(uint4 l) {
    if (l < 64) return 2;
    uint4 res = 0;
    for (uint4 target = 2.0 * l / MAX(avr_segm_on_vline, 8.0); target; target >>= 1)
      ++res;
    return MAX(res, 4);
  };

  // Apply insert/delete events at endpoint rank `end_rank` and update active set `L` via collection.
  // Additionally, if we are inside recursion (i.e., not top frame), "bubble up" newly inserted
  // segments and check them against staircases in higher recursion levels (`FindIntI`).
  //
  // Preconditions:
  // - `segments.SetCurStripe*(...)` is set consistently for `end_rank` in the caller.
  // - `self.L` contains the current active set just before applying events at `end_rank`.
  //
  // Postconditions:
  // - `segments.InsDel(...)` updates `self.L` / `self.L_size` according to endpoint event lists.
  template<class IntersectionFinder,class SegmentsColl>
  static void InsDel(IntersectionFinder &self, SegmentsColl& segments, uint4 end_rank, ProgramStackRec* stack_pos) {
    segments.InsDel(self.L_size, self.L, end_rank);
    if (stack_pos->isnot_top()) {
      // Find intersections for all inserted (first endpoints) and vertical segments
      // in staircase levels above (in the stack).
      for (auto p = segments.GetEndsListFirst(end_rank); p; p = segments.GetEndsListNext(p)) {
        auto sn = SegmentsColl::get_segm(*p);
        segments.SetCurSeg4Bubble(sn);
        self.FindIntI(segments, sn, stack_pos); // Find "internal" intersections against all parent ladders.
      }
    }

  }



  // Finds intersections of the current segment (as set in `segments`) with the ladder slice `[qb, qe)`,
  // using `l` as the "current" ladder step (result of a binary search by y-position).
  //
  // Sentinel usage:
  // - `*qb` and `*qe` are temporarily overwritten with sentinel values (if supported) to allow
  //   range-based intersection checks without explicit bounds checks.
  //
  // Line-segment optimization:
  // - For strict line segments, after the first intersection below is found, there cannot be
  //   further intersections above in this monotone sequence, so the upward scan is skipped.
  template <class SegmentsColl>
  static void FindInt(SegmentsColl& segments, uint4* const qb, uint4* const qe, uint4* l)
  {
    // Use range overloads with sentinel instead of single-element overloads.
    // Bottom and top sentinel are placed using CSentinel RAII.
    CSentinel bsentinel(segments, *qb);
    auto c = segments.FindCurSegIntDownWith(l, qb); // First, get intersections below.

    constexpr const bool line_seg = (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) == _Coll_flag_state::state_true);

    if (line_seg && (c != l))
      return; // If found and segment is a line, there cannot be any more intersections above.

    CSentinel tsentinel(segments, *qe, true);
    segments.FindCurSegIntUpWith(l + 1, qe); // Get intersections above.
  };

  // Simplified sequential "search in strip" for non-line segments:
  // - Repeatedly partitions `L` into segments that do not intersect the current ladder `Q`
  //   (moved into `Q`) and segments that intersect it (kept in `L` for next iteration).
  // - Uses `Q_tail` indices stored at the end of `Q` to drive upward checks efficiently.
  //
  // Preconditions:
  // - `_L[0..L_size)` is the active set for the stripe, sorted at the left boundary.
  // - `_Q` points to a valid location where the ladder can be built (`_Q[0]` is used as sentinel cell).
  //
  // Postconditions:
  // - All intersections inside the stripe are reported via `segments.FindCurSegInt*With` calls.
  // - On exit, `_Q` contains the moved segments (ladder), and `_L` contains remaining segments (if any).
  template<class SegmentsColl>
  void SearchInStripNonLineSeg(SegmentsColl& segments, uint4* _L, uint4* _Q)//simplified version for SearchInStrip
  {
    auto last_L = _L + L_size;
    auto first_L = _L + 1;
    auto* _Q_pos = _Q + 1;
    auto T_pos = GetQTail();
    do
    {
      auto new_L_pos = _L;
      *_Q_pos = _L[0];
      if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::needs_SetCurSegCutBE_at_start) == _Coll_flag_state::state_true)
        segments.SetCurSegCutBE(_L[0]);
      auto Q_tail = T_pos;
      auto cur_L = first_L;
      for (CSentinel sentinel(segments, *_Q); cur_L < last_L; ++cur_L)
      {
        segments.SetCurSegCutBE(*cur_L);

        if (_Q_pos == segments.FindCurSegIntDownWith(_Q_pos, _Q))
          *++_Q_pos = *cur_L;
        else
        {
          *new_L_pos++ = *cur_L;
          *--Q_tail = _Q_pos - _Q;
        }
      }
      if (new_L_pos == _L)// L is empty
        return;
      *--Q_tail = _Q_pos - _Q;//place a guard for the loop below
      ++_Q; ++_Q_pos;
      if constexpr (has_sentinels<SegmentsColl>)
        *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel

      Q_tail = T_pos;
      auto cur_Q = _Q + *--Q_tail;
      for (cur_L = _L; cur_Q != _Q_pos; ++cur_L)
      {
        segments.SetCurSegCutBE(*cur_L);
        segments.FindCurSegIntUpWith(cur_Q, _Q_pos);
        cur_Q = _Q + *--Q_tail;
      }
      _Q = _Q_pos - 1;
      last_L = new_L_pos;
    } while (last_L > first_L);//L contains >1 segment
    *++_Q = _L[0];
  };

  // "Search in strip" for line segments:
  // Performs insertion-sort-like processing from top to bottom.
  // For line segments, "intersects inside stripe" acts as a comparison that swaps order at the
  // right boundary. Sorting by this comparison both detects intersections and produces the order.
  //
  // Preconditions:
  // - `L_` contains the active set, sorted by y at the left stripe boundary.
  //
  // Postconditions:
  // - Intersections inside the stripe are reported.
  // - `L_` becomes sorted as required for the next step.
  template<class SegmentsColl>
  void SearchInStripLineSeg(SegmentsColl& segments, uint4* L_) {
    auto last = L_ + L_size;
    if constexpr (has_sentinels<SegmentsColl>)
      *last = segments.get_sentinel(true);
    decltype(last) cur = last - 1, next;
    if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::needs_SetCurSegCutBE_at_start) == _Coll_flag_state::state_true)
      segments.SetCurSegCutBE(*cur);//we don't need SetCurSegCutBE in all cases, only for Y caching collections
    while (cur != L_) {// insertion sort from last to first
      next = cur--;
      segments.SetCurSegCutBE(*cur);
      auto no_int = segments.FindCurSegIntUpWith(next, last);
      std::rotate(cur, next, no_int);
    }
  }

  // "Internal" intersection checks for a newly inserted segment `sn` against all parent ladders.
  //
  // The algorithm walks stack frames towards the top and, at each frame, performs:
  // - binary search in the level ladder `[Qb, Qe)` to locate the step around current segment y,
  // - then `FindInt` to scan down/up with sentinels.
  //
  // Preconditions:
  // - `segments.SetCurSeg4Bubble(sn)` was called (the collection knows which segment is current).
  // - `stack_pos` points to the current recursion frame in which insertion happened.
  //
  // Postconditions:
  // - Intersections between `sn` and ladder segments of all higher recursion levels are reported.
  template <class SegmentsColl>
  void FindIntI(SegmentsColl& segments, uint4 sn, ProgramStackRec* stack_pos) const
  {
    auto r_index = segments.GetSegR(sn);
    while (stack_pos->right_bound <= r_index)
      stack_pos = stack_pos->prev; // Go from bottom to top and find staircase to start.
    int4  m, len;
    uint4* Qb, * Ql, * Qe = Q + (stack_pos->Q_pos + 1);
    for (stack_pos = stack_pos->prev; stack_pos != nullptr; stack_pos = stack_pos->prev) {
      Qb = Ql = Q + stack_pos->Q_pos;
      len = Qe - Qb;
      while (len > 1) // binary search
      {
        m = len / 2; //
        if (segments.UnderCurPoint(Ql[m]))
          Ql += m;
        len -= m;
      }
      FindInt(segments, Qb, Qe, Ql);
      Qe = Qb + 1; // Move staircase bound to parent staircase.
    };
  };


  // Stripe processing for a narrow stripe (single x-step) or as a terminal stage.
  // Chooses the best specialized method depending on segment type.
  //
  // Preconditions:
  // - `segments.SetCurStripe(left,right)` is set by the caller (or will be inside `SISFindR` loop).
  // - `L[0..L_size)` is sorted at the stripe left boundary.
  //
  // Postconditions:
  // - All intersections inside the stripe are reported.
  // - `L` is sorted for the stripe right boundary.
  template<class SegmentsColl>
  void SearchInStrip(SegmentsColl &segments, uint4 qp)
  {
    if constexpr(SegmentsColl::get_coll_flag(_Coll_flags::line_segments)==_Coll_flag_state::state_true)
    {
      // For line segments we can do more efficient insertion sorting using intersection check as comparison.
      // If `s1 < s2` at the left bound and `s1` intersects `s2` inside the stripe, then `s1 > s2`
      // at the right bound of the stripe.
      // All segments are sorted at the left bound by precondition, so intersection is equivalent to
      // the comparison at the right bound. A byproduct of the sorting is intersection detection.
      SearchInStripLineSeg(segments, L);
    }
    else {
      // Use `R` as a temporary `Q`, performing a sequential simplified split.
      SearchInStripNonLineSeg(segments, L, R - 1);
      // At the end, all intersections are found, and all segments are moved from `L` to `Q`.
      // Move them back.
      std::swap(L, R);
      // Now segments are in `L` again, but in incorrect order; sort them.
      std::sort(L, L + L_size, [&](uint4 s1, uint4 s2) {return segments.RBelow(s1, s2); });
    }
  };

  // Simple iteration over endpoint ranks for narrow stripes.
  // For each unit stripe `[i-1,i]`: update stripe, detect intersections in stripe, then handle events at `i`.
  //
  // Postconditions:
  // - On exit, `segments` is set to stripe `[interval_right_rank-1, interval_right_rank]`.
  // - If `L_size > 1`, `SearchInStrip` has been applied for the final stripe.
  template <class SegmentsColl>
  void SISFindR(SegmentsColl& segments, uint4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    for (auto i = interval_left_rank + 1; i != interval_right_rank; ++i) {
      if (L_size > 1) {
        segments.SetCurStripe(i - 1, i);
        SearchInStrip(segments, ladder_start_index);
      }
      else
        segments.SetCurStripeRight(i);
      InsDel(*this,segments, i, stack_pos);
    }
    // The stripe right bound needs to be installed correctly even if `L_size <= 1`.
    segments.SetCurStripe(interval_right_rank - 1, interval_right_rank);
    if (L_size > 1) SearchInStrip(segments, ladder_start_index);
  }

  // Core recursive step without stripe-width checks (caller ensures stripe is "proper" for recursion).
  //
  // Steps:
  // 1) `Split` current active set `L` against ladder `Q` to build a new staircase (`Q_pos`).
  // 2) Decide whether to cut further (middle split) or re-run without further cutting (`dont_cut_stripe`).
  // 3) `Merge` ladder and remaining segments to restore consistent `L` ordering for parent.
  //
  // If `Split` produces empty ladder or `L` is empty, falls back to `MultipleCutting`.
  template <class SegmentsColl>
  void FindRNoChecks(SegmentsColl& segments, uint4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    uint4 Q_pos;
    // If `L` or `Q` is empty, cut the stripe into `2^divide_pow` sub-stripes.
    constexpr const uint4 divide_pow = 4;
    if ((L_size == 0)||(0 == (Q_pos = Split(segments, Q + ladder_start_index, interval_right_rank))))
      return MultipleCutting(segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos, divide_pow);
    // Normal processing otherwise.
    Q_pos += ladder_start_index;
    ProgramStackRec stack_rec(Q_pos, interval_right_rank, stack_pos);// Create new stack entry for bubbling up.
    if (dont_cut_stripe ) { // If found a lot of intersections, repeat FindR on the same interval.
      FindRNoChecks(segments, Q_pos, interval_left_rank, interval_right_rank, &stack_rec);
    }else{// Cut at the middle into two sub-stripes.
      uint4 middle = (interval_left_rank + interval_right_rank) / 2;
      FindR(segments, Q_pos, interval_left_rank, middle, &stack_rec);
      InsDel(*this,segments, middle, &stack_rec);
      FindR(segments, Q_pos, middle, interval_right_rank, &stack_rec);
      // Works without SetCurStripeLeft, but it simplifies segment collection class protocol.
      //segments.SetCurStripeLeft(interval_left_rank);
    }
    Merge(segments, interval_left_rank, ladder_start_index, Q_pos);
  };

  // Recursive driver with termination on narrow stripes.
  //
  // Preconditions:
  // - Caller has set a consistent active set `L` for `interval_left_rank`.
  template <class SegmentsColl>
  void FindR(SegmentsColl& segments, uint4 ladder_start_index, uint4 interval_left_rank, uint4 interval_right_rank, ProgramStackRec* stack_pos)
  {
    if (interval_right_rank < interval_left_rank + min_strip_width) {// End of recursion for narrow stripes.
      SISFindR(segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos);
    } else {// Normal step.
      segments.SetCurStripe(interval_left_rank, interval_right_rank);
      FindRNoChecks(segments, ladder_start_index, interval_left_rank, interval_right_rank, stack_pos);
    }
  }

  // Cuts the stripe into `2^divide_pow` sub-stripes and processes each sub-stripe sequentially.
  //
  // Purpose:
  // - Fallback strategy when `Split` cannot build a useful ladder (or active set is empty),
  //   or as an initial optimization for very wide stripes.
  template <class SegmentsColl>
  void MultipleCutting(SegmentsColl& segments, uint4 ladder_start_index, uint4 left_bound, uint4 interval_right_rank, ProgramStackRec* stack_pos, uint4 divide_pow)
  {
    if (interval_right_rank < left_bound + (min_strip_width << divide_pow)) // End of recursion for narrow stripes.
      return SISFindR(segments, ladder_start_index, left_bound, interval_right_rank, stack_pos);  // If stripe is narrow, just apply SISFindR.

    uint8 step = interval_right_rank - left_bound;
    uint8 rb = (uint8)left_bound << divide_pow;
    uint4 right_bound = (rb += step) >> divide_pow;
    segments.SetCurStripe(left_bound, right_bound);
    FindRNoChecks(segments, ladder_start_index, left_bound, right_bound, stack_pos);
    do {
      left_bound = right_bound;
      right_bound = (rb += step) >> divide_pow;
      InsDel(*this,segments, left_bound, stack_pos);
      segments.SetCurStripe(left_bound, right_bound);
      FindRNoChecks(segments, ladder_start_index, left_bound, right_bound, stack_pos);
    } while (right_bound != interval_right_rank);
  }

  // Merges the current ladder `[QB+1 .. QE]` with the remaining active set (stored in `L`/`R` buffers),
  // while also detecting some intersections needed for correctness when moving up recursion.
  //
  // High-level idea:
  // - Merge by y-order at the stripe right boundary, producing a new `L` (in `R` buffer).
  // - For segments that start inside current interval (GetSegL > LBoundIdx),
  //   emit required intersection checks against nearby ladder steps.
  //
  // Preconditions:
  // - `Q` contains a staircase built for the interval, starting at `QB` and ending at `QE`.
  // - Current `L` (aka `_R` here) contains segments not in the staircase, sorted by y at right boundary.
  //
  // Postconditions:
  // - `L` points to the merged buffer and `L_size` is updated.
  template<class SegmentsColl>
  void Merge(SegmentsColl& segments, uint4 LBoundIdx, uint4 QB, uint4 QE)
  {
    const auto Size = L_size;                // Number of non-ladder segments currently in the active set.
    uint4 cur_R_pos = 0, new_size = 0;       // `cur_R_pos` iterates old active set; `new_size` builds merged output size.
    auto cur_stair = Q + QB + 1;             // Iterator over ladder steps; `Q[QB]` is a bottom-sentinel/guard slot here.
    uint4* const __restrict _R = L;         // Old active set (non-ladder part) (read side of merge).
    uint4* const __restrict _L = R;         // Destination buffer for merged active set (write side of merge).
    auto cur_seg = _R[cur_R_pos];            // Current segment from `_R` being merged.
    auto const bot_Q = Q + QB;               // Pointer to bottom ladder cell used for bottom sentinel placement.
    auto const top_Q = Q + QE + 1;           // One-past-last ladder cell; may be overwritten with a top sentinel.
    CSentinel sentinel(segments, *bot_Q);    // Place bottom sentinel (restored on scope exit if supported).
    if constexpr (has_sentinels<SegmentsColl>)
      *top_Q = segments.get_sentinel(true);  // Place top sentinel to allow range scans up to `top_Q` safely.
    while ((cur_stair != top_Q) && (cur_R_pos < Size))
    {
      // Standard merge step: choose the next element by y-order at the stripe right boundary.
      // Invariant: `_R[0..Size)` and `Q[QB+1..QE]` are each sorted by y at the right boundary.
      if (segments.RBelow(cur_seg, *cur_stair))
      {
        // `cur_seg` should appear below the current stair step in the merged order.
        if (segments.GetSegL(cur_seg) > LBoundIdx)
        {
          // This segment started strictly inside the current recursion interval.
          // For correctness, it must be checked against neighboring ladder steps when joining levels.
          segments.SetCurSegCutEnd(cur_seg); // "Current segment" is set with cut-end semantics for this stripe.

          {
            // First, scan upward: find intersections with ladder steps above `cur_stair`.
            // Returns a pointer to the first ladder step that does NOT intersect (for some collections),
            // which enables a line-segment optimization below.
            auto c = segments.FindCurSegIntUpWith(cur_stair, top_Q);

            // For non-line segments, an intersection above does not imply anything about intersections below,
            // so we always also scan down.
            // For strict line segments, if an intersection above was found (`c != cur_stair`),
            // then (by monotonicity) scanning down from below `cur_stair` would be redundant.
            constexpr const bool non_line_seg =
              (SegmentsColl::get_coll_flag(_Coll_flags::line_segments) != _Coll_flag_state::state_true);
            if (non_line_seg || (c == cur_stair))
              segments.FindCurSegIntDownWith(cur_stair - 1, bot_Q); // Scan downward (range ends at bottom sentinel).
          };
        }

        _L[new_size++] = cur_seg;            // Emit `cur_seg` into merged output.
        cur_seg = _R[++cur_R_pos];           // Advance in `_R`.
      }
      else
      {
        // Next stair step should appear below `cur_seg` in the merged order.
        _L[new_size++] = *cur_stair++;       // Emit ladder step and advance ladder iterator.
      }
    }

    // If ladder is exhausted, append remaining `_R` segments.
    // (The opposite case—`_R` exhausted—is handled by the final `while (cur_stair != top_Q)` below.)
    for (auto last_Q = top_Q - 1; cur_R_pos < Size;)
    {
      if (segments.GetSegL(cur_seg) > LBoundIdx)
      {
        // For segments that start inside the interval, we still need at least a downward scan
        // against the (already fully included) ladder. `last_Q` anchors the scan from the top.
        segments.SetCurSegCutEnd(cur_seg);
        segments.FindCurSegIntDownWith(last_Q, bot_Q);
      }
      _L[new_size++] = cur_seg;              // Emit the segment.
      cur_seg = _R[++cur_R_pos];             // Advance in `_R`.
    }

    // If `_R` is exhausted, append remaining ladder steps.
    while (cur_stair != top_Q)
      _L[new_size++] = *cur_stair++;

    // Swap buffers: merged output becomes the new active set.
    L = _L; R = _R;
    L_size = new_size;                       // New size of active set at the stripe right boundary.
  };

  // Split step specialized for line segments.
  //
  // Builds a new ladder into `_Q[1..Q_size]` and keeps intersecting/non-covering segments in `L`.
  // Also records, in the tail of `Q`, ladder indices needed for the upward checks phase.
  //
  // Returns:
  // - `Q_size` (number of ladder steps), or 0 if no ladder could be built.
  //
  // Side effects:
  // - Updates `L_size`, writes into `R`, and sets `dont_cut_stripe` heuristic.
  template<class SegmentsColl>
  uint4 Split4LineSeg(SegmentsColl& segments, uint4* _Q, uint4 RBoundIdx)
  {
    auto R_pos = R;
    auto new_L_pos = L;
    auto Q_tail = GetQTail();
    assert(Q_tail > _Q + L_size);
    auto const last_L = L + L_size;
 //at first, place all lowest segments not covering current stripe to R (and implicitly L)
    while ((new_L_pos < last_L) && (segments.GetSegR(*new_L_pos) < RBoundIdx)) {
            *(R_pos++) = *(new_L_pos++);//place segment in  R
            //storing segment position in Q_tail
  //we have the empty ladder so current position is 0 and position above is 1
            *(--Q_tail) = 0;  //it should be position above i.e. 1. We save
            // one addition per loop by incrementing _Q later.
    }
    if (new_L_pos == last_L){
      return 0;
    }

    long long n_int = 0;
    auto  _Q_pos = _Q;
    //first segment covering current stripe we place to Q
    //it can't intersect any of the steps (stairs) because there are no stairs yet
    *++_Q_pos = *new_L_pos;
    if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::needs_SetCurSegCutBE_at_start) == _Coll_flag_state::state_true)
      segments.SetCurSegCutBE(*new_L_pos);//we don't need SetCurSegCutBE in all cases, only for Y caching collections
    auto cur = new_L_pos + 1;
    for (CSentinel sentinel(segments, *_Q); cur < last_L; ++cur) {
      auto cur_seg = *cur;
      segments.SetCurSegCutBE(cur_seg);
      auto k= _Q_pos - segments.FindCurSegIntDownWith(_Q_pos, _Q);
      if (k) {//segment  intersects upper ladder stair
        // add number of intersections with stairs to counter
        n_int += k;//increment found int number
        *(new_L_pos++) = cur_seg;//place segment in L
        continue;
      }

      //segment doesn't intersect the ladder stairs
      if (segments.GetSegR(cur_seg) >= RBoundIdx) {//segment is covering current stripe
        *++_Q_pos = cur_seg; //place segment in Q
        continue;
      }
      //segment is not covering current stripe
      *(R_pos++) = *(new_L_pos++) = cur_seg; //place segment in L and R
      //storing segment position in Q_tail
      *(--Q_tail) = _Q_pos - _Q;  //it should be _Q_pos+1. We save
      // one addition per loop by incrementing _Q later.
    }
    L_size = new_L_pos-L;
    uint4 Q_size=_Q_pos - _Q;
    *--Q_tail = Q_size;//place a guard for the loop below

    // important to start from stair above current segment; meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we increment _Q and `_Q[*Q_tail]` becomes stair above
    ++_Q_pos;

    if constexpr (has_sentinels<SegmentsColl>)//placing sentinel if collection supports
      *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel

    Q_tail = GetQTail() - 1;
    cur = R;
    for (auto Q_idx = *Q_tail; Q_idx != Q_size; Q_idx = *--Q_tail, ++cur)
    {
      segments.SetCurSegCutBE(*cur);//SetCurSegCutBeg
      segments.FindCurSegIntUpWith(_Q + Q_idx, _Q_pos);
    }
    dont_cut_stripe = n_int > L_size + cut_margin;
    return Q_size;
  };

  // Split step specialized for non-line segments.
  //
  // Logic parallels `Split4LineSeg` but uses a slightly different condition set.
  // Returns `Q_size` or 0 if ladder cannot be formed.
  template<class SegmentsColl>
  uint4 Split4NonLineSeg(SegmentsColl& segments, uint4* _Q, uint4 RBoundIdx)
  {
    auto  _Q_pos = _Q;
    long long n_int = 0;
    auto Q_tail = GetQTail();
    auto new_L_pos = L;
    CSentinel sentinel(segments, *_Q);
    for (auto L_pos = L, last_L = L + L_size; L_pos < last_L; ++L_pos) {
      auto cur_seg = *L_pos;
      segments.SetCurSegCutBE(cur_seg);
      auto cur_Q = segments.FindCurSegIntDownWith(_Q_pos, _Q);
      if ((cur_Q == _Q_pos) && (segments.GetSegR(cur_seg) >= RBoundIdx)) //segment doesn't intersect ladder stairs
        //and is covering current stripe and doesn't intersect ladder stairs
          *++_Q_pos = cur_seg;
      else {   //segment  intersects ladder stairs or is not covering
        n_int += _Q_pos - cur_Q; //inc int counter
        *(new_L_pos++)= cur_seg; //place segment in L
        //storing segment position in Q_tail
        *(--Q_tail) = _Q_pos -_Q;  //it should be _Q_pos+1. We save
        // one addition per loop by incrementing _Q later.
      }
    }
    if (_Q_pos == _Q) {
      return 0;
    }
    L_size = new_L_pos - L;
    uint4 Q_size = _Q_pos - _Q;
    *--Q_tail = Q_size;//place a guard for the loop below

    // important to start from stair above current segment; meanwhile _Q[*Q_tail] is stair below
    ++_Q;// so we increment _Q and `_Q[*Q_tail]` becomes stair above
    ++_Q_pos;
    if constexpr (has_sentinels<SegmentsColl>)
      *_Q_pos = segments.get_sentinel(true);//we don't need to restore _Q_pos just place sentinel

    Q_tail = GetQTail() - 1;
    auto L_pos = L;
    for (auto Q_idx =*Q_tail; Q_idx != Q_size; Q_idx = *--Q_tail, ++L_pos)
    {
      segments.SetCurSegCutBE(*L_pos);
      segments.FindCurSegIntUpWith(_Q + Q_idx, _Q_pos);
    }
    dont_cut_stripe = n_int > L_size + cut_margin;
    return Q_size;
  };


  // Split dispatcher: chooses a specialization based on collection flags.
  //
  // Inputs:
  // - `_Q` is the base pointer where the new ladder will be built.
  // - `RBoundIdx` is the stripe right boundary rank.
  //
  // Returns:
  // - ladder size (>0) on success; 0 triggers fallback cutting strategy.
  template<class SegmentsColl>
  uint4 Split(SegmentsColl& segments, uint4* _Q, uint4 RBoundIdx)
  {
    if constexpr (SegmentsColl::get_coll_flag(_Coll_flags::line_segments)==_Coll_flag_state::state_true)
      return Split4LineSeg(segments, _Q, RBoundIdx);
    return Split4NonLineSeg(segments, _Q, RBoundIdx);
  };



  // Copies statistics from master into this instance (used by worker threads).
  // Note: buffers are not shared; each thread allocates its own `L/R/Q`.
  void clone(CFastIntFinder* master)
  {
      nTotSegm = master->nTotSegm;
      LR_len = master->LR_len;
      len_of_Q = LR_len;
      avr_segm_on_vline = master->avr_segm_on_vline;
      nTotX = master->nTotX;
      clone_of = master;
  };

  void unclone()
  {
      if (clone_of != nullptr) {
          clone_of = nullptr;
      }
  };





  // Reconstructs `L` at a given endpoint rank by scanning all segments.
  //
  // Used when processing starts in the middle of the x-range (parallel workers),
  // because the active set is not available from previous `InsDel` steps.
  //
  // Preconditions:
  // - `L` has enough capacity for `LR_len` elements.
  //
  // Postconditions:
  // - `L[0..Size)` contains all segments covering `end_rank`, sorted by y at that x.
  template<class SegmentsColl>
  uint4 CalcLAt(SegmentsColl& segments, uint4 end_rank)
  {
      uint4 Size = 0;
      for (uint4 cur_seg = 0; cur_seg < nTotSegm; ++cur_seg)
          if ((segments.GetSegL(cur_seg) <= end_rank) && (segments.GetSegR(cur_seg) > end_rank))
              L[Size++] = cur_seg;
      segments.SortAt(end_rank, Size, L);
      return Size;
  };

};

