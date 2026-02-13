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

#include "utils.h"
#include "registrator.h"
#include <algorithm>
#include <cassert>
#include <concepts>


// CRTP base class for segment collections.
//
// Extracts boilerplate shared by CLine1SegmentCollection, CLine2SegmentCollection,
// CArcSegmentCollection and CGraphSegmentCollection.
//
// `Derived` must provide:
//   static bool   is_last(uint4 pt);
//   static uint4  get_segm(uint4 pt);
//   uint4  GetSegmNumb() const;
//   bool   UnderCurPoint(int4 s) const;
//   void   SetCurSegAndPoint(uint4 s);
//   void   SetCurPointAtBeg(uint4 s);
//   void   SetCurStripeLeft(uint4 left_rank);
//   bool   LBelow(int4 s1, int4 s2) const;
//   uint4  PrepareEndpointsSortedList(uint4* epoints);
//   void   InitDerived(uint4 n, void* c, IntersectionRegistrator* r);   -- set collection pointer, N, etc.
//   void   ResetDerived();                                               -- free derived-only resources
//   void   CloneDerived(Derived& src);                                   -- copy derived shared pointers
//   void   UncloneDerived();                                             -- null derived shared pointers

template<class Derived, class IntersectionRegistrator>
class CollectionBase {
public:

  // --- accessors for common arrays ---

  uint4 GetSegR(uint4 sn) const { return seg_R_rank[sn]; }
  uint4 GetSegL(uint4 sn) const { return seg_L_rank[sn]; }
  uint4 PointAtRank(uint4 rank) const { return ENDS[rank]; }

  // --- Prepare: build ENDS, seg_L_rank, seg_R_rank and compute statistics ---

  PrepareResult Prepare()
  {
    Reset();
    uint4 Nn = self().GetSegmNumb();
    if (Nn == 0) return {};
    uint4 NN = Nn << 1;
    ENDS = new uint4[NN];
    self().PrepareEndpointsSortedList(ENDS);
    seg_L_rank = new uint4[Nn];
    seg_R_rank = new uint4[Nn];
    uint4 max_segm_on_vline = 0, nsegm_on_vline = 0;
    double avr = 0;
    for (uint4 i = 0; i < NN; ++i) {
      if (Derived::is_last(ENDS[i])) {
        seg_R_rank[Derived::get_segm(ENDS[i])] = i;
        --nsegm_on_vline;
      }
      else {
        seg_L_rank[Derived::get_segm(ENDS[i])] = i;
        ++nsegm_on_vline;
        if (nsegm_on_vline > max_segm_on_vline) max_segm_on_vline = nsegm_on_vline;
      }
      avr += nsegm_on_vline;
    }
    avr /= (double)NN;
    return { NN, max_segm_on_vline, avr };
  }

  // --- Standard InsDel for non-degenerate (real-coordinate) collections ---

  template<class IntersectionFinder, class ProgramStackRec>
  void InsDel(IntersectionFinder& i_f, uint4& L_size, int4* L, uint4 end_rank, ProgramStackRec* stack_pos) {
    auto pt = PointAtRank(end_rank);
    auto sn = Derived::get_segm(pt);
    if (Derived::is_last(pt)) // if endpoint - remove
    {
      auto last = std::remove(L, L + L_size, sn);
      assert(last + 1 == L + L_size);
      --L_size;
    }
    else // if startpoint - insert
    {
      if (stack_pos->isnot_top()) {
        self().SetCurSegAndPoint(sn);
        i_f.FindIntI(self(), sn, stack_pos); // get internal intersections
      }
      else {
        self().SetCurPointAtBeg(sn);
      }
      int4 i = L_size;
      for (auto _L = L - 1; (i != 0) && (!self().UnderCurPoint(_L[i])); --i)
        L[i] = _L[i];
      L[i] = sn;
      ++L_size;
    }
  }

  // --- Reset / Registration helpers ---

  void Reset()
  {
    MY_FREE_ARR_MACRO(ENDS);
    MY_FREE_ARR_MACRO(seg_L_rank);
    MY_FREE_ARR_MACRO(seg_R_rank);
    self().ResetDerived();
  }

  void ResetRegistration()
  {
    if (factory) {
      if (registrator) registrator->Reset();
      factory->Reset();
    }
  }

  void CombineRegData()
  {
    if (factory) factory->combine_reg_data();
  }

  void IntersectionsFindingDone() {}

  // --- SortAt ---

  void SortAt(uint4 rank, uint4 n, int4* L)
  {
    self().SetCurStripeLeft(rank);
    std::sort(L, L + n, [this](int4 s1, int4 s2) { return self().LBelow(s1, s2); });
  }

  // --- Registrator access ---

  void SetRegistrator(IntersectionRegistrator* r) { registrator = r; }
  IntersectionRegistrator* GetRegistrator() { return registrator; }

  // --- InitClone ---

  void InitClone(uint4 n_threads) {
    if (factory) factory->InitClone(n_threads);
  }

  // --- clone / unclone ---

  void clone(Derived& src, IntersectionRegistrator* r)
  {
    clone_of = &src;
    ENDS = src.ENDS;
    seg_L_rank = src.seg_L_rank;
    seg_R_rank = src.seg_R_rank;
    factory = nullptr; // clones don't own factory
    self().CloneDerived(src);
    SetRegistrator(r);
  }

  void unclone()
  {
    if (registrator) registrator->Flash();
    if (clone_of == nullptr) return;
    ENDS = nullptr;
    seg_L_rank = nullptr;
    seg_R_rank = nullptr;
    clone_of = nullptr;
    self().UncloneDerived();
  }

protected:
  uint4* ENDS = nullptr;
  uint4* seg_L_rank = nullptr;
  uint4* seg_R_rank = nullptr;
  CRegistratorFactory<IntersectionRegistrator>* factory = nullptr;
  IntersectionRegistrator* registrator = nullptr;
  Derived* clone_of = nullptr;

  Derived& self() { return static_cast<Derived&>(*this); }
  const Derived& self() const { return static_cast<const Derived&>(*this); }
};

// =============================================================================
// Segment collection concept hierarchy.
//
// The concepts are layered so that each algorithm requires only what it actually
// uses. A collection author can check exactly which subset of the interface is
// needed for the algorithms they care about.
//
//   SegCollCore              -- shared by ALL algorithms
//     ├─ TrivialSegColl      -- + TrivCurSegIntWith              (CTrivialIntFinder)
//     ├─ SweepSegColl         -- + endpoint sorting, SSCurSegIntWith  (CSimpleSweepIntFinder)
//     └─ BalabanSegColl       -- + stripe/ordering/InsDel/Split/Merge (CFastIntFinder, COptimalIntFinder)
//          └─ OptimalSegColl -- + clone/registrator management   (parallel find_intersections)
//          └─ ParallelSegColl -- + clone/registrator management   (parallel find_intersections)
// =============================================================================

// ---------------------------------------------------------------------------
// Core: the absolute minimum every algorithm needs.
// ---------------------------------------------------------------------------
template<class C>
concept SegCollCore = requires(C c, C & cref, uint4 u, int4 i) {
  // Segment count
  { c.GetSegmNumb() } -> std::convertible_to<uint4>;
  // Set current segment for subsequent intersection queries
  { c.SetCurSeg(u) };
  // Post-processing hook (may be empty)
  { c.IntersectionsFindingDone() };
  // Compile-time flags
  { C::get_coll_flag(_Coll_flags::line_segments) } -> std::same_as<_Coll_flag_state>;
};

// ---------------------------------------------------------------------------
// Trivial O(N^2) algorithm: Core + pairwise check.
// Used by: CTrivialIntFinder
// ---------------------------------------------------------------------------
template<class C>
concept TrivialSegColl = SegCollCore<C> && requires(C c, int4 i) {
  // Check and register intersection of current segment with segment i
  { c.TrivCurSegIntWith(i) } -> std::convertible_to<bool>;
};

// ---------------------------------------------------------------------------
// Simple plane sweep: Core + endpoint enumeration + sweep check.
// Used by: CSimpleSweepIntFinder
// ---------------------------------------------------------------------------
template<class C>
concept SweepSegColl = SegCollCore<C> && requires(C c, uint4 u, uint4 * p) {
  // Endpoint encoding (static)
  { C::is_last(u) } -> std::convertible_to<bool>;
  { C::get_segm(u) } -> std::convertible_to<uint4>;
  // Build sorted endpoint list (caller-allocated buffer)
  { c.PrepareEndpointsSortedList(p) } -> std::convertible_to<uint4>;
  // Sweep-style intersection check with current segment
  { c.SSCurSegIntWith((int4)0) } -> std::convertible_to<bool>;
};

// ---------------------------------------------------------------------------
// Balaban-family algorithms (fast, optimal).
// Used by: CFastIntFinder, COptimalIntFinder
//
// This is the largest interface because these algorithms use stripes,
// staircases, Split/Merge, InsDel, and directional intersection queries.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Fast algorithm: Balaban + stripe-based intersection finding + InsDel + Split/Merge.
// Used by: CFastIntFinder::find_intersections(n_threads, segments)
// ---------------------------------------------------------------------------

template<class C>
concept BalabanSegColl = SegCollCore<C> && requires(
  C c, const C cc,
  uint4 u, int4 i, int4 * p, int4 * q
  ) {
  // Endpoint encoding (static)
    { C::is_last(u) } -> std::convertible_to<bool>;
    { C::get_segm(u) } -> std::convertible_to<uint4>;

    // Prepare sorted endpoints, ranks, statistics
    { c.Prepare() } -> std::same_as<PrepareResult>;

    // Stripe bounds
    { c.SetCurStripe(u, u) };
    { c.SetCurStripeLeft(u) };
    { c.SetCurStripeRight(u) };

    // Ordering at stripe boundaries
    { cc.LBelow(i, i) } -> std::convertible_to<bool>;
    { cc.RBelow(i, i) } -> std::convertible_to<bool>;

    // Current point operations (for InsDel and FindIntI)
    { c.SetCurPointAtBeg(u) };
    { c.SetCurSegAndPoint(u) };
    { cc.UnderCurPoint(i) } -> std::convertible_to<bool>;

    // Rank lookups
    { cc.GetSegR(u) } -> std::convertible_to<uint4>;
    { cc.GetSegL(u) } -> std::convertible_to<uint4>;
    { cc.PointAtRank(u) } -> std::convertible_to<uint4>;

    // Segment clipping to stripe
    { c.SetCurSegCutBE(u) };
    { c.SetCurSegCutBeg(u) };
    { c.SetCurSegCutEnd(u) };

    // Directional intersection finding (single-segment overloads)
    { c.FindCurSegIntDownWith(i) } -> std::convertible_to<bool>;
    { c.FindCurSegIntUpWith(i) } -> std::convertible_to<bool>;

    // Directional intersection finding (range overloads)
    { c.FindCurSegIntDownWith(p, q) };
    { c.FindCurSegIntUpWith(p, q) };


    // Sweep-style check (used in SearchInStrip via SISFindR)
    //{ c.SSCurSegIntWith(i) } -> std::convertible_to<bool>;

    // Sorting segments at a given endpoint rank
    { c.SortAt(u, u, p) };
};

// ---------------------------------------------------------------------------
// Optimal algorithm: Balaban + optimal intersection finding (no InsDel, no Split/Merge).
// Used by: CFastIntFinder::find_intersections(n_threads, segments)
// ---------------------------------------------------------------------------
template<class C>
concept OptimalSegColl = BalabanSegColl<C> && requires(C c, int4 i) {
  // Intersection check without registration 
  { c.IsIntersectsCurSegDown(i) } -> std::convertible_to<bool>;
  { c.IsIntersectsCurSegUp(i) } -> std::convertible_to<bool>;
};
// ---------------------------------------------------------------------------
// Parallel execution: Balaban + clone/registrator management.
// Used by: CFastIntFinder::find_intersections(n_threads, segments)
// ---------------------------------------------------------------------------
template<class C>
concept ParallelSegColl = BalabanSegColl<C> && requires(C c, C & cref, uint4 u) {
  // Clone construction: C(C&, uint4 thread_index)
  { C(cref, u) };
  // Prepare clone registrators
  { c.InitClone(u) };
  // Registrator access
  { c.GetRegistrator() };
  { c.SetRegistrator(c.GetRegistrator()) };
  // Post-parallel data merging
  { c.ResetRegistration() };
  { c.CombineRegData() };
};