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