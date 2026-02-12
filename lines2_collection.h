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
// NEW IMPLEMENTATION

#include "segments.h"
#include "collection_base.h"

template<class IntersectionRegistrator>
class CLine2SegmentCollection
  : public CollectionBase<CLine2SegmentCollection<IntersectionRegistrator>, IntersectionRegistrator>
{
  using Base = CollectionBase<CLine2SegmentCollection<IntersectionRegistrator>, IntersectionRegistrator>;
  friend Base;

  using Base::ENDS;
  using Base::seg_L_rank;
  using Base::seg_R_rank;
  using Base::factory;
  using Base::registrator;
  using Base::clone_of;

public:
  using this_T = CLine2SegmentCollection;
  static constexpr _Coll_flag_state get_coll_flag(_Coll_flags flag) {
    if (flag == _Coll_flags::line_segments)
      return _Coll_flag_state::state_true;
    return _Coll_flag_state::state_unimplemented;
  }

  using Base::Prepare;
  using Base::InsDel;
  using Base::Reset;
  using Base::ResetRegistration;
  using Base::CombineRegData;
  using Base::IntersectionsFindingDone;
  using Base::SortAt;
  using Base::InitClone;
  using Base::GetSegR;
  using Base::GetSegL;
  using Base::PointAtRank;
  using Base::SetRegistrator;
  using Base::GetRegistrator;

  // --- Endpoint encoding ---

  static bool is_last(uint4 pt) { return pt & 2; }
  static uint4 get_segm(uint4 pt) { return pt >> 2; }

  // --- Segment count ---

  uint4 GetSegmNumb() const { return N; }

  // --- Stripe bounds ---

  void SetCurStripe(uint4 left_rank, uint4 right_rank) {
    B = GetX(ENDS[left_rank]);
    E = GetX(ENDS[right_rank]);
  }

  void SetCurStripeLeft(uint4 left_rank) {
    B = GetX(ENDS[left_rank]);
  }

  void SetCurStripeRight(uint4 right_rank) {
    E = GetX(ENDS[right_rank]);
  }

  // --- Current point / segment ---

  void SetCurPointAtBeg(uint4 s) {
    cur_point = collection[s].BegPoint();
  }

  void SetCurSegAndPoint(uint4 s) {
    cur_point = collection[s].BegPoint();
    SetCurSeg(s);
    if constexpr ((IntersectionRegistrator::reg_type & _RegistrationType::point) == 0)
      active_end = cur_seg.EndPoint();
  }

  void SetCurSeg(uint4 s) {
    cur_seg_idx = s;
    cur_seg = collection[s];
  }

  void SetCurSegCutBE(uint4 s) {
    SetCurSeg(s);
    cur_seg.x2 = MIN(cur_seg.x2, E);
    if constexpr ((IntersectionRegistrator::reg_type & _RegistrationType::point) == 0)
      active_end = cur_seg.EndPoint();
    else
      cur_seg.x1 = MAX(cur_seg.x1, B);
  }

  void SetCurSegCutBeg(uint4 s) {
    SetCurSeg(s);
    cur_seg.x1 = MAX(cur_seg.x1, B);
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      active_end = cur_seg.EndPoint();
  }

  void SetCurSegCutEnd(uint4 s) {
    SetCurSeg(s);
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      active_end = cur_seg.BegPoint();
    else
      cur_seg.x2 = MIN(cur_seg.x2, E);
  }

  // --- Ordering predicates ---

  bool LBelow(int4 s_1, int4 s_2) const {
    auto s1 = collection + s_1;
    auto s2 = collection + s_2;
    return ((B * (s2->a - s1->a) + s2->b - s1->b) > 0);
  }

  bool RBelow(int4 s_1, int4 s_2) const {
    auto s1 = collection + s_1;
    auto s2 = collection + s_2;
    return ((E * (s2->a - s1->a) + s2->b - s1->b) > 0);
  }

  // --- Intersection detection & registration ---

  bool TrivCurSegIntWith(int4 s_) {
    auto& s1 = collection[s_];
    auto& s2 = cur_seg;
    auto x1 = MAX(s1.x1, s2.x1);
    auto x2 = MIN(s1.x2, s2.x2);
    if (x1 >= x2) return false;
    auto da = s2.a - s1.a;
    auto db = s1.b - s2.b;
    if ((da * x1 > db) ^ (da * x2 < db)) return false;
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0) {
      auto x = db / da;
      registrator->register_pair_and_point(cur_seg_idx, s_, TPlaneVect(x, s1.YAtX(x)));
    }
    else
      registrator->register_pair(cur_seg_idx, s_);
    return true;
  }

  bool SSCurSegIntWith(int4 s_) {
    auto& s1 = collection[s_];
    auto& s2 = cur_seg;
    auto da = s2.a - s1.a;
    auto db = s1.b - s2.b;
    if ((da * s2.x1 > db) ^ (da * MIN(s1.x2, s2.x2) < db)) return false;
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0) {
      auto x = db / da;
      registrator->register_pair_and_point(cur_seg_idx, s_, TPlaneVect{ x, s1.YAtX(x) });
    }
    else
      registrator->register_pair(cur_seg_idx, s_);
    return true;
  }

  template<bool do_register = true>
  bool FindIntWith(int4 s_) {
    auto& s1 = collection[s_];
    auto& s2 = cur_seg;
    auto da = s2.a - s1.a;
    auto db = s1.b - s2.b;
    if ((da * MAX(s1.x1, s2.x1) - db > 0) ^ (da * MIN(s1.x2, s2.x2) - db < 0)) return false;
    if constexpr (!do_register) return true;
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0) {
      REAL x = db / da;
      registrator->register_pair_and_point(cur_seg_idx, s_, TPlaneVect(x, s1.YAtX(x)));
    }
    else
      registrator->register_pair(cur_seg_idx, s_);
    return true;
  }

  bool FindCurSegIntDownWith(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return UnderActiveEnd(s_) ? false : (registrator->register_pair(cur_seg_idx, s_), true);
    return FindIntWith<true>(s_);
  }

  auto FindCurSegIntDownWith(int4* s_, int4* last) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0) {
      auto r = registrator;
      auto cs = cur_seg_idx;
      while ((THIS_HAS_SENTINELS || (s_ != last)) && !UnderActiveEnd(*s_)) {
        r->register_pair(cs, *s_);
        --s_;
      }
      return s_;
    }
    while ((THIS_HAS_SENTINELS || (s_ != last)) && FindIntWith(*s_))
      --s_;
    return s_;
  }

  auto FindCurSegIntUpWith(int4* s_, int4* last) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0) {
      auto r = registrator;
      auto cs = cur_seg_idx;
      while ((THIS_HAS_SENTINELS || (s_ != last)) && UnderActiveEnd(*s_)) {
        r->register_pair(cs, *s_);
        ++s_;
      }
      return s_;
    }
    while ((THIS_HAS_SENTINELS || (s_ != last)) && FindIntWith(*s_))
      ++s_;
    return s_;
  }

  bool FindCurSegIntUpWith(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return UnderActiveEnd(s_) ? registrator->register_pair(cur_seg_idx, s_), true : false;
    return FindIntWith(s_);
  }

  bool IsIntersectsCurSegDown(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return !UnderActiveEnd(s_);
    return FindIntWith<false>(s_);
  }

  bool IsIntersectsCurSegUp(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return UnderActiveEnd(s_);
    return FindIntWith<false>(s_);
  }

  bool UnderCurPoint(int4 s_) const { return collection[s_].a * cur_point.x + collection[s_].b < cur_point.y; }
  bool UnderActiveEnd(int4 s_) const { return collection[s_].a * active_end.x + collection[s_].b < active_end.y; }

  // --- Endpoint sorting ---

  uint4 PrepareEndpointsSortedList(uint4* epoints) {
    auto NN = N << 1;
    for (uint4 i = 0; i < NN; ++i) epoints[i] = i << 1;
    std::sort(epoints, epoints + NN,
      [x = ends](uint4 pt1, uint4 pt2) {
        return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
      }
    );
    return NN;
  }

  int4 get_sentinel(bool is_top_sentinel) {
    return N + is_top_sentinel;
  }

  // --- SVG ---

  void coll_to_SVG(chostream* SVG_stream) {
    if (!SVG_stream) return;
    int4 n = MIN(max_SVG_items, N);
    auto mmr = get_mmrect01(collection, n);
    *SVG_stream << "<svg height='100%' width='100%' viewBox='";
    *SVG_stream << mmr.ld.x << " " << mmr.ld.y << " "
      << mmr.rt.x - mmr.ld.x << " " << mmr.rt.y - mmr.ld.y << "' transform='scale(1, -1)'>\n";
    for (int4 i = 0; i < n; ++i) collection[i].write_SVG(i, SVG_stream);
  }

  // --- Construction / Destruction ---

  CLine2SegmentCollection() {}

  CLine2SegmentCollection(const CollectionOptions& co, void* c, IntersectionRegistrator* r) {
    InitDerived(co.n, c, r);
  }

  CLine2SegmentCollection(const CollectionOptions& co, void* c, CRegistratorFactory<IntersectionRegistrator>* f) {
    factory = f;
    factory->PrepareAlloc(co.n);
    InitDerived(co.n, c, factory->GetRegistrator(0));
  }

  CLine2SegmentCollection(CLine2SegmentCollection& coll, IntersectionRegistrator* r) {
    Base::clone(coll, r);
  }

  CLine2SegmentCollection(CLine2SegmentCollection& coll, uint4 thread_index) {
    if (coll.factory)
      Base::clone(coll, coll.factory->GetRegistrator(thread_index));
    else
      Base::clone(coll, coll.GetRegistrator());
  }

  ~CLine2SegmentCollection() {
    Base::unclone();
    Base::Reset();
  }

private:
  auto GetX(uint4 pt) const { return ends[pt]; }

  void InitDerived(uint4 n, void* c, IntersectionRegistrator* r) {
    N = n;
    collection = reinterpret_cast<TLineSegment2*>(c);
    static_assert(sizeof(TLineSegment2) == 4 * sizeof(REAL), "one segment should occupy 4 reals space");
    static_assert(offsetof(TLineSegment2, x1) == 0, "TLineSegment2.x1 should have 0 offset");
    static_assert(offsetof(TLineSegment2, x2) == 2 * sizeof(REAL), "TLineSegment2.x2 should have 2 reals offset");
    ends = reinterpret_cast<REAL*>(collection);
    Base::SetRegistrator(r);
  }

  void ResetDerived() {}

  void CloneDerived(CLine2SegmentCollection& src) {
    N = src.N;
    collection = src.collection;
    ends = src.ends;
  }

  void UncloneDerived() {
    collection = nullptr;
    ends = nullptr;
  }

  uint4 N = 0;
  REAL B, E;
  TPlaneVect cur_point, active_end;
  uint4 cur_seg_idx = 0xFFFFFFFF;
  TLineSegment2 cur_seg;
  TLineSegment2* collection = nullptr;
  REAL* ends = nullptr;
};
