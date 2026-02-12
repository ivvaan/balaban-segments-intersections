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
class CLine1SegmentCollection
  : public CollectionBase<CLine1SegmentCollection<IntersectionRegistrator>, IntersectionRegistrator>
{
  using Base = CollectionBase<CLine1SegmentCollection<IntersectionRegistrator>, IntersectionRegistrator>;
  friend Base;

  // --- fields accessible by Base via CloneDerived/UncloneDerived/InitDerived ---
  using Base::ENDS;
  using Base::seg_L_rank;
  using Base::seg_R_rank;
  using Base::factory;
  using Base::registrator;
  using Base::clone_of;

public:
  using this_T = CLine1SegmentCollection;
  static constexpr _Coll_flag_state get_coll_flag(_Coll_flags flag) {
    if (flag == _Coll_flags::line_segments)
      return _Coll_flag_state::state_true;
    if (flag == _Coll_flags::needs_SetCurSegCutBE_at_start)
      return _Coll_flag_state::state_true;
    return _Coll_flag_state::state_unimplemented;
  }

  // --- Prepare, InsDel, Reset, SortAt, etc. come from Base ---
  // Bring them into public scope:
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

  static bool is_last(uint4 pt) { return pt & 1; }
  static uint4 get_segm(uint4 pt) { return pt >> 1; }

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
    curB = cur_seg.org.x;
    curE = cur_seg.org.x + cur_seg.shift.x;
    is_chopped = false;
  }

  void SetCurSegCutBE(uint4 s) {
    SetCurSeg(s);
    curB = MAX(B, curB);
    curE = MIN(E, curE);
    if constexpr ((IntersectionRegistrator::reg_type & _RegistrationType::point) == 0) {
      active_end.x = curE;
      active_end.y = cur_seg.YAtX(curE);
      if (E == curE) {
        is_chopped = true;
        chopped_Y[s] = active_end.y;
      }
    }
  }

  void SetCurSegCutBeg(uint4 s) {
    SetCurSeg(s);
    curB = MAX(B, curB);
    if constexpr ((IntersectionRegistrator::reg_type & _RegistrationType::point) == 0)
      active_end = cur_seg.EndPoint();
  }

  void SetCurSegCutEnd(uint4 s) {
    SetCurSeg(s);
    curE = MIN(E, curE);
    if constexpr ((IntersectionRegistrator::reg_type & _RegistrationType::point) == 0)
      active_end = cur_seg.BegPoint();
  }

  // --- Ordering predicates ---

  bool LBelow(int4 s_1, int4 s_2) const {
    auto& s1 = collection[s_1];
    auto& s2 = collection[s_2];
    return s1.YAtX_Numerator(B) * s2.shift.x < s2.YAtX_Numerator(B) * s1.shift.x;
  }

  bool RBelow(int4 s_1, int4 s_2) const {
    auto& s1 = collection[s_1];
    auto& s2 = collection[s_2];
    return s1.YAtX_Numerator(E) * s2.shift.x < s2.YAtX_Numerator(E) * s1.shift.x;
  }

  // --- Intersection detection & registration ---

  template <bool do_register = true>
  bool FindIntWith(REAL x1, REAL x2, uint4 s_) {
    auto& s2 = cur_seg;
    auto& s1 = collection[s_];
    TPlaneVect delt = s2.org - s1.org;
    REAL prod = s1.shift % s2.shift;
    auto mul = s1.shift % delt;
    if ((mul > 0) ^ (mul + prod > 0)) {
      mul /= prod;
      auto xc = s2.org.x - mul * s2.shift.x;
      if (((xc <= x1) || (xc > x2))) return false;
      if (!do_register) return true;
      if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0)
        registrator->register_pair_and_point(cur_seg_idx, s_, TPlaneVect(xc, s2.org.y - mul * s2.shift.y));
      else
        registrator->register_pair(cur_seg_idx, s_);
      return true;
    }
    return false;
  }

  bool TrivCurSegIntWith(int4 s_) {
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    auto delt = s2.org - s1.org;
    if ((s1.shift.x < delt.x) || (s2.shift.x < -delt.x))
      return false;
    REAL prod = s1.shift % s2.shift, mul;
    if (((mul = s1.shift % delt) > 0) ^ (mul + prod > 0))
      if (((mul = delt % s2.shift) > 0) ^ (mul - prod > 0)) {
        if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0)
          registrator->register_pair_and_point(cur_seg_idx, s_, s1.org + (mul / prod) * s1.shift);
        else
          registrator->register_pair(cur_seg_idx, s_);
        return true;
      };
    return false;
  }

  bool SSCurSegIntWith(int4 s_) {
    auto& s1 = cur_seg;
    auto& s2 = collection[s_];
    auto delt = s2.org - s1.org;
    REAL prod = s1.shift % s2.shift, mul;
    if (((mul = s1.shift % delt) > 0) ^ (mul + prod > 0))
      if (((mul = delt % s2.shift) > 0) ^ (mul - prod > 0)) {
        if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) != 0)
          registrator->register_pair_and_point(cur_seg_idx, s_, s1.org + (mul / prod) * s1.shift);
        else
          registrator->register_pair(cur_seg_idx, s_);
        return true;
      };
    return false;
  }

  bool FindCurSegIntDownWith(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return UnderActiveEnd(s_) ? false : (registrator->register_pair(cur_seg_idx, s_), true);
    return FindIntWith(curB, curE, s_);
  }

  auto FindCurSegIntDownWith(int4* s_, int4* last) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0) {
      auto r = registrator;
      auto cs = cur_seg_idx;
      if (is_chopped) {
        REAL Y = active_end.y;
        while ((THIS_HAS_SENTINELS || (s_ != last)) && (Y < chopped_Y[*s_])) {
          r->register_pair(cs, *s_);
          --s_;
        }
      }
      else
        while ((THIS_HAS_SENTINELS || (s_ != last)) && !UnderActiveEnd(*s_)) {
          r->register_pair(cs, *s_);
          --s_;
        }
      return s_;
    }
    while ((THIS_HAS_SENTINELS || (s_ != last)) && FindIntWith(curB, curE, *s_))
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
    while ((THIS_HAS_SENTINELS || (s_ != last)) && FindIntWith(curB, curE, *s_))
      ++s_;
    return s_;
  }

  bool FindCurSegIntUpWith(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return UnderActiveEnd(s_) ? registrator->register_pair(cur_seg_idx, s_), true : false;
    return FindIntWith(curB, curE, s_);
  }

  bool IsIntersectsCurSegDown(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return !UnderActiveEnd(s_);
    return FindIntWith<false>(curB, curE, s_);
  }

  bool IsIntersectsCurSegUp(int4 s_) {
    if constexpr ((_RegistrationType::point & IntersectionRegistrator::reg_type) == 0)
      return UnderActiveEnd(s_);
    return FindIntWith<false>(curB, curE, s_);
  }

  bool UnderCurPoint(int4 s_) const { return collection[s_].under(cur_point); }
  bool UnderActiveEnd(int4 s_) const { return collection[s_].under(active_end); }

  // --- Endpoint sorting ---

  uint4 PrepareEndpointsSortedList(uint4* epoints) {
    auto NN = N * 2;
    for (uint4 i = 0; i < NN; ++i) epoints[i] = i * 2;

    for (uint4 i = 0; i < N; ++i)
      collection[i].shift.x += collection[i].org.x;

    static_assert(sizeof(TLineSegment1) == 4 * sizeof(REAL), "one segment should occupy 4 reals space");
    static_assert(offsetof(TLineSegment1, org.x) == 0, "TLineSegment1.org.x should have 0 offset");
    static_assert(offsetof(TLineSegment1, shift.x) == 2 * sizeof(REAL), "TLineSegment1.shift.x should have 2 reals offset");

    std::sort(epoints, epoints + NN,
      [x = reinterpret_cast<REAL*>(collection)](uint4 pt1, uint4 pt2) {
        return ((x[pt1] < x[pt2]) || ((x[pt1] == x[pt2]) && (pt1 < pt2)));
      }
    );
    for (uint4 i = 0; i < N; ++i)
      collection[i].shift.x -= collection[i].org.x;
    for (uint4 i = 0; i < NN; ++i) epoints[i] /= 2;
    return NN;
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

  CLine1SegmentCollection() {}

  CLine1SegmentCollection(const CollectionOptions& co, void* c, IntersectionRegistrator* r) {
    InitDerived(co.n, c, r);
    AllocChoppedY();
  }

  CLine1SegmentCollection(const CollectionOptions& co, void* c, CRegistratorFactory<IntersectionRegistrator>* f) {
    factory = f;
    factory->PrepareAlloc(co.n);
    InitDerived(co.n, c, factory->GetRegistrator(0));
    AllocChoppedY();
  }

  CLine1SegmentCollection(CLine1SegmentCollection& coll, IntersectionRegistrator* r) {
    Base::clone(coll, r);
    AllocChoppedY();
  }

  CLine1SegmentCollection(CLine1SegmentCollection& coll, uint4 thread_index) {
    if (coll.factory)
      Base::clone(coll, coll.factory->GetRegistrator(thread_index));
    else
      Base::clone(coll, coll.GetRegistrator());
    AllocChoppedY();
  }

  ~CLine1SegmentCollection() {
    Base::unclone();
    Base::Reset();
    MY_FREE_ARR_MACRO(chopped_Y);
  }

private:
  auto GetX(uint4 pt) const {
    return is_last(pt) ? collection[get_segm(pt)].ex() : collection[get_segm(pt)].bx();
  }

  // --- CRTP hooks called by Base ---

  void InitDerived(uint4 n, void* c, IntersectionRegistrator* r) {
    N = n;
    collection = reinterpret_cast<TLineSegment1*>(c);
    Base::SetRegistrator(r);
  }

  void ResetDerived() { /* nothing extra beyond ENDS/seg_L_rank/seg_R_rank */ }

  void CloneDerived(CLine1SegmentCollection& src) {
    N = src.N;
    collection = src.collection;
  }

  void UncloneDerived() {
    collection = nullptr;
  }

  void AllocChoppedY() {
    chopped_Y = new REAL[N + 2];
    if constexpr (THIS_HAS_SENTINELS) {
      chopped_Y[this->get_sentinel(false)] = std::numeric_limits<REAL>::lowest();
      chopped_Y[this->get_sentinel(true)] = std::numeric_limits<REAL>::max();
    }
  }

  TLineSegment1* collection = nullptr;
  REAL* chopped_Y = nullptr;
  TLineSegment1 cur_seg;
  uint4 N = 0;
  REAL B, E, curB, curE;
  TPlaneVect cur_point, active_end;
  uint4 cur_seg_idx = 0xFFFFFFFF;
  bool is_chopped = false;
};

