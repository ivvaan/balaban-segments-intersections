#ifndef REGISTRATOR_FOR_SEGMENT_INTERSECTION
#define REGISTRATOR_FOR_SEGMENT_INTERSECTION


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
#include <algorithm>
#define MY_FREE_ARR_MACRO(a) \
    if (a!=nullptr) {                 \
        delete[] a;          \
        a = nullptr;            \
    }



#include "utils.h"
template <class ipoint>
class JustCountingRegistrator
{
  class CIntersectionBuilder {
  public:
    inline CIntersectionBuilder* register_segments(uint4 s1, uint4 s2)
    {
#ifdef PRINT_SEG_AND_INT
      if (s1<s2)
        printf("alt int %i %i\n", s1, s2);
      else
        printf("alt int %i %i\n", s2, s1);

#endif 
      return this;
    };
    inline CIntersectionBuilder* register_point(ipoint &) { return this; };
    inline void end_registration() {};

  };
public:
  static constexpr uint4 reg_type = _RegistrationType::count;
  double counter = 0;
  inline CIntersectionBuilder* begin_registration(uint4 inters_numb) { counter += inters_numb; return nullptr; };
  void combine_reg_data(uint4 n_threads, JustCountingRegistrator *additional_reg_obj[])
  {
    for (int i = 0; i < n_threads - 1; i++)
      counter += additional_reg_obj[i]->counter;
  };
  inline void Alloc(uint4 _N) {};

  // to return some statistics about registered intersections;
  double get_stat(uint4 stat_type = 0) { return counter; };

  void write_SVG(uint4 alg, chostream* SVG_text) {};

};

template <class ipoint>
class PerSegmCountingRegistrator
{
  class CIntersectionBuilder {
  public:
    inline CIntersectionBuilder* register_segments(uint4 s1, uint4 s2)
    {
#ifdef PRINT_SEG_AND_INT
      if (s1<s2)
        printf("alt int %i %i\n", s1, s2);
      else
        printf("alt int %i %i\n", s2, s1);

#endif 
      ++segm_counters[s1];
      ++segm_counters[s2];
      return this;
    };
    inline CIntersectionBuilder* register_point(ipoint &) { return this; };
    inline void end_registration() {};

    inline void set_counters(uint4 *c) { segm_counters = c; };
  private:
    uint4 *segm_counters = nullptr;
  };
public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments;
  ~PerSegmCountingRegistrator() { if (segm_counters != nullptr) { delete[] segm_counters; segm_counters = nullptr; } };

  double counter = 0;
  inline CIntersectionBuilder* begin_registration(uint4 inters_numb) { counter += inters_numb; return &builder; };
  inline void Alloc(uint4 _N)
  {
    N = _N;
    segm_counters = new uint4[N];
    std::memset(segm_counters, 0, N * sizeof(*segm_counters));
    builder.set_counters(segm_counters);
  };
  void combine_reg_data(uint4 n_threads, PerSegmCountingRegistrator *to_add[])
  {
    for (uint4 r = 0; r < n_threads-1; ++r)
    {
      counter += to_add[r]->counter;
      uint4 *added = to_add[r]->segm_counters;
      for (uint4 i = 0; i < N; ++i)
        segm_counters[i] += added[i];
    }

  };

  // to return some statistics about registered intersections;
  double get_stat(uint4 stat_type = 0) 
  { 
    // if stat_type nonzero return maximal number of intersections per segment;
    if (stat_type== _Registrator::per_segm_reg_max_per_segm_stat)
      return *std::max_element(segm_counters, segm_counters+N);
    return counter; 
  };

  void write_SVG(uint4 alg, chostream* SVG_text) {};

private:
  uint4 N;
  uint4 *segm_counters = nullptr;

  CIntersectionBuilder builder;

};

template <class ipoint>
class PairAndPointRegistrator
{
  struct IntersectingPair {
    uint4 s1, s2, points_beg,points_end;
    void Init(uint4 s_1, uint4 s_2, PairAndPointRegistrator& r) {
      points_beg = r.get_free_pt_pos();
      if (s_2 < s_1) {
        s1 = s_2; s2 = s_1;
      }
      s1 = s_1; s2 = s_2;
    }
    void add_point(ipoint &p, PairAndPointRegistrator& r) {
      points_end = r.add_point(p);
    }
    void inc_pt_num(uint4 inc) {
      points_beg += inc;
      points_end += inc;
    }
  };

  class CIntersectionBuilder {
    IntersectingPair cur_pair;
  public:
    CIntersectionBuilder* register_segments(uint4 s1, uint4 s2)
    {
#ifdef PRINT_SEG_AND_INT
      if (s1 < s2)
        printf("alt int %i %i\n", s1, s2);
      else
        printf("alt int %i %i\n", s2, s1);

#endif
      cur_pair.Init(s1, s2, *r);
      return this;
    };

    CIntersectionBuilder* register_point(ipoint& p) {
      cur_pair.add_point(p,*r);
      return this;
    };

    void end_registration() {
      r->add_pair(cur_pair);
    };

    inline void set_registrator(PairAndPointRegistrator *reg) { r = reg; };
  private:
    PairAndPointRegistrator *r;
  };

public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments + _RegistrationType::point;

  PairAndPointRegistrator() {
    N_pairs = N_points = 0;
    pairs_size = points_size = 1024 * 1024;
    segm_pairs = new IntersectingPair[pairs_size];
    int_points = new ipoint[points_size];
    builder.set_registrator(this);
  };

  ~PairAndPointRegistrator() {
    if (segm_pairs != nullptr) { delete[] segm_pairs; segm_pairs = nullptr; }
    if (int_points != nullptr) { delete[] int_points; int_points = nullptr; }
    N_pairs = N_points = pairs_size = points_size = 0;
  };

  double counter = 0;

  CIntersectionBuilder* begin_registration(uint4 inters_numb){
    counter += inters_numb; 
    if (N_pairs >= pairs_size)realloc_pairs();
    if (N_points + inters_numb > points_size)realloc_points();
    return &builder;
  };

  void add_pair(IntersectingPair& p) {
    segm_pairs[N_pairs++] = p;
  };

  uint4 add_point(ipoint& p) {
    int_points[N_points]=p;
    return N_points++;
  };

  uint4 get_free_pt_pos() { return N_points; };

  void realloc_points(uint4 sz=0) {
    auto ns = sz?sz:points_size * 2;
    auto sp = new ipoint[ns];
    std::move(int_points, int_points+ points_size,sp);
    delete[] int_points;
    int_points = sp;
    points_size = ns;
  };

  void realloc_pairs(uint4 sz = 0) {
    auto ns = sz?sz:pairs_size * 2;
    auto sp = new IntersectingPair[ns];
    memcpy(sp, segm_pairs, points_size * sizeof(*segm_pairs));
    delete[] segm_pairs;
    segm_pairs = sp;
    pairs_size = ns;
  };

  void inc_pt_num(uint4 inc) {
    for (uint4 i = 0; i < N_pairs; i++)
      segm_pairs[i].inc_pt_num(inc);
  };

  void Alloc(uint4) {};

  void combine_reg_data(uint4 n_threads, PairAndPointRegistrator* to_add[])
  {
    auto tot_pairs = N_pairs;
    auto tot_points = N_points;

    for (uint4 r = 0; r < n_threads - 1; ++r) {
      counter += to_add[r]->counter;
      tot_pairs += to_add[r]->N_pairs;
      tot_points += to_add[r]->N_points;
    }
    if (tot_pairs >= pairs_size)realloc_pairs(tot_pairs);
    if (tot_points >= points_size)realloc_points(tot_points);
    for (uint4 r = 0; r < n_threads - 1; ++r)
    {
      PairAndPointRegistrator* added = to_add[r];
      std::move(added->int_points, added->int_points + added->N_points, int_points + N_points);
      added->inc_pt_num(N_points);
      memcpy(segm_pairs + N_pairs, added->segm_pairs, added->N_pairs * sizeof(*segm_pairs));
      N_points += added->N_points;
      N_pairs += added->N_pairs;
    }

  };

  // to return some statistics about registered intersections;
  double get_stat(uint4 stat_type = 0)
  {
    return counter;
  };

  void write_SVG(uint4 alg, chostream* SVG_text) {
    if (SVG_text) {
      for (uint4 i = 0; i < min(N_points, max_SVG_items); ++i) {
        *SVG_text << "<circle id='int" << i;
        *SVG_text << "' cx='" << int_points[i].getX();
        *SVG_text << "' cy='" << int_points[i].getY();
        switch (alg)
        {
        case _Algorithm::triv:*SVG_text << "' class='triv'/>\n"; break;
        case _Algorithm::simple_sweep:*SVG_text << "' class='ssw'/>\n"; break;
        case _Algorithm::fast:*SVG_text << "' class='fast'/>\n"; break;
        case _Algorithm::optimal:*SVG_text << "' class='optimal'/>\n"; break;
        case _Algorithm::fast_parallel:*SVG_text << "' class='parallel'/>\n"; break;
        case _Algorithm::mem_save:*SVG_text << "' class='mem_save'/>\n"; break;
        default:
          *SVG_text << "' class='algorithm'/>\n";
          break;
        }
      }
    };
  };

private:
  uint4 N_pairs,N_points,pairs_size,points_size;
  IntersectingPair* segm_pairs;
  ipoint* int_points;

  CIntersectionBuilder builder;

};


using SimpleCounter=JustCountingRegistrator<TPlaneVect> ;
using PerSegmCounter=PerSegmCountingRegistrator<TPlaneVect> ;
using TrueRegistrator = PairAndPointRegistrator<TPlaneVect> ;

#endif // !REGISTRATOR_FOR_SEGMENT_INTERSECTION