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
#include <tuple>
#include <vector>
#include "utils.h"

enum _RegistrationType
{
  count = 1,
  segments = 2,
  point = 4,
  count_and_segments = 3,
  full = 7
};

template <class ipoint>
class JustCountingRegistrator
{
  double counter = 0;
public:
  static constexpr uint4 reg_type = _RegistrationType::count;
  void register_pair(uint4 s1, uint4 s2) noexcept {
#ifdef PRINT_SEG_AND_INT
    if (s1 < s2)
      printf("alt int %i %i\n", s1, s2);
    else
      printf("alt int %i %i\n", s2, s1);

#endif 

    ++counter; 
  };
  void register_pair_and_point(uint4 s1, uint4 s2, const ipoint& p) { register_pair(s1, s2); };
 
  void combine_reg_data(uint4 n_threads, JustCountingRegistrator* additional_reg_obj[])
  {
    for (int i = 0; i < n_threads - 1; i++)
      counter += additional_reg_obj[i]->counter;
  };
  void Alloc(uint4 _N) {};

  // to return some statistics about registered intersections;
  double get_stat(uint4 stat_type = 0) { return counter; };

  void write_SVG(uint4 alg, chostream* SVG_text) {};

};

template <class ipoint>
class PerSegmCountingRegistrator
{
  uint4 N;
  uint4* segm_counters = nullptr;
  double counter = 0;

public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments;
  ~PerSegmCountingRegistrator() { if (segm_counters != nullptr) { delete[] segm_counters; segm_counters = nullptr; } };


  void register_pair(uint4 s1, uint4 s2) noexcept {
    ++counter;
    ++segm_counters[s1];
    ++segm_counters[s2];
  };
  void register_pair_and_point(uint4 s1, uint4 s2, const ipoint& p) { register_pair(s1, s2); };

  void Alloc(uint4 _N)
  {
    N = _N;
    segm_counters = new uint4[N];
    std::memset(segm_counters, 0, N * sizeof(*segm_counters));
  };
  void combine_reg_data(uint4 n_threads, PerSegmCountingRegistrator* to_add[])
  {
    for (uint4 r = 0; r < n_threads - 1; ++r)
    {
      counter += to_add[r]->counter;
      uint4* added = to_add[r]->segm_counters;
      for (uint4 i = 0; i < N; ++i)
        segm_counters[i] += added[i];
    }

  };

  // to return some statistics about registered intersections;
  double get_stat(uint4 stat_type = 0)
  {
    // if stat_type nonzero return maximal number of intersections per segment;
    if (stat_type == _Registrator::per_segm_reg_max_per_segm_stat)
      return *std::max_element(segm_counters, segm_counters + N);
    return counter;
  };

  void write_SVG(uint4 alg, chostream* SVG_text) {};

};

template <class ipoint>
class PairAndPointRegistrator
{
  using int_rec = std::tuple<uint4, uint4, ipoint>;
  std::vector<int_rec>  intersections;
  double counter = 0;

public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments + _RegistrationType::point;

  PairAndPointRegistrator() {
  };

  ~PairAndPointRegistrator() {
  };

 

  void register_pair(uint4 s1, uint4 s2) noexcept {
    if (counter == 0)
      intersections.reserve(1024 * 1024);
    ++counter;
  };
  void register_pair_and_point(uint4 s1, uint4 s2, const ipoint& p) {
    if (counter == 0)
      intersections.reserve(1024 * 1024);
    ++counter;
    if (s1 < s2)
      intersections.emplace_back(s1, s2, p);
    else
      intersections.emplace_back(s2, s1, p);
  };

  void Alloc(uint4) {};

  void combine_reg_data(uint4 n_threads, PairAndPointRegistrator* to_add[])
  {
    auto tot = intersections.size();

    for (uint4 r = 0; r < n_threads - 1; ++r) {
      counter += to_add[r]->counter;
      tot += to_add[r]->intersections.size();
    }
    intersections.reserve(tot);
    for (uint4 r = 0; r < n_threads - 1; ++r)
    {
      auto &added = to_add[r]->intersections;
      std::move(added.begin(), added.end(), std::back_inserter(intersections));
//      intersections.insert(intersections.end(), added.begin(), added.end());
    }

  };

  // to return some statistics about registered intersections;
  double get_stat(uint4 stat_type = 0)
  {
    return counter;
  };

  void write_SVG(uint4 alg, chostream* SVG_text) {
    if (SVG_text) {
      for (uint4 i = 0; i < MIN(intersections.size(), max_SVG_points); ++i) {
        uint4 s1 = std::get<0>(intersections[i]);
        uint4 s2 = std::get<1>(intersections[i]);
        auto& p = std::get<2>(intersections[i]);
        *SVG_text << "<circle cx='" << p.getX();
        *SVG_text << "' cy='" << p.getY();
        
        switch (alg)
        {
        case _Algorithm::triv:*SVG_text << "' class='triv' id='t" << s1<<"_"<<s2<<"' />\n"; break;
        case _Algorithm::simple_sweep:*SVG_text << "' class='ssw' id='s" << s1<<"_"<<s2<<"' />\n"; break;
        case _Algorithm::fast:*SVG_text << "' class='fast' id='f" << s1<<"_"<<s2<<"' />\n"; break;
        case _Algorithm::optimal:*SVG_text << "' class='optimal' id='o" << s1<<"_"<<s2<<"' />\n"; break;
        case _Algorithm::fast_parallel:*SVG_text << "' class='parallel' id='p" << s1<<"_"<<s2<<"' />\n"; break;
        case _Algorithm::mem_save:*SVG_text << "' class='mem_save' id='m" << s1<<"_"<<s2<<"' />\n"; break;
        default:
          *SVG_text << "' class='algorithm' id='a" << s1<<"_"<<s2<<"' />\n";
          break;
        }
      }
    };
  }


};


using SimpleCounter=JustCountingRegistrator<TPlaneVect> ;
using PerSegmCounter=PerSegmCountingRegistrator<TPlaneVect> ;
using TrueRegistrator = PairAndPointRegistrator<TPlaneVect> ;

#endif // !REGISTRATOR_FOR_SEGMENT_INTERSECTION