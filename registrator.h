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
#include "utils.h"




#define MY_FREE_ARR_MACRO(a) \
    if (a!=nullptr) {                 \
        delete[] a;          \
        a = nullptr;            \
    }




template <class ipoint>
class JustCountingRegistrator
{
public:
  static constexpr uint4 reg_type = _RegistrationType::count;
  double counter = 0;
  void register_pair(uint4 s1, uint4 s2) { 
#ifdef PRINT_SEG_AND_INT
    if (s1 < s2)
      printf("alt int %i %i\n", s1, s2);
    else
      printf("alt int %i %i\n", s2, s1);

#endif 

    ++counter; 
  };
  void register_pair_and_point(uint4 s1, uint4 s2, ipoint& p) { register_pair(s1, s2); };
  void combine_reg_data(uint4 n_threads, JustCountingRegistrator* additional_reg_obj[])
  {
    for (int i = 0; i < n_threads - 1; i++)
      counter += additional_reg_obj[i]->counter;
  };
  void Alloc(uint4 _N) {};

  // to return some statistics about registered intersections;
  double get_stat(uint4 stat_type = 0) { return counter; };
};

template <class ipoint>
class PerSegmCountingRegistrator
{
public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments;
  ~PerSegmCountingRegistrator() { if (segm_counters != nullptr) { delete[] segm_counters; segm_counters = nullptr; } };

  double counter = 0;

  void register_pair(uint4 s1, uint4 s2) {
    ++counter;
    ++segm_counters[s1];
    ++segm_counters[s2];
  };
  void register_pair_and_point(uint4 s1, uint4 s2, ipoint& p) { register_pair(s1, s2); };

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

private:
  uint4 N;
  uint4* segm_counters = nullptr;
};

using SimpleCounter = JustCountingRegistrator<TPlaneVect>;
using PerSegmCounter = PerSegmCountingRegistrator<TPlaneVect>;

#endif // !REGISTRATOR_FOR_SEGMENT_INTERSECTION