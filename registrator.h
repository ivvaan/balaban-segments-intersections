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

/*
#ifdef PRINT_SEG_AND_INT
    if (s1 < s2)
      printf("alt int %i %i\n", s1, s2);
    else
      printf("alt int %i %i\n", s2, s1);

#endif */



template <class ipoint>
class JustCountingRegistrator
{
public:
  static constexpr uint4 reg_type = _RegistrationType::count;
  double counter = 0;
  void register_pair(uint4 s1, uint4 s2) { ++counter; };
  void register_pair_and_point(uint4 s1, uint4 s2,ipoint &&p) { register_pair(s1, s2); };
  void combine_reg_data(uint4 n_threads, JustCountingRegistrator *additional_reg_obj[])
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
public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments;
  ~PerSegmCountingRegistrator() { if (segm_counters != nullptr) { delete[] segm_counters; segm_counters = nullptr; } };

  double counter = 0;

  void register_pair(uint4 s1, uint4 s2) {
    ++counter;
    ++segm_counters[s1];
    ++segm_counters[s2];
  };
  void register_pair_and_point(uint4 s1, uint4 s2, ipoint&& p) { register_pair(s1,s2); };

  void Alloc(uint4 _N)
  {
    N = _N;
    segm_counters = new uint4[N];
    std::memset(segm_counters, 0, N * sizeof(*segm_counters));
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
};

template <class ipoint>
class PairAndPointRegistrator
{
  struct IntersectingPair {
    uint4 s1, s2;
    ipoint p;
    void set(uint4 s_1, uint4 s_2, ipoint &&pp) {
      p = pp;
      if (s_2 < s_1) {
        s1 = s_2; s2 = s_1;
        return;
      }
      s1 = s_1; s2 = s_2;
    }
   
  };


public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments + _RegistrationType::point;
  double counter = 0;

  PairAndPointRegistrator() {
    N_pairs = N_pairs = 0;
    pairs_size =  1024 * 1024;
    segm_pairs = new IntersectingPair[pairs_size];
  };

  ~PairAndPointRegistrator() {
    if (segm_pairs != nullptr) { delete[] segm_pairs; segm_pairs = nullptr; }
    N_pairs = pairs_size = 0;
  };

 

  void register_pair(uint4 s1, uint4 s2) {
    ++counter;
  };
  void register_pair_and_point(uint4 s1, uint4 s2, ipoint&& p) {
    if (N_pairs >= pairs_size)realloc_pairs();
    segm_pairs[N_pairs++].set(s1, s2, p);
  };

  void realloc_pairs(uint4 sz=0) {
    auto ns = sz?sz:pairs_size * 2;
    auto sp = new IntersectingPair[ns];
    std::move(segm_pairs, segm_pairs+ pairs_size,sp);
    delete[] segm_pairs;
    segm_pairs = sp;
    pairs_size = ns;
  };

  void Alloc(uint4) {};

  void combine_reg_data(uint4 n_threads, PairAndPointRegistrator* to_add[])
  {
    auto tot_pairs = N_pairs;

    for (uint4 r = 0; r < n_threads - 1; ++r) {
      counter += to_add[r]->counter;
      tot_pairs += to_add[r]->N_pairs;
    }
    if (tot_pairs >= pairs_size)realloc_pairs(tot_pairs);
    for (uint4 r = 0; r < n_threads - 1; ++r)
    {
      PairAndPointRegistrator* added = to_add[r];
      std::move(added->segm_pairs, added->segm_pairs + added->N_pairs, segm_pairs + N_pairs);
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
      for (uint4 i = 0; i < MIN(N_pairs, max_SVG_points); ++i) {
        *SVG_text << "<circle id='int" << i;
        *SVG_text << "' cx='" << segm_pairs[i].p.getX();
        *SVG_text << "' cy='" << segm_pairs[i].p.getY();
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
  uint4 N_pairs,pairs_size;
  IntersectingPair* segm_pairs;

};


using SimpleCounter=JustCountingRegistrator<TPlaneVect> ;
using PerSegmCounter=PerSegmCountingRegistrator<TPlaneVect> ;
using TrueRegistrator = PairAndPointRegistrator<TPlaneVect> ;

#endif // !REGISTRATOR_FOR_SEGMENT_INTERSECTION