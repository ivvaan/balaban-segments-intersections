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
#include <algorithm>
#include <tuple>
#include <vector>
#include <cstdint>
#include <cstdio>
#include "utils.h"

enum _RegistrationType
{
  count = 1,
  segments = 2,
  point = 4,
  count_and_segments = 3,
  full = 7
};

using int_rec = std::tuple<uint4, uint4>;

constexpr const uint4 cache_line_size = std::hardware_constructive_interference_size;


//registrar factory can create different types of registrators
template <class Registrar>
class CRegistratorFactory
{
  constexpr static const uint4 reg_margin = (cache_line_size + sizeof(Registrar) - 1) / sizeof(Registrar);// for reg objects to be in different CPU cache lines
  //constexpr static const uint4 initial_size = 4;
  std::vector<Registrar> reg_objects;
  std::vector<Registrar*> reg_pointers;
  uint4 nSegm = 0;

  std::unique_ptr<Registrar> main_reg_owner;

public:
  CRegistratorFactory()
    : main_reg_owner(std::make_unique<Registrar>())
  {}

  void Reset()
  {
    reg_objects.clear();
    reg_pointers.clear();
    //nSegm = 0;
  }

  void InitClone(uint4 n_threads)
  {
      reg_objects.resize((n_threads-1) * reg_margin);
      reg_pointers.resize(n_threads-1);
  };

  void PrepareAlloc(uint4 _nSegm)
  {
    nSegm = _nSegm;
    //main_reg_owner->Alloc(nSegm);
  };


  Registrar* GetRegistrator(uint4 thread_index) 
  {
    Registrar *ptr;
    if (thread_index) {
      ptr = &reg_objects.at((thread_index - 1) * reg_margin);
      reg_pointers.at(thread_index - 1) = ptr;
    }
    else
      ptr = main_reg_owner.get();
    ptr->Alloc(nSegm);
    return ptr;
  };

  void remove_null() 
  {
    reg_pointers.erase(
      std::remove(reg_pointers.begin(), reg_pointers.end(), nullptr),
      reg_pointers.end());
  };

  void combine_reg_data() {
    remove_null();
    uint4 n_threads = reg_pointers.size()+1;
    if (n_threads > 1) {
      main_reg_owner->combine_reg_data(n_threads, reg_pointers.data());
    }
  };
  uint8 get_stat(uint4 stat_type = 0) {
    return  main_reg_owner->get_stat(stat_type);
  };

  void write_SVG(uint4 alg, chostream* SVG_text) {
    main_reg_owner->write_SVG(alg, SVG_text);
  }

};

template <uint4 r_type = 0>
class JustCountingRegistrator
{
  uint8 counter = 0;
public:
  static constexpr uint4 reg_type = r_type | _RegistrationType::count;
  //static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments + _RegistrationType::point;
  void register_pair(uint4 s1, uint4 s2) noexcept {
#ifdef PRINT_SEG_AND_INT
    if (s1 > s2)
      std::swap(s1, s2);
    printf("alt int %i %i\n", s1, s2);

#endif 

    ++counter; 
  };

  template <class ipoint>
  void register_pair_and_point(uint4 s1, uint4 s2, const ipoint& p) { register_pair(s1, s2); };
 
  void combine_reg_data(uint4 n_threads, JustCountingRegistrator* additional_reg_obj[])
  {
    for (int i = 0; i < n_threads - 1; i++){
      counter += additional_reg_obj[i]->counter;
    }
  };
  void Alloc(uint4 _N) {};

  void Reset() { counter = 0; };

  void Flash() {};

  // to return some statistics about registered intersections;
  auto get_stat(uint4 stat_type = 0) { return counter; };

  void write_SVG(uint4 alg, chostream* SVG_text) {};

};

class PerSegmCountingRegistrator
{
  uint4 N=0;
  uint4* segm_counters = nullptr;
  uint8 counter = 0;

public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments;
  ~PerSegmCountingRegistrator() { if (segm_counters != nullptr) { delete[] segm_counters; segm_counters = nullptr; } };


  void register_pair(uint4 s1, uint4 s2) noexcept {
    ++counter;
    ++segm_counters[s1];
    ++segm_counters[s2];
  };

  template <class ipoint>
  void register_pair_and_point(uint4 s1, uint4 s2, const ipoint& p) { register_pair(s1, s2); };

  void Alloc(uint4 _N)
  {
    if(!segm_counters){
      N = _N;
      segm_counters = new uint4[N];
      std::fill_n(segm_counters,N,0);
    }
  };

  void Reset() 
  { 
    counter = 0; 
    if(segm_counters)
      std::fill_n(segm_counters,N,0);
  };

  void Flash() {};

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
  auto get_stat(uint4 stat_type = 0)
  {
    // if stat_type nonzero return maximal number of intersections per segment;
    if (stat_type == _Registrator::per_segm_reg_max_per_segm_stat)
      return (uint8)*std::max_element(segm_counters, segm_counters + N);
    return counter;
  };

  void write_SVG(uint4 alg, chostream* SVG_text) {};

};

class PairRegistrator
{
public:

  std::vector<int_rec>  intersections;
  uint8 counter = 0;

  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments;

  PairRegistrator() {
    intersections.reserve(1024 * 1024);
  };
  PairRegistrator(const  PairRegistrator&) = default;


  ~PairRegistrator() {
  };

  void Reset() {
    intersections.clear();
    counter = 0;
  };

  void Flash() {};

  void register_pair(uint4 s1, uint4 s2) {
    ++counter;
#ifdef NDEBUG
    if (s1 < s2)
      intersections.emplace_back(s1, s2);
    else
      intersections.emplace_back(s2, s1);
#else
    if (s1 > s2)std::swap(s1, s2);
    intersections.emplace_back(s1, s2);
#endif
  };

  void Alloc(uint4) {};

  // to return some statistics about registered intersections;
  auto get_stat(uint4 stat_type = 0)
  {
    return counter;
  };

  // Sort stored pairs in-place (lexicographic on tuple).
  void sort_inplace()
  {
    std::sort(intersections.begin(), intersections.end());
  }

  // Compute multiset difference: fill a_minus_b with elements in A\B (with multiplicities),
  // and b_minus_a with elements in B\A (with multiplicities).
  // The function does NOT modify the original registrators.
  static void multiset_difference(const PairRegistrator &A, const PairRegistrator &B,
                                  std::vector<int_rec> &a_minus_b,
                                  std::vector<int_rec> &b_minus_a)
  {
    a_minus_b.clear();
    b_minus_a.clear();

    // make copies and sort (we don't mutate originals)
    std::vector<int_rec> a = A.intersections;
    std::vector<int_rec> b = B.intersections;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    size_t i = 0, j = 0;
    const size_t na = a.size();
    const size_t nb = b.size();

    while (i < na && j < nb)
    {
      if (a[i] == b[j])
      {
        // count multiplicities for this key in both
        size_t ia = i + 1;
        while (ia < na && a[ia] == a[i]) ++ia;
        size_t jb = j + 1;
        while (jb < nb && b[jb] == b[j]) ++jb;

        size_t ca = ia - i;
        size_t cb = jb - j;
        if (ca > cb)
        {
          // append (ca - cb) copies of a[i]
          for (size_t t = 0; t < ca - cb; ++t)
            a_minus_b.push_back(a[i]);
        }
        else if (cb > ca)
        {
          for (size_t t = 0; t < cb - ca; ++t)
            b_minus_a.push_back(b[j]);
        }
        i = ia;
        j = jb;
      }
      else if (a[i] < b[j])
      {
        // all occurrences of a[i] are only in A
        size_t ia = i + 1;
        while (ia < na && a[ia] == a[i]) ++ia;
        size_t ca = ia - i;
        for (size_t t = 0; t < ca; ++t)
          a_minus_b.push_back(a[i]);
        i = ia;
      }
      else // b[j] < a[i]
      {
        size_t jb = j + 1;
        while (jb < nb && b[jb] == b[j]) ++jb;
        size_t cb = jb - j;
        for (size_t t = 0; t < cb; ++t)
          b_minus_a.push_back(b[j]);
        j = jb;
      }
    }

    // remaining in A
    while (i < na)
    {
      size_t ia = i + 1;
      while (ia < na && a[ia] == a[i]) ++ia;
      size_t ca = ia - i;
      for (size_t t = 0; t < ca; ++t)
        a_minus_b.push_back(a[i]);
      i = ia;
    }
    // remaining in B
    while (j < nb)
    {
      size_t jb = j + 1;
      while (jb < nb && b[jb] == b[j]) ++jb;
      size_t cb = jb - j;
      for (size_t t = 0; t < cb; ++t)
        b_minus_a.push_back(b[j]);
      j = jb;
    }
  }

 
};

template <class ipoint>
class PairAndPointRegistrator
{
  using int_rec = std::tuple<uint4, uint4, ipoint>;
  std::vector<int_rec>  intersections;
  uint8 counter = 0;

public:
  static constexpr uint4 reg_type = _RegistrationType::count + _RegistrationType::segments + _RegistrationType::point;

  PairAndPointRegistrator() {
    intersections.reserve(1024 * 1024);
  };

  ~PairAndPointRegistrator() {
  };
  void Reset() {
    intersections.clear();
    counter = 0;
  };
  void Flash() {};

  void register_pair(uint4 s1, uint4 s2) noexcept {
    ++counter;
  };
  void register_pair_and_point(uint4 s1, uint4 s2, const ipoint& p) {
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
  auto get_stat(uint4 stat_type = 0)  {
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
        default:
          *SVG_text << "' class='algorithm' id='a" << s1<<"_"<<s2<<"' />\n";
          break;
        }
      }
    };
  }

};


using SimpleCounter = JustCountingRegistrator<>;
using SimpleCounter2 = JustCountingRegistrator<_RegistrationType::point>;
using PerSegmCounter = PerSegmCountingRegistrator;
using TrueRegistrator = PairAndPointRegistrator<TPlaneVect> ;