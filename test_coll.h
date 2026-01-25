#ifndef TEST_COLL_H
#define TEST_COLL_H
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

#include "utils.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include "math.h"
#include <cstddef>
#include <ostream>

struct SegmentsAndRegOptions
{
  int4 n = 10000;
  int4 seg_type = 2;
  int4 distr_type = 0;
  REAL distr_param = 1.0;
  int4 range_for_int_seg = full_int_range;
  uint4 reg_stat = 2;
};

class SegCollection
{
public:
  SegCollection() = default;
  SegCollection(int4 seg_type, PSeg segs) noexcept;

  SegCollection(const SegCollection&) = delete;
  SegCollection& operator=(const SegCollection&) = delete;

  SegCollection(SegCollection&& other) noexcept;
  SegCollection& operator=(SegCollection&& other) noexcept;

  ~SegCollection();

  int4 seg_type() const noexcept { return seg_type_; }
  PSeg get() const noexcept { return segs_; }
  explicit operator bool() const noexcept { return segs_ != nullptr; }
  operator PSeg() const noexcept { return segs_; }

  void reset() noexcept;

private:
  int4 seg_type_ = -1;
  PSeg segs_ = nullptr;
};

SegCollection create_test_collection(const SegmentsAndRegOptions& opt, CRandomValueGen& random_gen);

using find_intersections_func = uint8(*)(const SegmentsAndRegOptions& opt, PSeg segs, int4 alg);
find_intersections_func get_find_intersections_func(uint4 reg_type);

void write_SVG(std::ostream* svg_stream, const SegmentsAndRegOptions& opt, PSeg segs, int4 algs);

//******************************************************************************************************

#endif