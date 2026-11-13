#pragma once
/*
*
*      Copyright (c)  2011-2026  Ivan Balaban
*      ivvaan@gmail.com
*
This file is part of Seg_int library.

SPDX-License-Identifier: BSL-1.0
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
https://www.boost.org/LICENSE_1_0.txt)
*/

#include "utils.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include "math.h"
#include <cstddef>
#include <ostream>

constexpr uint4 n_threads_default = 8;

struct SegmentsAndRegOptions:CollectionOptions
{

  uint4 seg_type = 2;
  uint4 distr_type = 0;
  REAL distr_param = 1.0;
  uint4 reg_stat = 2;
  uint4 n_threads = n_threads_default;
};

class SegCollection
{
public:
  SegCollection() = default;
  SegCollection(uint4 seg_type, PSeg segs) noexcept;

  SegCollection(const SegCollection&) = delete;
  SegCollection& operator=(const SegCollection&) = delete;

  SegCollection(SegCollection&& other) noexcept;
  SegCollection& operator=(SegCollection&& other) noexcept;

  ~SegCollection();

  uint4 seg_type() const noexcept { return seg_type_; }
  PSeg get() const noexcept { return segs_; }
  explicit operator bool() const noexcept { return segs_ != nullptr; }
  operator PSeg() const noexcept { return segs_; }

  void reset() noexcept;

private:
  uint4 seg_type_ = 0;
  PSeg segs_ = nullptr;
};

SegCollection create_test_collection(const SegmentsAndRegOptions& opt, CRandomValueGen& random_gen);

using find_intersections_func = uint8(*)(const SegmentsAndRegOptions& opt, PSeg segs, uint4 alg);
find_intersections_func get_find_intersections_func(uint4 reg_type);

void write_SVG(std::ostream* svg_stream, const SegmentsAndRegOptions& opt, PSeg segs, uint4 algs);

//******************************************************************************************************
