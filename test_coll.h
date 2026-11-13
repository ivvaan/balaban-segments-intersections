#ifndef TEST_COLL_H
#define TEST_COLL_H
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

PSeg create_test_collection(int4 seg_type,int4 n,int4 distr,REAL par, CRandomValueGen &random_gen, PSeg **seg_ptr_coll_ptr=nullptr);
void  delete_test_collection(int4 seg_type,PSeg,PSeg* );
double find_intersections(int4 seg_type, int4 SN, PSeg* colls, int4 alg, double* counters, bool dont_need_ip = false);

typedef uint8(*find_intersections_func)(int4, int4, PSeg, int4, uint4);
find_intersections_func get_find_intersections_func(uint4 reg_type);

void write_SVG(chostream* SVG_stream, int4 seg_type, int4 n, PSeg segs, int4 algs, uint4 stat);

//******************************************************************************************************

#endif