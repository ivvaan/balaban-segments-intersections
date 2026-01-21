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

PSeg create_test_collection(int4 seg_type,int4 n,int4 distr,REAL par, CRandomValueGen &random_gen, PSeg **seg_ptr_coll_ptr=nullptr);
void  delete_test_collection(int4 seg_type,PSeg,PSeg* );
double find_intersections(int4 seg_type, int4 SN, PSeg* colls, int4 alg, double* counters, bool dont_need_ip = false);

typedef uint8(*find_intersections_func)(int4, int4, PSeg, int4, uint4,int4);
find_intersections_func get_find_intersections_func(uint4 reg_type);

void write_SVG(chostream* SVG_stream, int4 seg_type, int4 n, PSeg segs, int4 algs, uint4 stat);

//******************************************************************************************************

#endif