#ifndef TEST_COLL_H
#define TEST_COLL_H
/*
*
*      Copyright (c)  2011  Ivan Balaban 
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
#include "intersection_finder.h"

enum _Algorithm
  {
  triv=1,
  simple_sweep=2,
  fast=4,
  optimal=8,
  fast_parallel=16
  };

enum _Segment
  {
  line1=0,line2,arc
  };
enum _Distribution
  {
  random=0,mixed,small,param_defined
  };


PSeg* create_test_collection(int4 seg_type,int4 n,int4 distr,REAL par);
void  delete_test_collection(int4 seg_type,PSeg* );
double  find_intersections(int4 seg_type,int4 SN,PSeg *colls,int4 alg,double *counters);

#endif