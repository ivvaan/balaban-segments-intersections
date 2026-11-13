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

#include <random>


#ifdef PRINT_SEG_AND_INT
bool print_at_lineseg1_init = true;
#endif

namespace {
#ifdef NDEBUG
  std::random_device rd_gen;
  std::default_random_engine drd_gen;
#else
  std::default_random_engine rd_gen;
  std::default_random_engine drd_gen;
#endif // NDEBUG

  double stduniform1()
  {
    static std::uniform_real_distribution<double> dist{ 0.0, 1.0 };
    return dist(rd_gen);
  }

  double stduniform2()
  {
    static std::uniform_real_distribution<double> dist{ 0.0, 1.0 };
    return dist(drd_gen);
  }

}


CRandomValueGen::CRandomValueGen()
{
  stduniform = stduniform2;
};

CRandomValueGen::CRandomValueGen(unsigned seed)
{
  SetSeeed(seed);
};

//zero seed means use of random gen
//nonzero - use of pseudo random setting seed
void CRandomValueGen::SetSeeed(unsigned seed)
{
  stduniform = seed ? drd_gen.seed(seed), stduniform2 : stduniform1;
};

 

double CRandomValueGen::GetRandomDouble()
  {
  return stduniform();
  }

bool CRandomValueGen::RandomChoose()
  {
  double q,r= stduniform();
  do {q = stduniform();} while (q == r);
  return q>r;
  }


