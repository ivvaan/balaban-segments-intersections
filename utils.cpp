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

#include <random>
#include <intrin.h>

#ifdef PRINT_SEG_AND_INT
bool print_at_lineseg1_init =  true;
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
  SetSeed(seed);
};

//zero seed means use of random gen
//nonzero - use of pseudo random setting seed
void CRandomValueGen::SetSeed(unsigned seed)
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

int8 operator%(const TIntegerVect& v1, const TIntegerVect& v2) {
  return ((int8)v1.x) * ((int8)v2.y) - ((int8)v2.x) * ((int8)v1.y);
}

int8 operator*(const TIntegerVect& v1, const TIntegerVect& v2) {
  return((int8)v1.x) * ((int8)v2.x) + ((int8)v2.y) * ((int8)v1.y);
}

TIntegerVect operator+(const TIntegerVect& v1, const TIntegerVect& v2) {
  return { v1.x + v2.x,v1.y + v2.y };
}

TIntegerVect operator-(const TIntegerVect& v1, const TIntegerVect& v2) {
  return { v1.x - v2.x,v1.y - v2.y };
}

TIntegerVect operator*(int4 r, const TIntegerVect& t) {
  return { t.x * r, t.y * r };
}

TIntegerVect operator*(const TIntegerVect& t, int4 r) {
  return { t.x * r, t.y * r };
}

std::strong_ordering comp_prod(int64_t a, int64_t b, int64_t c, int64_t d) {
  int64_t hi1, hi2;
  uint64_t lo1 = _mul128(a, b, &hi1); // product1 = (hi1, lo1)
  uint64_t lo2 = _mul128(c, d, &hi2); // product2 = (hi2, lo2)

  if (hi1 != hi2)
    return hi1 <=> hi2;

  return lo1 <=> lo2;
}

// some string
