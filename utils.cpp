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


