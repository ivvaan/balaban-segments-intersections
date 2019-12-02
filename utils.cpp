/*
*
*      Copyright (c)  2011-2017  Ivan Balaban 
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

#define PRIME1 33461
#define PRIME2 23833
#define PRIME3 28657
#define MUL1 1867
#define MUL2 2833
#define MUL3 7657

#include <random>
#ifdef _DEBUG
double randm()
{
    static int4 f1 = MUL1;
    static int4 f2 = MUL2;
    static int4 f3 = MUL3;
    f1 = (f1 * MUL1) % PRIME1;
    f2 = (f2 * MUL2) % PRIME2;
    f3 = (f3 * MUL3) % PRIME3;
    return ((double)f1 + (double)f2 * 200.0 + (double)f3 * 30000.0) / 100000000.0;
};
#else 
double randm()
{
    static   std::random_device generator;
    //static   std::default_random_engine generator;
    static   std::uniform_real_distribution<double> distribution(0.0, 1.0);
 return distribution(generator);
}
#endif

 

double CRandomValueGen::GetRandomDouble()
  {
  return randm();
  }

bool CRandomValueGen::RandomChoose()
  {
  double q,r=randm();
  do {q = randm();} while (q == r);
  return q>r;
  }


