#pragma once
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
// NEW IMPLEMENTATION

#include "utils.h"
template <class ipoint>
class JustCountingRegistrator
{
  class CIntersectionBuilder {
  public:
    inline CIntersectionBuilder* register_segments(uint4 s1, uint4 s2)
    {
#ifdef PRINT_SEG_AND_INT
      if (s1<s2)
        printf("alt int %i %i\n", s1, s2);
      else
        printf("alt int %i %i\n", s2, s1);

#endif 
      return this;
    };
    inline CIntersectionBuilder* register_point(ipoint &) { return this; };
    inline void end_registration() {};

  };
public:
  static const _RegistrationType reg_type = count;
  double counter = 0;
  inline CIntersectionBuilder* begin_registration(uint4 inters_numb) { counter += inters_numb; return nullptr; };


};

typedef JustCountingRegistrator<TPlaneVect> SimpleCounter;


