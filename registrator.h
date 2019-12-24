#pragma once
#include "utils.h"
template <class ipoint>
class JustCountingRegistrator
{
public:
  static const _RegistrationType reg_type = count;
  double counter = 0;
  inline void begin_registration(uint4 inters_numb) { counter += inters_numb; };
  inline void register_segments(uint4 s1, uint4 s2)
  {
#ifdef PRINT_SEG_AND_INT
    if (s1<s2)
      printf("alt int %i %i\n", s1, s2);
    else
      printf("alt int %i %i\n", s2, s1);

#endif 
  };
  inline void register_points(ipoint*) {};
  inline void end_registration() {};


};

typedef JustCountingRegistrator<TPlaneVect> SimpleCounter;


