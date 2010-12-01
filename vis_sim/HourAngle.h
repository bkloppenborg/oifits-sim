#include "Matrix.h"

#ifndef HOUR_ANGLE_H
#define HOUR_ANGLE_H

#include <list>
class Source;
class Observation;
class Array;

// class defining the hour angles of the observations
class HourAngle
{
  public:
    HourAngle(Row < double >HA);
    HourAngle(double HAstart, double HAend, double HAstep);
    HourAngle(const char *HourAngle_file);
    HourAngle(Source & target, Array & array, std::list<Observation, std::allocator<Observation> >& observations);
    
    int Nha;
    Row<double> HA;
};

#endif
