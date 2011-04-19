/// \file UVPoint.h
/// Header file for the UVPoint class.

#ifndef UVPOINT_H
#define UVPOINT_H

#include <string>
#include <sstream>

using namespace std;

/// \class UVPoint UVPoint.h
/// \brief A class representing a UV point
class UVPoint
{
  public:
    double u;
    double v;
    double w;                   // rarely used 
    
    UVPoint();
    void    Scale(double wavenumber);
    string  HashString(void);

    UVPoint operator+(UVPoint other);
};

#endif // UVPOINT_H
