/// \file Baseline.h
/// Header file for the Baseline class.

// Forward class declarations:
class UVPoint;
class Array;

/// \class Baseline simulator.h
/// \brief A class representing a baseline.
class Baseline
{
  public:
    double x;                   // east

    double y;                   // north
    
    double z;                   // up

    Baseline(Array & array, int index_station1, int index_station2);

    double Geometric_OPD(double hour_angle, double source_declination,
                         double array_latitude);
    UVPoint UVcoords(double time, double source_declination,
                     double array_latitude, double wavenumber);
};
