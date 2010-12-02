/// \file Baseline.h
/// Header file for the Baseline class.

// Forward class declarations:
class UVPoint;
class Array;
class Station;

/// \class Baseline simulator.h
/// \brief A class representing a baseline.
class Baseline
{
  private:
    double NEU[3];  // Array for (North, East, Up)
    double xyz[3];

  public:
    Baseline(void);
    Baseline(double array_lat, Station & station1, Station & station2);

    /// \todo Rewrite this function to work with the new class definition.
    //double Geometric_OPD(double hour_angle, double source_declination);
    
    UVPoint UVcoords(double hour_angle, double source_declination, double wavenumber);
};
