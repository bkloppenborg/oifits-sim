/// \file SpatialFilter.h
/// Header file for the Spatial Filter class

/// \class Spatialfilter
#include "Matrix.h"

Matrix<double> CircularMask(double diameter, int screensize);

class SpatialFilter
{
  public:
    SpatialFilter(int type, double pixellation, double aperture_number, double pinhole_diameter);
    virtual ~ SpatialFilter();

  public:
    /// \brief A flag to indicate the type of filtering
    /// \details
    /// 0: No spatial filtering,
    /// 1: simple pinhole,
    /// 2: FFT pinhole,
    /// 3: Fiber
    int filteringtype;          
    double pixellation;         // beam pixellation before and after filtering

    double aperture_number;     // f/D ratio before filtering

    double pinhole_diameter;

    // double beamwidth;//in meters
    Matrix < double >mask;      // unity disk function

    Matrix < Complex > pinhole; // stores the pinhole
    void ComputePinhole(int size, int padded_size);
};
