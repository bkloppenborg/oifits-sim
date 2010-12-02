/// \file Source.h
/// Header file for the source class

#include <complex>
#include <string>

#include "Matrix.h"

// Forward Declaration of classes
class Baseline;
class UVPoint;
class Array;

/// \class Source simulator.h 
/// \brief A class to represent the source object.
// AS 2010-06-23
// added new memebers to the class: name of the source, band,
// background_magnitude, sky_background_aperture
class Source
{
  public:
    Source(const char *filename, char band, double magnitude,
           double temperature, double background_magnitude,
           double sky_background_aperture, double declination, 
           double right_ascension, double pixellation);
    Source(const char *Source_file);

    double flux;                // integrated flux over lambda (should be
    // Jansky)
    double backgroundflux;      // flux of the background sky

    std::string source_name;
    char band;

    double magnitude;

    double temperature;

    double declination;
    double right_ascension;

    double background_magnitude;

    double sky_background_aperture;

    Matrix < double >source_image;

    double source_pixellation;

    Complex GetVis(UVPoint uv);
    
    /// \todo Rewrite this function to be compliant with the new Observation class.
//    Complex GetVis(double wavenumber, double time, Array & array, int index_station1, int index_station2);
    double Spectrum(double wavenumber);

    double BackgroundSpectrum(double wavenumber);

  private:
    void InitFluxes(char band, double background_magnitude, double sky_background_aperture);
};
