/// \file SpectralMode.h
/// Header file for the Spectral Mode class

#include "Matrix.h"

// Header files for other libraries
extern "C" {
    #include "exchange.h"
}


/// \class SpectralMode simulator.h
/// \brief A class to represent the spectral mode(s) of the instrument.
class SpectralMode
{
  public:
    SpectralMode(string filename, string comment_chars);

    SpectralMode(double minvalue, double maxvalue, int spectralchannels);
    SpectralMode(Row < double >wav);


  public:
    double minnumber;

    double maxnumber;

    int nchannels;
    int type;
    string insname; // The name of the instrument

    Row < double >mean_wavenumber;

    Row < double >delta_wavenumber;

    Row < double >mean_wavelength;

    Row < double >delta_wavelength;

    // AS 2010-06-24
    // added a median wavelength to class SpectralMode on which the
    // integration time computed by TimeInt will be based
    double median_wavelength;


    virtual ~ SpectralMode();
  private:
    Row < double >wavenumber;

    Row < double >wavelength;
    
    oi_wavelength GetOIWavelength(void);
};
