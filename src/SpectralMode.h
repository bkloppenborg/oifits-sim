/// \file SpectralMode.h
/// Header file for the Spectral Mode class

#ifndef SPECTRALMODE_H_
#define SPECTRALMODE_H_

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
	SpectralMode();

    SpectralMode(double minvalue, double maxvalue, int spectralchannels);
    SpectralMode(Row < double >wav);

    void ImportFile(string filename, string comment_chars);

  public:
    double minnumber;

    double maxnumber;

    int nchannels;
    int type;
    string combiner; // The name of the instrument
    string spec_mode; // the name of the spectral mode.

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

  public:    
    void ImportFile(string filename, string combiner, string comment_chars);
    oi_wavelength GetOIWavelength(void);
    
    vector<double> GetWavelengths(void);
    vector<double> GetWavenumbers(void);
};

#endif /* SPECTRALMODE_H_ */
