/// \file FTPedretti.h

#include "Matrix.h"

class SpectralMode;
class Modulator;

class FTPedretti
{
  public:
    void GDEstimate();

    FTPedretti(int n_coh, int ndelta_coh, int n_incoh, SpectralMode * spectr,
               Modulator * modulator);
    void StoreOneFrame(Row < int >&detectorframe, double time);

    Matrix < double >frames;

  private:
    int framecounter;

    int nchannels;              // number of channels

    int n_coh;                  // number of frames for coherent integration

    int ndelta_coh;             // number of frames between two sucessive
    // coherent integrations
    int n_incoh;                // number of frames for incoherent integration

    double previous_time;       // used in frame acquisition routine

    SpectralMode *spectr;

    Modulator *modulator;

    Row < double >modulator_position;   // [ framenumber ]

    Row < double >deltamod;     // [ framenumber ]

    Row < double >frametime;    // [ framenumber ]

    Row < double >frameduration;

    // Specific to Pedretti's algorithm
    Row < Complex > ntilde;
    Row < Complex > X;
    Complex X_avg;
    double delta_m12;

    void DFT();

    void ComputeCrossSpectra();

    void AverageSpectrum();

};
