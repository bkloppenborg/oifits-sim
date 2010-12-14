/// \file FTBasden.h

#include "Matrix.h"

class Modulator;
class SpectralMode;

class FTBasden                  // simple fringe tracking algorithm - inline
                                // version
{
  public:
    void GDEstimate();

    FTBasden(int n_coh, int ndelta_coh, int n_incoh, SpectralMode * spectr,
             Modulator * modulator, int w1type, int w2type);
    void StoreOneFrame(Row < int >&detectorframe, double time);

    Matrix < double >frames;

    double OPD_estimate;

  private:
    int framecounter;

    int nchannels;              // number of channels

    int ntrials;                // number of trials of the group delay

    double scalefactor;         // scaling factor for the group delay list

    int n_coh;                  // number of frames for coherent integration

    int ndelta_coh;             // number of frames between two sucessive
    // coherent integrations
    int n_incoh;                // number of frames for incoherent integration

    double a;                   // constant for the autoregressive filter

    double previous_time;       // used in frame acquisition routine

    SpectralMode *spectr;

    Modulator *modulator;

    Row < double >modulator_position;   // [ framenumber ]

    Row < double >deltamod;     // [ framenumber ]

    Row < double >meantime;     // [ framenumber ]

    Row < double >frameduration;

    // Row<double> wavenumber; // [ channels ], from spectralmode class
    Row < Complex > F1;         // [ channels ] complex fringe amplitude at
    // given wavenumber
    Row < Complex > F2;         // [ gd_index ] total fringe amplitude as a
    // function of the group delay
    Row < double >F3;           // [ gd_index ] Autoregressive filter of the
    // sq modulus of F2
    Row < double >gd;           // list of the candidate group delays


    double W1(int pixel);

    double W2(int pixel);

    int w1type;

    int w2type;

    void ComputeF1();

    void ComputeF2();

    void ComputeF3();
};

