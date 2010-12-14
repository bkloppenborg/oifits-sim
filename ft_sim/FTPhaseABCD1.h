/// \file FTPhaseAVCD1.h

#include "Matrix.h"

class SpectralMode;
class Modulator;

class FTPhaseABCD1
{
  public:
    FTPhaseABCD1(int nbins, int frames_per_bin, SpectralMode * spectr,
                 Modulator * modulator);
    virtual ~ FTPhaseABCD1();
  public:

    void StoreOneFrame(Row < int >&detectorframe, double time);

    void OPD_Estimate();

    void OPD_Search();

    double Unwrap(double wrapped_phase, int channel);

    Matrix < double > frames;

    Row < double >previous_phases;

    double opd_estimate;

    double delayline_command;

    // Locking/searching variables
    int current_mode;           // 0: search mode, 1: locking , 2: locked

    double last_locked_position;        // Last working OPD, to which we
    // revert in case of problem
    Row < double >SNR;          // list of the latest SNR values (actually
    // squared SNR)
    int SNR_index;              // index of the current SNR point - wraps
    // after SNR_samples samples.
    int SNR_samples;            // number of SNR samples in the average SNR

    double timethreshold_locking;       // wait this time in locking mode to
    // average the SNR
    double locking_time;        // time at the start of the latest locking
    // attempt
    double SNR_avg;             // current average of the SNR over a number of 
                                // 
    // 
    // samples = SNR_samples
    double snrthreshold_search; // SNR under which we revert to search mode

    double snrthreshold_locking;        // SNR over which we attempt locking
    // mode
    double snrthreshold_locked; // SNR over which we switch to locked mode,
    // after a while in locking mode
    int search_direction;

    double search_OPD_max;

    double search_OPD_steplength;

  private:
    int nbins;

    int frames_per_bin;

    int framecounter;

    int nchannels;              // number of channels

    double previous_time;       // used in frame acquisition routine

    double previous_opd;        // used for unwrapping

    SpectralMode *spectr;

    Modulator *modulator;

    Row < double >modulator_position;   // [ framenumber ]

    Row < double >deltamod;     // [ framenumber ]

    Row < double >frametime;    // [ framenumber ]

    Row < double >frameduration;

};

