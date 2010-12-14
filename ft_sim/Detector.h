/// \file Detector.h
/// Header file for the Detector class.

#include "Matrix.h"

class Combiner;
class SpectralMode;
class Source;

/// \class Detector Detector.h
class Detector
{
  public:
    Detector(int nspectral, int nreads, int ntimesteps, int nwavesteps,
             double read_noise, double quantum_efficiency, double dark_current,
             double ADC_gain);
    double read_noise;          // in e-
    double quantum_efficiency;  // in ph/e-
    double dark_current;        // in e-/s
    double ADC_gain;            // in e-/ADU

    // double instrument_visibility; //instrument losses
    // double instrument_throughput;
    int nreads;
    int ntimesteps;
    int nwavesteps;
    int nspectral;
    int nspatial;

    // storage for the interference pattern (double valued)
    Row<double> flux;         
    // storage for the actual discretized image (integer)
    Row<int> frame;
    // storage for the phase of the fringes
    Row<double> phase;

    void OneReadP(double current_time, Combiner * combiner, SpectralMode * Wav, Source * Source);
    void FastOneReadP(double current_time, Combiner * combiner, SpectralMode * Wav, Source * Source);
    void Addnoise(double integration_time);
    void Clear();

    double previous_time;
    double current_time;

    virtual ~ Detector();
    
  private:
    double _dummy;
};
