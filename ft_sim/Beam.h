/// \file Beam.h
/// Header file for the Beam class

#include "Matrix.h"

class SpatialFilter;
class Modulator;
class Array;
class ZernikeLib;
class DelayLine;

/// \class Beam
class Beam
{
  public:
    double diameter;            // in meters

    int gridsize;               // in pixels

    double previous_time;

    Matrix < double >amplitude;

    Matrix < double >phase;     // the high order phase = phase - piston

    double piston;              // the piston part of the phase

    double delay;               // the total mechanical delay (modulators +
    // delaylines)
    Matrix < Complex > pupil;   // the complex pupil storing A*exp(i*phi)
    Matrix < double >mask;      // unity disk function

    Array *array;               // pointer to the relevant array

    int station_index;          // index of the corresponding station

    DelayLine *delayline;       // pointer to the corresponding delay line

    Modulator *modulator;       // pointer to the corresponding modulator

    SpatialFilter *spatialfilter;

    Beam(Array * array, int station_index, DelayLine * delayline,
         Modulator * modulator, SpatialFilter * spatialfilter, int gridsize,
         int AOorder, ZernikeLib * zernlib);
    void Update(double time, double wavenumber);

    int AOorder;                // removes norder Zernike modes, 0 = no AO

    ZernikeLib *zernlib;

    virtual ~ Beam();
  private:
    void ComputePhase(double time);

    void ComputeAmplitude(double time);

    void ComputePhasor(double wavenumber);

    void FilterBeam(double wavenumber);
};
