//============================================================================
// Name        : modulator.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

// The modulator module simulates the behavior of the piezo modulators
// (ramping, ringing effects). The modulator position changes several
// times during the acquisition of a data array.

// Modes:
// - Instantaneous position update
// - Second-order filter in a later phase
// One modulator instance per beam

#include "Modulator.h"

#include <cmath>

Modulator::Modulator( int modulationtype , double period , double amplitude , double dephasing )
{
  //cout <<"Instantiating modulator"<<endl;
  //cout <<" Amplitude "<< amplitude<<endl;
  this->period = period;
  this->dephasing = dephasing;
  this->amplitude = amplitude;
  this->modulationtype = modulationtype;
}

double Modulator::Square( double period , double amplitude , double dephasing , double time )
{
  double delay = (double) ((int) floor(2. * (time - dephasing) / period) % 2) * amplitude - amplitude / 2.;
  return delay;
}

double Modulator::Ramp( double amplitude , double dephasing , double time )
{
  double delay = (time - dephasing) * amplitude;
  return delay;
}

double Modulator::Triangle( double period , double amplitude , double dephasing , double time )
{
  // get time modulo (period/2)
  double t = (time - dephasing) / period - floor((time - dephasing) / period);
  double delay = fabs(t - .5) * amplitude * 2. - amplitude / 2.;
  return delay;
}

double Modulator::Sawtooth( double period , double amplitude , double dephasing , double time )
{
  // get time modulo (period/2)
  double t = (time - dephasing) / period - floor((time - dephasing) / period);
  double delay = (t - .5) * amplitude;
  return delay;
}

double Modulator::NoModulation( double time )
{
  //No modulation : delay = 0.0 + noise
  double delay = 0.0 + 0.0 * time;
  return delay;
}

double Modulator::GetDelay( double time )
{
  double delay = 0.0;
  switch (this->modulationtype)
  {
    case 1:
    { // Square modulation
      delay = Square(this->period, this->amplitude, this->dephasing, time);
      break;
    }

    case 2:
    { // Triangle modulation
      delay = Triangle(this->period, this->amplitude, this->dephasing, time);
      break;
    }
    case 3:
    { // Ramp
      delay = Ramp(this->amplitude, this->dephasing, time);
      break;
    }

    case 4:
    { // Sawtooth
      delay = Sawtooth(this->period, this->amplitude, this->dephasing, time);
      break;
    }
    default:
    { // No modulation
      delay = NoModulation(time);
    }
  }

  return delay;
}

double Modulator::GetDifferential( double time )
{
  double diff = 0.0;
  switch (this->modulationtype)
  {
    case 1:
    { // Square modulation
      diff = 0; // check if this is ok
      break;
    }

    case 2:
    { // Triangle modulation
      double t = (time - dephasing) / period - floor((time - dephasing) / period);
      if (t >= .5)
        diff = amplitude * 2.;
      else
        diff = -amplitude * 2.;
      break;
    }
    case 3:
    { // Ramp
      diff = this->amplitude;
      break;
    }

    case 4:
    { // Sawtooth
      diff = this->amplitude; // warning - discontinuities should be treated apart
      break;
    }

    default:
    { // No modulation
      diff = 0;
    }

  }

  return diff;
}

