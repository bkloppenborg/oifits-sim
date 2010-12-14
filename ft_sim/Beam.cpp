//============================================================================
// Name        : beam.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

// Beam module
// Connects stations and delay lines.

#include <iostream>

#include "Beam.h"
#include "ZernikeLib.h"
#include "DelayLine.h"
#include "Modulator.h"
#include "Array.h"
#include "Station.h"
#include "SpatialFilter.h"
#include "AtmosphereLayer.h"

using std::cout;

Beam::Beam(Array * array, int station_index, DelayLine * delayline,
         Modulator * modulator, SpatialFilter * spatialfilter, int gridsize,
         int AOorder, ZernikeLib * zernlib)
{
  // Init constants and link to Zernike modes
  this->array = array;
  this->station_index = station_index;
  this->delayline = delayline;
  this->modulator = modulator;
  this->spatialfilter = spatialfilter;
  this->gridsize = gridsize;
  this->AOorder = AOorder;
  this->zernlib = zernlib;
  // Initialize matrix sizes and the disk function for this station
  amplitude.setsize(gridsize, gridsize);
  phase.setsize(gridsize, gridsize);
  pupil.setsize(gridsize, gridsize);
  mask.setsize(gridsize, gridsize);
  for (int ii = 0; ii < gridsize; ii++)
    for (int jj = 0; jj < gridsize; jj++)
      mask[ ii ][ jj ] = zernlib->modes[ 0 ][ ii ][ jj ] / zernlib->modes[ 0 ][ gridsize / 2 ][ gridsize / 2 ];

  previous_time = -1.;
  ComputePhase(0.0); // init phase screen
  ComputeAmplitude(0.0); // init amplitude screen
}

void Beam::ComputePhase(double time)
// Get the atmosphere phase and apply the AO correction (at lambda0)
// Add the modulator and delaylines delays
{
  // Apply the mask to the phase screen
  phase = this->array->GetStation(station_index).layer->phase_screen * mask ;

  // computes the AO + tip/tilt correction
  if (AOorder > 0)
  {
    for (int ii = 1; ii <= AOorder; ii++)
    {
      phase -= zernlib->modes[ ii ] * scalprod(phase, zernlib->modes[ ii ]);
    }
  }

  // Now update the piston at lambda0
  // No high order here, but if need be the phase could be modified to introduce some
  // Note: we do not use phase += piston*zernike[0] as the OPD alone is needed in later computations

  double delay_modulator = this->modulator->GetDelay(time);
  double delay_delaylines = this->delayline->GetCurrentPosition(time);

  delay = delay_modulator + delay_delaylines;
}

void Beam::ComputeAmplitude(double time)
{
    // TODO : add strehl if wanted
    // Apply normalized mask (= zernlib->modes[0] ) to the amplitude, as well as station gain
    // Wavelength dependency could be introduced

    Station sta = this->array->GetStation(station_index);
      
    if ( sta.layer->scintillation_diameter > 0 )
        amplitude = sta.layer->amplitude_screen * zernlib->modes[ 0 ] * sta.gain;
    else
        amplitude = zernlib->modes[ 0 ] * sta.gain;

}

void Beam::ComputePhasor(double wavenumber) // compute the complex high-order pupil at a given wavelength
{
    for (int ix = 0; ix < gridsize; ix++)
        for (int iy = 0; iy < gridsize; iy++)
            pupil[ ix ][ iy ] = polar(amplitude[ ix ][ iy ], this->array->GetStation(station_index).layer->lambda0 * wavenumber * phase[ ix ][ iy ]);
}

void Beam::Update(double time, double wavenumber)
{
    if (time != previous_time)
    {
        // Update the atmosphere
        this->array->GetStation(station_index).layer->Update(time);
        ComputePhase(time);
        ComputeAmplitude(time);
        previous_time = time;
    }
    ComputePhasor(wavenumber);
    FilterBeam(wavenumber);
}


Beam::~Beam()
{
}

void Beam::FilterBeam(double wavenumber)
{
  switch (spatialfilter->filteringtype)
  {
    case 0: // No filtering at all, so do nothing
    {
      //cout << "No filtering\n";
      break;
    }

    case 1: //Pinhole matched to size of airy disk in FFT of incoming beam
    {
      //set padding
      int padded_size = 1024;
      //spatial filtering
      if (spatialfilter->pinhole.GetRows() < 1) // if pinhole does not exist yet..
        spatialfilter->ComputePinhole(gridsize, padded_size); // ...create it

      Matrix<Complex> padded_pupil(padded_size, padded_size);
      Matrix<Complex> beam_spectrum(padded_size, padded_size);
      ZeroPad(pupil, padded_pupil);
      FFT(padded_pupil, beam_spectrum, 1, 0);
      beam_spectrum *= spatialfilter->pinhole;
      FFT(beam_spectrum, padded_pupil, -1, 0);
      ZeroUnPad(pupil, padded_pupil);

      //recompute phase and amplitude (+remask ?) in case they are used in other modules
      for (int ii = 0; ii < gridsize; ii++)
      {
        for (int jj = 0; jj < gridsize; jj++)
        {
          amplitude[ ii ][ jj ] = abs(pupil[ ii ][ jj ]);
          phase[ ii ][ jj ] = arg(pupil[ ii ][ jj ]);
        }
      }

      break;
    }

    default: // just in case...
    {
      cout << "Wrong spatial filter selection !\n";
      break;
    }
  }
}
