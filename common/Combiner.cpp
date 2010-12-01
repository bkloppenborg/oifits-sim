//============================================================================
// Name        : beamcombiner.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

// Beam Combiner Module

// Built-in includes
#include <cassert>
#include <cstdlib>
#include <cstdarg>

// Custom includes
#include "Combiner.h"
#include "Simulator.h"
#include "Beam.h"
#include "DelayLine.h"

Combiner::Combiner( double delay_offset , int Nbeams , Beam* beam , ... )
{
  assert(Nbeams > 1);
  this->Nbeams = Nbeams;
  flux.setsize(Nbeams);
  fringe.setsize(Nbeams, Nbeams);
  fringephase.setsize(Nbeams, Nbeams);
  OPD.setsize(Nbeams, Nbeams);
  this->beamlist.setsize(Nbeams);
  this->delay_offset = delay_offset;

  va_list beampointer;
  // First element
  va_start( beampointer , beam);
  this->beamlist[ 0 ] = beam;

  for (int ii = 1; ii < Nbeams; ii++)
  {
    this->beamlist[ ii ] = va_arg( beampointer , Beam* );
  }
  va_end( beampointer );
}

void Combiner::Update( double time , double wavenumber )
{
  // Update all beams and compute fluxes
  for (int i = 0; i < Nbeams; i++)
  {
    beamlist[ i ]->Update(time, wavenumber);
    flux[ i ] = norm2(beamlist[ i ]->pupil);
  }

  Complex phasor_atm;
  // Update atmospheric visibilities on all baselines
  // WARNING fringe matrix is not filled outside [i][j] with j<i
  for (int i = 0; i < Nbeams; i++)
  {
    for (int j = 0; j < i; j++)
    {
      // phasor_atm = total(pi.pj*) = complex visibility for atmosphere + instrument, but does not include the delays
      phasor_atm = scalconj(beamlist[ i ]->pupil, beamlist[ j ]->pupil);
      //phasor_atm = polar(1., 0.0);

      // Complex fringe visibility
      fringe[ i ][ j ] = phasor_atm * polar(1., 2. * PI * wavenumber * (beamlist[ i ]->delay - beamlist[ j ]->delay + delay_offset));

      // Compute the OPD (for the fast approximation computation)
      OPD[ i ][ j ] = arg(phasor_atm) / (2. * PI * wavenumber) + beamlist[ i ]->delay - beamlist[ j ]->delay + delay_offset; // total OPD

      // DIAGNOSTICS
      diagnostic_atm = arg(phasor_atm) / (2. * PI * wavenumber); // WARNING: potential wrapping problemm as OPD_max = lambda/2
      diagnostic_DL = beamlist[ i ]->delayline-> GetCurrentPosition(time) - beamlist[ j ]->delayline-> GetCurrentPosition(time);
  }
}

}

Combiner::~Combiner( )
{
  //delete beamlist;
}

