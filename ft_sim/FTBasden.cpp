//============================================================================
// Name        : FTBasden.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

// This module implements the group-delay tracking algorithm described in:
// Improvements for group delay fringe tracking, A.G. Basden and D.F. Buscher, Mnras, vol 357, pp 656-668

#include <iostream>

#include "FTBasden.h"
#include "Simulator.h"
#include "SpectralMode.h"
#include "Modulator.h"

using std::cout;

// Uses the matrix class and slow DFTs. Processes n_coh images at once then flush them
// Only one baseline
// inline will use FFTW and double* pointers


double FTBasden::W1( int pixel )
{
  double value;
  switch (w1type)
  {
    case 0:
    {
      value = 1.;
      break;
    }

    case 1:
    {
      // Welch window on [ 0 ... n_coh - 1 ]
      value = 1. - pow(double(2 * pixel - n_coh + 1) / double(n_coh - 1), 2);
      break;
    }
    default:
    {
      value = 1.;
      break;
    }
  }

  return value;

}

double FTBasden::W2( int pixel )
{
  double value;
  switch (w2type)
  {
    case 0:
    {
      value = 1.;
      break;
    }

    case 1:
    {
      // Welch window on [ 0 ... nchannels - 1 ]
      value = 1. - pow(double(2 * pixel - nchannels + 1) / double(nchannels - 1), 2);
      break;
    }
    default:
    {
      value = 1.;
      break;
    }
  }

  return value;
}

FTBasden::FTBasden( int n_coh , int ndelta_coh , int n_incoh , SpectralMode* spectr , Modulator* modulator , int w1type , int w2type )
{
  // Init classic variables
  this->nchannels = spectr->nchannels;
  this->spectr = spectr;
  this->modulator = modulator;
  this->n_coh = n_coh; // number of frames in a coherent integration
  this->ndelta_coh = ndelta_coh;
  this->n_incoh = n_incoh;
  this->a = 1. - exp(-(double) ndelta_coh / (double) n_incoh);
  this->w1type = w1type;
  this->w2type = w2type;

  // Init acquisition variables
  framecounter = 0;
  previous_time = 0.0; // assumes simulation starts a t(0)=0
  meantime.setsize(n_coh);
  frameduration.setsize(n_coh);
  modulator_position.setsize(n_coh);
  deltamod.setsize(n_coh);
  frames.setsize(nchannels, n_coh);

  // Specific to Basden's algo
  F1.setsize(nchannels);

  ntrials = 4000;
  scalefactor = 1e-3;

  F2.setsize(ntrials);
  F3.setsize(ntrials);
  gd.setsize(ntrials);

  F3 = 0.;

  for (int p = 0; p < ntrials; p++)
    gd[ p ] = double(p - ntrials / 2) * scalefactor / (spectr->maxnumber - spectr->minnumber);
  cout << "FT step resolution: " << scalefactor / (spectr->maxnumber - spectr->minnumber) * 1e6 << " microns." << "   Range: +/-" << .5
      * (gd[ ntrials - 1 ] - gd[ 0 ]) * 1e6 << " microns." << "\n";
  cout << "Coherent integration uses: " << n_coh << " frames" << " with " << ndelta_coh << " frames between successive integrations."
      << "\n";
  cout << "Incoherent integration uses: " << n_incoh << " frames (decay constant: " << a << " )." << "\n";
  cout << "Maximum tolerable error: " << 1. / (spectr->maxnumber - spectr->minnumber) * 1e6 << "\n";
}

void FTBasden::StoreOneFrame( Row<int>& detectorframe , double time )
{
  // acquire fringes, detect potential problems (bad frames), convert to double precision

  for (int channel = 0; channel < nchannels; channel++)
    frames[ channel ][ framecounter ] = (double) (detectorframe[ channel ]);

  // Storing data about the current frame
  frameduration[ framecounter ] = time - previous_time;
  deltamod[ framecounter ] = modulator->GetDelay(time) - modulator->GetDelay(previous_time); // path done by modulator during the frame integration
  meantime[ framecounter ] = (time + previous_time) / 2.; // mean time in this frame
  modulator_position[ framecounter ] = modulator->GetDelay(meantime[ framecounter ]);

  previous_time = time;

  if (framecounter == (n_coh - 1))
  {
    framecounter = 0; // reset frame counter for the next batch of coherent integration
    cout << "Estimation begins" << endl;
    GDEstimate(); // launch group delay estimation on n_coh frames
    cout << "Estimation done" << endl;
  }
  else
    framecounter++;

}

void FTBasden::ComputeF1( )
{
  double local_derivative = 0.;

  for (int channel = 0; channel < nchannels; channel++)
  {
    F1[ channel ] = ZEROCOMP;

    for (int framenumber = 0; framenumber < n_coh; framenumber++)
    {
      // Can be improved by taking the real derivative
      local_derivative = frameduration[ framenumber ] / deltamod[ framenumber ];

      F1[ channel ] += W1(framenumber) * exp(-2. * PI * I * spectr->mean_wavenumber[ channel ] * modulator_position[ framenumber ])
          * local_derivative * frames[ channel ][ framenumber ];

    }

  }

}

void FTBasden::ComputeF2( )
{
  double local_derivative, chromatic_term;

  for (int gd_index = 0; gd_index < ntrials; gd_index++)
  {

    F2[ gd_index ] = ZEROCOMP;

    for (int channel = 0; channel < nchannels; channel++)
    {
      // Pixel/wavelength relation
      local_derivative = 1. / spectr->delta_wavenumber[ channel ];
      chromatic_term = 1.; // to correct for non linear effects (none yet)

      F2[ gd_index ] += W2(channel) * exp(-2. * PI * I * spectr->mean_wavenumber[ channel ] * gd[ gd_index ]) * chromatic_term
          * local_derivative * F1[ channel ];

    }
  }

}

void FTBasden::ComputeF3( )
{

  for (int gd_index = 0; gd_index < ntrials; gd_index++)
  {
    F3[ gd_index ] *= 1. - a;
    F3[ gd_index ] += a * norm(F2[ gd_index ]);
  }

}

void FTBasden::GDEstimate( ) // get the peak of the squared sum
{

  ComputeF1();

  ComputeF2();

  ComputeF3();

  //ofstream file;
  //file.open( "frames.txt" );
  //file << frames;
  //file.close();
  //getchar();
  //file.open( "mod.txt");
  //file << modulator_position ;
  //file.close();

  //file.open( "F3.txt" );
  //file << F3;
  //file.close();
  //file.open( "gd.txt" );
  //file << gd;
  //file.close();
  //getchar();

  double pow;
  double powmax = 0.;
  int gd_estim_index = 0;

  for (int gd_index = 0; gd_index < ntrials; gd_index++)
  {
    pow = F3[ gd_index ];
    if (pow > powmax)
    {
      powmax = pow;
      gd_estim_index = gd_index;
    }

  }

  // Determine centroid
  double centroid = 0., sum = 0., gd_centroid;
  //file << F3;
  //file.close();
  //getchar();

  if ((gd_estim_index > 10) && (gd_estim_index < ntrials - 10))
  {
    // Centroid method possible
    for (int ii = 0; ii <= 20; ii++)
    {
      sum += F3[ gd_estim_index - 10 + ii ];
      centroid += double(gd_estim_index - 10 + ii) * F3[ gd_estim_index - 10 + ii ];
    }
    centroid /= sum;
    // Interpolate
    gd_centroid = (1 - centroid + floor(centroid)) * gd[ int(floor(centroid)) ] + (centroid - floor(centroid)) * gd[ int(floor(centroid))
        + 1 ];
  }
  else
    gd_centroid = gd[ gd_estim_index ];

  cout << "Estimated Group delay: " << gd_centroid << " with maximum power: " << powmax << " at index= " << int(floor(centroid)) << "\n";

  OPD_estimate = gd_centroid;
}
