//============================================================================
// Name        : FTPhaseABCD1.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

// This modules implements an ABCD algorithm

#include <iostream>

#include "FTPhaseABCD1.h"
#include "Simulator.h"
#include "SpectralMode.h"
#include "Modulator.h"

using std::cout;

FTPhaseABCD1::FTPhaseABCD1( int nbins , int frames_per_bin , SpectralMode* spectr , Modulator* modulator )
{
  // Init dependencies
  this->spectr = spectr;
  this->modulator = modulator;
  this->nbins = nbins;
  this->frames_per_bin = frames_per_bin;
  nchannels = spectr->nchannels;

  // Init acquisition variables
  framecounter = 0;
  previous_time = 0.0; // assumes simulation starts a t(0)=0
  frametime.setsize(nbins * frames_per_bin);
  frameduration.setsize(nbins * frames_per_bin);
  modulator_position.setsize(nbins * frames_per_bin);
  deltamod.setsize(nbins * frames_per_bin);
  frames.setsize(nchannels, nbins * frames_per_bin);
  previous_phases.setsize(nchannels); // unwrapper initialization

  // initialization of variables for the search, locking and locked modes
  //current_mode = 0; // start in search mode
  current_mode = 2;
  last_locked_position = 0.0; // start from 0 OPD
  SNR_index = 0; // begin to fill SNR values
  SNR_samples = 5; // number of SNR samples in the average SNR
  SNR.setsize(SNR_samples);// storage for previous/current/future SNR points
  timethreshold_locking = 0.1; // wait this time in locking mode to average the SNR
  locking_time = 0.0; // time at the start of the latest locking attempt
  SNR_avg = 0.0; // current average of the SNR
  snrthreshold_search = 10.; // (squared) SNR under which we revert to search mode
  snrthreshold_locking = 36.; // (squared) SNR over which we attempt locking mode
  snrthreshold_locked = 20; // (squared) SNR over which we switch to locked mode, after timethreshold_locking elapsed in locking mode
  search_direction = 1;
  search_OPD_max = 20. / spectr->mean_wavenumber[ nchannels - 1 ];
  search_OPD_steplength = 1. / spectr->mean_wavenumber[ nchannels - 1 ];
}

void FTPhaseABCD1::StoreOneFrame( Row<int>& detectorframe , double time )
{
  // acquire fringes, detect potential problems (bad frames), convert to double precision
  for (int channel = 0; channel < nchannels; channel++)
    frames[ channel ][ framecounter ] = double(detectorframe[ channel ]);
  // Storing data about the current frame
  frametime[ framecounter ] = (time + previous_time) / 2.; // mean time in this frame
  frameduration[ framecounter ] = time - previous_time;
  modulator_position[ framecounter ] = modulator->GetDelay(frametime[ framecounter ]);
  deltamod[ framecounter ] = modulator_position[ framecounter ] - modulator->GetDelay(previous_time);
  previous_time = time;

  if (framecounter == (nbins * frames_per_bin - 1))
  {
    framecounter = 0; // reset frame counter for the next batch of ABCD1

    OPD_Estimate(); // launch OPD estimation

    // Now that the OPD and SNR are estimated, we can decide on the strategy:
    // search mode 0
    // locking mode 1
    // locked mode 2

    switch (current_mode)
    {
      case 0: // algorithm currently searching for fringes
      {
        if (SNR[ SNR_index ] > snrthreshold_locking) // latest SNR sufficient to attempt locking
        {
          //reinitialize search engine
          search_OPD_max = 20. / spectr->mean_wavenumber[ nchannels - 1 ];
          search_direction = 1;

          // change into locking mode
          cout << "*************************************************************** Going into locking mode\n";
          current_mode = 1;
          locking_time = time;
          delayline_command = opd_estimate;

        }
        else // no luck yet, continue search
        {
          OPD_Search();
        }

        break;
      }

      case 1: // algorithm is currently locking
      {
        if ((locking_time - time) > timethreshold_locking) // enough time elapsed in locking time
        {
          if (SNR_avg > snrthreshold_locked) // average SNR is now sufficient to lock fully
          {
            delayline_command = opd_estimate;
            current_mode = 2;
            cout << "*************************************************************** Entering locked mode\n";
          }
          else // SNR still too low, locking failed, reverting to search
          {
            cout << "*************************************************************** SNR insufficient, reverting to search mode\n";
            current_mode = 0;
          }
        }
        else // not enough time elapsed, we stay in locking mode
        {
          delayline_command = opd_estimate;
        }

        break;
      }

      case 2: // algorithm is locked
      {
        if (SNR_avg > snrthreshold_search) // algorithm has not lost lock
        {
          last_locked_position = opd_estimate;
          delayline_command = opd_estimate;
        }
        else // loss of lock, revert to search
        {
          cout << "*************************************************************** Loss of lock, reverting to search mode\n";
          current_mode = 0;

        }

        break;
      }

    }

  }
  else
    framecounter++;
}

void FTPhaseABCD1::OPD_Search( )
{
  cout<<"*************************************************************** Scanning...\n";
  // Search the true OPD around the last locked position
  delayline_command = last_locked_position + search_direction * search_OPD_steplength;
  if ( abs( delayline_command ) > abs(search_OPD_max) ) // if the search exceeds the limit
  {
    search_OPD_max *= 2.; // increase search range
    search_direction = -search_direction; // and change direction
  }

}

void FTPhaseABCD1::OPD_Estimate( )
{
  // Note: this is for 4 bins
  double A, B, C, D, opd_monochannel, phi_monochannel;
  double X, Y, R, E;
  double N, Vsq;
  opd_estimate = 0.;

  for (int ichannel = 0; ichannel < nchannels; ichannel++)
  {
    R = PI / 2. * spectr->mean_wavenumber[ ichannel ] / spectr->mean_wavenumber[ 0 ];
    E = sin(R) / (cos(R) + 1.); // Note: R = PI does not happen inside J, H, K
    if (modulator_position[ 3 ] < modulator_position[ 0 ]) // Descending slope
      E = -E;

    A = B = C = D = 0.;
    for (int ii = 0; ii < frames_per_bin; ii++)
    {
      A += frames[ ichannel ][ ii ];
      B += frames[ ichannel ][ frames_per_bin + ii ];
      C += frames[ ichannel ][ 2 * frames_per_bin + ii ];
      D += frames[ ichannel ][ 3 * frames_per_bin + ii ];
    }

    X = A - C; // Note: this is not proportional to cos( phi )
    Y = B - D;
    if (X != Y)
      phi_monochannel = atan(E * (Y + X) / (Y - X));
    else
      phi_monochannel = PI / 2.;

    N = A + B + C + D;
    Vsq = PI * PI / 2. * (X * X + Y * Y) / (N * N);

    SNR[ SNR_index ] = 2. * (X * X + Y * Y) / N;
    SNR_avg = mean(SNR);

    // cout << "Waveband: " << spectr->mean_wavenumber[ ichannel ] << " SNR: " << SNR[ SNR_index ] << " SNR_avg: " << SNR_avg;

    SNR_index++;
    if (SNR_index > (SNR_samples - 1))
      SNR_index = 0;

    phi_monochannel = Unwrap(phi_monochannel, ichannel);
    opd_monochannel = phi_monochannel / (2. * PI * spectr->mean_wavenumber[ ichannel ]);

//    cout << " N: " << N << " V2: " << sqrt(Vsq) << " OPD: " << opd_monochannel << "\n";

    opd_estimate += opd_monochannel / double(nchannels);
  }

  // Unwrapper/filter using previous estimates - insert Kalman here
  //opd_estimate =  Unwrap( opd_estimate );

  cout << " OPD estimate " << opd_estimate << " and other possibilities: " << opd_estimate - 0.5 / spectr->mean_wavenumber[ 0 ] << " "
      << opd_estimate + 0.5 / spectr->mean_wavenumber[ 0 ] << "\n";
}

double FTPhaseABCD1::Unwrap( double wrapped_phase , int channel )
{
  // Simple unwrapper for the moment
  if (int(2. * (wrapped_phase - previous_phases[ channel ]) / PI ) > 0)
  {
    cout << "************************************* UNWRAPPED ************************\n";
    cin.get();
  }
  double unwrapped_phase = wrapped_phase + 2. * PI * double(int(2. * (wrapped_phase - previous_phases[ channel ]) / PI ));
  // Update phase info
  previous_phases[ channel ] = wrapped_phase;
  return unwrapped_phase;
}

FTPhaseABCD1::~FTPhaseABCD1( )
{
  // NTD
}

