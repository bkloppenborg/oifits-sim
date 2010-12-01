//============================================================================
// Name        : detector.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

// Detector Module, simulates data frames
// Includes sampling, dispersion, noise.

#include "Detector.h"
#include "Simulator.h"
#include "Combiner.h"
#include "SpectralMode.h"
#include "Source.h"
#include "Beam.h"
#include "random.h"
#include "Array.h"
#include "Modulator.h"

Detector::Detector( int nspectral , int nreads , int ntimesteps , int nwavesteps , double read_noise , double quantum_efficiency ,
    double dark_current , double ADC_gain )
{
  this->read_noise = read_noise;
  this->quantum_efficiency = quantum_efficiency;
  this->dark_current = dark_current;
  this->nspectral = nspectral;
  this->nspatial = nspatial;
  this->ADC_gain = ADC_gain;
  this->nreads = nreads; // number of reads per integration time
  this->ntimesteps = ntimesteps; // number of time increments integrated over per integration time
  this->nwavesteps = nwavesteps; //  number of wavenumber increments integrated over per channel
  this->previous_time = 0.0;
  flux.setsize(nspectral);
  frame.setsize(nspectral);
  phase.setsize(nspectral);
}

void Detector::Clear( )
{
  for (int ii = 0; ii < flux.size(); ii++)
    flux[ ii ] = 0.;
  for (int ii = 0; ii < frame.size(); ii++)
    frame[ ii ] = 0;
}

// Detector for a pupil-plane combination scheme - fast approximation
void Detector::FastOneReadP( double current_time , Combiner* combiner , SpectralMode* Wav , Source* Source )
{
  // Note: temporal and spatial correlation transfer factors of the atmosphere to take into account
  // Check what is included in our atmosphere model and what is not !
  // temporal_coherence_atmosphere_fluctuations  =  sqrt( 1.79 * t0( lambda[j] )/delta_t[i] ); // from David's formula ( 8 ) in 1988 paper
  // spatial_coherence_atmosphere = ... see Brummelaar eq ( 44 )
  // coherence_transfer_factor =  strehl ratio, see Brummelaar, include into OA ?

  double collecting_area = PI * pow(combiner->beamlist[ 0 ]->array->diameter[0] / 2., 2); // surface of the primary mirror - can insert this into the loop if variable
  double instrument_throughput = 0.1324 / 0.65; // global throughput, does not include quantum efficiency
  //double instrumental_visibility = 0.3; // global factor for the moment, would depend on wavelength
  double instrumental_visibility = 1.0;
  double source_flux;
  double background_flux;
  double time;
  Complex source_vis; // variable storing the complex source visibilityone the same while visiting other major city's. Of all the restaurants we've done, this was the most expensive. The quality of service was impressive. It amazed me to see 5 to
  this->current_time = current_time;

  // Clear the CCD before integration
  Clear();

  double total_integration_time = current_time - previous_time;
  double timestep = total_integration_time / double(ntimesteps);
  double flux_norm = 0.5 * collecting_area * instrument_throughput * timestep; // factor 0.5 <- using only one output of the combiner

  for (int itime = 0; itime < ntimesteps; itime++)
  {
    time = previous_time + double(itime + 1) * total_integration_time / double(ntimesteps + 1);

    for (int iwavenumber = 0; iwavenumber < Wav->nchannels; iwavenumber++) // loop on spectral channels / pixels
    {
      // Update combiner for current time and/or wavelength (this updates all beams)
      combiner->Update(time, Wav->mean_wavenumber[ iwavenumber ]);

      // Compute the flux normalisation due to the source, in the current waveband
      source_flux = flux_norm * Source->Spectrum(Wav->mean_wavenumber[ iwavenumber ]) * Wav->delta_wavenumber[ iwavenumber ];
      background_flux = flux_norm * Source->BackgroundSpectrum(Wav->mean_wavenumber[ iwavenumber ]) * Wav->delta_wavenumber[ iwavenumber ];

      // Add background flux
      flux[ iwavenumber ] += background_flux;

      // Add beam intensities < p1^2 > + < p2^2 > + ... + < pn^2 >
      // TODO: compute and use fringes[ i ][ i ] instead of combiner->flux
      for (int nbeam = 0; nbeam < combiner->Nbeams; nbeam++)
        flux[ iwavenumber ] += source_flux * combiner->flux[ nbeam ];

      // Now add the fringes

      for (int nbeam = 1; nbeam < combiner->Nbeams; nbeam++) // loop on beams
      {
        for (int mbeam = 0; mbeam < nbeam; mbeam++)
        {
          // Get source visibility for the current baseline
          source_vis = Source->GetVis(Wav->mean_wavenumber[ iwavenumber ], time, *combiner->beamlist[ 0 ]->array, combiner->beamlist[ nbeam ]->station_index,
              combiner->beamlist[ mbeam ]->station_index );

          // phase[ iwavenumber ] = combiner->fringephase[ nbeam ][ mbeam ] + arg( source_vis ) ;

          // Now we add the coherence terms <p2.p1*> + <p3.p2*> + <p3.p1*> + <pn.pn-1*> to the flux

          // Hypothesis for the temporal smearing term: linear modulation on short timescale
          // the integration of A.cos(g(t)+p) over t0-Dt/2 to t0+Dt/2 gives I = 2A/[dg(t0)/dt] * cos[ ( g(t0+Dt/2) + g(t0-Dt/2) )/2 + p ] * sin[ ( g(t0+Dt/2) - g(t0-Dt/2) )/2 ]
          // When the timescale is small, g(t) is solely determined by the modulators (atmosphere + delaylines are frozen)
          // and can be approximated by a linear function: g(t) = alpha * t + beta
          // giving I = A cos[ g(t0) + p ] *  ( Dt * sinc[ alpha * Dt / 2] )
          flux[ iwavenumber ] += source_flux // flux normalisation factor
              * instrumental_visibility // visibility losses
              * sinc(PI * combiner->OPD[ nbeam ][ mbeam ] * Wav->delta_wavenumber[ iwavenumber ]) // spectral bandpass effect (top-hat)
              * sinc(combiner->beamlist[ nbeam ]->modulator->GetDifferential(time) * timestep / 2.) // temporal smearing term
              * 2. * real(combiner->fringe[ nbeam ][ mbeam ] * source_vis);
        }// end of beam loop 1

      } // end of beam loop 2


    } // end of channel loop
  } // end of time loop

  // Update detector timer
  previous_time = current_time;

  // Now we discretize the acquisition + add noise
  Addnoise(total_integration_time);
}

// Detector for a pupil-plane combination scheme -- loop version
void Detector::OneReadP( double current_time , Combiner* combiner , SpectralMode* Wav , Source* Source )
{
  // Note: temporal and spatial correlation transfer factors of the atmosphere to take into account
  // Check what is included in our atmosphere model and what is not !
  // temporal_coherence_atmosphere_fluctuations  =  sqrt( 1.79 * t0( lambda[j] )/delta_t[i] ); // from David's formula ( 8 ) in 1988 paper
  // spatial_coherence_atmosphere = ... see Brummelaar eq ( 44 )
  // coherence_transfer_factor =  strehl ratio, see Brummelaar, include into OA ?

  double collecting_area = PI * pow(combiner->beamlist[ 0 ]->array->diameter[0] / 2., 2); // surface of the primary mirror - can insert this into the loop if variable
  double instrument_throughput = 0.204; // global throughput, does not include quantum efficiency
  //double instrumental_visibility = 0.3; // global factor for the moment, would depend on wavelength
  double instrumental_visibility = 1.0; // 0.3 default global factor for the moment, would depend on wavelength
  double source_flux;
  double background_flux;
  double time;

  Complex source_vis; // variable storing the complex source visibilityone the same while visiting other major city's. Of all the restaurants we've done, this was the most expensive. The quality of service was impressive. It amazed me to see 5 to
  this->current_time = current_time;

  // Clear the CCD before integration
  Clear();

  double total_integration_time = current_time - previous_time;
  double timestep = total_integration_time / double(ntimesteps);
  double flux_norm = 0.5 * collecting_area * instrument_throughput * timestep; // factor 0.5 <- using only one output of the combiner

  double wavenumber; // stores the current wavenumber
  double delta_wavenumber; // stores the current waveband width

  for (int itime = 0; itime < ntimesteps; itime++)
  {
    time = previous_time + double(itime + 1) * total_integration_time / double(ntimesteps + 1);

    for (int ichannel = 0; ichannel < Wav->nchannels; ichannel++) // loop on channels
    {

      delta_wavenumber = Wav->delta_wavenumber[ ichannel ] / double(nwavesteps);

      for (int iwavenumber = 0; iwavenumber < nwavesteps; iwavenumber++) // loop on wavenumber inside each channel
      {
        wavenumber = Wav-> mean_wavenumber[ ichannel ] - Wav->delta_wavenumber[ ichannel ] / 2 // starting wavenumber in the waveband
            + double(iwavenumber + 1) * Wav->delta_wavenumber[ ichannel ] / double(nwavesteps + 1); // increment for equispaced sampling

        // Update combiner for current time and/or wavelength (this updates all beams)
        combiner->Update(time, wavenumber);

        // Compute the flux normalisation due to the source, in the current waveband
        source_flux = flux_norm * Source->Spectrum(wavenumber) * delta_wavenumber;
        background_flux = flux_norm * Source->BackgroundSpectrum(wavenumber) * delta_wavenumber;

        // Add background flux
        flux[ ichannel ] += background_flux;

        // Add beam intensities < p1^2 > + < p2^2 > + ... + < pn^2 >
        // TODO: compute and use fringes[ i ][ i ] instead of combiner->flux
        for (int nbeam = 0; nbeam < combiner->Nbeams; nbeam++)
        {
          flux[ ichannel ] += source_flux * combiner->flux[ nbeam ];
        }
        //        cout << "Incoherent Flux " << flux[ ichannel ] << endl;
        // Now add the fringes
        for (int nbeam = 1; nbeam < combiner->Nbeams; nbeam++) // loop on beams
        {
          for (int mbeam = 0; mbeam < nbeam; mbeam++)
          {
            // Get source visibility for the current baseline
            source_vis = Source->GetVis(wavenumber, time, *combiner->beamlist[ 0 ]->array, combiner->beamlist[ nbeam ]->station_index, combiner->beamlist[ mbeam ]->station_index);

            // Now we add the coherence terms <p2.p1*> + <p3.p2*> + <p3.p1*> + <pn.pn-1*> to the flux
            // Hypothesis for the temporal smearing term: linear modulation on short timescale
            // the integration of A.cos(g(t)+p) over t0-Dt/2 to t0+Dt/2 gives I = 2A/[dg(t0)/dt] * cos[ ( g(t0+Dt/2) + g(t0-Dt/2) )/2 + p ] * sin[ ( g(t0+Dt/2) - g(t0-Dt/2) )/2 ]
            // When the timescale is small, g(t) is solely determined by the modulators (atmosphere + delaylines are frozen)
            // and can be approximated by a linear function: g(t) = alpha * t + beta
            // giving I = A cos[ g(t0) + p ] *  ( Dt * sinc[ alpha * Dt / 2] )
            flux[ ichannel ] += source_flux // flux normalisation factor
                * instrumental_visibility // visibility losses
                * sinc(combiner->beamlist[ nbeam ]->modulator->GetDifferential(time) * timestep / 2.) // temporal smearing term
                * 2. * real(combiner->fringe[ nbeam ][ mbeam ] * source_vis);
          } // end of beam loop 1
        } // end of beam loop 2
      } // end of wavenumber integration loop
    } // end of channel loop
  } // end of time loop

  // Update detector timer
  previous_time = current_time;

  // Now we discretize the acquisition + add noise
  Addnoise(total_integration_time);
}

void Detector::Addnoise( double integration_time )
{
  // Poisson noise, read noise, dark current
  for (int ix = 0; ix < nspectral; ix++)
  {
    this->frame[ ix ] = (int) trunc((Ranpoiss(random_number_seed, quantum_efficiency * flux[ ix ] + dark_current * integration_time) + read_noise
      * Rangauss(random_number_seed)) / ADC_gain);

    //Uncomment this for no noise
    //frame[ ix ] = (int) trunc((quantum_efficiency * flux[ ix ] + dark_current * integration_time) / ADC_gain);
  }

  // To implement: CCD linearity + saturation
  //for(int ix = 0; ix < nspectral; ix++)
  //	frame[ ix ] = response( frame[ ix ] );
}

Detector::~Detector( )
{
  //
}

