//============================================================================
// Name        : AtmosphereLayer.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron and David Buscher
//============================================================================

/// \file AtmosphereLayer.cpp

// The Atmosphere module simulates the OPD variation as a function of
// time and location accross the aperture of each telescope.  An OPD
// mask is created for the entire interferometer array, then
// translated accross the apertures at the speed of the wind.

#include "AtmosphereLayer.h"
#include "Matrix.h"
#include "random.h"

void AtmosphereLayer::Update( double time )
{
  if ((time != previous_time) || (time == 0.0))
  {
    double dx = time * wind_speed * pixellation * cos(wind_direction * PI / 180.);
    double dy = time * wind_speed * pixellation * sin(wind_direction * PI / 180.);

    // Interpolate the screen
    Interp2d(phase_generator, phase_screen, dx, dy, 1.0, 0.0);

    if (scintillation_diameter > 0)
      Interp2d(amplitude_generator, amplitude_screen, dx, dy, 1.0, 0.0);

    // Note: those screens are square matrices. They will be masked in the beam module
    previous_time = time;
  }
}

/// \brief Initialization for the Atmosphere Layer
///
/// Set up parameters for generating a simulated atmospheric phase screen.
/// Assumes Kolmogorov spectrum of atmospheric disturbances
/// [c.f. Tatarski 1961,1971], s.t. the phase structure function
/// is given by
//
/// D(r) = < [ phi( r' ) - phi( r' + r ) ]**2 > = 6.88 * ( r / r0 )**5/3
/// where r0 is the Fried parameter.

void AtmosphereLayer::InitScreen( )
{
  double r, rsq, nyqsq;
  int nfft = generatorsize;
  int nyq = nfft / 2;
  nyqsq = nyq * nyq;
  Matrix<double> amp(nyq + 1, nyq + 1);
  Matrix<Complex> spectrum(nyq + 1, nfft);

  ///////////////////////////////////
  //
  // Phase screen
  //
  //////////////////////////////////


  phase_generator.setsize(nfft, nfft);
  phase_screen.setsize(screensize, screensize);

  // Initialize the spectrum with Gaussian random numbers

  for (int ix = 0; ix < nyq + 1; ix++)
  {
    for (int iy = 0; iy < nfft; iy++)
    {
      spectrum[ ix ][ iy ].real() = Rangauss(random_number_seed);
      spectrum[ ix ][ iy ].imag() = Rangauss(random_number_seed);
    }
  }

  // Calculate the amplitude filter array
  double generator_width = (double) nfft / pixellation; // size of the generator in meters
  double C = sqrt(0.0229 * pow(generator_width / r0, 5.0 / 3.0) / 2.);

  for (int ix = 0; ix < nyq + 1; ix++)
  {
    for (int iy = 0; iy < nyq + 1; iy++)
    {
      rsq = (double) (ix * ix + iy * iy);
      if (rsq > nyqsq || rsq == 0.0)
      {
        amp[ ix ][ iy ] = 0.0;
      }
      else
      {
        if (this->type == 0)
          amp[ ix ][ iy ] = C * pow(rsq, -11.0 / 12.0);
        if (this->type == 1)
          amp[ ix ][ iy ] = C * pow(rsq + generator_width * generator_width / (outerscale * outerscale), -11.0 / 12.0);
        if (this->type == 2)
          amp[ ix ][ iy ] = C * pow(rsq + generator_width * sqrt(rsq) / outerscale, -11.0 / 12.0);
        if ((this->type > 2) || (this->type < 0))
          amp[ ix ][ iy ] = 0.;

      }
    }
  }

  // Then filter the positive quarter-plane
  for (int iy = 0; iy < nyq + 1; iy++)
  {
    for (int ix = 0; ix < nyq + 1; ix++)
    {
      spectrum[ ix ][ iy ] *= amp[ ix ][ iy ];
    }
  }

  // Filter negative quarter-plane, reusing the amplitudes from the positive part
  for (int iy = nyq + 1; iy < nfft; iy++)
  {
    for (int ix = 0; ix < nyq + 1; ix++)
    {
      spectrum[ ix ][ iy ] *= amp[ ix ][ nfft - iy ];
    }
  }

  // Fix up the on-axis points
  for (int iy = 1; iy < nyq; iy++)
    spectrum[ 0 ][ nfft - iy ] = conj(spectrum[ 0 ][ iy ]);
  spectrum[ 0 ][ 0 ] = 0.0; // need to check this

  // FFTW3, Half complex to Real transformation
  // input : the spectrum
  // output : the generator
  // Matrix <-> fftw array and C++ complex <-> fftw_complex conversions are done here

  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nyq + 1) * nfft);
  double* phase_out = (double *) fftw_malloc(sizeof(double) * nfft * nfft);
  for (int iy = 0; iy < nfft; iy++)
  {
    for (int ix = 0; ix < nyq + 1; ix++)
    {
      in[ ix + (nyq + 1) * iy ][ 0 ] = spectrum[ ix ][ iy ].real();
      in[ ix + (nyq + 1) * iy ][ 1 ] = spectrum[ ix ][ iy ].imag();
    }
  }

  fftw_plan phs_p = fftw_plan_dft_c2r_2d(nfft, nfft, in, phase_out, FFTW_ESTIMATE);
  fftw_execute(phs_p);
  fftw_destroy_plan(phs_p);

  // Correct for the fact that FFTW is unnormalized
  for (int iy = 0; iy < nfft; iy++)
  {
    for (int ix = 0; ix < nfft; ix++)
    {
      phase_generator[ ix ][ iy ] = phase_out[ ix + nfft * iy ];
    }
  }

  fftw_free(in);
  fftw_free(phase_out);

  ///////////////////////////////////////////////////
  //
  // Amplitude/ scintillation (if required)
  //
  ///////////////////////////////////////////////////

  if (scintillation_diameter > 0)
  {
    amplitude_generator.setsize(nfft, nfft);
    amplitude_screen.setsize(screensize, screensize);

    // Initialize the spectrum with Gaussian random numbers

    for (int ix = 0; ix < nyq + 1; ix++)
    {
      for (int iy = 0; iy < nfft; iy++)
      {
        spectrum[ ix ][ iy ].real() = Rangauss(random_number_seed);
        spectrum[ ix ][ iy ].imag() = Rangauss(random_number_seed);
      }
    }

    // scintillation expressions adapted from G. Sedmak, JOSA 2004
    // "Implementation of fast Fourier-transform-based simulations of
    // extra-large atmospheric phase and scintillation screens"
    // chosen value S0 = 1

    for (int ix = 0; ix < nyq + 1; ix++)
    {
      for (int iy = 0; iy < nyq + 1; iy++)
      {
        rsq = (double) (ix * ix + iy * iy);
        if (rsq > nyqsq || rsq == 0.0)
        {
          amp[ ix ][ iy ] = 0.0;
        }
        else
        {
          rsq = rsq / generator_width;
          r = sqrt(rsq);
          C = 1.54 / (lambda0 * lambda0) * pow(r, -11. / 3.) * pow(2. * j1(PI * scintillation_diameter * r) / (PI * scintillation_diameter * r), 2.0);
          double h1 = 2250.;
          double h2 = 10000.;
          double delta_h1 = 3600.;
          double delta_h2 = 10000.;
          amp[ ix ][ iy ] = sqrt(C * (Cn2(h1) * delta_h1 * pow(sin(PI * lambda0 * h1 * r * r), 2.0) + Cn2(h2) * delta_h2 * pow(sin(PI
              * lambda0 * h2 * r * r), 2.0))) / generator_width;
        }
      }
    }

    // Then filter the positive quarter-plane
    for (int iy = 0; iy < nyq + 1; iy++)
    {
      for (int ix = 0; ix < nyq + 1; ix++)
      {
        spectrum[ ix ][ iy ] *= amp[ ix ][ iy ];
      }
    }

    // Fix up the on-axis points
    for (int iy = 1; iy < nyq; iy++)
      spectrum[ 0 ][ nfft - iy ] = conj(spectrum[ 0 ][ iy ]);
    spectrum[ 0 ][ 0 ] = 0.0; // need to check this

    // FFTW3, Half complex to Real transformation
    // input : the spectrum
    // output : the generator
    // Matrix <-> fftw array and C++ complex <-> fftw_complex conversions are done here

    double* amplitude_out = (double *) fftw_malloc(sizeof(double) * nfft * nfft);
    for (int iy = 0; iy < nfft; iy++)
    {
      for (int ix = 0; ix < nyq + 1; ix++)
      {
        in[ ix + (nyq + 1) * iy ][ 0 ] = spectrum[ ix ][ iy ].real();
        in[ ix + (nyq + 1) * iy ][ 1 ] = spectrum[ ix ][ iy ].imag();
      }
    }

    fftw_plan amp_p = fftw_plan_dft_c2r_2d(nfft, nfft, in, amplitude_out, FFTW_ESTIMATE);
    fftw_execute(amp_p);
    fftw_destroy_plan(amp_p);

    // Correct for the fact that FFTW is unnormalized
    for (int iy = 0; iy < nfft; iy++)
    {
      for (int ix = 0; ix < nfft; ix++)
      {
        amplitude_generator[ ix ][ iy ] = amplitude_out[ ix + nfft * iy ];
      }
    }

    fftw_free(in);
    fftw_free(amplitude_out);
  }

  // Finally, setup both the phase and amplitude screens using the newly created generators
  Update(0.0);

}

double AtmosphereLayer::Cn2( double h )
{
  double c0 = 1.027e-3;
  double h0 = 4632.;
  double h1 = 1000.;
  double h2 = 1500.;
  double C = c0 * c0 * pow(r0, -5.0 / 3.0) * pow(2. * PI / lambda0, -2.) * (pow(h / h0, 10.) * exp(-h / h1) + exp(-h / h2));
  return C;
}

AtmosphereLayer::AtmosphereLayer( int type , int generatorsize , int screensize , double lambda0 , double r0 , double pixellation ,
    double windspeed , double winddirection , double outerscale , double scintillation_diameter )
{
  // Create a Layer for wavelength = lambda0s
  this->type = type; // 0: Kolmogorov 1: Von Karman 2: Greenwood-Tarazano
  this->generatorsize = generatorsize;
  this->screensize = screensize;
  this->lambda0 = lambda0;
  this->r0 = r0;
  this->pixellation = pixellation;
  this->wind_speed = windspeed;
  this->wind_direction = winddirection;
  this->previous_time = 0.0;
  this->outerscale = outerscale;
  this->scintillation_diameter = scintillation_diameter;
  // Memory allocation for the screen and generator is done in Initscreen
  InitScreen();
}

AtmosphereLayer::~AtmosphereLayer( )
{
  // nothing to do
}
