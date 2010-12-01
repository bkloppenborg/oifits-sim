
#include <iostream>

#include "SpatialFilter.h"
#include "DelayLine.h"

using std::cout;

SpatialFilter::SpatialFilter(int filteringtype, double pixellation,
    double aperture_number, double pinhole_diameter)
{
  this-> filteringtype = filteringtype;
  this-> pixellation = pixellation;
  this-> aperture_number = aperture_number;
  this-> pinhole_diameter = pinhole_diameter;
}

SpatialFilter::~SpatialFilter()
{
}

void SpatialFilter::ComputePinhole(int size, int padded_size)//add in zero padding thing
{
  cout << "Initializing the pinhole\n";
  this->pinhole.setsize(padded_size, padded_size);
  Matrix<double> mask(padded_size, padded_size);
  //double pupil_pixellation = beamwidth/size; // m/pix
  //double focal_pixellation = padded_size / pupil_pixellation; // pix/m size of padded_sizeded screen / pixellation of pupil
  //double airy_pattern_in_pix = 2.0;
  //double pinhole_radius_m = focal_pixellation * airy_pattern_in_pix ; // 1*airy pattern is cutting off all wings
  double pinhole_radius_pix = 2.44 * padded_size / size;//multiply this by some factor alpha to take alpha airy pattern radii
  mask = CircularMask(1.0 * pinhole_radius_pix, padded_size);

  // Move DC component from center to (0,0)
  for (int iy = 0; iy < padded_size / 2; iy++)
  {
    for (int ix = 0; ix < padded_size / 2; ix++)
      pinhole[ ix ][ iy ] = Complex(mask[ padded_size / 2 + ix ][ padded_size
          / 2 + iy ], 0.);

    for (int ix = padded_size / 2; ix < padded_size; ix++)
      pinhole[ ix ][ iy ] = Complex(mask[ ix - padded_size / 2 ][ iy
          + padded_size / 2 ], 0.);
  }

  for (int iy = padded_size / 2; iy < padded_size; iy++)
  {

    for (int ix = 0; ix < padded_size / 2; ix++)
      pinhole[ ix ][ iy ] = Complex(mask[ ix + padded_size / 2 ][ iy
          - padded_size / 2 ], 0.);

    for (int ix = padded_size / 2; ix < padded_size; ix++)
      pinhole[ ix ][ iy ] = Complex(mask[ ix - padded_size / 2 ][ iy
          - padded_size / 2 ], 0.);

  }

}

Matrix<double> CircularMask(double diameter, int screensize)
{
  Matrix<double> disk(screensize, screensize);
  double rmax = diameter * diameter;
  double rsq;
  for (int ix = 0; ix < screensize; ix++)
  {
    for (int iy = 0; iy < screensize; iy++)
    {
      rsq = (double) (ix * ix + iy * iy);
      if (rsq <= rmax)
        disk[ ix ][ iy ] = 1.0;
      else
        disk[ ix ][ iy ] = 0.0;

    }
  }
  return disk;
}

/*
 void SpatialFilter::ComputePinhole(int size) {
 this->pinhole.setsize(size, size);
 Matrix<Complex> pinhole_spectrum(size, size);

 // The pinhole is computed as the Fourier Transform of an Airy pattern amplitude, 2*J1(r)/r
 double r;
 for (int i = 0; i < size; i++) {
 for (int j = 0; j < size; j++) {
 r = 0.01 * sqrt(double(i - size / 2) * double(i - size / 2)
 + double(j - size / 2) * (double(j - size / 2)));

 if (r > 0.)
 pinhole_spectrum[i][j] = Complex(2. * j1(r) / r, 0.);
 else
 pinhole_spectrum[i][j] = Complex(1., 0.);
 }
 }
 cout << "Pinhole spectrum norm: " << norm2(pinhole_spectrum) << "\n";
 double norm_pinhole = sqrt(norm2(pinhole_spectrum));
 FFT(pinhole_spectrum, pinhole, 1, 0);

 // Transform into complex (for future multiplications) and normalize
 for (int i = 0; i < size; i++)
 for (int j = 0; j < size; j++)
 pinhole[i][j] = Complex(pinhole[i][j].real(), 0.); // creal ? abs ?

 cout << "Pinhole norm: " << norm2(pinhole) << "\n";

 }
 */

/*
 Fiber::Fiber(int screensize, double diameter, double modesize=0.9)
 {
 double norm = sqrt(2. / pi) / ( modesize * diameter / 2.0 ) ;
 double scale=1.0 / pow ( modesize * diameter / 2.0 , 2) ;
 // Normalisation so that integral power over mode is unity
 mode.setsize( screensize, screensize );
 mask.setsize( screensize, screensize );
 mask = CircularMask( diameter, screensize );

 for(int ix=0 ; ix < screensize ; ix++)
 for(int iy=0 ; iy < screensize ; iy++)
 mode[ ix ][ iy ] =  norm * exp( - scale * (double) (  ix * ix + iy * iy ) ); // need to check the constant term

 this->opdmode = mode * mask ;
 this->opdmode /=  total( this->opdmode ) ;
 }

 Complex Fiber::Couple(Beam& Beam1)
 {
 // Returns complex amplitude coupled into fiber
 return scalprod( Beam1.pupil , mode );
 }

 double Fiber::Opd(Beam& Beam1)
 {
 // Returns mode-weighted average OPD over aperture
 return scalprod( Beam1.phase , opdmode );
 }

 //
 // class FiberCombiner
 // {
 // fiber = Fiber(screenSize,diameter,modeSize);
 // amplitudes
 // opd
 // amplitude
 // phase
 // flux
 //
 // void Couple(beams)
 // {
 // amplitudes = fiber.Couple(beams.pupils)
 // opd =  fiber.Opd(beams.screens)
 // }
 // void Combine(beams)
 // {
 // amplitudes = fiber.Couple( beams.pupils )
 // opd = fiber.Opd( beams.screens )
 // amplitude = abs( amplitudes[0] * conjugate(amplitudes[1]) )
 // phase     = arg( amplitudes[0] * conjugate(amplitudes[1]) )
 // fringe    = complex(amplitude , phase - opd[0] + opd[1] )
 // flux      = total( abs2(amplitudes) )
 // }
 //
 // }

 */
