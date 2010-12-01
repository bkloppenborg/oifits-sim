/// \file PowerSpectrum.h
/// Definition for the the PowerSpectrum class.


/// PowerSpectrum class
class PowerSpectrum
{
  public:

    // power spectrum without noise
    Matrix < double >pow;		
    // standard deviation of power spectrum
    Matrix < double >vis2err;	

    // power spectrum with noise
    Matrix < double >vis2data;	

    // u coordinate in the UV plane
    Matrix < double >u;		
    // v coordinate in the UV plane
    Matrix < double >v;		

    // index of telescope 1
    Matrix < int >t1;			
    // index of telescope 2
    Matrix < int >t2;			
    // incoherent integration time for all points
    double int_time;			
};
