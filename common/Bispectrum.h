/// \file Bispectrum.h

#include "Matrix.h"

// Bispectrum class
class Bispectrum
{
  public:

    // < triple amplitude without error
    Matrix < double >trueAmp;	
    // < standard deviation of  triple amplitude
    Matrix < double >t3amperr;	
    // < triple amplitude with  noise
    Matrix < double >t3amp;	

    // < closure phase without noise
    Matrix < double >truePhi;	
    // < standard deviation of closure phase
    Matrix < double >t3phierr;	
    // < closure phase with noise
    Matrix < double >t3phi;	
    // < u coordinate in the UV plane of 1st baseline
    Matrix < double >u1;		
    // < v coordinate in the UV plane of 1st baseline
    Matrix < double >v1;		
    // < u coordinate in the UV plane of 2nd baseline
    Matrix < double >u2;		
    // < v coordinate in the UV plane of 2nd baseline
    Matrix < double >v2;		

    // < index of telescope 1
    Matrix < int >t1;			
    // < index of telescope 2
    Matrix < int >t2;			
    // < index of telescope 3
    Matrix < int >t3;			
    // < incoherent integration time for all points
    double int_time;			
};

