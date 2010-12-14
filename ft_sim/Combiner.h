/// \file Combiner.h
/// Header file for the combiner class

#include "Matrix.h"

class Beam;

/// \class Combiner
class Combiner
{
  public:
    int Nbeams;

    Row < double >flux;

    Matrix < Complex > fringe;  // stores the complex atmospheric fringes only
    Matrix < double >OPD;       // stores the phase of the fringes due to
    // instrument + atmosphere - used by the
    // approximated spectral integration
    Matrix < double >fringephase;       // stores the phase of the fringes due 
                                        // 
    // 
    // to instrument + atmosphere
    Row < Beam * >beamlist;
    double delay_offset;        // use for testing various OPDs

    double diagnostic_atm;

    double diagnostic_DL;

    Combiner(double delay_offset, int Nbeams, Beam * beam, ...);

    void Update(double time, double wavenumber);

    virtual ~ Combiner();
};
