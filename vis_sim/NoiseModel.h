#include "Matrix.h"

class SpectralMode;
class Source;
class Array;
class Instrument;

double Strehl(double wl, double r0, double D);

double TimeInt(double wl, double r0, double v);

// AS 2010-06-24
// added a new dependency of PhCount on median_wl
double PhCount(Source & target, Array & s, Instrument & inst,
   double median_wl, double wl, double Dlambda);

// AS 2010-06-22
// eliminated the incoherent integration time from arguments of the function,
// since it is now contained in the instrument class
double VarPow1(double Pow, SpectralMode & spec, int specIndex,
   Source & target, Array & s, Instrument & inst);
double VarUnbiasedPow1(double Pow, SpectralMode & spec, int specIndex,
   Source & target, Array & s, Instrument & inst);

Matrix < double > VarCloPhase(SpectralMode & spec,
   Matrix < Complex > &visibility, Source & target,
   Array & s, Instrument & inst, int Nbs, int Nha);
