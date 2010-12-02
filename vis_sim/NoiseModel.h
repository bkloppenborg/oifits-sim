#include "Matrix.h"

class SpectralMode;
class Source;
class Array;
class Instrument;
class Observation;

double Strehl(double wl, double r0, double D);

double TimeInt(double wl, double r0, double v);

double PhCount(Source & target, Observation & obs, Instrument & inst, double median_wl, double wl, double Dlambda);
double CoherentFlux1(double Ni, double Pow, Observation & obs, Instrument & inst);
double VarCoherentFlux1(double Ni, double Pow, Observation & obs, Instrument & inst);
double VarPow1(double Pow, SpectralMode & spec, int specIndex, Source & target, Observation & obs, Instrument & inst);
double VarUnbiasedPow1(double Pow, SpectralMode & spec, int specIndex, Source & target, Observation & obs, Instrument & inst);
Matrix < double > VarCloPhase(SpectralMode & spec, Matrix < Complex > &visibility, Source & target, Observation & obs, Instrument & inst, int Nbs, int Nha);


