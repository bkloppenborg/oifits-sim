
/// \file Observation.h

/// \class Observation
/// Stores information about an observation.

class Observation
{
    // Probably shouldn't be public, but it seems to be the convention in this set of code.
  public:
    double mMJD;
    double mTime;
    
  public:
    Observation();
    Observation(double mjd);
    Observation(double mjd, double time);
    
    double GetJD();
    double GetHA(double targ_ra, double targ_dec, double array_lat, double array_long);
    double GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long);
    double GetSiderealTime(double jd_high, double jd_low, double ee);
     
};

bool SameObservation(Observation & A, Observation & B);
