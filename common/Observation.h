
/// \file Observation.h

/// \class Observation
/// Stores information about an observation.

#include <vector>
class Baselines;

typedef std::tr1::unordered_map<std::string, string> BLNameHash;

class Observation
{
    // Probably shouldn't be public, but it seems to be the convention in this set of code.
  public:
    double mJD;         // The Full Julian Date of the observation (including time)
    double mHA;         // The hour angle of the observation
    bool   mComputeHA;  // A boolean flag to indicate that the hour angle has not been computed.
    
    vector<Baseline> baselines;

  private:
    vector<Baseline> FindBaselines(string telescopes, string exclude_baselines);
    double GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long);
    double GetSiderealTime(double jd_high, double jd_low, double ee);
    
  public:
  // Constructors
    Observation();
    Observation(Array & array, double hour_angle, string telescopes, string exclude_baselines);
    Observation(Array * array, double mjd, double time, string telescopes, string exclude_baselines);
    
  // Other Methods
    double GetHA(double targ_ra, double array_lat, double array_long);
    bool friend SameObservation(Observation & A, Observation & B);
    
     
};
