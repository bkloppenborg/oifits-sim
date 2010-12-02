
/// \file Observation.h

/// \class Observation
/// Stores information about an observation.

#include <vector>
#include <string>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

#include "Array.h"
#include "Baseline.h"
#include "Station.h"

using namespace std;
typedef std::tr1::unordered_map<string, string> BLNameHash;

class Observation
{
    // Probably shouldn't be public, but it seems to be the convention in this set of code.
  public:
    Array * mArray;
    double mJD;         // The Full Julian Date of the observation (including time)
    double mHA;         // The hour angle of the observation
    bool   mComputeHA;  // A boolean flag to indicate that the hour angle has not been computed.
    
    vector<Baseline>    mBaselines;
    vector<Station>     mStations;

  private:
    vector<Baseline> FindBaselines(vector<Station> stations, string exclude_baselines);
    vector<Station> FindStations(string telescopes);
    double  GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long);
    double  GetSiderealTime(double jd_high, double jd_low, double ee);
    
  public:
  // Constructors
    Observation(Array * array, double hour_angle, string telescopes, string exclude_baselines);
    Observation(Array * array, double mjd, double time, string telescopes, string exclude_baselines);
    
  // Other Methods
    double  GetHA(double targ_ra);
    friend  bool SameObservation(Observation & A, Observation & B);
    
    int         GetNumStations(void);
    Station &   GetStation(int sta_index);
    
};
