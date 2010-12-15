
/// \file Observation.h

/// \class Observation
/// Stores information about an observation.

#include <vector>
#include <string>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

#include "Array.h"
#include "Station.h"
#include "Baseline.h"
#include "Triplet.h"

using namespace std;
typedef std::tr1::unordered_map<string, string> BLNameHash;
typedef std::tr1::unordered_map<string, string> TNameHash;

enum ObsFileType
{
    HOUR_ANGLE,
    DESCRIPTIVE,
    OIFITS
};

class Observation
{
    // Probably shouldn't be public, but it seems to be the convention in this set of code.
  private:
    Array * mArray;
    double mJD;         // The Full Julian Date of the observation (including time)
    double mHA;         // The hour angle of the observation
    bool   mComputeHA;  // A boolean flag to indicate that the hour angle has not been computed.
    
    // Note, these are vectors of pointers.  It is intended that this class only destroy the
    // pointers and NOT the objects themselves.  
    /// \bug Smart pointers would be better here.
    vector<Station*>     mStations;
    vector<Baseline*>    mBaselines;
    vector<Triplet*>     mTriplets;

  private:
    vector<Station*>     FindStations(string telescopes);
    vector<Baseline*>    FindBaselines(vector<Station*> stations, string exclude_baselines);
    vector<Triplet*>     FindTriplets(vector<Station*> stations, string exclude_baselines);  

    double  GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long);
    double  GetSiderealTime(double jd_high, double jd_low, double ee);
    
  public:
  // Constructors
    Observation(Array * array, vector<Station*> stations, string exclude_baselines);
    Observation(Array * array, double hour_angle, string telescopes, string exclude_baselines);
    Observation(Array * array, double mjd, double time, string telescopes, string exclude_baselines);
    
  public:
  // Methods for reading in data files with observations
    static vector <Observation> ReadObservations(Array * array, string filename, string comment_chars, ObsFileType file_type);
    static vector <Observation> ReadObservation_HA(Array * array, string filename, string comment_chars);
    static vector <Observation> ReadObservation_Descriptive(Array * array, string filename, string comment_chars);
    static vector <Observation> ReadObservation_OIFITS(Array * array, string filename);
    
  // Methods for simulating data.
    oi_vis2 GetVis2(string ins_name, Source & source, vector<double> & wavenumbers);
    oi_t3   GetT3(string ins_name, Source & source, vector<double> & wavenumbers);
    
  // Other Methods
    double  GetHA(double targ_ra);
    friend  bool SameObservation(Observation & A, Observation & B);
    
    int         GetNumStations(void);
    Station *   GetStation(int sta_index);
    bool        HasTriplets(void);
};

