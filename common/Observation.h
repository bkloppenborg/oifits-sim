
/// \file Observation.h

/// \class Observation
/// Stores information about an observation.

#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <vector>
#include <string>
#include <stdexcept>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

extern "C" {
    #include "exchange.h"
    #include "oifile.h"
}

#include "Array.h"
#include "Station.h"
#include "Baseline.h"
#include "Triplet.h"
#include "Source.h"

using namespace std;
typedef std::tr1::unordered_map<string, string> BLNameHash;
typedef std::tr1::unordered_map<string, string> TNameHash;

enum ObsType
{
    HOUR_ANGLE,
    DESCRIPTIVE,
    OIFITS
};


class Observation
{
  protected:
    Array * mArray;
    double  mJD;
    
    ObsType mObsType;
    
    vector<Station*>    mStations;
    vector<Baseline*>   mBaselines;
    vector<Triplet*>    mTriplets;
    
    vector<Station*>     FindStations(string telescopes);
    vector<Baseline*>    FindBaselines(vector<Station*> stations, string exclude_baselines);
    vector<Triplet*>     FindTriplets(vector<Station*> stations, string exclude_baselines);  

  public:
    Observation(void);
    ~Observation(void);
    
    static vector <Observation*> ReadObservations(Array * array, string filename, string comment_chars, ObsType file_type);
  
    virtual oi_vis2 GetVis2(string ins_name, Source & source, vector<double> & wavenumbers) {};
    virtual oi_t3   GetT3(string ins_name, Source & source, vector<double> & wavenumbers) {};

    int         GetNumStations(void);
    Station *   GetStation(int sta_index);
    bool        HasTriplets(void);
    
    ObsType     GetObsType(void);
  
};

#endif // OBSERVATION_H
