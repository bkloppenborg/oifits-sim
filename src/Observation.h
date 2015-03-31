
/// \file Observation.h

/// \class Observation
/// Stores information about an observation.

#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <iostream>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

extern "C" {
    #include "exchange.h"
    #include "oifile.h"
	#include "random.h"
}

class Array;
class Target;
class Baseline;
class NoiseModel;
class Combiner;
class SpectralMode;
class UVPoint;
class Station;
class Triplet;
class Quadruplet;

using namespace std;
typedef std::tr1::unordered_map<string, string> BLNameHash;
typedef std::tr1::unordered_map<string, string> TNameHash;

enum ObsType
{
    HOUR_ANGLE = 0,
    DESCRIPTIVE = 1,
    OIFITS = 2
};


class Observation
{
  protected:
    Array * mArray;
    double  mJD;
    bool    mbHasTriplets;
    bool    mbHasQuadruplets;
    
    ObsType mObsType;
    
    vector<Station*>       mStations;
    vector<Baseline*>      mBaselines;
    vector<Triplet*>       mTriplets;
    vector<Quadruplet*>    mQuadruplets;
    
    vector<Station*>       FindStations(string telescopes);
    vector<Baseline*>      FindBaselines(vector<Station*> stations, string exclude_baselines);
    vector<Triplet*>       FindTriplets(vector<Station*> stations, string exclude_baselines);  
    vector<Quadruplet*>    FindQuadruplets(vector<Station*> stations, string exclude_baselines);  


  public:
    Observation(void);
    ~Observation(void);
    
    //static vector <Observation*> ReadObservations(Array * array, string filename, string comment_chars);

   static vector <Observation*> ParseCommandLine(Array * array, char *argv[], int i, int argc, string comment_chars);

  private:
    static vector <Observation*> ParseCommandLineObs(Array * array, char *argv[], int i, int argc);
    static vector <Observation*> ImportFile(Array * array, string filename, string comment_chars);
  
  public:

    virtual oi_vis2 GetVis2(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed) =0 ;
    virtual oi_t3   GetT3(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed) =0 ;
    virtual oi_t4   GetT4(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed) =0;

    int         GetNumStations(void);
    Station *   GetStation(int sta_index);
    bool HasTriplets(void);
    bool HasQuadruplets(void);

    ObsType     GetObsType(void);
  
};

#endif // OBSERVATION_H
