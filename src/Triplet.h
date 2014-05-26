
#ifndef TRIPLET_H
#define TRIPLET_H

#include <string>
#include <complex>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

#include "Baseline.h"
#include "Station.h"
#include "UVPoint.h"

using namespace std;

typedef std::tr1::unordered_map<std::string, complex<double> > T3Hash;
typedef std::tr1::unordered_map<std::string, double> T3ErrHash;

class Station;
class Array;
class Target;

/// \todo The hash lookup code in this class is in common with the Baseline code, factor this.

class Triplet
{
  private:
    Station     *   mStations[3];
    Baseline    *   mBaselines[3]; 
    string name;
    
    T3Hash     mT3Values; // Stores computed bispectrum values
    T3ErrHash  mT3Errors; // Stores computed/stored bispectrum error values.

  public:
    Triplet();
    Triplet(Array * array, Station * station1, Station * station2, Station * station3);

  private:    
//    complex<double> ComputeT3Error(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    complex<double> ComputeT3(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    string  GetHashKey(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    
  public:
    string  GetName(void);
    bool    ContainsBaseline(string bl_name);
    
    complex<double> GetT3(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    complex<double> GetT3(Target & target, double hour_angle, double wavenumber);
    
    double  GetT3Amp(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    double  GetT3AmpErr(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    double  GetT3Phi(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    double  GetT3PhiErr(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);

    // Wrapper functions for the UV-coordinate taking versions.
    double  GetT3Amp(Target & target, double hour_angle, double wavenumber);
    double  GetT3AmpErr(Target & target, double hour_angle, double wavenumber);
    double  GetT3Phi(Target & target, double hour_angle, double wavenumber);
    double  GetT3PhiErr(Target & target, double hour_angle, double wavenumber);

        
    int     GetStationID(int station_num);
    
    Baseline * GetBaseline(int baseline_num);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Triplet*> TripletHash;

vector<Triplet> ComputeTriplets(Array * array, vector<Station> & stations);
TripletHash ComputeTripletHash(vector<Triplet> & triplets);

#endif // TRIPLET_H
