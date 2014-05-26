
#ifndef QUADRUPLET_H
#define QUADRUPLET_H

#include <string>
#include <complex>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

#include "Baseline.h"
#include "Station.h"
#include "UVPoint.h"

using namespace std;

typedef std::tr1::unordered_map<std::string, complex<double> > T4Hash;
typedef std::tr1::unordered_map<std::string, double> T4ErrHash;

class Station;
class Array;
class Target;

/// \todo The hash lookup code in this class is in common with the Baseline code, factor this.

class Quadruplet
{
  private:
    Station     *   mStations[4];
    Baseline    *   mBaselines[4]; 
    string name;
    
    T4Hash     mT4Values; // Stores computed bispectrum values
    T4ErrHash  mT4Errors; // Stores computed/stored bispectrum error values.

  public:
    Quadruplet();
    Quadruplet(Array * array, Station * station1, Station * station2, Station * station3, Station * station4);

  private:    
//    complex<double> ComputeT4Error(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac);
    complex<double> ComputeT4(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc);
    string  GetHashKey(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc);
    
  public:
    string  GetName(void);
    bool    ContainsBaseline(string bl_name);
    
    complex<double> GetT4(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc);
    complex<double> GetT4(Target & target, double hour_angle, double wavenumber);
    
    double  GetT4Amp(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc);
    double  GetT4AmpErr(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc);
    double  GetT4Phi(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc);
    double  GetT4PhiErr(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc);

    // Wrapper functions for the UV-coordinate taking versions.
    double  GetT4Amp(Target & target, double hour_angle, double wavenumber);
    double  GetT4AmpErr(Target & target, double hour_angle, double wavenumber);
    double  GetT4Phi(Target & target, double hour_angle, double wavenumber);
    double  GetT4PhiErr(Target & target, double hour_angle, double wavenumber);

        
    int     GetStationID(int station_num);
    
    Baseline * GetBaseline(int baseline_num);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Quadruplet*> QuadrupletHash;

vector<Quadruplet> ComputeQuadruplets(Array * array, vector<Station> & stations);
QuadrupletHash ComputeQuadrupletHash(vector<Quadruplet> & quadruplets);

#endif // QUADRUPLET_H
