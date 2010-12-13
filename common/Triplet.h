
#ifndef TRIPLET_H
#define TRIPLET_H

#include <string>
#include <complex>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

#include "Baseline.h"
#include "Station.h"

using namespace std;

typedef std::tr1::unordered_map<std::string, complex<double> > BisHash;

class Station;
class Array;

/// \todo The hash lookup code in this class is in common with the Baseline code, factor this.

class Triplet
{
  private:
    Station     mStations[3];
    Baseline    mBaselines[3]; 
    string name;
    
    BisHash mBisValues; // Stores computed bispectrum values
    BisHash mBisErrors; // Stores computed/stored bispectrum error values.

  public:
    Triplet();
    Triplet(Array * array, Station & station1, Station & station2, Station & station3);

  private:    
    complex<double> ComputeBisError(Source & source, double hour_angle, double wavenumber);
    complex<double> ComputeBispectra(Source & source, double hour_angle, double wavenumber);
    
    string  GetHashKey(Source & source, double hour_angle, double wavenumber);
    
  public:
    string  GetName(void);
    bool    ContainsBaseline(string bl_name);
    complex<double> GetBispectra(Source & source, double hour_angle, double wavenumber);
    
    complex<double> GetBisError(Source & source, double hour_angle, double wavenumber);
    void    SetBisError(Source & source, double hour_angle, double wavenumber);
    int     GetStationID(int station_num);
    
    Baseline & GetBaseline(int baseline_num);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Triplet> TripletHash;

vector<Triplet> ComputeTriplets(Array * array, vector<Station> stations);
TripletHash ComputeTripletHash(vector<Triplet> triplets);

#endif // TRIPLET_H
