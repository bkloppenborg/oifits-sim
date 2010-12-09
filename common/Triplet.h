
#ifndef TRIPLET_H
#define TRIPLET_H

#include <string>
#include <complex>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

#include "Baseline.h"
#include "Station.h"

using namespace std;


class Station;
class Array;

class Triplet
{
  private:
    Station     mStations[3];
    Baseline    mBaselines[3]; 
    string name;

  public:
    Triplet();
    Triplet(Array * array, Station & station1, Station & station2, Station & station3);
    
  public:
    string  GetName(void);
    bool    ContainsBaseline(string bl_name);
    complex<double> GetBispectra(Source & source, double hour_angle, double wavenumber);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Triplet> TripletHash;

vector<Triplet> ComputeTriplets(Array * array, vector<Station> stations);
TripletHash ComputeTripletHash(vector<Triplet> triplets);

#endif // TRIPLET_H
