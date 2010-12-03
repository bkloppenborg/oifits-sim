/// \file Baseline.h
/// Header file for the Baseline class.

#ifndef BASELINE_H
#define BASELINE_H

// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>
#include <string>
#include <vector>
using namespace std;

// Forward class declarations:
class UVPoint;
class Array;
class Station;

/// \class Baseline simulator.h
/// \brief A class representing a baseline.
class Baseline
{
  private:
    double xyz[3];  // True XYZ coordinates
    
    string name;
    int indicies[2];

  public:
    Baseline(void);
    Baseline(Array * array, int station1, int station2);
    Baseline(Station & station1, Station & station2);
    
    /// \todo Rewrite this function to work with the new class definition.
    //double Geometric_OPD(double hour_angle, double source_declination);
    
    UVPoint UVcoords(double hour_angle, double source_declination, double wavenumber);
    string  GetName(void);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Baseline> BaselineHash;

vector<Baseline> ComputeBaselines(vector<Station> stations);
BaselineHash ComputeBaselineHash(vector<Baseline> baselines);


#endif // #ifndef BASELINE_H
