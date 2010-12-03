/// \class Array array.h
/// \brief A class to represent the array (i.e. the interferometer).

#ifndef ARRAY_H
#define ARRAY_H
#include <string>
#include <vector>


#include "Matrix.h"
#include "Baseline.h"
#include "Station.h"
using namespace std;


// Forward declarations
class Station;
class Baseline;

// A quick struct for the hash_map below.
struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

class Array
{
    /// \todo Make all datamembers private unless absolutely necessary.
    /// \todo Compute Station indicies as stations are read in.
  private:
    /// latitude of the array in radians
    double latitude;
    /// longitude of the array in radians
    double longitude;
    /// altitude of the array in meters
    double altitude;
    
    /// The name of the array
    std::string arrayname;
    
    // Vector types to store the stations and baselines for this array.
    vector<Station> stations;
    vector<Baseline> baselines;
    
  private:
    BaselineHash    bl_hash;
    StationHash     sta_hash;

  public:
    // Constructor
    Array(string filename, string comment_chars);
    Array(std::string arrayname, double latitude, double longitude, double altitude, int nstations, Station * stations);

  public:
    Baseline &  GetBaseline(string baseline_name);
    Station &   GetStation(int station_index);
    Station &   GetStation(string sta_name);
    double      GetLatitude(void);
    double      GetLongitude(void);
    double      GetAltitude(void);
    int         GetNumStations(void);
    string      GetArrayName(void);
    
    vector<Station> GetAllStations(void);
};

#endif // #ifndef ARRAY_H
