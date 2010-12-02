/// \class Array array.h
/// \brief A class to represent the array (i.e. the interferometer).

#include <string>
#include <vector>
// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

#include "Matrix.h"
#include "Baseline.h"

// Forward declarations
class Station;
class Baseline;

typedef std::tr1::unordered_map<std::string, Baseline> BaselineHash;

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
    BaselineHash bl_hash;


  public:
    // Constructor
    Array(const char * Array_file);
    Array(std::string arrayname, double latitude, double longitude, double altitude, int nstations, Station * stations);
    
  private:
    void        ComputeBaselines(void);

  public:
    Baseline &  GetBaseline(string baseline_name);
    Station &   GetStation(int station_index);
    double      GetLatitude(void);
    double      GetLongitude(void);
    double      GetAltitude(void);
};

