/// \class Array array.h
/// \brief A class to represent the array (i.e. the interferometer).

#ifndef ARRAY_H
#define ARRAY_H
#include <cstring>
#include <vector>

#include "Station.h"
#include "Baseline.h"
#include "Triplet.h"
#include "Quadruplet.h"

using namespace std;

// Header files for other libraries
extern "C" {
    #include "exchange.h"
}

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
    /// speed of the wind in km/s
    double wind_speed;
    /// r0 value in meters
    double r0;
    /// The array throughput:
    double throughput;
    
    /// The name of the array
    std::string arrayname;
    
    // Vector types to store the stations and baselines for this array.
    vector<Station>     stations;
    vector<Baseline>    baselines;
    vector<Triplet>     triplets;
    vector<Quadruplet>  quadruplets;
    
  private:
    BaselineHash       bl_hash;
    StationHash        sta_hash;
    TripletHash        tri_hash;
    QuadrupletHash     quad_hash;

  public:
    // Constructor
    Array();

    ~Array();

  public:

    void ImportFile(string filename, string comment_chars);
    
    void ParseOptions(char *argv[], int i, int argc);
    
    /// \todo The Get* functions that deal with hashes should throw an exception if the key is not found.
    Station     *   GetStation(int station_index);
    Station     *   GetStation(string sta_name);
    Baseline    *   GetBaseline(string baseline_name);
    Baseline    *   GetBaseline(int sta1, int sta2);

    Quadruplet  *   GetQuadruplet(string quadruplet_name);
    Quadruplet  *   GetQuadruplet(int sta1, int sta2, int sta3, int sta4);

    Triplet     *   GetTriplet(string triplet_name);
    Triplet     *   GetTriplet(int sta1, int sta2, int sta3);    

    
    
    double      GetLatitude(void);
    double      GetLongitude(void);
    double      GetAltitude(void);
    int         GetNumStations(void);
    string      GetArrayName(void);
    double	Get_r0(void);
    double	GetWindSpeed(void);
    double	GetThroughput(void);
    
    vector<Station*> GetAllStations(void);
    string	GetAllStationNames(void);
    
    oi_array    GetOIArray(void);
};

#endif // #ifndef ARRAY_H
