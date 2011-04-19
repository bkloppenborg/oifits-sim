/// \file Station.h
/// Header file for the Station Class

#ifndef STATION_H
#define STATION_H

// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>
#include <string>
using namespace std;

//#include "Simulator.h"

class AtmosphereLayer;

/// \class Station Station.h
/// \brief A class to represent a station (i.e. telescope) in an
// interferometer
class Station
{
  private:
    string  sta_name;
    int     sta_index;

  public:
    // AS 2010-06-18
    double NEU[3];
    double xyz[3];

    AtmosphereLayer * layer;

    double diameter;            // effective diameter of the station
    double gain;                // telescope gain (=1.0 by default, 0.0 to shut down a station)
    
  public:
    Station();
    Station(double array_lat, string station_name, int sta_index, bool ABC_is_XYZ, double A, double B, double C, double gain, double diameter);
    Station(double array_lat, string station_name, int sta_index, bool ABC_is_XYZ, double A, double B, double C, AtmosphereLayer * atm, double gain, double diameter);

  private:
    void    ComputeXYZ(double phi);
    
  public:
    string  GetName(void);
    int     GetIndex(void);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Station*> StationHash;
StationHash ComputeStationHash(vector<Station> & stations);

#endif // STATION_H
