// The station module

#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <vector>

#include "Station.h"

class AtmosphereLayer;

Station::Station()
{
    this->sta_name = "None";
    this->NEU[0] = 0;
    this->NEU[1] = 0;
    this->NEU[2] = 0;
    this->gain = 0;
    this->diameter = 0;  
    this->sta_index = 0;
}

Station::Station(double array_lat, string station_name, int sta_index, double North, double East, double Up, double gain, double diameter)
{
    this->sta_name = station_name;
    this->sta_index = sta_index;
    this->NEU[0] = North;
    this->NEU[1] = East;
    this->NEU[2] = Up;
    this->layer = NULL;
    this->gain = gain;
    this->diameter = diameter;
    
    this->ComputeXYZ(array_lat);   
}

Station::Station(double array_lat, string station_name, int sta_index, double North, double East, double Up, AtmosphereLayer * atm, double gain, double diameter)
{
    this->sta_name = station_name;
    this->sta_index = sta_index;
    this->NEU[0] = North;
    this->NEU[1] = East;
    this->NEU[2] = Up;
    this->layer = atm;
    this->gain = gain;
    this->diameter = diameter;
    
    this->ComputeXYZ(array_lat);
}

void    Station::ComputeXYZ(double phi)
{
    // Convert (North, East, Up) to (x,y,z) as defined by the APIS++ standards
    // see:    http://aips2.nrao.edu/docs/glossary
    // for more information.
    // the xyz array is ordered (x, y, z)
    this->xyz[0] = -sin(phi) * NEU[0] + cos(phi) * NEU[2];
    this->xyz[1] = NEU[1];
    this->xyz[2] = cos(phi) * NEU[0] + sin(phi) * NEU[2];
}

string  Station::GetName(void)
{
    return this->sta_name;
}

int     Station::GetIndex(void)
{
    return this->sta_index;
}


////////////////////////////////////////////////////////////////////
// Non Class Functions Below
////////////////////////////////////////////////////////////////////
StationHash ComputeStationHash(vector<Station> stations)
{
    StationHash hash;
    
    for(int i = 0; i < stations.size(); i++)
    {
        hash.insert(StationHash::value_type(stations[i].GetName(), stations[i]) );
    }
    
    return hash;
}


