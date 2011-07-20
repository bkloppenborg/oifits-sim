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

/// Constructs a Station object.
/// ABC_is_XYZ permits specification of the XYZ coordinates if true.  If false input the (A=North, B=East, C=Up) coords
Station::Station(double array_lat, string station_name, int sta_index, bool ABC_is_XYZ, double A, double B, double C, double gain, double diameter)
{
    this->sta_name = station_name;
    this->sta_index = sta_index;

    this->layer = NULL;
    this->gain = gain;
    this->diameter = diameter;
    
    if(ABC_is_XYZ)
    {
        this->xyz[0] = A;
        this->xyz[1] = B;
        this->xyz[2] = C;   
    }
    else
    {
        this->NEU[0] = A;
        this->NEU[1] = B;
        this->NEU[2] = C;    
        this->ComputeXYZ(array_lat);  
    } 
}

Station::Station(double array_lat, string station_name, int sta_index, bool ABC_is_XYZ, double A, double B, double C, AtmosphereLayer * atm, double gain, double diameter)
{
    this->sta_name = station_name;
    this->sta_index = sta_index;
    this->layer = atm;
    this->gain = gain;
    this->diameter = diameter;
    
    if(ABC_is_XYZ)
    {
        this->xyz[0] = A;
        this->xyz[1] = B;
        this->xyz[2] = C;   
    }
    else
    {
        this->NEU[0] = A;
        this->NEU[1] = B;
        this->NEU[2] = C;    
        this->ComputeXYZ(array_lat);  
    } 
}

void    Station::ComputeXYZ(double phi)
{
    // Convert (North, East, Up) to (x,y,z) as defined by the APIS++ standards
    // see:    http://aips2.nrao.edu/docs/glossary
    // for more information.
    
    // Parameters:
    //  phi is the array latitude
    // Output:
    // xyz (class variable) is ordered (x, y, z)
    
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
StationHash ComputeStationHash(vector<Station> & stations)
{
    StationHash hash;
    
    for(int i = 0; i < stations.size(); i++)
    {
        hash.insert(StationHash::value_type(stations[i].GetName(), &stations[i]) );
    }
    
    return hash;
}


