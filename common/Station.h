/// \file Station.h
/// Header file for the Station Class

class AtmosphereLayer;

/// \class Station Station.h
/// \brief A class to represent a station (i.e. telescope) in an
// interferometer
class Station
{

  public:
    // AS 2010-06-18
    std::string staname;        // Name of the station
    double north;               // Onfloor coordinate in meters
    double east;                // Onfloor coordinate in meters
    double up;                  // Above/below floor coordinate in meters

    AtmosphereLayer * layer;

    double diameter;            // effective diameter of the station
    double gain;                // telescope gain (=1.0 by default, 0.0 to shut down a station)

  public:
    Station(string station_name, double North, double East, double Up, double gain, double diameter);
    
    virtual ~ Station();
    
    string GetName();
};
