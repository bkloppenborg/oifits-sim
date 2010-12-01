/// \file Station.h
/// Header file for the Station Class

class AtmosphereLayer;

/// \class Station Station.h
/// \brief A class to represent a station (i.e. telescope) in an
// interferometer
class Station
{
  public:
    Station(double coordx, double coordy, double coordz, AtmosphereLayer * atm, double gain, double diameter);

  public:
    // AS 2010-06-18
    std::string staname;        // Name of the station
    double x;                   // Onfloor coordinate in meters
    double y;                   // Onfloor coordinate in meters
    double z;                   // Above/below floor coordinate in meters

    AtmosphereLayer *layer;

    double diameter;            // effective diameter of the station
    double gain;                // telescope gain (=1.0 by default, 0.0 to
    // shut down a station)

    // Removed JSY 2009-10-05 as it introduces an unwanted dependency on the
    // AtmosphereLayer implementation
    // void Update( double current_time );
    virtual ~ Station();
};
