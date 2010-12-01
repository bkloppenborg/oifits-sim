/// \class Array array.h
/// \brief A class to represent the array (i.e. the interferometer).

#include <string>
#include <vector>

#include "Matrix.h"

// Forward declarations
class Station;
class AtmosphereLayer;

class Array
{
  public:
    // AS 2010-06-18 
    // added parameters to the array class (array name, longitude, altitude
    // and names of the #include <cstdarg>stations)
    
    /// The name of the array
    std::string arrayname;
    /// The number of stations in the array.
    int nstations;
    /// vector containing the station names
    std::vector < std::string > staname;
    /// latitude of the array
    double latitude;
    /// longitude of the array
    double longitude;
    /// altitude of the array
    double altitude;
    /// Onfloor station coordinates in meters
    Row < double >x;
    /// Onfloor station coordinates in meters
    Row < double >y;
    /// Above/below floor station coordinates in meters
    Row < double >z;

    /// Effective diameter of the primary mirror of the stations.
    Row < double >diameter;

    /// \brief Station Gains
    ///
    /// Set to 1.0 by default, 0.0 to shut down a station.  
    Row < double >gain;
    /// Atmosphere layers
    Row < AtmosphereLayer * >layer;
    Array(const char *Array_file);

    /// \deprecated An old constructor for the array that is not used anymore 
    // 
    // (AS 2010-06-22)
    Array(double latitude, Row < double >&x, Row < double >&y, Row < double >&z,
          Row < double >&diameter, Row < double >gain,
          Row < AtmosphereLayer * >layer);
    // 
    Array(std::string arrayname, double latitude, double longitude, double altitude, int nstations, Station * station, ...);    // 
    void Update(double current_time);
};

