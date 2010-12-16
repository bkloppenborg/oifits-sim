/// \file Baseline.h
/// Header file for the Baseline class.

#ifndef BASELINE_H
#define BASELINE_H

// Include headers for the unordered map.  Note, this may need to be just <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>
#include <string>
#include <vector>
#include <complex>

// Header files for other libraries
extern "C" {
    #include "exchange.h"
}

using namespace std;

typedef std::tr1::unordered_map<std::string, complex<double> > VisHash;
typedef std::tr1::unordered_map<std::string, double > Vis2Hash;

// Forward class declarations:
class UVPoint;
class Array;
class Station;
class Source;

/// \class Baseline simulator.h
/// \brief A class representing a baseline.
class Baseline
{
  private:
    double xyz[3];  // True XYZ coordinates
    
    string name;
    int indicies[2];
    VisHash     mVisValues; // Stores computed visibility values
    Vis2Hash    mVis2Errors; // Stores computed/stored visibility squared error values.

  public:
    Baseline(void);
    Baseline(Array * array, int station1, int station2);
    Baseline(Station * station1, Station * station2);
    
    /// \todo Rewrite this function to work with the new class definition.
    //double Geometric_OPD(double hour_angle, double source_declination);
    
  private:
    complex<double> ComputeVisibility(Source & source, UVPoint uv);
    double  ComputeVis2Error(Source & source, UVPoint uv);
    
    string  GetHashKey(Source & source, UVPoint uv);
    
  public:
    string  GetName(void);
    complex<double> GetVisibility(Source & source, double hour_angle, double wavenumber);
    complex<double> GetVisibility(Source & source, UVPoint uv_coords);   
    
    double  GetVis2Error(Source & source, double hour_angle, double wavenumber);
    double  GetVis2Error(Source & source, UVPoint uv_coords);    
    
    double  GetVis2(Source & source, double hour_angle, double wavenumber);
    double  GetVis2(Source & source, UVPoint uv_coords);
    
    
    double  GetVis2Err(Source & source, double hour_angle, double wavenumber);   
    
    void    SetVis2Error(Source & source, double hour_angle, double wavenumber, double vis2error);
    void    SetVis2Error(Source & source, UVPoint uv, double vis2error);
    
    UVPoint UVcoords(double hour_angle, double source_declination);
//    oi_vis2_record GetVis2Record(Source & source, double hour_angle, vector<double> wavenumbers);
    
    int GetStationID(int station_num);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Baseline*> BaselineHash;

vector<Baseline> ComputeBaselines(vector<Station> & stations);
BaselineHash ComputeBaselineHash(vector<Baseline> & baselines);


#endif // #ifndef BASELINE_H
