/// \file Baseline.cpp
/// Implements the required functions for the baseline class.

#include <cmath>

#include "Baseline.h"
#include "UVPoint.h"
#include "Array.h"
#include "Simulator.h"
#include "Station.h"

Baseline::Baseline()
{
    this->xyz[0] = this->xyz[1] = this->xyz[2] = 0;
    this->name = "";
    this->indicies[0] = 0;
    this->indicies[1] = 0;
}

/// \todo Stop using NEU coordinates here, switch to station.x, station.y, and station.z.
Baseline::Baseline(Station & station1, Station & station2)
{
    // Now calculate the xyz coordinates:
    this->xyz[0] = station1.xyz[0] - station2.xyz[0];
    this->xyz[1] = station1.xyz[1] - station2.xyz[1];
    this->xyz[2] = station1.xyz[2] - station2.xyz[2];
        
    this->name = station1.GetName() + "-" + station2.GetName();
    this->indicies[0] = station1.GetIndex();
    this->indicies[1] = station2.GetIndex();
}

UVPoint Baseline::UVcoords(double hour_angle, double source_declination, double wavenumber)
{
    // First convert all values into radians (they should be degrees or decimal hours (of time) before now.
    double h = hour_angle * PI;
    double delta = source_declination;

    // Now compute the UV coordinates, again according to the APIS++ standards.
    UVPoint uv;
    uv.u = sin(h) * xyz[0]                    + cos(h) * xyz[1];
    uv.v = -sin(delta) * cos(h) * xyz[0]      + sin(delta) * sin(h) * xyz[1]    + cos(delta) * xyz[2];
    uv.w = cos(delta) * cos(h)*xyz[0]         - cos(delta) * sin(h) * xyz[1]    + sin(delta) * xyz[2];
    
    // Scale by the wavenumber (1/wavelength)
    /// \note The coordinates (u,v) calculated here appear to be (-u, -v) with respect to CHARA+MIRC data.
    /// Given that the UV plane is symmetric, this doesn't matter.
    uv.u *= wavenumber;
    uv.v *= wavenumber;
    uv.w *= wavenumber;
    
    return uv;
}

string Baseline::GetName(void)
{
    return this->name;
}

////////////////////////////////////////////////////////////////////
// Non Class Functions Below
////////////////////////////////////////////////////////////////////

/// Computes all possible baselines formed by the specified stations.
vector<Baseline> ComputeBaselines(vector<Station> stations)
{
    int num_stations = stations.size();
    
    vector<Baseline> baselines;
    
    // Now compute all of the baselines and make a hash table for each baseline value
    for(int i = 0; i < num_stations; i++)
    {
        for(int j = i+1; j < num_stations; j++)
        {
            // Create a new baseline, append it to the list of baselines
            baselines.push_back(Baseline(stations[i], stations[j]));
        }
    }
    
    return baselines;
}

/// Computes a (baseline_name, baseline_object) hash table.
BaselineHash ComputeBaselineHash(vector<Baseline> baselines)
{
    BaselineHash hash;
    
    for(int i = 0; i < baselines.size(); i++)
    {
        hash.insert(BaselineHash::value_type(baselines[i].GetName(), baselines[i]) );
    }
    
    return hash;
}


/// \todo Rewrite this function to work with the new class definition.
//double Baseline::Geometric_OPD(double hour_angle, double source_declination,
//                               double array_latitude)
//{
//    double trad, drad, lrad;

//    trad = hour_angle * PI / (12.0 * 3600.);
//    drad = source_declination * PI / 180.0;
//    lrad = array_latitude * PI / 180.0;
//    Row < double >XYZ(3);

//    /// \todo This needs to be updated to use the z-coordinate of the telescopes.  See
//    /// functions above for more information.
//    XYZ[0] = -this->y * sin(lrad);
//    XYZ[1] = this->x;
//    XYZ[2] = this->y * cos(lrad);

//    double sin_alt, altitude, azimuth;

//    // Compute the altitude
//    sin_alt = cos(lrad) * cos(drad) * cos(trad) + sin(lrad) * sin(drad);
//    altitude = atan(sin_alt / sqrt(1 - sin_alt * sin_alt));

//    // Compute the azimuth
//    azimuth = atan(sin(trad) * cos(drad) / (cos(lrad) * sin(drad) - cos(trad)
//                                            * cos(drad) * sin(lrad)));

//    // Compute the unit vector
//    Row < double >s(3);

//    s[0] = cos(altitude) * cos(azimuth);
//    s[1] = cos(altitude) * sin(azimuth);
//    s[2] = sin(altitude);

//    // Return the OPD for the baseline
//    return s[0] * XYZ[0] + s[1] * XYZ[1] + s[2] * XYZ[2];
//}
