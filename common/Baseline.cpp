/// \file Baseline.cpp
/// Implements the required functions for the baseline class.

#include "Baseline.h"
#include "UVPoint.h"
#include "Array.h"
#include "Simulator.h"
#include "Station.h"

Baseline::Baseline()
{
    this->NEU[0] = this->NEU[1] = this->NEU[2] = 0;
    this->xyz[0] = this->xyz[1] = this->xyz[2] = 0;
}

Baseline::Baseline(double array_lat, Station & station1, Station & station2)
{
    // First calculate the vector offsets between the scopes in (North, East, Up) coordinates
    this->NEU[0] = station1.north - station2.north;
    this->NEU[1] = station1.east - station2.east;
    this->NEU[2] = station1.up - station2.up;
    
    double phi = array_lat * PI / 180.0;
    // Convert (North, East, Up) to (x,y,z) as defined by the APIS++ standards
    // see:    http://aips2.nrao.edu/docs/glossary
    // for more information.
    // the xyz array is ordered (x, y, z)
    this->xyz[0] = -sin(phi) * NEU[0] + cos(phi) * NEU[2];
    this->xyz[1] = NEU[1];
    this->xyz[2] = cos(phi) * NEU[0] + sin(phi) * NEU[2];
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
