/// \file Baseline.cpp
/// Implements the required functions for the baseline class.

#include "Baseline.h"
#include "UVPoint.h"
#include "Array.h"
#include "Simulator.h"

Baseline::Baseline(Array & array, int index_station1, int index_station2)
{
    this->x = array.x[index_station1] - array.x[index_station2];
    this->y = array.y[index_station1] - array.y[index_station2];
    this->z = array.z[index_station1] - array.z[index_station2];
}

UVPoint Baseline::UVcoords(double time, double source_declination,
                           double array_latitude, double wavenumber)
{
    // IN : baseline , OUT : corresponding uv point
    // B.x is East coordinate
    // B.y is North coordinate
    // float H : hour angle, converted into radians from time in seconds
    // float dec: declination of the source (in degrees)
    // float lat: latitute of the array (in degrees)
    // float wav: wavelength (meters)
    // Conversion degrees-> radians or seconds -> radians

    // First convert all values into radians (they should be degrees or seconds (of time) before now
    double phi = array_latitude * PI / 180.0;
    double h = time * PI / (12.0 * 3600.);
    double delta = source_declination * PI / 180.0;

    UVPoint uv;
    
    // First array is in (north, east, up) coordinates
    Row < double >A(3);
    A[0] = this->y; // North
    A[1] = this->x; // East
    A[2] = this->z; // Up (elevation)
    
    // Convert this to (x,y,z) as defined by the APIS++ standards
    // see:    http://aips2.nrao.edu/docs/glossary
    // for more information.
    Row < double>B(3);
    B[0] = -sin(phi) * A[0] + cos(phi) * A[2];    // X
    B[1] = A[1];                                  // Y
    B[2] = cos(phi) * A[0] + sin(phi) * A[2];     // Z

    
    // Now compute the UV coordinates, again according to the APIS++ standards.
    uv.u = sin(h) * B[0]                    + cos(h) * B[1];
    uv.v = -sin(delta) * cos(h) * B[0]      + sin(delta) * sin(h) * B[1]    + cos(delta) * B[2];
    uv.w = cos(delta) * cos(h)*B[0]         - cos(delta) * sin(h) * B[1]    + sin(delta) * B[2];
    
    // Scale by the wavenumber (1/wavelength)
    /// \note The coordinates (u,v) calculated here appear to be (-u, -v) with respect to CHARA+MIRC data.
    /// Given that the UV plane is symmetric, this doesn't matter.
    uv.u *= wavenumber;
    uv.v *= wavenumber;
    uv.w *= wavenumber;
    
    return uv;
}

double Baseline::Geometric_OPD(double hour_angle, double source_declination,
                               double array_latitude)
{
    double trad, drad, lrad;

    trad = hour_angle * PI / (12.0 * 3600.);
    drad = source_declination * PI / 180.0;
    lrad = array_latitude * PI / 180.0;
    Row < double >XYZ(3);

    /// \todo This needs to be updated to use the z-coordinate of the telescopes.  See
    /// functions above for more information.
    XYZ[0] = -this->y * sin(lrad);
    XYZ[1] = this->x;
    XYZ[2] = this->y * cos(lrad);

    double sin_alt, altitude, azimuth;

    // Compute the altitude
    sin_alt = cos(lrad) * cos(drad) * cos(trad) + sin(lrad) * sin(drad);
    altitude = atan(sin_alt / sqrt(1 - sin_alt * sin_alt));

    // Compute the azimuth
    azimuth = atan(sin(trad) * cos(drad) / (cos(lrad) * sin(drad) - cos(trad)
                                            * cos(drad) * sin(lrad)));

    // Compute the unit vector
    Row < double >s(3);

    s[0] = cos(altitude) * cos(azimuth);
    s[1] = cos(altitude) * sin(azimuth);
    s[2] = sin(altitude);

    // Return the OPD for the baseline
    return s[0] * XYZ[0] + s[1] * XYZ[1] + s[2] * XYZ[2];
}
