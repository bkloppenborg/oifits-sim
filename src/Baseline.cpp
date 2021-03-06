/// \file Baseline.cpp
/// Implements the required functions for the baseline class.

#include <cmath>

#include "Baseline.h"
#include "UVPoint.h"
#include "Array.h"
#include "Common.h"
#include "Station.h"
#include "Target.h"

Baseline::Baseline()
{
    this->xyz[0] = this->xyz[1] = this->xyz[2] = 0;
    this->name = "";
    this->indicies[0] = 0;
    this->indicies[1] = 0;
}

/// Computes the XYZ position of a baseline composed of two stations
Baseline::Baseline(Station * station1, Station * station2)
{
    // Now calculate the xyz coordinates:
    this->xyz[0] = station1->xyz[0] - station2->xyz[0];
    this->xyz[1] = station1->xyz[1] - station2->xyz[1];
    this->xyz[2] = station1->xyz[2] - station2->xyz[2];
        
    this->name = station1->GetName() + "-" + station2->GetName();
    this->indicies[0] = station1->GetIndex();
    this->indicies[1] = station2->GetIndex();
}

/// Computes the UV coordinates for the given information.
///     hour_angle : The hour angle (decimal hours)
///     target_dec : The declination of the target (radians)
///     wavenumber : Wavenumber (1/meter)
UVPoint Baseline::UVcoords(double hour_angle, double target_declination)
{
    // First convert all values into radians (they should be degrees or decimal hours (of time) before now.
    double h = hour_angle * PI / 12;
    double delta = target_declination;

    // Now compute the UV coordinates, again according to the APIS++ standards.
    UVPoint uv = UVPoint();
    uv.u = sin(h) * xyz[0]                    + cos(h) * xyz[1];
    uv.v = -sin(delta) * cos(h) * xyz[0]      + sin(delta) * sin(h) * xyz[1]    + cos(delta) * xyz[2];
    uv.w = cos(delta) * cos(h)*xyz[0]         - cos(delta) * sin(h) * xyz[1]    + sin(delta) * xyz[2];
    
    /// \note The coordinates (u,v) calculated here appear to be (-u, -v) with respect to CHARA+MIRC data.
    /// Given that the UV plane is symmetric, this doesn't matter.
    
    return uv;
}

string Baseline::GetName(void)
{
    return this->name;
}


/// Computes the complex visibility of the target at the specified UV point.
/// \todo Switch the hash lookup over to a multihash instead of using strings?
complex<double> Baseline::GetVisibility(Target & target, double hour_angle, double wavenumber)
{
    // First look up the UV coordinates
    UVPoint uv = UVcoords(hour_angle, target.declination);
    uv.Scale(wavenumber);    
    return GetVisibility(target, uv);
}

/// Returns the complex visibility of the given target
complex<double> Baseline::GetVisibility(Target & target, UVPoint uv)
{
    string hash_key = GetHashKey(target, uv);
    complex <double> visibility(0.0, 0.0);
    
    if(mVisValues.find(hash_key) != mVisValues.end())
    {
        visibility = mVisValues[hash_key];  
    }
    else
    {
        visibility = ComputeVisibility(target, uv);
        mVisValues[hash_key] = visibility;
    }
    
    return visibility;
}


// Computes a hash key from the target, hour angle, and wavenumber.
string  Baseline::GetHashKey(Target & target, UVPoint uv)
{
    /// \todo It may be necessary for the doubles coming into this function to be cast into some 
    /// finite floating point format.
    
    /// \todo This function is in common with the Triplet class, need to factor this code.
    
    std::ostringstream sstream;
    sstream << target.GetName() << '-' << uv.HashString();
    std::string str = sstream.str();
    return str;
}

/// Computes the error in the visibility
/// \todo Implement computations of error in visibility.  This function returns an error of zero for now.
double Baseline::ComputeVis2Error(Target & target, UVPoint uv)
{
    return 0.001;
}

/// Returns the error in the visibility assoicated with this target, hour angle, and wavenumber.
double Baseline::GetVis2Error(Target & target, double hour_angle, double wavenumber)
{
     // First look up the UV coordinates
    UVPoint uv = UVcoords(hour_angle, target.declination);
    uv.Scale(wavenumber);
    return GetVis2Error(target, uv);
}

double Baseline::GetVis2Error(Target & target, UVPoint uv)
{
    string hash_key = GetHashKey(target, uv);
    double vis_error = 0.0;
    
    // First try looking up the value in the hash table:
    if(mVis2Errors.find(hash_key) != mVis2Errors.end())
    {
        vis_error = mVis2Errors[hash_key];
    }
    else
    {
        vis_error = ComputeVis2Error(target, uv);
        mVis2Errors[hash_key] = vis_error;
    }
    
    return vis_error;
}

/// Sets the visibility error for the given target, hour angle and wavenumber.
/// Note: It is intended that this function is used to set the error from existing OIFITS data files.
void    Baseline::SetVis2Error(Target & target, double hour_angle, double wavenumber, double vis2error)
{
     UVPoint uv = UVcoords(hour_angle, target.declination);
     uv.Scale(wavenumber);
     SetVis2Error(target, uv, vis2error);
}

/// Sets the visibility squared error based upon the target and UV coordiantes.
void    Baseline::SetVis2Error(Target & target, UVPoint uv, double vis2error)
{
    string hash_key = GetHashKey(target, uv);
    mVis2Errors[hash_key] = vis2error;
}

/// Computes the visibility at the specified UV point.
/// Note, uv should already be scaled by wavenumber.
complex<double> Baseline::ComputeVisibility(Target & target, UVPoint uv)
{
    complex <double> visibility(0.0, 0.0);

    // Point target
    if (target.image.GetCols() == 0)
    {
        visibility = 1.0;
    }
    else    // A resolved object
    {      
        int nx = target.image.GetRows();
        int ny = target.image.GetCols();

        /// \todo This calculation could be farmed out to a GPU very easily.
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < ny; jj++)
            {
	      visibility += target.image[ii][jj] * polar(1., -2.0 * PI * target.pixellation * milliarcsec * (uv.u * (double)(ii - nx / 2 ) + uv.v * (double)(jj - ny / 2)));
            }
        }
    }

    return visibility;
}

/// Returns an OIFITS oi_vis2_record for the given target, hour angle and wavenumber
//oi_vis2_record Baseline::GetVis2Record(Target & target, double hour_angle, <double> wavenumbers)
//{
//    // init local vars
//    complex<double> vis;
//    complex<double> vis_err;
//    double vis2;
//    double vis2_err;

//    // Finally assemble the oi_vis2_record:
//    oi_vis2_record record;
//    record.target_id = target.GetTargetID();
//    
//    /// \bug The time and MJD recorded here are nonsense values.
//    record.time = 0;
//    record.mjd = 0.0;
//    
//    /// \bug Integration time is set to 1-second by default regardless of instrument setting.
//    record.int_time = 1;
//    
//    // Get the UV point at the mid-point of the data.
//    UVPoint uv = UVcoords(hour_angle, target.declination, wavenumbers[wavenumbers.size()/2]);
//    
//    record.ucoord = uv.u;
//    record.vcoord = uv.v;
//    record.sta_index[0] = this->indicies[0];
//    record.sta_index[1] = this->indicies[1];

//    // Now write out the 
//    for(unsigned int i = 0; i < wavenumbers.size(); i++)
//    {
//        // First get the visibility and error.
//        vis = GetVisibility(target, hour_angle, wavenumbers[i]);
//        vis_err = GetVisError(target, hour_angle, wavenumbers[i]);
//        
//        // Now compute the powerspectra for this point:
//        vis2 = norm(vis);
//        vis2_err = norm(vis_err);

//        record.vis2data[i] = vis2;
//        record.vis2err[i] = vis2_err;
//        record.flag[i] = FALSE;
//    }

//    return record;
//}

// Returns the station ID of stations 0 or 1
int Baseline::GetStationID(int num)
{
    /// \todo This should probably query the telesocopes directly.
    return this->indicies[num];
}

double  Baseline::GetVis2(Target & target, double hour_angle, double wavenumber)
{
    complex<double> vis = this->GetVisibility(target, hour_angle, wavenumber);
    
    return norm(vis);
}

double Baseline::GetVis2(Target & target, UVPoint uv)
{
    complex<double> vis = this->GetVisibility(target, uv);
    return norm(vis);   
}

////////////////////////////////////////////////////////////////////
// Non Class Functions Below
////////////////////////////////////////////////////////////////////

/// Computes all possible baselines formed by the specified stations.
vector<Baseline> ComputeBaselines(vector<Station> & stations)
{
    int num_stations = stations.size();
    
    vector<Baseline> baselines;
    
    // Now compute all of the baselines and make a hash table for each baseline value
    for(int i = 0; i < num_stations; i++)
    {
        for(int j = i+1; j < num_stations; j++)
        {
            // Create a new baseline, append it to the list of baselines
            printf("Creating baseline %i %i. \n", stations[i].GetIndex(), stations[j].GetIndex());
            
            baselines.push_back(Baseline(&stations[i], &stations[j]));
        }
    }
    
    return baselines;
}

/// Computes a (baseline_name, baseline_object) hash table.
BaselineHash ComputeBaselineHash(vector<Baseline> & baselines)
{
    BaselineHash hash;
    
    for(unsigned int i = 0; i < baselines.size(); i++)
    {
        hash.insert(BaselineHash::value_type(baselines[i].GetName(), &baselines[i]) );
    }
    
    return hash;
}


/// \todo Rewrite this function to work with the new class definition.
//double Baseline::Geometric_OPD(double hour_angle, double target_declination,
//                               double array_latitude)
//{
//    double trad, drad, lrad;

//    trad = hour_angle * PI / (12.0 * 3600.);
//    drad = target_declination * PI / 180.0;
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
