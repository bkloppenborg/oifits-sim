
/// \file Observation.cpp

#include <cmath>
#include <cstdio>

#include "Observation.h"
#include "Simulator.h"

/// The minimum timespan
const double SAMEOBS = 5 / (60 * 24); // 5 minutes in days.

inline bool SameObservation(Observation & A, Observation & B)
{
    // Because observations can be specified by either JDs or hour angles, we need to allow for both cases here.
    if(A.mComputeHA && B.mComputeHA)
    {
        // Both use Julain Dates
        if(fabs(A.mJD - B.mJD) < SAMEOBS)
            return true;
    }
    else if(!A.mComputeHA && !B.mComputeHA)
    {
        // Both use hour angles:
        if(fabs(A.mHA - B.mHA) < SAMEOBS)
            return true;
    }
    
    // Uh oh, mixed observation types... this is bad.   
    /// \todo Re-enable this exception.
    //throw std::runtime_error("Warning: Mixed Observation Types (JD and Hour Angles) Detected!");
    
    return false;
}

// Construct an Observation object from the hour angle, and included/excluded telescopes.
Observation::Observation(Array * array, double hour_angle, string telescopes, string exclude_baselines)
{
    this->mArray = array;
    this->mHA = hour_angle;
    this->mComputeHA = false;
    
    // Now Make the baselines
    this->mBaselines = this->FindBaselines(telescopes, exclude_baselines);
}

/// Construct an Observation object from the MJD, time, and included/excluded telescopes.
Observation::Observation(Array * array, double MJD, double time, string telescopes, string exclude_baselines)
{
    this->mArray = array;
    this->mHA = 0;
    // Compute the (full) Julian date from the Modified Julian Date (MJD)
    this->mJD = MJD + time / 24 + 2400000.5;
    this->mComputeHA = true;
    
    // Now Make the baselines
    this->mBaselines = this->FindBaselines(telescopes, exclude_baselines);    
}

// Finds the baselines specified in the Array object.
vector<Baseline> Observation::FindBaselines(string telescopes, string exclude_baselines)
{
    vector<Baseline> baselines;
    vector<string> station_names;
    vector<string> excluded_baselines;
    string str;
    string sta1_name;
    string sta2_name;
    string bl_name;
    
    BLNameHash bl_names;
    
    // First extract the names of the telescopes from the CSV string, "telescopes"
    StringSplit(telescopes, ",", station_names);
    StripWhitespace(station_names);
    
    // Now compute all of the baselines and make a hash table for each baseline value
    int num_stations = station_names.size();
    for(int i = 0; i < num_stations; i++)
    {
        sta1_name = station_names[i];
        
        for(int j = i; j < num_stations; j++)
        {
            string sta2_name = station_names[j];
            bl_name = sta1_name + "-" + sta2_name;
            
            bl_names.insert( BLNameHash::value_type(bl_name, bl_name) );
        }   
    }
    
    // Now parse out the baselines that should be excluded.
    StringSplit(exclude_baselines, ",", excluded_baselines);
    StripWhitespace(excluded_baselines);
    
    // Remove the excluded baseline names from the hash
    int num_excluded = excluded_baselines.size();
    for(int i = 0; i < num_excluded; i++)
    {
        bl_names.erase(bl_names[excluded_baselines[i]]);
    }
    
    // Now query for all of the remaining baselines and add them to the output "baselines"
    for (BLNameHash::iterator it = bl_names.begin(); it != bl_names.end(); ++it)
    {
        bl_name = it->second;
        baselines.push_back( this->mArray->GetBaseline(bl_name) );
    }
    
    return baselines;
}

double Observation::GetHA(double targ_ra)
{
    // If the hour angle is already specified, just reutrn it directly.
    if(!this->mComputeHA)
        return this->mHA;
    
    double lst = this->GetLocalSiderealTime(this->mJD, 0, 0, this->mArray->GetLongitude());
    return lst - targ_ra * 24.0/360;
}

double Observation::GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long)
{
    double lst = this->GetSiderealTime(jd_high, jd_low, ee) + array_long * 24 / 360;
    
	if(lst < 0)
        lst += 24;
        
    return lst;
}

// Get the Sidereal time as a double given the High and Low Julian Date
double  Observation::GetSiderealTime(double jd_high, double jd_low, double ee)
{
	// Code slightly modified from NOVAS to use the current SystemTime
	// Naval Observatory Vector Astrometry Subroutines (C Language Version 2.0)
	const double T0 = 2451545.00000000;
	double t_hi = 0;
	double t_lo = 0;
	double t = 0;
	double t2 = 0;
	double t3 = 0;
	double st = 0;
	double gst = 0;

	t_hi = (jd_high -  T0) / 36525.0;
	t_lo = jd_low / 36525.0;
	t = t_hi + t_lo;
	t2 = t * t;
	t3 = t2 * t;

	st =  ee - 6.2e-6 * t3 + 0.093104 * t2 + 67310.54841 + 8640184.812866 * t_lo  + 3155760000.0 * t_lo + 8640184.812866 * t_hi + 3155760000.0 * t_hi;

	gst = fmod ((st / 3600.0), 24.0);

	if (gst < 0.0)
		gst += 24.0;

	return gst;
}
