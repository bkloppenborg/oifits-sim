
/// \file Observation.cpp

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>

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

Observation::Observation(Array * array, vector<Station> stations, string exclude_baselines)
{
    this->mArray = array;
    this->mStations = stations;
    this->mBaselines = this->FindBaselines(mStations, exclude_baselines);
}

// Construct an Observation object from the hour angle, and included/excluded telescopes.
Observation::Observation(Array * array, double hour_angle, string telescopes, string exclude_baselines)
{
    this->mArray = array;
    this->mHA = hour_angle;
    this->mComputeHA = false;
    
    // Now Make the baselines
    this->mStations = this->FindStations(telescopes);
    this->mBaselines = this->FindBaselines(mStations, exclude_baselines);  
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
    this->mStations = this->FindStations(telescopes);
    this->mBaselines = this->FindBaselines(mStations, exclude_baselines);    
}

/// Queries the Array object for the scopes found in the comma separated string "telescopes"
vector<Station> Observation::FindStations(string telescopes)
{
    vector<Station> stations;
    vector<string> station_names;

    // First extract the names of the telescopes from the CSV string, "telescopes"
    StringSplit(telescopes, ",", station_names);
    StripWhitespace(station_names);
    
    // Now query the array for the station objects.
    int num_stations = station_names.size();
    for(int i = 0; i < num_stations; i++)
    {
        stations.push_back( this->mArray->GetStation(station_names[i]) );
    }
    
    return stations;
}

// Finds the baselines specified in the Array object.
vector<Baseline> Observation::FindBaselines(vector <Station> stations, string exclude_baselines)
{
    vector<Baseline> baselines;
    vector<string> excluded_baselines;
    string str;
    string sta1_name;
    string sta2_name;
    string bl_name;
    
    BLNameHash bl_names;
    
    // Now compute all of the baselines and make a hash table for each baseline value
    int num_stations = stations.size();
    for(int i = 0; i < num_stations; i++)
    {
        sta1_name = stations[i].GetName();
        
        for(int j = i; j < num_stations; j++)
        {
            string sta2_name = stations[j].GetName();
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


int Observation::GetNumStations(void)
{
    return this->mStations.size();
}
Station & Observation::GetStation(int sta_index)
{
    return this->mStations[sta_index];
}

/// Reads in an properly formatted observation file in a formats defined by file_type:
///     0: A list of hour angles (in decimal hours)
///     1: A descriptive list of the observation (see ReadObservation_Descriptive() for more info)
vector<Observation> Observation::ReadObservations(Array * array, string filename, string comment_chars, int file_type)
{
    if(file_type == 0)
        return ReadObservation_HA(array, filename, comment_chars);
    else if (file_type == 1)
        return ReadObservation_Descriptive(array, filename, comment_chars);
    else
    {
        /// \exception runtime_error Invalid Observation File format Specifier
        throw std::runtime_error("Invalid Observation File format Specifier.");
    }
}

/// Reads in a file that consists of lines of hour angles with or without comments.
vector <Observation> Observation::ReadObservation_HA(Array * array, string filename, string comment_chars)
{
    vector<Observation> observations;
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Observation Definition File");
    
    // Now parse the file, make the observations
    double ha;
    for(int i = 0; i < lines.size(); i++)
    {
        if (!isdigit(lines[i][0]))
        {
            throw std::runtime_error("Non numeric character found in observation entry " + i);
        }

        ha = atof(lines[i].c_str());
        
        // Make a new observation with all of the stations included.
        observations.push_back( Observation(array, array->GetAllStations(), "") );
    }
    
}

/// Reads in a series of formatted lines that consist of keywords:
///     hour_angle = DECIMAL NUMBER
///     telescopes = STATION_NAMES
///     exclude = EXCLUDED_BASELINES
/// DECIMAL NUMBER is a decimal number denoting the hour angle of the observation.
/// STATION_NAMES is a CSV list of stations included in this observation in the following format:
///     S1, S2, ... , E1
/// white space is automatically stripped between commas.
/// EXCLUDED_BASELINES is a CSV list of excluded baselines consisting of two station names separated
/// by a hyphen, e.g.:
///     S1-S2, S2-E2
/// with the lower station index (as defined in the array definition file) always appearing first.
/// white space is automatically stripped between commas.  EXCLUDED_BASELINES may be blank.
vector <Observation> Observation::ReadObservation_Descriptive(Array * array, string filename, string comment_chars)
{
    vector <Observation> observations;

    vector < string > lines = ReadFile(comment_chars, filename, "Cannot Open Observation Definition File");
    
    
    return observations;
}
