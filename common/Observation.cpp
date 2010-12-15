
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
#include "Source.h"
#include "UVPoint.h"

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

Observation::Observation(Array * array, vector<Station*> stations, string exclude_baselines)
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
    this->mTriplets = this->FindTriplets(mStations, exclude_baselines);
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
    this->mTriplets = this->FindTriplets(mStations, exclude_baselines);  
}

/// Queries the Array object for the scopes found in the comma separated string "telescopes"
vector<Station*> Observation::FindStations(string telescopes)
{
    vector<Station*> stations;
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
vector<Baseline*> Observation::FindBaselines(vector <Station*> stations, string exclude_baselines)
{
    vector<Baseline*> baselines;
    vector<string> excluded_baselines;
    string str;
    string sta1_name;
    string sta2_name;
    string bl_name;
    
    BLNameHash bl_names;
    
    // Now compute all of the baselines and make a hash table for each baseline value
    int num_stations = stations.size();
    int i, j;
    for(i = 0; i < num_stations; i++)
    {
        sta1_name = stations[i]->GetName();
        
        for(j = i + 1; j < num_stations; j++)
        {
            string sta2_name = stations[j]->GetName();
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

vector<Triplet*>    Observation::FindTriplets(vector<Station*> stations, string exclude_baselines)
{
    vector<Triplet*> triplets;
    vector<string> excluded_baselines;
    string sta1_name;
    string sta2_name;
    string sta3_name;
    string tri_name;
    string bl_name;
    
    TNameHash tri_names;
    
    // Start by first computing all of the names of all of the possible triplets
    int num_stations = stations.size();
    int i, j, k;
    for(i = 0; i < num_stations - 2; i++)
    {
        sta1_name = stations[i]->GetName();
        
        for(j = i + 1; j < num_stations - 1; j++)
        {
            sta2_name = stations[j]->GetName();
           
            for(k = j + 1; k < num_stations; k++)
            {
                sta3_name = stations[k]->GetName();
                tri_name = sta1_name + "-" + sta2_name + "-" + sta3_name;
                
                // Insert the triplet name into the hash.
                tri_names.insert( TNameHash::value_type(tri_name, tri_name) );
            }
        }
    }
    
    // Now parse out the baselines that should be excluded.
    StringSplit(exclude_baselines, ",", excluded_baselines);
    StripWhitespace(excluded_baselines);
    
    // Now query for all of the remaining baselines and add them to the output "baselines"
    for (TNameHash::iterator it = tri_names.begin(); it != tri_names.end(); ++it)
    {
        tri_name = it->second;
        
        for(unsigned int j = 0; j < excluded_baselines.size(); j++)
        {
            bl_name = excluded_baselines[j];
            
            if(this->mArray->GetTriplet(tri_name)->ContainsBaseline(bl_name))
            {   
                // The excluded baseline is contained in this triplet, remove it from the hash
                tri_names.erase(it);
                
                // Useful text if you want to see what triplets are removed.
                //printf("Excluding triplet %s due to baseline %s \n", tri_name.c_str(), bl_name.c_str());
                
                // Now that we have found a match, invalidate the loop over baselines (index j)
                // so that we inspect the next triplet.
                j = excluded_baselines.size();
            }
        
        }
    } 
    
    // Now grab the triplets from the array and put them into the triplets variable.
    for (TNameHash::iterator it = tri_names.begin(); it != tri_names.end(); ++it)
    {
        tri_name = it->second;  
        triplets.push_back( this->mArray->GetTriplet(tri_name) );
    }  
    
    return triplets;
}   

double Observation::GetHA(double targ_ra)
{
    // If the hour angle is already specified, just reutrn it directly.
    if(!this->mComputeHA)
        return this->mHA;
    
    double lst = this->GetLocalSiderealTime(this->mJD, 0, 0, this->mArray->GetLongitude());
    return lst - targ_ra * 12.0/PI;
}

double Observation::GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long)
{
    double lst = this->GetSiderealTime(jd_high, jd_low, ee) + array_long * 12 / PI;
    
	if(lst < 0)
        lst += 24;
        
    return lst;
}

// Get the Sidereal time as a double given the High and Low Julian Date
double  Observation::GetSiderealTime(double jd_high, double jd_low, double ee)
{
	// Code slightly modified from Naval Observatory Vector Astrometry Subroutines 
	// (C Language Version 2.0)
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
Station * Observation::GetStation(int sta_index)
{
    return this->mStations[sta_index];
}

/// Creates an OIFITS oi_vis2 compliant entry for this observation.
oi_vis2 Observation::GetVis2(string ins_name, Source & source, vector<double> & wavenumbers)
{
    // init local vars
    oi_vis2 vis2;
    int npow = this->mBaselines.size();
    int nwave = wavenumbers.size();
    int iwave = 0;
    string arrname = this->mArray->GetArrayName();
    double wavenumber = 0;
    
    double ra = source.right_ascension;
    double dec = source.declination;
    UVPoint uv;
    
    // Allocate room for the vis2 records and their data.
	vis2.record = (oi_vis2_record *) malloc(npow * sizeof(oi_vis2_record));
	for (int i = 0; i < npow; i++)
	{
		vis2.record[i].vis2data = (double *) malloc(nwave * sizeof(double));
		vis2.record[i].vis2err = (double *) malloc(nwave * sizeof(double));
		vis2.record[i].flag = (char *) malloc(nwave * sizeof(char));
	}
	
	vis2.revision = 1;
	/// \bug The observation date is set to all zeros by default.  
	/// This is to ensure the user knows this is simulated data, but may not be compliant
	/// with the OIFITS format, or good "note taking"
	strncpy(vis2.date_obs, "0000-00-00", 11);
	strncpy(vis2.arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(vis2.insname, ins_name.c_str(), FLEN_VALUE);
	vis2.numrec = npow;
	vis2.nwave = nwave;
	for (int i = 0; i < npow; i++)
	{
		vis2.record[i].target_id = source.GetTargetID();
		/// \bug The time is set to zero by default 
		vis2.record[i].time = 0.0;
		vis2.record[i].mjd = this->mJD;
		/// \bug The integration time is set to 10 seconds by default.
		vis2.record[i].int_time = 10;
		
		// Compute the UV coordinates and record the station positions:
		uv = this->mBaselines[i]->UVcoords(this->GetHA(ra), dec);
		vis2.record[i].ucoord = uv.u;
		vis2.record[i].vcoord = uv.v;
		vis2.record[i].sta_index[0] = this->mBaselines[i]->GetStationID(0);
		vis2.record[i].sta_index[1] = this->mBaselines[i]->GetStationID(1);
		
		// Now compute the individual visibilities and uncertainties
		for (iwave = 0; iwave < nwave; iwave++)
		{
		    // look up the present wavenumber, and then find the data
		    wavenumber = wavenumbers[iwave];
			vis2.record[i].vis2data[iwave] = this->mBaselines[i]->GetVis2(source, mHA, wavenumber);
			vis2.record[i].vis2err[iwave] = this->mBaselines[i]->GetVis2Err(source, mHA, wavenumber);
			vis2.record[i].flag[iwave] = FALSE;
		}
	}
	
	return vis2;
}

/// Create an OIFITS-compliant t3 table for this observation.
oi_t3 Observation::GetT3(string ins_name, Source & source, vector<double> & wavenumbers)
{
    oi_t3 t3;
    int nTriplets = this->mTriplets.size();
    int nwave = wavenumbers.size();
    string arrname = this->mArray->GetArrayName();
    complex<double> bis;
    complex<double> bis_err;
    double wavenumber = 0.0;
    
    UVPoint uv_AB;
    UVPoint uv_BC;
    
	t3.record = (oi_t3_record *) malloc(nTriplets * sizeof(oi_t3_record));
	for (int i = 0; i < nTriplets; i++)
	{
		t3.record[i].t3amp = (double *) malloc(nwave * sizeof(double));
		t3.record[i].t3amperr = (double *) malloc(nwave * sizeof(double));
		t3.record[i].t3phi = (double *) malloc(nwave * sizeof(double));
		t3.record[i].t3phierr = (double *) malloc(nwave * sizeof(double));
		t3.record[i].flag = (char *) malloc(nwave * sizeof(char));
	}
	t3.revision = 1;
	/// \bug Observation date is set to 0000-00-00 by default
	strncpy(t3.date_obs, "0000-00-00", FLEN_VALUE);
	strncpy(t3.arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(t3.insname, ins_name.c_str(), FLEN_VALUE);
	t3.numrec = nTriplets;
	t3.nwave = nwave;
	
	// Now copy the data into t3 records:
	for (int i = 0; i < nTriplets; i++)
	{
	   
		t3.record[i].target_id = source.GetTargetID();
		/// \bug Time is set to zero by default.
		t3.record[i].time = 0.0;
		t3.record[i].mjd = this->mJD;
		/// \bug Integration time set to 10 seconds by default.
		t3.record[i].int_time = 10;
		
		// Get the UV coordinates for the AB and BC baselines
		uv_AB = mTriplets[i]->GetBaseline(0)->UVcoords(this->mHA, source.declination);
		uv_BC = mTriplets[i]->GetBaseline(1)->UVcoords(this->mHA, source.declination);
		
		t3.record[i].u1coord = uv_AB.u;
		t3.record[i].v1coord = uv_AB.v;
		t3.record[i].u2coord = uv_BC.u;
		t3.record[i].v2coord = uv_BC.v;
		t3.record[i].sta_index[0] = mTriplets[i]->GetStationID(0);
		t3.record[i].sta_index[1] = mTriplets[i]->GetStationID(1);
		t3.record[i].sta_index[2] = mTriplets[i]->GetStationID(2);
		
		for(int j = 0; j < nwave; j++)
		{
		    bis = mTriplets[i]->GetBispectra(source, this->mHA, wavenumbers[j]);
		    bis_err = mTriplets[i]->GetBisError(source, this->mHA, wavenumbers[j]);
			t3.record[i].t3amp[j] = abs(bis);
			t3.record[i].t3phi[j] = arg(bis) * 180 / PI;
			t3.record[i].t3amperr[j] = abs(bis_err);
			t3.record[i].t3phierr[j] = arg(bis_err) * 180 / PI;
			t3.record[i].flag[j] = FALSE;
		}
		
	}
		
	return t3;
}



/// Reads in an properly formatted observation file in a formats defined by file_type:
///     0: A list of hour angles (in decimal hours)
///     1: A descriptive list of the observation (see ReadObservation_Descriptive() for more info)
vector<Observation> Observation::ReadObservations(Array * array, string filename, string comment_chars, ObsFileType file_type)
{
    if(file_type == HOUR_ANGLE)
        return ReadObservation_HA(array, filename, comment_chars);
    else if (file_type == DESCRIPTIVE)
        return ReadObservation_Descriptive(array, filename, comment_chars);
    else if(file_type == OIFITS)
        return ReadObservation_OIFITS(array, filename);
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
    for(unsigned int i = 0; i < lines.size(); i++)
    {
        if (!isdigit(lines[i][0]))
        {
            throw std::runtime_error("Non numeric character found in observation entry " + i);
        }

        ha = atof(lines[i].c_str());
        
        // Make a new observation with all of the stations included.
        observations.push_back( Observation(array, array->GetAllStations(), "") );
    }
    
    return observations;
}

/// Reads in a series of formatted lines that consist of keywords:
///     hour_angle = DECIMAL NUMBER
///     stations = STATION_NAMES
///     exclude = EXCLUDED_BASELINES
/// DECIMAL NUMBER is a decimal number denoting the hour angle of the observation.
/// STATION_NAMES is a CSV list of stations included in this observation in the following format:
///     S1, S2, ... , E1
/// white space is automatically stripped between commas.
/// EXCLUDED_BASELINES is a CSV list of excluded baselines consisting of two station names separated
/// by a hyphen, e.g.:
///     S1-S2, S2-E2
/// with the lower station index (as defined in the array definition file) always appearing first.
/// white space is automatically stripped between commas.  EXCLUDED_BASELINES may be blank or omitted.
vector <Observation> Observation::ReadObservation_Descriptive(Array * array, string filename, string comment_chars)
{
    // init local vars
    vector <Observation> observations;
    vector<string> results;
    
    // local vars for an observation
    double hour_angle = 0;
    string stations = "";
    string excluded_baselines = "";
    
    // Flags to keep track of whether or not we have the corresponding data
    bool bStations = false;
    bool bHourAngle = false;
    bool bExclude = false;

    // Get the non-comment, non-blank lines
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Observation Definition File");
    
    unsigned int j = 1;  // A counter for the number of observations.  Only used in error messages below.
    // Iterate over the lines in the file
    for(unsigned int i = 0; i < lines.size(); i++)
    {
        // First split the line and strip out white space
        results.clear();
        StringSplit(lines[i], "=", results);
        StripWhitespace(results);
        
        if(results[0] == "hour_angle")
        {
            // We allow for the "exclude" parameter to be omitted.  If we have already found
            // a "hour_angle" flag and a "stations" flag we need to write out an observation block
            if(bHourAngle && bStations)
            {
                // Set the excluded stations to be blank, mark it as found so we will write out an observation
                excluded_baselines = "";
                bExclude = true;
                
                // Step back one in "i" so at we may revisit this line
                i--;
            }
            else if(bHourAngle && !(bStations))
            {
                /// \todo Perhaps issue a warning, but keep going?
                // We've found another an additional hour angle, but no stations.
                // The file is not formatted correctly
                throw std::runtime_error("Found an 'hour_angle' specifier without stations in observation block " + j);
            }
            else
            {
                try
                {
                    hour_angle = atof(results[1].c_str());
                    bHourAngle = true;
                }
                catch (const std::exception&)
                {
                    throw std::runtime_error("Hour angle specification is not numeric in observation block " + j);
                }
            }
        }
        else if(results[0] == "stations")
        {

            if(results[1].size() < 2)
                throw std::runtime_error("Two few stations specified in observation entry " + j);
            
            stations = results[1];
            bStations = true;
        }
        else if(results[0] == "exclude")
        {
            excluded_baselines = results[1];
            bExclude = true;
        }
        else
        {
            // This is an unknown type, thrown an exception.
            throw std::runtime_error("Unknown parameter specified in observation block " + j);
        }
        
        if(bHourAngle && bStations && (bExclude || i == lines.size() - 1))
        {
            // Push the observation on to the list
            observations.push_back( Observation(array, hour_angle, stations, excluded_baselines) );
            
            // Increment the observation block counter, j
            j++;
            
            // Reset the station, hour angle, and exclude flags.
            bStations = false;
            bHourAngle = false;
            bExclude = false;
        }
    }
    
    return observations;
}


/// Reads an OIFITS data file and creates a series of observation based upon the data.
/// Note, presently this function only supports ONE array, combiner, spectral mode per OIFITS file.
vector <Observation> Observation::ReadObservation_OIFITS(Array * array, string filename)
{
    vector<Observation> observations;
    
    
    return observations;
}

bool    Observation::HasTriplets(void)
{
    if(this->mTriplets.size() > 0)
        return true;
        
    return false;
}
