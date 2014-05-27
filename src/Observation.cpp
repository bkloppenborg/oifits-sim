
/// \file Observation.cpp

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>

#include "Observation.h"
#include "Obs_HA.h"
//#include "Obs_OIFITS.h"

#include "Common.h"
#include "Target.h"
#include "UVPoint.h"
#include "Array.h"


/// Default constructor.
Observation::Observation()
{
    mbHasTriplets = false;
    mbHasQuadruplets = false;
}

Observation::~Observation()
{

}

/// Queries the Array object for the scopes found in the comma separated string "telescopes"
vector<Station*> Observation::FindStations(string telescopes)
{
    vector<Station*> stations;
    vector<string> station_names;

    // First extract the names of the telescopes from the CSV string, "telescopes"
    station_names = SplitString(telescopes, ',');
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
    excluded_baselines = SplitString(exclude_baselines, ',');
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
    excluded_baselines = SplitString(exclude_baselines, ',');
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


vector<Quadruplet*>    Observation::FindQuadruplets(vector<Station*> stations, string exclude_baselines)
{
    vector<Quadruplet*> quadruplets;
    vector<string> excluded_baselines;
    string sta1_name;
    string sta2_name;
    string sta3_name;
    string sta4_name;
    string quad_name;
    string bl_name;
    
    TNameHash quad_names;
    
    // Start by first computing all of the names of all of the possible quadruplets
    int num_stations = stations.size();
    int i, j, k, l;

    for(i = 0; i < num_stations - 3; i++)
    {
        sta1_name = stations[i]->GetName();
        
        for(j = i + 1; j < num_stations - 2; j++)
        {
            sta2_name = stations[j]->GetName();
           
            for(k = j + 1; k < num_stations - 1 ; k++)
            {
                sta3_name = stations[k]->GetName();

				for(l = k + 1; l < num_stations; l++)
				{
					sta4_name = stations[l]->GetName();

					quad_name = sta1_name + "-" + sta2_name + "-" + sta3_name + "-" + sta4_name;

					// Insert the quadruplet name into the hash.
					quad_names.insert( TNameHash::value_type(quad_name, quad_name) );
				}
            }
        }
    }
    
    // Now parse out the baselines that should be excluded.
    excluded_baselines = SplitString(exclude_baselines, ',');
    StripWhitespace(excluded_baselines);
    
    // Now query for all of the remaining baselines and add them to the output "baselines"
    for (TNameHash::iterator it = quad_names.begin(); it != quad_names.end(); ++it)
    {
        quad_name = it->second;
        
        for(unsigned int j = 0; j < excluded_baselines.size(); j++)
        {
            bl_name = excluded_baselines[j];
            
            if(this->mArray->GetQuadruplet(quad_name)->ContainsBaseline(bl_name))
            {   
                // The excluded baseline is contained in this quadruplet, remove it from the hash
                quad_names.erase(it);
                
                // Useful text if you want to see what quadruplets are removed.
                //printf("Excluding quadruplet %s due to baseline %s \n", quad_name.c_str(), bl_name.c_str());
                
                // Now that we have found a match, invalidate the loop over baselines (index j)
                // so that we inspect the next quadruplet.
                j = excluded_baselines.size();
            }
        
        }
    } 
    
    // Now grab the quadruplets from the array and put them into the quadruplets variable.
    for (TNameHash::iterator it = quad_names.begin(); it != quad_names.end(); ++it)
    {
        quad_name = it->second;  
        quadruplets.push_back( this->mArray->GetQuadruplet(quad_name) );
    }  
    
    return quadruplets;
}   






int Observation::GetNumStations(void)
{
    return this->mStations.size();
}

Station * Observation::GetStation(int sta_index)
{
    return this->mStations[sta_index];
}

ObsType Observation::GetObsType(void)
{
    return this->mObsType;
}

vector <Observation*> Observation::ParseCommandLine(Array * array, char *argv[], int i, int argc, string comment_chars)
{
	// Parsing the command line is a little funny as we permit some observations to
	// be defined on the command line OR we can get a file.

	// This first case is for a single observation defined on the command line:
	if(strcmp(argv[i+1], "--start") == 0 || strcmp(argv[i+1], "--end") == 0 || strcmp(argv[i+1], "--every") == 0 || strcmp(argv[i+1], "--T") == 0)
	{
		// The user is specifying things on the command line.  Parse those
		return ParseCommandLineObs(array, argv, i, argc);
	}
	else
	{
		// we should be expecting a filename
		return ImportFile(array, string(argv[i+1]), comment_chars);
	}
}

vector <Observation*> Observation::ParseCommandLineObs(Array * array, char *argv[], int i, int argc)
{
	double start;
	double end;
	double every;
	string telescopes;

	bool have_all[3]; // for start, end, and every
	for(int j = 0; j < 3; j++)
		have_all[j] = false;


	for(int j = i+1; j < argc; j++)
	{
		// TODO: Add in some way to jump out of this loop before we hit argc otherwise
		// we're doing to multiply parse the observations.

		if ((strcmp(argv[j], "--start") == 0) && (j < argc - 1))
		{
			try
			{
				start = atof(argv[j+1]);
				have_all[0] = true;
			}
			catch(...)
			{
    			throw std::runtime_error("Observation start time on command line isn't a number!");
			}
		}

		if ((strcmp(argv[j], "--end") == 0) && (j < argc - 1))
		{
			try
			{
				end = atof(argv[j+1]);
				have_all[1] = true;
			}
			catch(...)
			{
    			throw std::runtime_error("Observation end time on command line isn't a number!");
			}
		}

		if ((strcmp(argv[j], "--every") == 0) && (j < argc - 1))
		{
			try
			{
				every = atof(argv[j+1]);
				have_all[2] = true;
			}
			catch(...)
			{
    			throw std::runtime_error("Observation interval on command line isn't a number!");
			}
		}

		if ((strcmp(argv[j], "--T") == 0) && (j < argc - 1))
		{
			try
			{
				telescopes = string(argv[j+1]);
			}
			catch(...)
			{
    			throw std::runtime_error("Telescope definition on command line isn't a string!");
			}
		}
	}

	if(have_all[0] && have_all[1] && have_all[2])
	{
		return Obs_HA::MakeObservations(array, start, end, every, telescopes);
	}
	else
	{
		throw std::runtime_error("Not all observation parameters were specified on the command line.  You need --start, --stop, and --every at a minimum.");
	}

}

/// Reads in an properly formatted observation file in a formats defined by file_type:
///     0: A list of hour angles (in decimal hours)
///     1: A descriptive list of the observation (see ReadObservation_Descriptive() for more info)
vector<Observation*> Observation::ImportFile(Array * array, string filename, string comment_chars)
{
	// First determine the type of observation and fork it off to the
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Observation Definition File");
	vector <string> results;

	int obs_type = -1;
	unsigned int i = 0;

	for(i = 0; i < lines.size(); i++)
	{
		// Clear out the results, split the string and strip whitespace
        results.clear();
        results = SplitString(lines[i], '=');
        StripWhitespace(results);

        if(results[0] == "type")
        {
        	try
        	{
        		obs_type = atoi(results[1].c_str());
        		break;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid observation type field in observation file.");
        	}
        }
	}

	// This function only reads in two types of observations
    if(obs_type == HOUR_ANGLE)
        return Obs_HA::ReadObservation_HA(array, lines, i);
    else if (obs_type == DESCRIPTIVE)
        return Obs_HA::ReadObservation_Descriptive(array, lines, i);
    else
    {
        /// \exception runtime_error Invalid Observation File format Specifier
        throw std::runtime_error("Invalid Observation type field in observation file.");
    }
}

// A simple function to see if the observation has triplets.
bool    Observation::HasTriplets(void)
{
    return this->mbHasTriplets;
}

// A simple function to see if the observation has quadruplets.
bool    Observation::HasQuadruplets(void)
{
    return this->mbHasQuadruplets;
}
