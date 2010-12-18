
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
#include "Obs_OIFITS.h"

#include "Simulator.h"
#include "Source.h"
#include "UVPoint.h"


/// Default constructor.
Observation::Observation()
{

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


/// Reads in an properly formatted observation file in a formats defined by file_type:
///     0: A list of hour angles (in decimal hours)
///     1: A descriptive list of the observation (see ReadObservation_Descriptive() for more info)
vector<Observation*> Observation::ReadObservations(Array * array, string filename, string comment_chars, ObsType file_type)
{
    if(file_type == HOUR_ANGLE)
        return Obs_HA::ReadObservation_HA(array, filename, comment_chars);
    else if (file_type == DESCRIPTIVE)
        return Obs_HA::ReadObservation_Descriptive(array, filename, comment_chars);
    else if(file_type == OIFITS)
        return Obs_OIFITS::ReadObservation_OIFITS(array, filename);
    else
    {
        /// \exception runtime_error Invalid Observation File format Specifier
        throw std::runtime_error("Invalid Observation File format Specifier.");
    }
}

// A simple function to see if the observation has triplets.
bool    Observation::HasTriplets(void)
{
    if(this->mTriplets.size() > 0)
        return true;
        
    return false;
}
