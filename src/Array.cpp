// AS 2010-06-22
// modified this constructor to work with the new parameters added to the
// Array class on 2010-06-18

// First define all of the built-in includes:
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <iostream>
#include <fstream>
#include <stdexcept>


// Now custom includes
#include "Array.h"
#include "Station.h"
#include "Baseline.h"
//#include "Simulator.h"  // Required for constants, like PI.
#include "Triplet.h"
#include "Quadruplet.h"
#include "textio.hpp"
#include "Common.h"

Array::Array()
{

}

// Destructor for the array object.
Array::~Array()
{

}

// Create and initialize an array from data given in a text file.
void Array::ImportFile(string filename, string comment_chars)
{
	// Some local variables
	int max_params = 8; // This is the number of non-telescope parameters found in the array definition file.
	int n_params = 0; // the number of parameters read in from the file

    // Local varaibles for creating Station objects.
    string staname;
    int sta_index;
    double North;
    double East;
    double Up;
    double gain;
    double diameter;
    int coord_sys;
    bool xyz_coords = false;

    int line_number;

    // We permit shortcuts in the name of the arrays, this lets us do
    // $ oifits-sim -a CHARA
    // instead of
    // $ oifits-sim -a /home/.../etc/CHARA.txt
    // which I think is a nice feature.  To do this, we first need to see if a file with the
    // name, filename, really exists.
    // TODO: Use a global variable, read in from a configuration file, to specify "../etc" below.
    if(!FileExists(filename))
    {
    	filename = "../etc/" + filename + ".txt";
    }

    // stores non-blank, non-comment lines
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Array Definition File");
	vector <string> results;

	for(line_number = 0; line_number < max_params; line_number++)
	{
		// Clear out the results, split the string and strip whitespace
        results.clear();
        results = SplitString(lines[line_number], '=');
        StripWhitespace(results);


        if(results[0] == "name")
        {
        	try
        	{
        		this->arrayname = string(results[1]);
        		n_params += 1;
        		cout << "Loading Array " << this->arrayname << endl;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid array name");
        	}
        }

        if(results[0] == "lat")
        {
        	try
        	{
        		this->latitude = atof(results[1].c_str());

        		// Convert to radians:
                this->latitude *= PI / 180;
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid array latitude");
        	}
        }

        if(results[0] == "lon")
        {
        	try
        	{
        		this->longitude = atof(results[1].c_str());

        		// Convert to radians:
                this->longitude *= PI / 180;
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid array longitude");
        	}
        }

        if(results[0] == "alt")
        {
        	try
        	{
        		this->altitude = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid array altitude");
        	}
        }

        if(results[0] == "coord")
        {
        	try
        	{
                coord_sys = atoi(results[1].c_str());
                if(coord_sys == 0)
                {
                    cout << "Array Coordinate system specified in XYZ." << endl;
                    xyz_coords = true;
                }
                else
                    cout << "Array Coordinate system specified in (North, East, Up)." << endl;

        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid Coordinate definition specification for telescope positions.");
        	}
        }

        if(results[0] == "throughput")
        {
        	try
        	{
        		this->throughput = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid array throughput");
        	}
        }

        if(results[0] == "wind_speed")
        {
        	try
        	{
        		// Import the wind speed.  In the input parameter file it is in m/s
        		// all variables are stored in MKS units, so no conversion is necessary.
        		this->wind_speed = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid array average wind speed");
        	}
        }

        if(results[0] == "r0")
        {
        	try
        	{
        		this->r0 = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid average r0 for the array");
        	}
        }
	}

	// Do some quick error checking on the file format.  First, make sure we've got all of the parameters:
	if(n_params != max_params)
		throw std::runtime_error("Parameters are missing from the array definition file.");

	// Now double-check that there are at least two telescopes in this array
	// This is represented by at least two lines remaining in the file.
	if(lines.size() - line_number < 2)
		throw std::runtime_error("At least two telescopes are required for an array!  Check configuration file.");

    // Now iterate through the telescope definitions.  There are four entries above
    // this point in the data file so we start at 4.
    for (unsigned int i = line_number; i < lines.size(); i++)
    {
    	// First tokenize the line
    	vector < string > tokens = Tokenize(lines[i]);

    	// And now attempt to read in the line
    	try
    	{
			staname = tokens[0];
			sta_index = atoi(tokens[1].c_str());
			North = atof(tokens[2].c_str());
			East = atof(tokens[3].c_str());
			Up = atof(tokens[4].c_str());
			gain = atof(tokens[5].c_str());
			diameter = atof(tokens[6].c_str());
    	}
    	catch(...)
    	{
    		throw std::runtime_error("Error in telecope definitions in array configuration file.");
    	}
        
        // Push this station on to the list of stations for this array.
        this->stations.push_back( Station(this->latitude, staname, sta_index, xyz_coords, North, East, Up, gain, diameter) );
    }
    
    // Compute a hash table for the stations:
    this->sta_hash = ComputeStationHash(stations);
    
    // Now compute all possible baselines from these stations.
    this->baselines = ComputeBaselines(stations);
    this->bl_hash = ComputeBaselineHash(this->baselines);
    this->triplets = ComputeTriplets(this, stations);
    this->tri_hash = ComputeTripletHash(this->triplets);
    this->quadruplets = ComputeQuadruplets(this, stations);
    this->quad_hash = ComputeQuadrupletHash(this->quadruplets);

}

// This function permits a few target parameters to be overridden via. command line arguments
void Array::ParseOptions(char *argv[], int i, int argc)
{
	// Permit the r0 and wind_speed parameters to be overridden:
	for(int j = i; j < argc; j++)
	{
		// TODO: Figure out how to break out if we find an option with  "^-[a-z]"
		// otherwise we'll parse the command line arguments several times

		if ((strcmp(argv[j], "--r0") == 0) && (j < argc - 1))
		{
			try
			{
				this->r0 = atof(argv[j+1]);
				printf("r0 overridden from default value in array config file.  Now: %f m\n", this->r0);
			}
			catch(...)
			{
				throw std::runtime_error("Invalid r0 parameter override on Command Line");
			}
		}

		if ((strcmp(argv[j], "--wind_speed") == 0) && (j < argc - 1))
		{
			try
			{
				this->wind_speed = atof(argv[j+1]);
				printf("wind_speed overridden from default value in array config file.  Now: %f km/s\n", this->wind_speed);
			}
			catch(...)
			{
				throw std::runtime_error("Invalid wind speed parameter override on Command Line");
			}
		}
	}
}

// Returns a reference to a baseline object.
Baseline *  Array::GetBaseline(string baseline_name)
{
    Baseline * baseline = this->bl_hash[baseline_name];
    
    if(this->bl_hash.find(baseline_name) == this->bl_hash.end())
        throw std::runtime_error("Request for non-existant baseline, " + baseline_name + "!");
        
    return this->bl_hash[baseline_name];
}

Baseline *  Array::GetBaseline(int sta1, int sta2)
{
    // First enforce sta1 < sta2
    Sort(sta1, sta2);

    string sta1name = GetStation(sta1)->GetName();
    string sta2name = GetStation(sta2)->GetName();
    
    return GetBaseline(sta1name + "-" + sta2name);
}

Station *   Array::GetStation(int station_index)
{
    return &this->stations[station_index];
}

Station *   Array::GetStation(string sta_name)
{
    return this->sta_hash[sta_name];
}

double      Array::GetLatitude(void)
{
    return this->latitude;
}

double      Array::GetLongitude(void)
{
    return this->longitude;
}

double  Array::GetAltitude(void)
{
    return this->altitude;
}

int     Array::GetNumStations(void)
{
    return this->stations.size();
}

string  Array::GetArrayName(void)
{
    return this->arrayname;
}

vector<Station*> Array::GetAllStations(void)
{
    vector<Station*> temp;
    
    for(unsigned int i = 0; i < this->stations.size(); i++)
        temp.push_back(&this->stations[i]);
    
    return temp;
}

Triplet *   Array::GetTriplet(string triplet_name)
{
    return this->tri_hash[triplet_name];
}

Triplet *   Array::GetTriplet(int sta1, int sta2, int sta3)
{
    // Enforce sta1 < sta2 < sta3
    Sort(sta1, sta2, sta3);

    string sta1name = GetStation(sta1)->GetName();
    string sta2name = GetStation(sta2)->GetName();
    string sta3name = GetStation(sta3)->GetName();
    
    return GetTriplet(sta1name + '-' + sta2name + '-' + sta3name);
}

Quadruplet *   Array::GetQuadruplet(string quadruplet_name)
{
    return this->quad_hash[quadruplet_name];
}

Quadruplet *   Array::GetQuadruplet(int sta1, int sta2, int sta3, int sta4)
{
    // Enforce sta1 < sta2 < sta3 < sta4
    Sort(sta1, sta2, sta3, sta4);

    string sta1name = GetStation(sta1)->GetName();
    string sta2name = GetStation(sta2)->GetName();
    string sta3name = GetStation(sta3)->GetName();
    string sta4name = GetStation(sta4)->GetName();
    
    return GetQuadruplet(sta1name + '-' + sta2name + '-' + sta3name +  '-' + sta4name);
}




/// Converts this object to a valid OIFITSLIB OI_ARRAY struct.
oi_array    Array::GetOIArray(void)
{	
    /// \todo This could be claned up to use datamembers instead of function calls.
    
    // init local vars:
    oi_array array;
    int nstations = GetNumStations();
    string arrname = GetArrayName();
	double GeocLat, GeocRadius;
    
    // Allocate room for the station information
    array.elem = (element *) malloc(nstations * sizeof(element));
    
    // Set the revision number, array name, and coordiantes.
	array.revision = 1;
	strncpy(array.arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(array.frame, "GEOCENTRIC", FLEN_VALUE);

	wgs84_to_geoc(latitude, altitude, &GeocLat, &GeocRadius);
	array.arrayx = GeocRadius * cos(GeocLat) * cos(this->longitude);
	array.arrayy = GeocRadius * cos(GeocLat) * sin(this->longitude);
	array.arrayz = GeocRadius * sin(GeocLat);

	array.nelement = nstations;
	Station * station;
	string station_name;
	for (int i = 0; i < nstations; i++)
	{
	    station = this->GetStation(i); 
	    station_name = station->GetName();
	    station_name.resize(16);
	    strncpy(array.elem[i].tel_name, station->GetName().c_str(), 16);
	    strncpy(array.elem[i].sta_name, station->GetName().c_str(), 16);
	    array.elem[i].sta_index = station->GetIndex();
	    array.elem[i].diameter = station->diameter;
	    array.elem[i].staxyz[0] = station->xyz[0];
	    array.elem[i].staxyz[1] = station->xyz[1];
	    array.elem[i].staxyz[2] = station->xyz[2];
	}
	
	return array;
}

string	Array::GetAllStationNames(void)
{
	string names = "";
	string name = "";
	int num_stations = this->stations.size();
	for(unsigned int i = 0; i < this->stations.size(); i++)
	{
		name = stations[i].GetName();
		if(i == num_stations - 1)
			names += name;
		else
			names += name + ",";
	}

	return names;
}

double		Array::Get_r0(void)
{
	return this->r0;
}

double		Array::GetWindSpeed(void)
{
	return this->wind_speed;
}

double		Array::GetThroughput(void)
{
	return this->throughput;
}
