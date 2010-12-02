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
#include "Simulator.h"  // Required for constants, like PI.

/// \todo fix this function to comply with the class definition.
//Array::Array(std::string arrayname, double center_latitude,
//             double center_longitude, double altitude, int nstations,
//             Station * station, ...)
//{                               // create an Array from a number of stations
//    assert(nstations > 1);
//    this->arrayname = arrayname;
//    this->nstations = nstations;
//    this->latitude = center_latitude;
//    this->longitude = center_longitude;
//    this->altitude = altitude;
//    this->x.setsize(nstations);
//    this->y.setsize(nstations);
//    this->z.setsize(nstations);
//    this->diameter.setsize(nstations);
//    this->gain.setsize(nstations);
//    this->layer.setsize(nstations);
//    Station *pstation = NULL;

//    va_list stationpointer;

//    for (int ii = 0; ii < nstations; ii++)
//    {
//        if (ii == 0)
//        {
//            va_start(stationpointer, station);
//            pstation = station;
//        }
//        else
//        {
//            pstation = va_arg(stationpointer, Station *);
//        }
//        x[ii] = pstation->x;
//        y[ii] = pstation->y;
//        z[ii] = pstation->z;
//        layer[ii] = pstation->layer;
//        diameter[ii] = pstation->diameter;
//        gain[ii] = pstation->gain;
//    }
//    va_end(stationpointer);
//}

// Create and initialize an array from data given in a text file.
Array::Array(const char *Array_file)
{                               
    /// \todo It would be wise to parse these values from the file in key->value pairs
    const string comments("\\/#~$&??%");
    
    // Local varaibles for creating Station objects.
    string staname;
    double North;
    double East;
    double Up;
    double gain;
    double diameter;

    // stores non-blank, non-comment lines
    vector < string > lines;    
    ifstream fil(Array_file);

    if (fil.is_open())
    {
        string line;

        while (!fil.eof())
        {
            getline(fil, line);
            while ((line.size() == 0 || comments.find(line[0]) != string::npos)
                   && !fil.eof())
            {
                getline(fil, line);
            }
            if (!fil.eof())
                lines.push_back(line);
        }
        fil.close();
    }
    else
    {
        /// \exception runtime_error Error opening array file
        /// The array file could not be located.  It is likely that the user just specified an invalid path.
        throw std::runtime_error("Error opening array file");
    }

    // Before we start parsing the input, allocate the array and baseline vectors:

    // read in the name of the telescope array
    if (!isalpha(lines[0][0]))
    {
        throw std::runtime_error("Invalid array name");
    }
    else
    {
        arrayname = lines[0];
        cout << "Array name: " << arrayname << endl;
    }

    // read in the latitude of the telescope array
    if (!(isdigit(lines[1][0]) || lines[2][0] == '+' || lines[1][0] == '-'))
    {
        throw std::runtime_error("Invalid latitude");
    }
    else
    {
        double lat = atof(lines[1].c_str());
        cout << "Array latitude (deg): " << lat << endl;
        
        // Convert to radians
        this->latitude = lat * PI / 180;
    }

    // read in the longitude of the telescope array
    if (!(isdigit(lines[2][0]) || lines[2][0] == '+' || lines[2][0] == '-'))
    {
        throw std::runtime_error("Invalid longitude");
    }
    else
    {
        double lon = atof(lines[2].c_str());
        cout << "Array longitude (deg): " << longitude << endl;
        
        this->longitude = lon * PI / 180;
    }

    // read in the altitude of the telescope array
    if (!isdigit(lines[3][0]))
    {
        throw std::runtime_error("Invalid altitude");
    }
    else
    {
        altitude = atof(lines[3].c_str());
        cout << "Array altitude (m): " << altitude << endl;
    }

    for (int i = 0; i < lines.size(); i++)
    {

        istringstream lineStream(lines[i + 4]);

        vector < string > tokens;
        while (lineStream)
        {
            string item;

            lineStream >> item;
            tokens.push_back(item);
        }

        staname = tokens[0];
        North = atof(tokens[1].c_str());
        East = atof(tokens[2].c_str());
        Up = atof(tokens[3].c_str());
        gain = atof(tokens[4].c_str());
        diameter = atof(tokens[5].c_str());
        
        // Push this station on to the list of stations for this array.
        this->stations.push_back(Station(staname, North, East, Up, gain, diameter));
    }
    
    // Now compute all possible baselines from these stations.
    this->ComputeBaselines();
}

void Array::ComputeBaselines(void)
{
    int num_stations = this->stations.size();
    string sta1_name;
    string sta2_name;
    string bl_name;
    
    // Now compute all of the baselines and make a hash table for each baseline value
    for(int i = 0; i < num_stations; i++)
    {
        sta1_name = this->stations[i].GetName();
        
        for(int j = i; j < num_stations; j++)
        {
            sta2_name = this->stations[j].GetName();
            
            // Create a name for the baseline in the form STATION1-STATION2
            bl_name = sta1_name + "-" + sta2_name;
            
            // Create a new baseline, append it to the list of baselines
            this->baselines.push_back(Baseline(this->latitude, this->stations[i], this->stations[j]));
            
            // Create a new hash entry for this baseline.  We use this to look things up later.
            this->bl_hash.insert( BaselineHash::value_type(bl_name, this->baselines.back()) );
        }
    }
}

// Returns a reference to a baseline object.
Baseline & Array::GetBaseline(string baseline_name)
{
    return this->bl_hash[baseline_name.c_str()];
}

Station & Array::GetStation(int station_index)
{
    return this->stations[station_index];
}

double  Array::GetLatitude(void)
{
    return this->latitude;
}

double  Array::GetLongitude(void)
{
    return this->longitude;
}

double  Array::GetAltitude(void)
{
    return this->altitude;
}

