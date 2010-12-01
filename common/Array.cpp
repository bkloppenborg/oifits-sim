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

Array::Array(std::string arrayname, double center_latitude,
             double center_longitude, double altitude, int nstations,
             Station * station, ...)
{                               // create an Array from a number of stations
    assert(nstations > 1);
    this->arrayname = arrayname;
    this->nstations = nstations;
    this->latitude = center_latitude;
    this->longitude = center_longitude;
    this->altitude = altitude;
    this->x.setsize(nstations);
    this->y.setsize(nstations);
    this->z.setsize(nstations);
    this->diameter.setsize(nstations);
    this->gain.setsize(nstations);
    this->layer.setsize(nstations);
    Station *pstation = NULL;

    va_list stationpointer;

    for (int ii = 0; ii < nstations; ii++)
    {
        if (ii == 0)
        {
            va_start(stationpointer, station);
            pstation = station;
        }
        else
        {
            pstation = va_arg(stationpointer, Station *);
        }
        x[ii] = pstation->x;
        y[ii] = pstation->y;
        z[ii] = pstation->z;
        layer[ii] = pstation->layer;
        diameter[ii] = pstation->diameter;
        gain[ii] = pstation->gain;
    }
    va_end(stationpointer);
}

Array::Array(const char *Array_file)
{                               
    /// \todo It would be wise to parse these values from the file in key->value pairs
// create an Array from data read from txt
                                // file
    const string comments("\\/#~$&??%");

    vector < string > lines;    // stores non-blank, non-comment lines

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

    nstations = lines.size() - 4;       // NB line count is now correct
    x.setsize(nstations);
    y.setsize(nstations);
    z.setsize(nstations);
    gain.setsize(nstations);
    diameter.setsize(nstations);

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
        latitude = atof(lines[1].c_str());
        cout << "Array latitude: " << latitude << endl;
    }

    // read in the longitude of the telescope array
    if (!(isdigit(lines[2][0]) || lines[2][0] == '+' || lines[2][0] == '-'))
    {
        throw std::runtime_error("Invalid longitude");
    }
    else
    {
        longitude = atof(lines[2].c_str());
        cout << "Array longitude: " << longitude << endl;
    }

    // read in the altitude of the telescope array
    if (!isdigit(lines[3][0]))
    {
        throw std::runtime_error("Invalid altitude");
    }
    else
    {
        altitude = atof(lines[3].c_str());
        cout << "Array altitude: " << altitude << endl;
    }

    for (int i = 0; i < nstations; i++)
    {

        istringstream lineStream(lines[i + 4]);

        vector < string > tokens;
        while (lineStream)
        {
            string item;

            lineStream >> item;
            tokens.push_back(item);
        }

        staname.push_back(tokens[0]);
        x[i] = atof(tokens[1].c_str());
        y[i] = atof(tokens[2].c_str());
        z[i] = atof(tokens[3].c_str());
        gain[i] = atof(tokens[4].c_str());
        diameter[i] = atof(tokens[5].c_str());
    }

}
