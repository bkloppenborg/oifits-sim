#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <list>

#include "HourAngle.h"
#include "Source.h"
#include "Array.h"
#include "Observation.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;

HourAngle::HourAngle(double HAstart, double HAend, double HAstep)
{
	this->Nha = (int) ((HAend - HAstart + HAstep) / HAstep);
	HA.setsize(Nha);
	for (int i = 0; i < Nha; i++)
	{
		this->HA[i] = HAstart + i * HAstep;
	}
}

// Derive the hour angle list from a list of observations:
HourAngle::HourAngle(Source & target, Array & array, std::list<Observation, std::allocator<Observation> >& observations)
{
    // Set the number of hour angles
    int num_obs = observations.size();
    double ra = target.right_ascension;
    double dec = target.declination;
    double array_lat = array.latitude;
    double array_lon = array.longitude;
    double hour_angle;
    
    this->Nha = num_obs;
    this->HA.setsize(num_obs);

    // Calculate the hour angles from the observations:
    for(int i = 0; i < num_obs; i++)
    {        
        Observation obs = observations.front();
        hour_angle = obs.GetHA(ra, dec, array_lat, array_lon);
        this->HA[i] = hour_angle * 3600;
        //printf("JD %f ObsHA %f \n", obs.GetJD(), this->HA[i]/3600);
        observations.pop_front();
    }
}

HourAngle::HourAngle(Row < double >HA)
{
	this->HA = HA;				// defining the HA vector
	this->Nha = HA.size();
}

HourAngle::HourAngle(const char *HourAngle_file)
{
	const string comments("\\/#~$&£%");
	vector < string > lines;	// stores non-blank, non-comment lines

	ifstream fil(HourAngle_file);
	if (fil.is_open())
	{
		string line;
		while (!fil.eof())
		{
			getline(fil, line);
			while ((line.size() == 0
				  || comments.find(line[0]) != string::npos) && !fil.eof())
			{
				getline(fil, line);
			}
			if (!fil.eof())
				lines.push_back(line);
		}
		fil.close();
	}

	this->Nha = lines.size();
	this->HA.setsize(this->Nha);

	for (int i = 0; i < this->Nha; i++)
	{
		if (!
		   (isdigit(lines[i][0]) || lines[i][0] == '+'
			  || lines[i][0] == '-'))
			throw std::runtime_error("Invalid hour angle");
		else
		this->HA[i] = atof(lines[i].c_str());
	}
}
