#include "Combiner.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include "ReadTextFile.h"

using std::vector;
using std::string;
using std::ifstream;
using std::cout;
using std::endl;

Combiner::Combiner()
{
	this->name = "";
}

Combiner::~Combiner()
{

}

void Combiner::ImportFile(string filename, string comment_chars)
{
	// Some local variables
	int max_params = 15; // This is the number of non-telescope parameters found in the array definition file.
	int n_params = 0; // the number of parameters read in from the file

    int line_number;

    // We permit shortcuts in the name of the combiners.  See if the absolute path was given, otherwise
    // generate a new filename from the shortcut:
    if(!FileExists(filename))
    {
    	filename = "../etc/" + filename + ".txt";
    }

    // stores non-blank, non-comment lines
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Combiner Definition File");
	vector <string> results;

	for(line_number = 0; line_number < max_params; line_number++)
	{
		// Clear out the results, split the string and strip whitespace
        results.clear();
        StringSplit(lines[line_number], "=", results);
        StripWhitespace(results);


        if(results[0] == "name")
        {
        	try
        	{
        		this->name = string(results[1]);
        		n_params += 1;
        		cout << "Using Combiner " << this->name << endl;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid combiner name");
        	}
        }

        if(results[0] == "int_trans")
        {
        	try
        	{
        		this->int_trans = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid Internal Transmission Value for the combiner");
        	}
        }

        if(results[0] == "n_pix_fringe")
        {
        	try
        	{
        		this->n_pix_fringe = atoi(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid number of fringe pixels for the combiner");
        	}
        }

        if(results[0] == "n_pix_photometry")
        {
        	try
        	{
        		this->n_pix_photometry = atoi(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid number of photometry pixels for the combiner");
        	}
        }

        if(results[0] == "flux_frac_photometry")
        {
        	try
        	{
        		this->flux_frac_photometry = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid photometry flux fraction for the combiner");
        	}
        }

        if(results[0] == "flux_frac_fringes")
        {
        	try
        	{
        		this->flux_frac_fringes = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid fringe flux fraction for the combiner");
        	}
        }

        if(results[0] == "throughput_photometry")
        {
        	try
        	{
        		this->throughput_photometry = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid photometry throughput for the combiner");
        	}
        }

        if(results[0] == "throughput_fringes")
        {
        	try
        	{
        		this->throughput_fringes = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid photometry throughput for the combiner");
        	}
        }

        if(results[0] == "n_splits")
        {
        	try
        	{
        		this->n_splits = atoi(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid number splits specified for the combiner");
        	}
        }

        if(results[0] == "vis")
        {
        	try
        	{
        		this->visibility = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid combiner visibility value");
        	}
        }

        if(results[0] == "read_noise")
        {
        	try
        	{
        		this->read_noise = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid read noise for combiner");
        	}
        }

        if(results[0] == "quantum_efficiency")
        {
        	try
        	{
        		this->quantum_efficiency = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid QE combiner");
        	}
        }

        if(results[0] == "v2_cal_err")
        {
        	try
        	{
        		this->v2_frac_cal_err = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid V2 error for combiner");
        	}
        }

        if(results[0] == "phase_cal_err")
        {
        	try
        	{
        		this->phase_cal_err = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid phase error for combiner");
        	}
        }

        if(results[0] == "incoh_int_time")
        {
        	try
        	{
        		this->incoh_int_time = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid incoherent integraion time for combiner");
        	}
        }
	}

	if (n_params != max_params)
		throw std::runtime_error("Missing parameters in combiner file!");
}

string Combiner::GetName(void)
{
	return this->name;
}

double Combiner::GetThroughput(void)
{
	return this->int_trans;
}

int 	Combiner::GetNumPhotometricPixels()
{
	return this->n_pix_photometry;
}
int		Combiner::GetNumFringePixels()
{
	return this->n_pix_fringe;
}

double	Combiner::GetPhotometryFluxFrac()
{
	return this->flux_frac_photometry;
}
double	Combiner::GetFringeFluxFrac()
{
	return this->flux_frac_fringes;
}

double 	Combiner::GetPhotometryThroughput()
{
	return this->throughput_photometry;
}
double 	Combiner::GetFringeThroughput()
{
	return this->throughput_fringes;
}

double 	Combiner::GetNumSplits()
{
	return this->n_splits;
}
