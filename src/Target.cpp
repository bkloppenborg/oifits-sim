/*
 * Target.cpp
 *
 *  Created on: Apr 8, 2011
 *      Author: bkloppenborg
 */

#include "Target.h"
#include "Common.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <locale>

#include "textio.hpp"

Target::Target()
{
	this->name = "Uninitialized";
	this->pixellation = 1;
	this->mag = 1;
	this->band = 'A';
	this->temperature = 0;
	this->background = 0;
	this->background_ap = 0;
	this->declination = 0;
	this->right_ascension = 0;
}

Target::~Target()
{

}

void Target::ImportFile(string filename, string comment_chars)
{
	vector <string> lines = ReadFile(filename, comment_chars, "Cannot Open Target Definition File");
	vector <string> results;

	int max_params = 9;
	int n_params = 0;

	for(unsigned int i = 0; i < lines.size(); i++)
	{
		// Clear out the results, split the string and strip whitespace
        results.clear();
        StringSplit(lines[i], "=", results);
        StripWhitespace(results);

        // Now we parse the file by keyword
        if(results[0] == "name")
        {
        	this->name = string(results[1]);
        	n_params += 1;
        }

        if(results[0] == "res")
        {
        	try
        	{
        		this->pixellation = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid Pixellation Definition in Target File");
        	}
        }

        if(results[0] == "mag")
        {
        	try
        	{
        		this->mag = atof(results[1].c_str());
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid Magnitude Definition in Target File");
        	}
        }

        if(results[0] == "band")
        {
        	try
        	{
        		this->band = results[1][0];
				n_params += 1;
			}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid Band Definition in Target File");
        	}
        }

        if(results[0] == "temp")
        {
        	try
        	{
        		this->temperature = atof(results[1].c_str());
				n_params += 1;
			}
        	catch(...)
			{
        		throw std::runtime_error("Invalid Temperature Definition in Target File");
			}
        }

        if(results[0] == "bkgrnd")
        {
        	try
        	{
        		this->background = atof(results[1].c_str());
				n_params += 1;
			}
        	catch(...)
			{
        		throw std::runtime_error("Invalid Background Definition in Target File");
			}
        }

        if(results[0] == "bkgrnd_ap")
        {
        	try
        	{
        		this->background_ap = atof(results[1].c_str());
				n_params += 1;
			}
        	catch(...)
			{
        		throw std::runtime_error("Invalid Background Definition in Target File");
			}
        }

        if(results[0] == "src_dec")
        {
        	try
        	{
        		this->declination = atof(results[1].c_str()) * PI / 180;
				n_params += 1;
			}
        	catch(...)
			{
        		throw std::runtime_error("Invalid source declination Target File");
			}
        }

        if(results[0] == "src_ra")
        {
        	try
        	{
        		this->right_ascension = atof(results[1].c_str()) * PI / 180;
				n_params += 1;
			}
        	catch(...)
			{
        		throw std::runtime_error("Invalid source right ascention Target File");
			}
        }

	}

	if(n_params != max_params)
		throw std::runtime_error("The Target file is missing some parameters.  Please check the file.");

	// Now, init the flux.
	this->InitFlux();
}

// This function permits a few target parameters to be overridden via. command line arguments
void Target::ParseOptions(char *argv[], int i, int argc)
{
	// Permit a few things to be overridden via. command line:
	for(int j = i; j < argc; j++)
	{
		// TODO: Figure out how to break out if we find an option with  "^-[a-z]"
		// otherwise we'll parse the command line arguments several times.

		if ((strcmp(argv[j], "--res") == 0) && (j < argc - 1))
		{
			try
			{
				this->pixellation = atof(argv[j+1]);
				printf("Pixellation/Resolution overridden from default value target config file.  Now: %f\n", this->pixellation);
			}
			catch(...)
			{
				throw std::runtime_error("Invalid Pixellation/Resolution parameter override on Command Line");
			}
		}
		if ((strcmp(argv[j], "--mag") == 0) && (j < argc - 1))
		{
			try
			{
				this->mag = atof(argv[j+1]);
				printf("Magnitude overridden from default value target config file.  Now: %f\n", this->mag);

				// Now, init the flux.
				this->InitFlux();
			}
			catch(...)
			{
				throw std::runtime_error("Magnitude parameter override on Command Line");
			}
		}
	}
}

void Target::ParseFileOptions(char *argv[], int i, int argc)
{
	// Right now the image and file use the same parser
	ParseOptions(argv, i, argc);
}

void Target::ParseImageOptions(char *argv[], int i, int argc)
{
	// Right now the image and file use the same parser
	ParseOptions(argv, i, argc);
}

void Target::SetImage(string filename)
{
    if (filename.size() > 0)
    {
    	// Read the file into a matrix:
        fits2mat(filename.c_str(), this->image);

        // Normalize it, just in case it wasn't normalized to begin with.
        this->image /= total(this->image);
    }
}

string Target::GetName()
{
	return this->name;
}

/// Returns the target ID
/// NOTE: It is assume a single target will be simulated so this value is fixed to one, but kept
/// as a function for future expansion of this program.
int     Target::GetTargetID(void)
{
    return 1;
}

// Computes the flux at the specified wavelength
double Target::GetTargNPhotons(double wavelength, double bandwidth, double area, double time)
{
	return this->flux * pow(wavelength, -2.0) * bandwidth * area * time;
}

double Target::GetBackNPhotons(double wavelength, double bandwidth, double area, double time)
{
	return this->flux_bg * pow(wavelength, -2.0) * bandwidth * area * time;
}

// Compute the flux of the source in nph.m.m-2.s-1
void Target::InitFlux()
{
	// Multiply by waveband width (in wavenumber unit = m-1 ), collecting surface, and integration time
	// to get the number of photons
	// Using M0 fluxes for Vega in Jansky from Tokunaga, A. T. Vacca, W. D. (2005)
	// see http://adsabs.harvard.edu/abs/2005PASP..117..421T


	double flux_m0 = 0.0;
	switch(band)
	{
		case 'V':
			flux_m0 = 3630 * 1e-26 * 0.5446e-6 / 6.626068e-34;
			break;
		case 'J':
			flux_m0 = 1560* 1e-26 * 1.250e-6 / 6.626068e-34;
			break;
		case 'H':
			flux_m0 = 1040* 1e-26 * 1.644e-6 / 6.626068e-34;
			break;
		case 'K':
			flux_m0 = 645* 1e-26 * 2.198e-6 / 6.626068e-34;
			break;

		default:
			flux_m0 = 0;
	}

	// Flux density for the chosen magnitude
	this->flux = flux_m0 * pow(10., -0.4 * this->mag);

	// Background noise calculations
	// The background magnitude is in magnitude/ sq arcsec and assumed given in the same band as the source
	// The sky background aperture is in arcseconds on sky (MROI = 0.5 )
	this->flux_bg = pow(this->background_ap, 2) * flux_m0 * pow(10.0, -0.4 * this->background);

//	double SNR = flux / backgroundflux;

//	cout << "Flux = " << this->flux << " photons per m, per second." << endl;
//	cout << "Background flux = " << this->backgroundflux
//			<< " photons per m, per second." << endl;
//	cout << "Source SNR: " << SNR << endl;
}

/// Returns an OIFITS target struct representing this source.
oi_target Target::GetOITarget(void)
{
    // Init local vars:
	target * targ = (target*) malloc(1 * sizeof(target));
    // Default values following
    string veltyp_str = "LSR     ";
    string veldef_str = "OPTICAL ";
    string spectyp_str = "Unknown         ";
	string targ_name = "Simulated " + this->GetName();

	// TARGETS
	/*
	 * use malloc as free_oi_target() uses free:
	 */
	targ->target_id = 1;
	strncpy(targ->target, targ_name.c_str(), targ_name.length());
	targ->raep0 = this->right_ascension;
	targ->decep0 = this->declination;
	targ->equinox = 2000.0;
	targ->ra_err = 0.0;
	targ->dec_err = 0.0;
	targ->sysvel = 0.0;
	strncpy(targ->veltyp, veltyp_str.c_str(), veltyp_str.size());
	strncpy(targ->veldef, veldef_str.c_str(), veldef_str.size());
	targ->pmra = 0.0;
	targ->pmdec = 0.0;
	targ->pmra_err = 0.0;
	targ->pmdec_err = 0.0;
	targ->parallax = 0.0;
	targ->para_err = 0.0;

	/// \bug Spectral type output here always has some junk after it.  Not sure why though.
	strncpy(targ->spectyp, spectyp_str.c_str(), spectyp_str.size());

	oi_target oi_targ;
	oi_targ.revision = 1;
	oi_targ.ntarget = 1;
	oi_targ.targ = targ;


	return oi_targ;
}
