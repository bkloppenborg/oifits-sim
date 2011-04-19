/*
 * oifits-sim.cpp
 *
 *  Created on: Apr 8, 2011
 *      Author: bkloppenborg
 */
#include "oifits-sim.h"

#include <string>
#include <iostream> 
#include <cstdio>
#include <stdexcept>

#include "Target.h"
#include "Array.h"
#include "Combiner.h"
#include "SpectralMode.h"
#include "Observation.h"
#include "Obs_HA.h"
#include "Obs_OIFITS.h"
#include "NoiseModel.h"
#include "NoiseModel_Tatulli2006.h"

extern "C" {
	#include "random.h"
}
using namespace std;

// Prints out help describing the options on the command line
void PrintHelp()
{
	string usage = "The OIFITS Simulator\n"
	"Usage: \n"
	" oifits-sim arguments \n"
	" \n"
	"Arguments: \n"
	" -h     Prints this message \n"
	" -t     Target definition file \n"
	" -i     Input image \n"
	" -o     Output OIFITS file \n"
	" \n"
	"There are two Simulation Options: \n"
	"(1) from an existing OIFITS file, copying uncertainties from the real data: \n"
	" -a     The array used.\n"
	" -m     Spectral Mode for the combiner (combiner must be specified before -m).\n"
	" -d     Input OIFITS file \n"
	" \n"
	"or (2) from a list of observations at an array with a specific combiner in\n"
	"which the noise is estimated as described in the documentation: \n"
	" -a     The array used.\n"
	" -c     The combiner.\n"
	" -m     Spectral Mode for the combiner (combiner must be specified before -m).\n"
	" -obs   Observation Definition File \n"
	" \n"
	"Note: Some of the parameters also have double-dash overrides that \n"
	"immediately follow other parameters that change default/averaged values \n"
	"that are found in configuration files.  See the documentation for details.";
	
	cout << usage << "\n";
}

// The main routine.  Basically just used to parse out some parameters before handing
// things off to other functions.
int main(int argc, char *argv[])
{
	// Define characters that are used as comments:
	// TODO: Update the documentation to reflect these characters are all valid comment chars.
    const string comment_chars("\\#~$&£%");

    // Define a few variables for the simulation.
	Target target;
	string output_fname;

	Array array;
	Combiner combiner;
	SpectralMode spec_mode;
	vector <Observation*> observations;

	// TODO: For now we only have one noise model, so we load it by default
	NoiseModel * noisemodel = new NoiseModel_Tatulli2006();

	int n_params = 0;


    if(argc == 1)
        PrintHelp();

	for (int i = 1; i < argc; i++)
	{
		// First see if the user is requesting help:
		if(strcmp(argv[i], "-h") == 0)
		{
			PrintHelp();
		}	
	
		// We need to know some information about the target:
		if ((strcmp(argv[i], "-t") == 0) && (i < argc - 1))
		{
			target.ImportFile(string(argv[i+1]), comment_chars);
			target.ParseFileOptions(argv, i+2, argc);
			n_params += 1;
		}		
			
		// And the image file we're going to use for the simulation:
		if ((strcmp(argv[i], "-i") == 0) && (i < argc - 1))
		{
			target.SetImage(string(argv[i + 1]));
			target.ParseImageOptions(argv, i+2, argc);
			n_params += 1;
		}

		// And where we're going to write the output:
		if ((strcmp(argv[i], "-o") == 0) && (i < argc - 1))
		{
			try
			{
				output_fname = string(argv[i + 1]);
				n_params += 1;
			}
			catch(...)
			{
				throw std::runtime_error("Invalid Output File Definition");
			}
		}

		//
		// After this there are two cases:
		//

		// ################
		// (1) we use an existing OIFITS data file.  Really this has everything
		// need to run the simulation, but I don't want to write a module to parse
		// the OIFITS array and wavelength tables.  So, we'll just require
		// that the array and spectral mode be specified too.
		if ((strcmp(argv[i], "-d") == 0) && (i < argc - 1))
		{
			observations = Obs_OIFITS::ReadObservation_OIFITS(string(argv[i+1]));
			n_params += 1;
		}

		// ################
		// or (2) we are simulating something from scratch in which case we need
		// the array
		if ((strcmp(argv[i], "-a") == 0) && (i < argc - 1))
		{
			array.ImportFile(string(argv[i+1]), comment_chars);
			array.ParseOptions(argv, i, argc);
			n_params += 1;
		}

		// the combiner
		if ((strcmp(argv[i], "-c") == 0) && (i < argc - 1))
		{
			combiner.ImportFile(string(argv[i+1]), comment_chars);
			n_params += 1;
		}

		// the spectral mode.  Note, the combiner must be specified first so we can check that
		// the two are indeed to be used with eachother.
		if ((strcmp(argv[i], "-m") == 0) && (i < argc - 1))
		{
			if(combiner.name == "")
			{
				cout << "The combiner, -c, must be specified before the spectral mode, -m.\n";
				exit(0);
			}

			spec_mode.ImportFile(string(argv[i+1]), combiner.GetName(), comment_chars);
			n_params += 1;

		}

		// Now for observations
		if ((strcmp(argv[i], "-obs") == 0) && (i < argc - 1))
		{
			// First ensure that the array has been defined.  If not, quit.
			if(array.GetArrayName() == "")
			{
				cout << "The array, -a, must be specified before the observations, -obs.\n";
				exit(0);
			}

			// Observations are a little funny.  We permit both files and command line options.
			// Just pass things off to the observation class so it can decide what to do:
			observations = Observation::ParseCommandLine(&array, argv, i, argc, comment_chars);
			n_params += 1;
		}
	}

	// Check that all parameters were specified
	// TODO: Better checking here should be implemented
	if(n_params == 7)
		run_sim(&target, &array, &combiner, &spec_mode, noisemodel, observations, output_fname);
	else
		cout << "Something is missing on the command line, quitting!" << endl;

	// Clean up memory
	delete noisemodel;

	return 0;
}

void run_sim(Target * target, Array * array, Combiner * combiner, SpectralMode * spec, NoiseModel * noisemodel, vector<Observation*> observations, string output_filename)
{
	// Pull up the random number generator.
	static Rand_t random_seed;

	// See if the output file aready exists:
	// TODO: this could be bad/unexpected behavior.
	if(FileExists(output_filename))
	{
		cout << "Warning, output file aready exists!  Overwriting..." << endl;
		output_filename = "!" + output_filename;
	}

    // Open up the OIFITS file.
	fitsfile *fptr;
	int status = 0;
	fits_create_file(&fptr, output_filename.c_str(), &status);
	if (status)
	{
		fits_report_error(stderr, status);
		return;
	}

    oi_array oi_arr = array->GetOIArray();
    write_oi_array(fptr, oi_arr, 1, &status);
    oi_target oi_targ = target->GetOITarget();
    write_oi_target(fptr, oi_targ, &status);
    oi_wavelength oi_wave = spec->GetOIWavelength();
    write_oi_wavelength(fptr, oi_wave, 1, &status);

    vector<double> wavenumbers = spec->GetWavenumbers();

    // Now compute the vis2 records and t3s:
    oi_vis2 vis2table;
    oi_t3 t3table;
    Observation * observation;
    int n_observations = observations.size();

	for(unsigned int i = n_observations; i > 0; i--)
	{
		observation = observations.back();
		// First look up the type of observation
		ObsType type = observation->GetObsType();

		// Do a dymamic cast to get the subclass object back
		if(type == HOUR_ANGLE || type == DESCRIPTIVE)
		{
			Obs_HA * observation = dynamic_cast<Obs_HA *>(observation);
		}
		else    //(type == OIFITS)
		{
			Obs_OIFITS * observation = dynamic_cast<Obs_OIFITS *>(observation);
		}

		//printf("Simulating Observation at HA %f \n", observation->GetHA(target->right_ascension));

		vis2table = observation->GetVis2(array, combiner, spec, target, noisemodel, random_seed);
		write_oi_vis2(fptr, vis2table, 1, &status);

		if(observation->HasTriplets())
		{
			t3table = observation->GetT3(array, combiner, spec, target, noisemodel, random_seed);
			write_oi_t3(fptr, t3table, 1, &status);
		}

		// All done with this observation object.  Pop it off the vector and free memory.
		observations.pop_back();
		delete observation;

		cout << "Completed Observation " << (n_observations - i + 1) << endl;
	}

	if (status)
	{
		fits_delete_file(fptr, &status);
		return;
	}
	else
	{
		fits_close_file(fptr, &status);
	}
}
