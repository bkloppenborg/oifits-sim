#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "Simulator.h"
#include "oifits_sim.h"
#include "VisSimParams.h"
#include "random.h"

#include "userinterface.h"

using std::string;


/// \todo Determine where a forward declaration will suffice, rather than including header files

int main(int argc, char *argv[])
{
	if (argc > 1)
	{
		VisSimParams p = VisSimParams();

		for (int i = 1; i < argc; i++)
		{
			if ((strcmp(argv[i], "-d") == 0) && (i < argc - 1))
			{
				p.input_oifits_filename = argv[i + 1];
				p.bFromOIFITSFile = true;
				cout << "Input OIFITS data file is  " << p.input_oifits_filename << endl;
			}
			if ((strcmp(argv[i], "-t") == 0) && (i < argc - 1))
			{
				p.target_filename = argv[i + 1];
				cout << "Target is " << p.target_filename << endl;
			}
			if ((strcmp(argv[i], "-a") == 0) && (i < argc - 1))
			{
				p.array_filename = argv[i + 1];
				cout << "Array is " << p.array_filename << endl;
			}
			if ((strcmp(argv[i], "-i") == 0) && (i < argc - 1))
			{
				p.instrument_filename = argv[i + 1];
				cout << "Instrument is " << p.instrument_filename << endl;
			}
			if ((strcmp(argv[i], "-s") == 0) && (i < argc - 1))
			{
				p.SpectralMode_filename = argv[i + 1];
				cout << "SpectralMode is " << p.
				   SpectralMode_filename << endl;
			}
			if ((strcmp(argv[i], "-obs") == 0) && (i < argc - 1))
			{
				p.observation_filename = argv[i + 1];
				cout << "Observations are found in " << p.observation_filename << endl;
			}
			if ((strcmp(argv[i], "-o") == 0) && (i < argc - 1))
			{
				p.oifits_filename = argv[i + 1];
				cout << "The output file is " << p.oifits_filename << endl;
			}
		}
		if (!p.have_all_params())
		{
			cout << "Missing at least one of the required arguments. \n";
			exit(2);
		}
		else
		{
			run_sim(&p);
		}
	}
	else
	{
		gui_main(argc, argv);
	}
	exit(EXIT_SUCCESS);
}
