/// \file SpectralMode.cpp

#include <iostream>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "SpectralMode.h"
#include "Common.h"
#include "ReadTextFile.h"

using std::cout;
using std::string;

SpectralMode::SpectralMode()
{

}

SpectralMode::~SpectralMode( )
{
  //
}

void SpectralMode::ImportFile(string filename, string combiner, string comment_chars)
{
	// Some local variables
	int max_params = 2; // This is the number of non-telescope parameters found in the array definition file.
	int n_params = 0; // the number of parameters read in from the file

    int line_number;

    // Permit the usage of shortcut names.  Note
    // TODO: Use a global variable, read in from a configuration file, to specify "../etc" below.
    if(!FileExists(filename))
    {
    	filename = "../etc/" + combiner + "_" + filename + ".txt";
    }

    // stores non-blank, non-comment lines
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Spectral Mode Definition File");
	vector <string> results;

	for(line_number = 0; line_number < max_params; line_number++)
	{
		// Clear out the results, split the string and strip whitespace
        results.clear();
        StringSplit(lines[line_number], "=", results);
        StripWhitespace(results);


        if(results[0] == "combiner")
        {
        	try
        	{
        		this->combiner = string(results[1]);
        		n_params += 1;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Error in combiner name name");
        	}
        }

        if(results[0] == "mode")
        {
        	try
        	{
        		this->spec_mode = string(results[1]);
        		n_params += 1;
        		cout << "Combiner is in spectral mode: " << this->spec_mode << endl;
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Error in combiner's spectral mode.");
        	}
        }

	}

	// A little error checking
	if(n_params != max_params)
		throw std::runtime_error("Parameters are missing from the spectral mode definition file");

	// Now double-check that there are at least two telescopes in this array
	// This is represented by at least two lines remaining in the file.
	if(lines.size() - line_number < 1)
		throw std::runtime_error("At least one wavelength + bandpass is required for a spectral mode.  Check configuration file.");
    
    // Now, set the number of channels
    this->nchannels = lines.size() - line_number;
    this->mean_wavelength.setsize(this->nchannels);
    this->delta_wavelength.setsize(this->nchannels);
    this->mean_wavenumber.setsize(this->nchannels);
    this->delta_wavenumber.setsize(this->nchannels);
    
    // Parse the text file, pull out the wavelengths and compute the associated wavenumbers.
    double waven_low = 0;
    double waven_high = 0;
    for(int i = 0; i < this->nchannels; i++)
    {
    	// Tokenize the line:
        vector <string> tokens = Tokenize(lines[i + n_params]);

        mean_wavelength[i] = atof(tokens[0].c_str());
        delta_wavelength[i] = atof(tokens[1].c_str());
        
        waven_low = 1/(mean_wavelength[i] - delta_wavelength[i]/2);
        waven_high = 1/(mean_wavelength[i] + delta_wavelength[i]/2);
        
        mean_wavenumber[i] = (waven_high + waven_low)/2;
        delta_wavenumber[i] = (waven_high - waven_low);
        
    }
    
    // Lastly compute the median wavelength (added in original code by AS 2010-06-24)
    median_wavelength = (mean_wavelength[0] + mean_wavelength[nchannels - 1]) / 2.0;
}

/// Returns an oi_wavelength object
oi_wavelength SpectralMode::GetOIWavelength(void)
{
	oi_wavelength wave;
	wave.nwave = this->nchannels;
	wave.eff_wave = (float *) malloc(this->nchannels * sizeof(float));
	wave.eff_band = (float *) malloc(this->nchannels * sizeof(float));
	wave.revision = 1;
	for(int i=0; i<FLEN_VALUE;i++)
	   wave.insname[i]='\0';
	strncpy(wave.insname, this->spec_mode.c_str(), this->spec_mode.size());
	// Now copy the wavelengths over:
	for (int i = 0; i < wave.nwave; i++)
	{
		wave.eff_wave[i] = this->mean_wavelength[i];
		wave.eff_band[i] = this->delta_wavelength[i];
	}
	
	return wave;
}

/// Returns the wavelengths as a vector.
vector<double> SpectralMode::GetWavelengths(void)
{
    vector<double> wavelengths;
    
    for(int i = 0; i < this->nchannels; i++)
        wavelengths.push_back(1.0*this->mean_wavelength[i]);
        
    return wavelengths;
}

/// Returns the wavenumbers as a vector.
vector<double> SpectralMode::GetWavenumbers(void)
{
    vector<double> wavenumbers;
    
    for(int i = 0; i < this->nchannels; i++)
        wavenumbers.push_back(1.0 / this->mean_wavelength[i]);
        
    return wavenumbers;
}
