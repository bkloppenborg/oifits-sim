/// \file SpectralMode.cpp

#include <iostream>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "SpectralMode.h"

using std::cout;
using std::string;

/// Constructor for a SpectralMode object.  
/// Reads in space delineated text file of one or more lines in the following format:
///     mean_wavelength effective_bandwidth
/// Comment characters are specified by comment_chars
SpectralMode::SpectralMode(const char* SpectralMode_file, string comment_chars) 
{

    vector < string > lines;    // stores non-blank, non-comment lines

    ifstream fil(SpectralMode_file);
    if (fil.is_open()) 
    {
        string line;
        while (!fil.eof()) 
        {
            getline(fil, line);
            while ((line.size() == 0 || comment_chars.find(line[0]) != string::npos) && !fil.eof()) 
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
        throw std::runtime_error("Error opening spectral mode file");
    }
    
    // Now, set the number of channels
    this->nchannels = lines.size();
    this->mean_wavelength.setsize(this->nchannels);
    this->delta_wavelength.setsize(this->nchannels);
    this->mean_wavenumber.setsize(this->nchannels);
    this->delta_wavenumber.setsize(this->nchannels);
    
    // Parse the text file, pull out the wavelengths and compute the associated wavenumbers.
    double waven_low = 0;
    double waven_high = 0;
    for(int i = 0; i < this->nchannels; i++)
    {
        // Parse the line into two tokens:
        istringstream lineStream(lines[i]);
        vector <string> tokens;
        string item;
        while(lineStream >> item)
        {
            tokens.push_back(item);
        }
    
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

//// Generates the wavelength to pixel mapping
//SpectralMode::SpectralMode(int type, double minlambda , double maxlambda , int spectralchannels )
//{
//  nchannels = spectralchannels;
//  minnumber = 1. / maxlambda;
//  maxnumber = 1. / minlambda;
//  wavenumber.setsize(nchannels + 1);
//  wavelength.setsize(nchannels + 1);

//  mean_wavenumber.setsize(nchannels);
//  delta_wavenumber.setsize(nchannels);

//  mean_wavelength.setsize(nchannels);
//  delta_wavelength.setsize(nchannels);

//  switch (type)
//  // Note: lower wavenumber = lower channel number
//  {
//    case 0:
//    {
//      // Linear wavenumber (e.g. grism)
//      for (int j = 0; j < nchannels + 1; j++)
//      {
//        wavenumber[ j ] = minnumber + (maxnumber - minnumber) * double(j) / double(nchannels);
//        wavelength[ j ] = 1. / wavenumber[ j ];
//      }

//      break;
//    }

//    case 1:
//    {
//      // Linear wavelength
//      for (int j = 0; j < nchannels + 1; j++)
//      {
//        wavelength[ j ] = minlambda + (maxlambda - minlambda) * double(nchannels - j) / double(nchannels);
//        wavenumber[ j ] = 1. / wavelength[ j ];

//      }

//      break;
//    }

//    default:
//      cout << "Error in spectral mode selection " << endl;
//      break;
//  }

//  for (int j = 0; j < nchannels; j++)
//  {
//    mean_wavenumber[ j ] = (wavenumber[ j ] + wavenumber[ j + 1 ]) / 2.; // mean wavenumber in pixel j
//    delta_wavenumber[ j ] = wavenumber[ j + 1 ] - wavenumber[ j ]; // waveband for pixel j
//    mean_wavelength[ j ] = (wavelength[ j ] + wavelength[ j + 1 ]) / 2.;
//    delta_wavelength[ j ] = wavelength[ j + 1 ] - wavelength[ j ];
//  }
//}

//// Get directly the wavelength channels from a list, wavenumber increases with channel number
//SpectralMode::SpectralMode(int type, Row<double> wav )
//{
//  nchannels = wav.size() - 1 ;
//  wavenumber.setsize(nchannels + 1);
//  wavelength.setsize(nchannels + 1);

//  mean_wavenumber.setsize(nchannels);
//  delta_wavenumber.setsize(nchannels);

//  mean_wavelength.setsize(nchannels);
//  delta_wavelength.setsize(nchannels);

//  switch (type)
//  {
//    case 0: // the input list corresponds to wavenumbers
//    {
//      wavenumber = wav;
//      for (int j = 0; j < nchannels + 1; j++)
//      {
//        wavelength[ j ] = 1. / wav[ j ];
//      }
//      break;
//    }

//    case 1: // the input list corresponds to wavelengths
//    {
//      wavelength = wav;
//      for (int j = 0; j < nchannels + 1; j++)
//      {
//        wavenumber[ j ] = 1. / wav[ j ];
//      }

//      break;
//    }

//    default:
//      cout << "Error in spectral mode selection " << endl;
//      break;
//  }

//  for (int j = 0; j < nchannels; j++)
//  {
//    mean_wavenumber[ j ] = (wavenumber[ j ] + wavenumber[ j + 1 ]) / 2.; // mean wavenumber in pixel j
//    delta_wavenumber[ j ] = wavenumber[ j + 1 ] - wavenumber[ j ]; // waveband for pixel j
//    mean_wavelength[ j ] = (wavenumber[ j ] + wavenumber[ j + 1 ]) / 2.;
//    delta_wavelength[ j ] = wavelength[ j + 1 ] - wavelength[ j ];
//  }

//}

SpectralMode::~SpectralMode( )
{
  //
}
