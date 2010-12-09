// ============================================================================
// Name : source.cpp
// Version : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
// ============================================================================

// The source module simulates the visibility of sources according to
// spectra, structure and observing angle. It returns the visibility
// at the current time whenever it is called for the baselines defined
// by the station module.

#include <iostream>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "Simulator.h"
#include "Source.h"

//Complex Source::GetVis(double wavenumber, double hour_angle, Array & array, int index_station1, int index_station2)
//{
//    Complex visibility;

//    if (source_image.GetCols() == 0)
//    {
//        // point source, no need to calculate uv coords
//        visibility = 1.0;
//    }
//    else
//    {

//        // resolved object
//        // according to time and wavenumber, get the uv points, then the
//        // corresponding visibility
//        // to implement - update the visibility only if enough time elapsed
//        // between measurements
//        Baseline B = Baseline(array.GetStation(index_station1), array.GetStation(index_station2));
//        UVPoint uv = B.UVcoords(hour_angle, declination, wavenumber); // from 
//                                                                                // UVPoint 
//                                                                                // uv 
//                                                                                // is 
//                                                                                // now 
//                                                                                // = 
//                                                                                // B.UVcoords. 
//                                                                                // function 
//                                                                                // defined 
//                                                                                // in 
//                                                                                // stations.cpp
//        visibility = GetVis(uv);        // source visibility is equal to the
//                                        // equation above: visibility +=
//                                        // source_image[ ii ][ jj ] *
//                                        // polar(1., -2.0 * PI *
//                                        // source_pixellation * milliarcsec *
//                                        // (uv.u * (double) ii + uv.v *
//                                        // (double) jj));

//    }

//    return visibility;
//}

double Source::Spectrum(double wavenumber)
{
    // return the black-body intensity in photons
    // Spectral radiance W.m-2.s-1.sr-1.m-1
    // double L = (2 * h * c * c / pow( wavelength , 5) / ( exp( h * c / (
    // wavelength * k * temperature ) -1.0)
    // Number of photons = L*wavelength/(h*c), nph.m-2.s-1.sr-1.m-1
    // double nph = L*wavelength/(h*c);

    // flux * nph / nph_max, where nph_max gives the maximum of emission
    return (flux + 0. * wavenumber);
}

double Source::BackgroundSpectrum(double wavenumber)
{
    return (backgroundflux + 0. * wavenumber);
}

void Source::InitFluxes(char band, double background_magnitude,
                        double sky_background_aperture)
{
    // Flux density for a magnitude 0 star in nph.m.m-2.s-1
    // Multiply by waveband width (in wavenumber unit = m-1 ), collecting
    // surface, and integration time to get the number of photons
    // see e.g. Campins, Reike, & Lebovsky (1985) -- M0 fluxes in Jansky =
    // 1590, 1080, 667 for J, H , K
    double flux_m0 = 0.0;

    /// \todo Implement more than just JHK bands here.  Needs to work at optical wavelengths too.
    if (band == 'J')
        flux_m0 = 1.59e-23 / 6.626068e-34 * 1.235e-6;   // J band
    if (band == 'H')
        flux_m0 = 1.08e-23 / 6.626068e-34 * 1.635e-6;   // H band
    if (band == 'K')
        flux_m0 = 6.67e-24 / 6.626068e-34 * 2.159e-6;   // K band
    // Flux density for the chosen magnitude
    this->flux = flux_m0 * pow(10., -0.4 * magnitude);

    // Background noise calculations
    // The background magnitude is in magnitude/ sq arcsec and assumed given
    // in the same band as the source
    // The sky background aperture is in arcseconds on sky (MROI = 0.5 )
    this->backgroundflux = pow(sky_background_aperture, 2) * flux_m0 * pow(10.,-0.4 * background_magnitude);
    double SNR = flux / backgroundflux;

    cout << "Flux = " << this->flux << " photons per m, per second." << endl;
    cout << "Background flux = " << this->backgroundflux
        << " photons per m, per second." << endl;
    cout << "Source SNR: " << SNR << endl;
}

// AS 2010-06-23
// modified the constructor to compile with the new members of the source
// class, ie target file name, band, background_magnitude 
// and sky_background_aperture
// Source class is being called in simulator.cpp: 
// Source target = Source(NULL, 'H', magnitude, 2600., 14.7, 0.5, 0.0);
Source::Source(const char *filename, char band, double magnitude,
               double temperature, double background_magnitude,
               double sky_background_aperture, double declination, double right_ascension, double pixellation)
{
    this->source_name = filename;       // name of the file containing the
                                        // fake image of the source
    this->band = band;          // band in which the image was taken
    this->magnitude = magnitude;        // completely defines the star flux
    this->temperature = temperature;    // completely defines the spectrum
    this->declination = declination;    // declination of the source
    this->right_ascension = right_ascension;
    this->background_magnitude = background_magnitude;  // background
                                                        // magnitude of the
                                                        // sky
    this->sky_background_aperture = sky_background_aperture;    // sky
                                                                // background
                                                                // aperture
    InitFluxes(this->band, this->background_magnitude, this->sky_background_aperture);  // Call 
                                                                                        // flux 
                                                                                        // function 
                                                                                        // and 
                                                                                        // passes 
                                                                                        // values 
                                                                                        // defined 
                                                                                        // above
                                                                                        
    this->source_pixellation = pixellation;

    if (filename != NULL)       // if there is a file then
    {
        // extended source
        // Get the image and pixellation (keyword PIXEL should be in the FITS
        // file)
        fits2mat(filename, source_image);

        // Normalize it
        source_image /= total(source_image);
    }
}

// AS 2010-06-23
// modified the constructor to compile with the new members of the source
// class, ie target file name, band, background_magnitude 
// and sky_background_aperture
// Input values into class passed via source file
Source::Source(string filename, string comment_chars)
{
    vector < string > lines = ReadFile(filename, comment_chars, "Cannot Open Observation Definition File");

    if (lines.size() != 9)
        throw std::
            runtime_error("Invalid number of parameters in the source file");

    string fits_filename = lines[0];

    // AS 2010-06-24
    // added forced exit of the programme if the parameters are of the wrong
    // type
    // \todo Allow point source to be specified
    if (lines[0].find(".fits") == string::npos)
    {
        /// \exception runtime_error Invalid FITS image filename
        /// The file specified by the user evidently does not exist.
        throw std::runtime_error("Invalid FITS image filename");
    }
    else
    {
        source_name = lines[0];
    }
    
    // Pull out the pixellation constant:
    if(!isdigit(lines[1][0]))
    {
        throw std::runtime_error("Invalid Pixellation Definition.");
    }
    else
    {
        this->source_pixellation = atof(lines[1].c_str());
        cout << "Image pixellation = " << source_pixellation << " mas/pixel" << endl;
    }

    if (!(isalpha(lines[2][0])))
    {   
        /// \exception runtime_error Invalid source band
        /// The source band specified bythe user is not valid.
        throw std::runtime_error("Invalid source band");
    }
    else
    {
        band = lines[2][0];
    }

    if (!isdigit(lines[3][0]))
    {
        /// \exception runtime_error Invalid source magnitude
        /// The source magnitude was not defined
        throw std::runtime_error("Invalid source magnitude");
    }
    else
    {
        magnitude = atof(lines[3].c_str());
    }

    if (!isdigit(lines[4][0]))
    {
        /// \exception runtime_error Invalid source temperature
        /// 
        throw std::runtime_error("Invalid source temperature");
    }
    else
    {
        temperature = atof(lines[4].c_str());
    }

    if (!isdigit(lines[5][0]))
    {
        /// \exception runtime_error Invalid sky background magnitude
        /// 
        throw std::runtime_error("Invalid sky background magnitude");
    }
    else
    {
        background_magnitude = atof(lines[5].c_str());
    }

    if (!isdigit(lines[6][0]))
    {
        /// \exception runtime_error Invalid sky background aperture
        /// 
        throw std::runtime_error("Invalid sky background aperture");
    }
    else
    {
        sky_background_aperture = atof(lines[6].c_str());
    }

    if (!isdigit(lines[7][0]))
    {
        /// \exception runtime_error Invalid source declination
        /// 
        throw std::runtime_error("Invalid source declination");
    }
    else
    {
        declination = atof(lines[7].c_str());
    }
    if (!isdigit(lines[8][0]))
    {
        /// \exception runtime_error Invalid source declination
        /// 
        throw std::runtime_error("Invalid source right ascention");
    }
    else
    {
        right_ascension = atof(lines[8].c_str());
    }


    InitFluxes(band, background_magnitude, sky_background_aperture);    // calls 
                                                                        // flux 
                                                                        // equation

    if (filename.size() > 0)
    {
        // extended source
        // Get the image and pixellation (keyword PIXEL should be in the FITS
        // file)
        fits2mat(fits_filename.c_str(), source_image);

        // Normalize it
        source_image /= total(source_image);
    }
}

string  Source::GetName(void)
{
    return this->source_name;
}
