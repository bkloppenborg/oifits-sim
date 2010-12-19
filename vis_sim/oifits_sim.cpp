/*
 * simulator.cpp
 *
 *  Created on: 21 Jul 2009
 *      Author: Andra Stroe
 */

// Include built-in libraries
#include <cstdlib>
#include <string>

#include <iostream>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "oifits_sim.h"

// Header files for remaining classes.
#include "Array.h"
#include "Baseline.h"
//#include "Bispectrum.h"   // removed
#include "Instrument.h"
#include "NoiseModel.h"
#include "PowerSpectrum.h"
#include "random.h"
#include "Simulator.h"
#include "Source.h"
#include "SpectralMode.h"
#include "UVPoint.h"
#include "VisSimParams.h"
// Various includes for the observation types
#include "Observation.h"
#include "Obs_HA.h"
#include "Obs_OIFITS.h"


using std::cout;

// function that computes power spectra from the complex visibilities
//Matrix < double > Vis2Pow(Matrix < Complex > &visibility, int Nwav, int Npow)
//{
//	Matrix < double >pow(Nwav, Npow);

//	for (int i = 0; i < Nwav; i++)
//	{
//		for (int j = 0; j < Npow; j++)
//		{						// power spectrum as the
//			// norm of the visibility
//			pow[i][j] = norm(visibility[i][j]);
//		}
//	}
//	return (pow);
//}

// function that computes bispectra from the complex visibilities
// the values depend on: the visibility, the number of telescopes, the
// number of
// wavelengths, the number of bispectra and the number of hour angles
//Matrix < Complex > Vis2Bis(Matrix < Complex > &visibility, int N, int Nwav, int Nbs, int Nha)
//{
//	int i = 0;					// just a counter
//	// declaring the bispectrum under the variable name bis; it is a
//	// complex matrix
//	Matrix < Complex > bis(Nwav, Nbs);

//	for (int ii = 0; ii < Nwav; ii++)
//	{
//		i = 0;
//		for (int jj = 0; jj < Nha; jj++)
//		{
//			for (int j = 0; j < N - 2; j++)
//			{
//				for (int k = j + 1; k < N - 1; k++)
//				{
//					bis[ii][i] = visibility[ii][j + jj * N * (N - 1) / 2]
//					   * conj(visibility[ii][k + jj * N * (N - 1) / 2])
//					   * visibility[ii][i + (N - 1) * (jj + 1)];
//					i++;
//				}
//			}
//		}
//	}
//	return (bis);
//}

// function computing the noisy power spectra
// the values depend on a series of parameters: the wavelengths, the
// visibility,
// the target, the array, the instrument, the incoherent integration time,
// the number of power spectra, the
// telescope indexes defining the baselines and the u and v coordinates of
//// the baselines
//PowerSpectrum GenPower(SpectralMode & spec, Matrix < Complex > &visibility,
//   Source & target, Observation & obs, Instrument & inst, int Npow,
//   Matrix < int >&t1, Matrix < int >&t2, Matrix < double >&u,
//   Matrix < double >&v)
//{
//	PowerSpectrum Power;
//	Matrix < double >Pow(spec.nchannels, Npow);
//	Matrix < double >Err(spec.nchannels, Npow);
//	Matrix < double >noisyPow(spec.nchannels, Npow);

//	Pow = Vis2Pow(visibility, spec.nchannels, Npow);

//	for (int i = 0; i < spec.nchannels; i++)
//	{
//		for (int j = 0; j < Npow; j++)
//		{
//			Err[i][j] = sqrt(VarUnbiasedPow1(Pow[i][j], spec, i, target, obs, inst));
//			noisyPow[i][j] = Pow[i][j] + Err[i][j] * Rangauss(random_seed);
//		}
//	}

//	// assigning the values to different members of the class
//	// PowerSpectrum
//	// for the object Power
//	Power.pow = Vis2Pow(visibility, spec.nchannels, Npow);
//	Power.vis2err = Err;
//	Power.vis2data = noisyPow;
//	Power.u = u;
//	Power.v = v;
//	Power.t1 = t1;
//	Power.t2 = t2;
//	Power.int_time = inst.incoh_time;
//	/*
//	 * cout << "POWER SPECTRUM DATA: " << endl; cout << "T1 and T2 are: "
//	 * << endl << Power.t1 << endl << Power.t2 << endl; cout << "u and v
//	 * are: " << endl << Power.u << endl << Power.v << endl; cout << "The
//	 * power spectrum is: " << endl << Power.pow << endl; cout << "The
//	 * error of the power spectrum is: " << endl << Power.vis2err << endl;
//	 * cout << "The power spectrum with error is: " << endl <<
//	 * Power.vis2data << endl; 
//	 */
//	return (Power);
//}

// function that computes the triple amplitude with added noise and the
// closure
// phase with added noise
// the values depend on: the wavelengths, the visibility, the Target, the
// array,
// the instrument, the incoherent integration time,
// the number of triangles and sampled hour angles,
// the telescope indexes and the UV coordinates
//Bispectrum GenBispec(SpectralMode & spec, Matrix < Complex > &visibility,
//   Source & target, Observation & obs, Instrument & inst, int Nbs, int Nha,
//   Matrix < int >&t1, Matrix < int >&t2, Matrix < double >&u,
//   Matrix < double >&v)
//{
//	Bispectrum bispec;
//	int i = 0;

//    int nstations = obs.GetNumStations();

//	// matrix containing the U coordinate of the 1st baseline of the
//	// triangle
//	Matrix < double >u1(spec.nchannels, Nbs);
//	// matrix containing the V coordinate of the 1st baseline of the
//	// triangle
//	Matrix < double >v1(spec.nchannels, Nbs);
//	// matrix containing the U coordinate of the 2nd baseline of the
//	// triangle
//	Matrix < double >u2(spec.nchannels, Nbs);
//	// matrix containing the V coordinate of the 2nd baseline of the
//	// triangle
//	Matrix < double >v2(spec.nchannels, Nbs);
//	// matrix containing the indexes of the 1st telescope in the triangle
//	Matrix < int >tel1(spec.nchannels, Nbs);
//	// matrix containing the indexes of the 2nd telescope in the triangle
//	Matrix < int >tel2(spec.nchannels, Nbs);
//	// matrix containing the indexes of the 3rd telescope in the triangle
//	Matrix < int >tel3(spec.nchannels, Nbs);

//	// matrix containing the bispectrum
//	Matrix < Complex > bis(spec.nchannels, Nbs);

//	// matrix containing the triple amplitude
//	Matrix < double >T3true(spec.nchannels, Nbs);
//	// matrix containing the err of the triple amplitude
//	Matrix < double >ErrT3(spec.nchannels, Nbs);
//	// matrix containing the triple amplitude with added noise
//	Matrix < double >T3(spec.nchannels, Nbs);

//	// matrix containing the closure phase
//	Matrix < double >Phi(spec.nchannels, Nbs);
//	// matrix containing the variance of the closure phase
//	Matrix < double >VarPhi(spec.nchannels, Nbs);
//	// matrix containing the error of the closure phase
//	Matrix < double >ErrPhi(spec.nchannels, Nbs);
//	// matrix containing the closure phase with added noise
//	Matrix < double >PhiErr(spec.nchannels, Nbs);

//	bis = Vis2Bis(visibility, nstations, spec.nchannels, Nbs, Nha);

//	// separate bispectrum into amplitude and phase
//	for (int i = 0; i < spec.nchannels; i++)
//	{
//		for (int j = 0; j < Nbs; j++)
//		{
//			Phi[i][j] = arg(bis[i][j]);
//			T3true[i][j] = abs(bis[i][j]);
//		}
//	}

//	// add noise
//	VarPhi = VarCloPhase(spec, visibility, target, obs, inst, Nbs, Nha);
//	for (int i = 0; i < spec.nchannels; i++)
//	{
//		for (int j = 0; j < Nbs; j++)
//		{
//			ErrPhi[i][j] = sqrt(VarPhi[i][j]);
//			PhiErr[i][j] = Phi[i][j] + ErrPhi[i][j] * Rangauss(random_seed);
//			// assume circular noise cloud
//			ErrT3[i][j] = sqrt(T3true[i][j] * T3true[i][j] * VarPhi[i][j]);
//			T3[i][j] = T3true[i][j] + ErrT3[i][j] * Rangauss(random_seed);
//		}
//	}

//	for (int ii = 0; ii < spec.nchannels; ii++)
//	{
//		i = 0;
//		for (int jj = 0; jj < Nha; jj++)
//		{
//			for (int j = 0; j < nstations - 2; j++)
//			{
//				for (int k = j + 1; k < nstations - 1; k++)
//				{
//					u1[ii][i] =
//					   u[ii][j + jj * nstations * (nstations - 1) / 2];
//					v1[ii][i] =
//					   v[ii][j + jj * nstations * (nstations - 1) / 2];
//					u2[ii][i] = u[ii][i + (jj + 1) * (nstations - 1)];
//					v2[ii][i] = v[ii][i + (jj + 1) * (nstations - 1)];
//					tel1[ii][i] =
//					   t1[ii][j + jj * nstations * (nstations -
//						  1) / 2];
//					tel2[ii][i] =
//					   t2[ii][j + jj * nstations * (nstations -
//						  1) / 2];
//					tel3[ii][i] =
//					   t2[ii][k + jj * nstations * (nstations -
//						  1) / 2];
//					i++;
//				}
//			}
//		}
//	}
//	bispec.u1 = u1;
//	bispec.v1 = v1;
//	bispec.u2 = u2;
//	bispec.v2 = v2;
//	bispec.t1 = tel1;
//	bispec.t2 = tel2;
//	bispec.t3 = tel3;

//	bispec.trueAmp = T3true;
//	bispec.t3amperr = ErrT3;
//	bispec.t3amp = T3;

//	bispec.truePhi = Phi * (180 / PI);
//	bispec.t3phierr = ErrPhi * (180 / PI);
//	bispec.t3phi = PhiErr * (180 / PI);
//	bispec.int_time = inst.incoh_time;

//	/*
//	 * cout << "BISPECTRUM DATA: " << endl; cout << "U points and V points 
//	 * are:" << endl << bispec.u1 << endl << bispec.v1 << endl; cout <<
//	 * "Telescope 1 , 2 and 3 are:" << endl << bispec.t1 << endl <<
//	 * bispec.t2 << endl << bispec.t3 << endl; cout << "The triple
//	 * amplitude is: " << endl << bispec.t3true << endl; cout << "The
//	 * error of the triple amplitude is: " << endl << bispec.t3amperr <<
//	 * endl; cout << "The triple amplitude with error is: " << endl <<
//	 * bispec.t3amp << endl;
//	 * 
//	 * cout << "The closure phase is: " << endl << bispec.phi << endl;
//	 * cout << "The error of the closure phase is: " << endl <<
//	 * bispec.t3phierr << endl; cout << "The closure phase with error is:
//	 * " << endl << bispec.t3phi << endl; 
//	 */

//	return (bispec);
//}

// function that computes the visibilities
// as arguments: the wavelengths, the hour angles, the source, the array
// configuration, and the instrument
//void generate_oifits(SpectralMode & spec, HourAngle & ha, Source & target, Array & s, Instrument & inst, string oifits_filename)
//{
//	int Npow;					// number of sampled hour angles
//	int Nbs;					// number of bispectra
//    int nstations = s.GetNumStations();

//	// number of columns in the power spectra matrix
//	Npow = ha.Nha * nstations * (nstations - 1) / 2;
//	// number of columns in the bispectrum matrix
//	Nbs = ha.Nha * (nstations - 1) * (nstations - 2) / 2;

//	// matrix recording the visibilities at each time step
//	Matrix < Complex > visibility(spec.nchannels, Npow);
//	// matrix recording the u coordinate
//	Matrix < double >u(spec.nchannels, Npow);
//	// matrix recording the v coordinate
//	Matrix < double >v(spec.nchannels, Npow);
//	// matrix recording the index of the 1st telescope
//	Matrix < int >t1(spec.nchannels, Npow);
//	// matrix recording the index of the 2nd telescope
//	Matrix < int >t2(spec.nchannels, Npow);
//	// row vector recording the hour angles
//	Row < double >time(Npow);

//	int ii = 0;					// counter running through all the
//	// wavelengths
//	int jj = 0;					// counter running through all the HA

//	double wl = 0;				// temporary variable storing the current
//	// wavelength
//	double wn = 0;				// temporary variable storing the current
//	// wavenumber
//	double h = 0;				// temporary variable storing the current
//	// hour angle

//	UVPoint uv;

//	cout << "Integration time = " << 1e3 * TimeInt(spec.median_wavelength,
//	   inst.ast_seeing, inst.wind_speed) << "ms" << endl;
//	for (int l = 0; l < spec.nchannels; l++)	// loop over the
//		// wavelengths
//	{
//		wl = spec.mean_wavelength[l];
//		wn = spec.mean_wavenumber[l];
//		jj = 0;

//		// print photon count for this spectral channel for user to check
//		double Ni = PhCount(target, s, inst, spec.median_wavelength, wl, spec.delta_wavelength[l]);

//		cout << "Channel " << l << "  Wavelength " << wl <<
//		   " photons = " << Ni << "  (assumed Strehl " << Strehl(wl,
//		   inst.ast_seeing, s.diameter[0]) << ")" << endl;

//		// 
//		double snr = pow(inst.visibility * Ni / nstations,
//		   2.0) / pow((pow(Ni, 2.0) + 2.0 * pow(Ni,
//				 3.0) * pow(inst.visibility,
//				 2.0) / pow(nstations,
//				 2.0) + 2.0 * pow(inst.Npix,
//				 2.0) * pow(inst.read_noise,
//				 4.0)), 0.5);
//		cout << "  point source unbiased V^2 SNR = " << snr << endl;
//		// AS 2010-06-22
//		// modified the VarUnbiasedPow1 function by eliminating the
//		// incoherent integration time
//		cout << "  point source unbiased V^2 SNR = " <<
//		   pow(VarUnbiasedPow1(1.0, spec, l, target, s, inst),
//		   -0.50) << endl;
//		// AS 2010-06-22
//		// modified the VarPow1 function by eliminating the incoherent
//		// integration time
//		cout << "  point source biased   V^2 SNR = " <<
//		   pow(VarPow1(1.0, spec, l, target, s, inst), -0.50) << endl;
//		// 

//		for (int k = 0; k < ha.HA.size(); k++)	// for loop over all the
//			// hour angles
//		{
//			h = ha.HA[k];
//			// for loop over the 1st telescope defining the baseline
//			for (int i = 0; i < (s.nstations - 1); i++)
//			{
//				// for loop over the 2nd telescope defining the baseline
//				for (int j = i + 1; j < s.nstations; j++)
//				{
//					Baseline B(s, i, j);
//					uv = B.UVcoords(h, target.declination, wn);
//					visibility[ii][jj] = target.GetVis(uv);
//					u[ii][jj] = uv.u * wl;	// want u, v in meters
//					v[ii][jj] = uv.v * wl;
//					t1[ii][jj] = i;
//					t2[ii][jj] = j;
//					time[jj] = h;
//					jj++;		// counter running through all the columns
//				}
//			}
//		}
//		ii++;					// counter running through all the rows
//	}

//	cout << "Npow is:" << Npow << endl;
//	cout << "Nbs is:" << Nbs << endl;
//	PowerSpectrum Power = GenPower(spec, visibility, target, s, inst, Npow, t1, t2, u, v);

//	// AS 2010-06-21
//	// fixed :Bug: segfault if Nbs zero
//	if (Nbs == 0)
//	{
//		cout << "Only 2 telescopes, no bispectra computed." << endl;
//		write_oifits_file(target.declination, oifits_filename, s, spec,
//		   Power, Npow, NULL, Nbs, time);
//	}
//	else
//	{
//		Bispectrum Bispectrum =
//		   GenBispec(spec, visibility, target, s, inst,
//		   Nbs, ha.Nha, t1, t2, u, v);
//		write_oifits_file(target.declination, oifits_filename, s, spec,
//		   Power, Npow, &Bispectrum, Nbs, time);
//	}
//}

void run_sim(const VisSimParams * p)
{
    /// \todo Check that all parameters are specified here instead of in the gui_main and main functions.
    /// \todo Pass comment_chars to file reading functions rather than have it be specified differently everywhere.
    const string comment_chars("\\#~$&£%");

    // Now pull out the information that must be in a file:
    // We use pointers in this function for consistency across all variables.
    Instrument * inst = new Instrument(p->instrument_filename, comment_chars);
    
    /// \todo Read in the spectral dispersion information from the OIFITS file
	SpectralMode * spec = new SpectralMode(p->SpectralMode_filename, comment_chars);
	Source * target = new Source(p->target_filename, comment_chars);
	
	// Allocate pointers for the array and hour angle information:
	Array * array = new Array(p->array_filename, comment_chars);

    // Read in the observations.
    /// \todo read in the file format type, right now it's locked to the descriptive only format.
    /// \bug Note, when observations is destroyed it will not free memory occupied by the observation
    ///     objects.  We must do this explicitly in the code below.    
    vector<Observation*> observations = Observation::ReadObservations(array, p->input_oifits_filename, comment_chars, OIFITS);
    
    // Open up the OIFITS file.
	string filename = "!test.oifits";
	fitsfile *fptr;
	int status = 0;
	fits_create_file(&fptr, filename.c_str(), &status);
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
    
    /// \todo Use an iterator for this instead, much cleaner.
    /// For now iterate from the back of the vector to the front.
    for(unsigned int i = observations.size(); i > 0; i--)
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
        
        vis2table = observation->GetVis2(spec->insname, *target, wavenumbers);
        write_oi_vis2(fptr, vis2table, 1, &status);
        
        if(observation->HasTriplets())
        {
            t3table = observation->GetT3(spec->insname, *target, wavenumbers);
            write_oi_t3(fptr, t3table, 1, &status);
        }

        // All done with this observation object.  Pop it off the vector and free memory.
        observations.pop_back();
        delete observation;
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
    
    // Free up memory
    delete inst;
    delete spec;
    delete target;
    delete array;
//    delete ha;
}
