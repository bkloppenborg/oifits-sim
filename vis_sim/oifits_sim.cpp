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
#include <list>

#include "oifits_sim.h"

// Header files for remaining classes.
#include "Array.h"
#include "Baseline.h"
#include "Bispectrum.h"
#include "Instrument.h"
#include "NoiseModel.h"
#include "PowerSpectrum.h"
#include "random.h"
#include "Simulator.h"
#include "Source.h"
#include "SpectralMode.h"
#include "UVPoint.h"
#include "VisSimParams.h"
#include "Observation.h"


using std::cout;

// function that computes power spectra from the complex visibilities
Matrix < double > Vis2Pow(Matrix < Complex > &visibility, int Nwav, int Npow)
{
	Matrix < double >pow(Nwav, Npow);

	for (int i = 0; i < Nwav; i++)
	{
		for (int j = 0; j < Npow; j++)
		{						// power spectrum as the
			// norm of the visibility
			pow[i][j] = norm(visibility[i][j]);
		}
	}
	return (pow);
}

// function that computes bispectra from the complex visibilities
// the values depend on: the visibility, the number of telescopes, the
// number of
// wavelengths, the number of bispectra and the number of hour angles
Matrix < Complex > Vis2Bis(Matrix < Complex > &visibility, int N, int Nwav, int Nbs, int Nha)
{
	int i = 0;					// just a counter
	// declaring the bispectrum under the variable name bis; it is a
	// complex matrix
	Matrix < Complex > bis(Nwav, Nbs);

	for (int ii = 0; ii < Nwav; ii++)
	{
		i = 0;
		for (int jj = 0; jj < Nha; jj++)
		{
			for (int j = 0; j < N - 2; j++)
			{
				for (int k = j + 1; k < N - 1; k++)
				{
					bis[ii][i] = visibility[ii][j + jj * N * (N - 1) / 2]
					   * conj(visibility[ii][k + jj * N * (N - 1) / 2])
					   * visibility[ii][i + (N - 1) * (jj + 1)];
					i++;
				}
			}
		}
	}
	return (bis);
}

// function computing the noisy power spectra
// the values depend on a series of parameters: the wavelengths, the
// visibility,
// the target, the array, the instrument, the incoherent integration time,
// the number of power spectra, the
// telescope indexes defining the baselines and the u and v coordinates of
// the baselines
PowerSpectrum GenPower(SpectralMode & spec, Matrix < Complex > &visibility,
   Source & target, Observation & obs, Instrument & inst, int Npow,
   Matrix < int >&t1, Matrix < int >&t2, Matrix < double >&u,
   Matrix < double >&v)
{
	PowerSpectrum Power;
	Matrix < double >Pow(spec.nchannels, Npow);
	Matrix < double >Err(spec.nchannels, Npow);
	Matrix < double >noisyPow(spec.nchannels, Npow);

	Pow = Vis2Pow(visibility, spec.nchannels, Npow);

	for (int i = 0; i < spec.nchannels; i++)
	{
		for (int j = 0; j < Npow; j++)
		{
			Err[i][j] = sqrt(VarUnbiasedPow1(Pow[i][j], spec, i, target, obs, inst));
			noisyPow[i][j] = Pow[i][j] + Err[i][j] * Rangauss(random_seed);
		}
	}

	// assigning the values to different members of the class
	// PowerSpectrum
	// for the object Power
	Power.pow = Vis2Pow(visibility, spec.nchannels, Npow);
	Power.vis2err = Err;
	Power.vis2data = noisyPow;
	Power.u = u;
	Power.v = v;
	Power.t1 = t1;
	Power.t2 = t2;
	Power.int_time = inst.incoh_time;
	/*
	 * cout << "POWER SPECTRUM DATA: " << endl; cout << "T1 and T2 are: "
	 * << endl << Power.t1 << endl << Power.t2 << endl; cout << "u and v
	 * are: " << endl << Power.u << endl << Power.v << endl; cout << "The
	 * power spectrum is: " << endl << Power.pow << endl; cout << "The
	 * error of the power spectrum is: " << endl << Power.vis2err << endl;
	 * cout << "The power spectrum with error is: " << endl <<
	 * Power.vis2data << endl; 
	 */
	return (Power);
}

// function that computes the triple amplitude with added noise and the
// closure
// phase with added noise
// the values depend on: the wavelengths, the visibility, the Target, the
// array,
// the instrument, the incoherent integration time,
// the number of triangles and sampled hour angles,
// the telescope indexes and the UV coordinates
Bispectrum GenBispec(SpectralMode & spec, Matrix < Complex > &visibility,
   Source & target, Observation & obs, Instrument & inst, int Nbs, int Nha,
   Matrix < int >&t1, Matrix < int >&t2, Matrix < double >&u,
   Matrix < double >&v)
{
	Bispectrum bispec;
	int i = 0;

    int nstations = obs.GetNumStations();

	// matrix containing the U coordinate of the 1st baseline of the
	// triangle
	Matrix < double >u1(spec.nchannels, Nbs);
	// matrix containing the V coordinate of the 1st baseline of the
	// triangle
	Matrix < double >v1(spec.nchannels, Nbs);
	// matrix containing the U coordinate of the 2nd baseline of the
	// triangle
	Matrix < double >u2(spec.nchannels, Nbs);
	// matrix containing the V coordinate of the 2nd baseline of the
	// triangle
	Matrix < double >v2(spec.nchannels, Nbs);
	// matrix containing the indexes of the 1st telescope in the triangle
	Matrix < int >tel1(spec.nchannels, Nbs);
	// matrix containing the indexes of the 2nd telescope in the triangle
	Matrix < int >tel2(spec.nchannels, Nbs);
	// matrix containing the indexes of the 3rd telescope in the triangle
	Matrix < int >tel3(spec.nchannels, Nbs);

	// matrix containing the bispectrum
	Matrix < Complex > bis(spec.nchannels, Nbs);

	// matrix containing the triple amplitude
	Matrix < double >T3true(spec.nchannels, Nbs);
	// matrix containing the err of the triple amplitude
	Matrix < double >ErrT3(spec.nchannels, Nbs);
	// matrix containing the triple amplitude with added noise
	Matrix < double >T3(spec.nchannels, Nbs);

	// matrix containing the closure phase
	Matrix < double >Phi(spec.nchannels, Nbs);
	// matrix containing the variance of the closure phase
	Matrix < double >VarPhi(spec.nchannels, Nbs);
	// matrix containing the error of the closure phase
	Matrix < double >ErrPhi(spec.nchannels, Nbs);
	// matrix containing the closure phase with added noise
	Matrix < double >PhiErr(spec.nchannels, Nbs);

	bis = Vis2Bis(visibility, nstations, spec.nchannels, Nbs, Nha);

	// separate bispectrum into amplitude and phase
	for (int i = 0; i < spec.nchannels; i++)
	{
		for (int j = 0; j < Nbs; j++)
		{
			Phi[i][j] = arg(bis[i][j]);
			T3true[i][j] = abs(bis[i][j]);
		}
	}

	// add noise
	VarPhi = VarCloPhase(spec, visibility, target, obs, inst, Nbs, Nha);
	for (int i = 0; i < spec.nchannels; i++)
	{
		for (int j = 0; j < Nbs; j++)
		{
			ErrPhi[i][j] = sqrt(VarPhi[i][j]);
			PhiErr[i][j] = Phi[i][j] + ErrPhi[i][j] * Rangauss(random_seed);
			// assume circular noise cloud
			ErrT3[i][j] = sqrt(T3true[i][j] * T3true[i][j] * VarPhi[i][j]);
			T3[i][j] = T3true[i][j] + ErrT3[i][j] * Rangauss(random_seed);
		}
	}

	for (int ii = 0; ii < spec.nchannels; ii++)
	{
		i = 0;
		for (int jj = 0; jj < Nha; jj++)
		{
			for (int j = 0; j < nstations - 2; j++)
			{
				for (int k = j + 1; k < nstations - 1; k++)
				{
					u1[ii][i] =
					   u[ii][j + jj * nstations * (nstations - 1) / 2];
					v1[ii][i] =
					   v[ii][j + jj * nstations * (nstations - 1) / 2];
					u2[ii][i] = u[ii][i + (jj + 1) * (nstations - 1)];
					v2[ii][i] = v[ii][i + (jj + 1) * (nstations - 1)];
					tel1[ii][i] =
					   t1[ii][j + jj * nstations * (nstations -
						  1) / 2];
					tel2[ii][i] =
					   t2[ii][j + jj * nstations * (nstations -
						  1) / 2];
					tel3[ii][i] =
					   t2[ii][k + jj * nstations * (nstations -
						  1) / 2];
					i++;
				}
			}
		}
	}
	bispec.u1 = u1;
	bispec.v1 = v1;
	bispec.u2 = u2;
	bispec.v2 = v2;
	bispec.t1 = tel1;
	bispec.t2 = tel2;
	bispec.t3 = tel3;

	bispec.trueAmp = T3true;
	bispec.t3amperr = ErrT3;
	bispec.t3amp = T3;

	bispec.truePhi = Phi * (180 / PI);
	bispec.t3phierr = ErrPhi * (180 / PI);
	bispec.t3phi = PhiErr * (180 / PI);
	bispec.int_time = inst.incoh_time;

	/*
	 * cout << "BISPECTRUM DATA: " << endl; cout << "U points and V points 
	 * are:" << endl << bispec.u1 << endl << bispec.v1 << endl; cout <<
	 * "Telescope 1 , 2 and 3 are:" << endl << bispec.t1 << endl <<
	 * bispec.t2 << endl << bispec.t3 << endl; cout << "The triple
	 * amplitude is: " << endl << bispec.t3true << endl; cout << "The
	 * error of the triple amplitude is: " << endl << bispec.t3amperr <<
	 * endl; cout << "The triple amplitude with error is: " << endl <<
	 * bispec.t3amp << endl;
	 * 
	 * cout << "The closure phase is: " << endl << bispec.phi << endl;
	 * cout << "The error of the closure phase is: " << endl <<
	 * bispec.t3phierr << endl; cout << "The closure phase with error is:
	 * " << endl << bispec.t3phi << endl; 
	 */

	return (bispec);
}

// AS 2010-06-21 
// Added to avoid crash of the program in case there are only 2 
// telescopes
// AS 2010-06-24
// merged with the other similar function writing the OIFITS file without
// bispectrum
// now working with pointers for bispectrum, rather than objects
// function writing the data obtained into the OIFITS format
// arguments are: the declination of the source, the datafile name, the
// array,
// the wavelengths, the power spectrum, the number of power spectra, the
// bispectrum,
// the number of bispectra and the hour angles (time)
void write_oifits_file(float declination, string datafile, Array s,
   SpectralMode spec, PowerSpectrum Power, int Npow,
   Bispectrum * pBispec, int Nbs, Row < double >time)
{
/*
//	oi_array array;
//	oi_target targets;
//	oi_wavelength wave;
	oi_vis2 vis2;
	oi_t3 t3;
	fitsfile *fptr;
	char zerostring[FLEN_VALUE];
	char insname[FLEN_VALUE];
	char arrname[FLEN_VALUE];
	char target_filename[FLEN_VALUE];
	int status = 0;
	int i, j, iwav;
	double nulval;
	string filename = "!" + datafile;

	cout << "Generating OIFITS: " << filename << endl;
	// Set to strings
	strncpy(zerostring, " ", FLEN_VALUE);
	strncpy(insname, "Fake_Ins", FLEN_VALUE);
	strncpy(arrname, "Fake_Arr", FLEN_VALUE);
	strncpy(target_filename, "Fake_Targ", FLEN_VALUE);
	// FITS NULL
	nulval = 0.0;
	nulval /= nulval;


	// ARRAY
	// AS 2010-06-17
	// added the array table in the OIFITS file
	double nstations = s.GetNumStations();
	double longitude = s.GetLongitude();
	double latitude = s.GetLatitude();
	double altitude = s.GetAltitude();
	string arrayname = s.GetArrayName();
	

	// WAVE
//	wave.nwave = spec.nchannels;
//	wave.eff_wave = (float *) malloc(wave.nwave * sizeof(float));
//	wave.eff_band = (float *) malloc(wave.nwave * sizeof(float));
//	wave.revision = 1;
//	strncpy(wave.insname, insname, FLEN_VALUE);
//	for (iwav = 0; iwav < wave.nwave; iwav++)
//	{
//		wave.eff_wave[iwav] = spec.mean_wavelength[iwav];
//		wave.eff_band[iwav] = spec.delta_wavelength[iwav];
//	}

	// VIS2
//	vis2.record = (oi_vis2_record *) malloc(Npow * sizeof(oi_vis2_record));
//	for (j = 0; j < Npow; j++)
//	{
//		vis2.record[j].vis2data = (double *) malloc(wave.nwave * sizeof(double));
//		vis2.record[j].vis2err = (double *) malloc(wave.nwave * sizeof(double));
//		vis2.record[j].flag = (char *) malloc(wave.nwave * sizeof(char));
//	}
//	vis2.revision = 1;
//	strncpy(vis2.date_obs, "2009-08-06", FLEN_VALUE);
//	strncpy(vis2.arrname, arrname, FLEN_VALUE);
//	strncpy(vis2.insname, insname, FLEN_VALUE);
//	vis2.numrec = Npow;
//	vis2.nwave = wave.nwave;
//	for (j = 0; j < Npow; j++)
//	{
//		vis2.record[j].target_id = 1;
//		vis2.record[j].time = time[j];
//		vis2.record[j].mjd = 0.0;
//		vis2.record[j].int_time = Power.int_time;
//		vis2.record[j].ucoord = Power.u[0][j];
//		vis2.record[j].vcoord = Power.v[0][j];
//		vis2.record[j].sta_index[0] = Power.t1[0][j];
//		vis2.record[j].sta_index[1] = Power.t2[0][j];
//		for (iwav = 0; iwav < wave.nwave; iwav++)
//		{
//			vis2.record[j].vis2data[iwav] = Power.vis2data[iwav][j];
//			vis2.record[j].vis2err[iwav] = Power.vis2err[iwav][j];
//			vis2.record[j].flag[iwav] = FALSE;
//		}
//	}

	// T3
	if (pBispec != NULL)
	{
		t3.record = (oi_t3_record *) malloc(Nbs * sizeof(oi_t3_record));
		for (j = 0; j < Nbs; j++)
		{
			t3.record[j].t3amp = (double *) malloc(wave.nwave * sizeof(double));
			t3.record[j].t3amperr = (double *) malloc(wave.nwave * sizeof(double));
			t3.record[j].t3phi = (double *) malloc(wave.nwave * sizeof(double));
			t3.record[j].t3phierr = (double *) malloc(wave.nwave * sizeof(double));
			t3.record[j].flag = (char *) malloc(wave.nwave * sizeof(char));
		}
		t3.revision = 1;
		strncpy(t3.date_obs, "2009-04-15", FLEN_VALUE);
		strncpy(t3.arrname, arrname, FLEN_VALUE);
		strncpy(t3.insname, insname, FLEN_VALUE);
		t3.numrec = Nbs;
		t3.nwave = wave.nwave;
		for (j = 0; j < Nbs; j++)
		{
			t3.record[j].target_id = 1;
			t3.record[j].time = time[j];
			t3.record[j].mjd = 0.0;
			t3.record[j].int_time = pBispec->int_time;
			t3.record[j].u1coord = pBispec->u1[0][j];
			t3.record[j].v1coord = pBispec->v1[0][j];
			t3.record[j].u2coord = pBispec->u2[0][j];
			t3.record[j].v2coord = pBispec->v2[0][j];
			t3.record[j].sta_index[0] = pBispec->t1[0][j];
			t3.record[j].sta_index[1] = pBispec->t2[0][j];
			t3.record[j].sta_index[2] = pBispec->t3[0][j];
			for (iwav = 0; iwav < wave.nwave; iwav++)
			{
				t3.record[j].t3amp[iwav] = pBispec->t3amp[iwav][j];
				t3.record[j].t3phi[iwav] = pBispec->t3phi[iwav][j];
				t3.record[j].t3amperr[iwav] = pBispec->t3amperr[iwav][j];
				t3.record[j].t3phierr[iwav] = pBispec->t3phierr[iwav][j];
				t3.record[j].flag[iwav] = FALSE;
			}
		}
	}
	// Write OI-FITS file
	fits_create_file(&fptr, filename.c_str(), &status);
	if (status)
	{
		fits_report_error(stderr, status);
		return;
	}

	write_oi_array(fptr, array, 1, &status);
	write_oi_target(fptr, targets, &status);
	write_oi_vis2(fptr, vis2, 1, &status);
	if (pBispec != NULL)
	{
		write_oi_t3(fptr, t3, 1, &status);
	}
	write_oi_wavelength(fptr, wave, 1, &status);

	if (status)
	{
		fits_delete_file(fptr, &status);
		return;
	}
	else
	{
		fits_close_file(fptr, &status);
	}

	free_oi_array(&array);
	free_oi_target(&targets);
	free_oi_wavelength(&wave);
	free_oi_vis2(&vis2);
	if (pBispec != NULL)
	{
		free_oi_t3(&t3);
	}
	cout << "File written.\n";
	*/
}

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
	SpectralMode * spec = new SpectralMode(p->SpectralMode_filename, comment_chars);
	Source * target = new Source(p->target_filename, comment_chars);
	
	// Allocate pointers for the array and hour angle information:
	Array * array = new Array(p->array_filename, comment_chars);

    // Read in the observations.
    /// \todo read in the file format type, right now it's locked to the descriptive only format.
    vector<Observation> observations = Observation::ReadObservations(array, p->observation_filename, comment_chars, 1);
    
    // Start getting the OIFITS-formatted data:
    oi_array oi_arr = array->GetOIArray();    
    /// \bug The target
    oi_target oi_targ = target->GetOITarget();
    oi_wavelength oi_wave = spec->GetOIWavelength();
    
	// Write OI-FITS file
	string filename = "!test.oifits";
	fitsfile *fptr;
	int status = 0;
	fits_create_file(&fptr, filename.c_str(), &status);
	if (status)
	{
		fits_report_error(stderr, status);
		return;
	}

	write_oi_array(fptr, oi_arr, 1, &status);
	write_oi_target(fptr, oi_targ, &status);
    
    vector<double> wavenumbers = spec->GetWavenumbers();
    
    // Now compute the vis2 records and t3s:
    oi_vis2 vis2table;
    oi_t3 t3table;
    for(unsigned int i = 0; i < observations.size(); i++)
    {
        vis2table = observations[i].GetVis2(spec->insname, *target, wavenumbers);
        write_oi_vis2(fptr, vis2table, 1, &status);
        
        if(observations[i].HasTriplets())
        {
            t3table = observations[i].GetT3(spec->insname, *target, wavenumbers);
            write_oi_t3(fptr, t3table, 1, &status);
        }
        

    }
    


//	write_oi_vis2(fptr, vis2table, 1, &status);
//	if (pBispec != NULL)
//	{
//		write_oi_t3(fptr, t3, 1, &status);
//	}
	write_oi_wavelength(fptr, oi_wave, 1, &status);

	if (status)
	{
		fits_delete_file(fptr, &status);
		return;
	}
	else
	{
		fits_close_file(fptr, &status);
	}

	free_oi_array(&oi_arr);
	free_oi_target(&oi_targ);
	free_oi_wavelength(&oi_wave);
//	free_oi_vis2(&vis2);
//	if (pBispec != NULL)
//	{
//		free_oi_t3(&t3);
//	}
	cout << "File written.\n";    
    
    
    
//    get array info
//    get target info
//    get wavelength info
//    
//    for observation in observations
//        get oi_vis2
//        get oi_t3
//        write_data
//        
//    write oifile
        
	
//	// If we are using an OIFITS file as input, we need to call different constructors (hence the need for pointers here)
//	if(p->bFromOIFITSFile)
//    {
//        oi_fits input_oifits;
//        GList * entry;
//        oi_vis2 * pVis2 = NULL;
//        int status = 0;
//        list<Observation> observations;
//        
//        // Read in the OIFITS file using the library routine.
//        read_oi_fits(p->input_oifits_filename.c_str(), &input_oifits, &status);
//        
//        // Get some information about the data file
//        //print_oi_fits_summary(&input_oifits);
//        
//        // Read out the UV coordinates:
//        // First iterate over the vis2 tables:
//        entry = g_list_first(input_oifits.vis2List);
//        while(entry)
//        {
//            pVis2 = (oi_vis2*) entry->data;
//            oi_vis2_record* records = pVis2->record;
//              
//            // Now iterate over the records in the table
//            for(int j = 0; j < pVis2->numrec; j++)
//            {                                          
//                observations.push_back(Observation(records[j].mjd));
//            }
//            
//            entry = g_list_next(entry);
//        }
//        
//        // Now find the unique mjd,times in the data
//        //printf("List size: %i \n", observations.size());
//        observations.unique(SameObservation);
//        //printf("List size: %i \n", observations.size());
//        
//        // Now generate all list of hour angles from the specified observations:
//        ha = new HourAngle((*target), (*s), observations);
//        
//        // Now that we have what we need, free the memory used by the OIFITS file:
//        free_oi_fits(&input_oifits);
//        
//        //exit(0);
//    }
//    else
//    {
//	    // Now get the information that could come from an OIFITS file
//	    ha = new HourAngle(p->HourAngles_filename.c_str());
//	}
//	
//	
//	// Now call the generator.  We de-reference the pointers here to avoid rewriting subsequent functions.
//    generate_oifits((*spec), (*ha), (*target), (*s), (*inst), p->oifits_filename);
    
    // Free up memory
    delete inst;
    delete spec;
    delete target;
    delete array;
//    delete ha;
}
