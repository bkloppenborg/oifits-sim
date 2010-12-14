#include "NoiseModel.h"

#include "Simulator.h"
#include "SpectralMode.h"
#include "Instrument.h"
#include "Array.h"
#include "Source.h"
#include "Observation.h"

// AS 2010-06-24
// removed redundant function Flux0 and the limits of the channels as they 
// 
// were used in Flux0

/// function for setting the correct values of the Strehl ratio according
/// to the
/// wavelenght, the astronomical seeing at 500 nm and the station diameter
/// Assumes 1.4m telescopes and 1 arcsec seeing
/// AS 2010-06-22
/// modified the function to match equation given in the error budget
/// spreadsheet
double Strehl(double wl, double r0, double D)
{
	double strehl = 0.0;		// Strehl ratio
	double r0lambda = 0.0;		// astronomical seeing at current
	// wavelength
	r0lambda = r0 * pow(wl / 5e-7, 6.0 / 5.0);
	strehl = exp(-0.134 * pow(D / r0lambda, 5.0 / 3.0));
	return (strehl);
}

/// function for setting the correct values of the integration time
/// according to
/// the mean wavelength, the Fried parameter at 500 nm and the seeing wind
/// speed 
// AS 2010-06-23
// modified the function to match equation given in the error budget
// spreadsheet
// now takes as an argument the median wavelength
double TimeInt(double median_wl, double r0, double v)
{
	double ti = 0.0;			// integration time
	double r0lambda = 0.0;		// Fried parameter at current
	// wavelength
	r0lambda = r0 * pow(median_wl / 5e-7, 6.0 / 5.0);

	// formula taken from Bridget O'Donovan's Phd Theis
	// see equation (1.23) 
	ti = 2 * 0.314 * r0lambda / v;
	return (ti);
}


double PhCount(Source & target, Observation & obs, Instrument & inst, double median_wl, double wl, double Dlambda)
{
	double Ni;

	Ni = target.flux * pow(wl, -2.0) * TimeInt(median_wl, inst.ast_seeing, inst.wind_speed) *
	   obs.GetNumStations() * PI * pow((obs.GetStation(0)->diameter), 2.0) * 
	   Dlambda * inst.GlobalThroughput() * Strehl(wl, inst.ast_seeing, obs.GetStation(0)->diameter) / 4;
	return (Ni);
}


double CoherentFlux1(double Ni, double Pow, Observation & obs, Instrument & inst)
{
	double Vsq = Pow * inst.visibility * inst.visibility;
	double CFlux = pow(Ni, 2.0) * Vsq / pow(obs.GetNumStations(),
	   2.0) + Ni + inst.Npix * pow(inst.read_noise, 2.0);
	return CFlux;
}

/// compute the variance of the squared coherent flux
double VarCoherentFlux1(double Ni, double Pow, Observation & obs, Instrument & inst)
{
	double Vsq = Pow * inst.visibility * inst.visibility;
    int nstations = obs.GetNumStations();

	// the components are computed separately for readability

	// photon noise term
	double VarCFlux_ph =
	   2 * pow(Ni, 3.0) * Vsq / pow(nstations, 2.0) + 4 * pow(Ni,
	   2.0) * Vsq / pow(nstations, 2.0) + pow(Ni, 2.0) + Ni;

	// detector noise terms
	double VarCFlux_det = pow(inst.Npix, 2.0) * pow(inst.read_noise,
	   4.0) + 3 * inst.Npix * pow(inst.read_noise, 4.0);

	// coupled terms
	double VarCFlux_coupled = 2 * inst.Npix * pow(inst.read_noise,
	   2.0) * Ni + 2 * inst.Npix * pow(inst.read_noise, 2.0) * pow(Ni,
	   2.0) * Vsq / pow(nstations, 2.0);

	// the total variance
	return VarCFlux_ph + VarCFlux_det + VarCFlux_coupled;
}

/// compute the variance of a power spectrum point
// AS 2010-06-22
// eliminated the incoherent integration time which is now contained in
// the instrument class
double VarPow1(double Pow, SpectralMode & spec, int specIndex, Source & target, Observation & obs, Instrument & inst)
{
	int n_coh = inst.incoh_time / TimeInt(spec.median_wavelength,
	   inst.ast_seeing,
	   inst.wind_speed);
	double Ni = PhCount(target, obs, inst, spec.median_wavelength, spec.mean_wavelength[specIndex], 
	    spec.delta_wavelength[specIndex]);
	double Flux = CoherentFlux1(Ni, Pow, obs, inst);
	double VarC = VarCoherentFlux1(Ni, Pow, obs, inst);

	// JSY 2009-10-21: prev expression for VarPow_coh was:
	// pow(Pow * inst.Vinst * inst.Vinst, 2.0) * VarC / pow(Flux, 2.0);
	double VarPow_coh = pow(Pow, 2.0) * VarC / pow(Flux, 2.0);
	return VarPow_coh / n_coh + pow(Pow * inst.vsq_frac_cal_err, 2.0);
}

/// compute the variance of a power spectrum point, assuming unbiased
/// estimator
// added JSY 24/03/2010 for James Gordon
// AS 2010-06-22
// eliminated the incoherent integration time which is now contained in
// the instrument class
//double VarUnbiasedPow1(double Pow, SpectralMode & spec, int specIndex, Source & target, Array & s, Instrument & inst)
//{
//	int n_coh = inst.incoh_time / TimeInt(spec.median_wavelength,
//	   inst.ast_seeing,
//	   inst.wind_speed);
//	double Ni = PhCount(target, s, inst, spec.median_wavelength,
//	   spec.mean_wavelength[specIndex],
//	   spec.delta_wavelength[specIndex]);
//	double VarPow_coh =
//	   (pow(Ni, 2.0) +
//	   2.0 * pow(Ni, 3.0) * pow(inst.visibility, 2.0) / pow(s.nstations,
//		  2.0) + 2.0 * pow(inst.Npix, 2.0) * pow(inst.read_noise,
//		  4.0)) / (pow(inst.visibility * Ni / s.nstations, 4.0));
//	return VarPow_coh / n_coh + pow(Pow * inst.vsq_frac_cal_err, 2.0);
//}

double VarUnbiasedPow1(double Pow, SpectralMode & spec, int specIndex, Source & target, Observation & obs, Instrument & inst)
{
    int nstations = obs.GetNumStations();

	int n_coh = inst.incoh_time / TimeInt(spec.median_wavelength,
	   inst.ast_seeing,
	   inst.wind_speed);
	double Ni = PhCount(target, obs, inst, spec.median_wavelength,
	   spec.mean_wavelength[specIndex],
	   spec.delta_wavelength[specIndex]);
	double VarPow_coh =
	   (pow(Ni, 2.0) +
	   2.0 * pow(Ni, 3.0) * pow(inst.visibility, 2.0) / pow(nstations,
		  2.0) + 2.0 * pow(inst.Npix, 2.0) * pow(inst.read_noise,
		  4.0)) / (pow(inst.visibility * Ni / nstations, 4.0));
	return VarPow_coh / n_coh + pow(Pow * inst.vsq_frac_cal_err, 2.0);
}

/// \TODO: VarCloPhase1() etc.

// Function that computes the variance of the closure phase
// the values depend on: the wavelengths, the visibility, the target, the
// array,
// the instrument, the incoherent integration time,
// the number of bispectra, and the number of hour angles
// Adds contribution from instrument calibration error
// AS 2010-06-22
// eliminated the incoherent integration time which is now contained in
// the instrument class
Matrix < double > VarCloPhase(SpectralMode & spec, Matrix < Complex > &visibility,
   Source & target, Observation & obs, Instrument & inst, int Nbs, int Nha)
{
	int i, n_coh;
	double Ni;					// variable containing the photon count

	Matrix < double >VarPhi(spec.nchannels, Nbs);	// declaring the
	// variance of the 
	// phase
	// components of the variance of the closure phase
	double VarPhi_coh, VarPhi_cal;
	double VarPhi_ph,
	   VarPhi_ph1, VarPhi_ph2, VarPhi_ph3, VarPhi_ph4, VarPhi_ph5;
	double VarPhi_det;
	double VarPhi_det1, VarPhi_det2, VarPhi_det3;

	Matrix < Complex > bis(spec.nchannels, Nbs);	// matrix
	// containing the
	// bispectra

	// some local variables to hold the values of the matrices
	// employed to make the code more readable
	double b = 0, q = 0, p, x = 0, y = 0, z = 0;
	int N = obs.GetNumStations();
	double SigDet = inst.read_noise;
	int Npix = inst.Npix;

	VarPhi_cal = pow((inst.clo_cal_err_deg * PI / 180), 2.0);
	for (int ii = 0; ii < spec.nchannels; ii++)
	{
		i = 0;
		n_coh = inst.incoh_time / TimeInt(spec.median_wavelength, inst.ast_seeing, inst.wind_speed);
		Ni = PhCount(target, obs, inst, spec.median_wavelength, spec.mean_wavelength[ii], spec.delta_wavelength[ii]);
		p = N / Ni;
		for (int jj = 0; jj < Nha; jj++)
		{
			for (int j = 0; j < N - 2; j++)
			{
				for (int k = j + 1; k < N - 1; k++)
				{
					bis[ii][i] =
					   visibility[ii][j +
					   jj * N * (N - 1) / 2]
					   * conj(visibility[ii]
					   [k + jj * N * (N - 1) / 2])
					   * visibility[ii][i + (N - 1) * (jj + 1)];
					x = norm(visibility[ii][j + jj * N * (N - 1) / 2]);
					y = norm(visibility[ii][k + jj * N * (N - 1) / 2]);
					z = norm(visibility[ii][i + (N - 1) * (jj + 1)]);
					b = abs(bis[ii][i]);
					q = 1 / (2 * norm(bis[ii][i]));

					// some partial results, written like this for
					// readability
					VarPhi_ph1 = (N * N * N - 2 * b) * p * p * p * q;
					VarPhi_ph2 = N * N * (x + y + z) * p * p * q;
					VarPhi_ph3 =
					   -(x * x + y * y + z * z +
					   2 * x * y + 2 * x * z + 2 * y * z) * p * p * q;
					VarPhi_ph4 = (x * y + y * z + z * x) * N * p * q;
					VarPhi_ph5 = -(x + y + z) * 2 * b * p * q;
					// photon noise
					VarPhi_ph =
					   VarPhi_ph1 + VarPhi_ph2 +
					   VarPhi_ph3 + VarPhi_ph4 + VarPhi_ph5;

					// some partial results, written like this for
					// readability
					VarPhi_det1 =
					   (Npix * Npix * Npix +
					   3 * Npix * Npix) * pow(p, 6.0) * pow(SigDet, 6.0);
					VarPhi_det2 =
					   (x + y + z) * (3 * Npix +
					   Npix * Npix) * pow(SigDet, 4.0) * pow(p, 4.0) * q;
					VarPhi_det3 =
					   (x * y + x * z +
					   y * z) * Npix * SigDet * SigDet * p * p * q;
					// detector noise
					VarPhi_det = VarPhi_det1 + VarPhi_det2 + VarPhi_det3;

					// total variance per coherent integration
					VarPhi_coh = VarPhi_ph + VarPhi_det;

					// total variance after incoherent integration and
					// calibration
					VarPhi[ii][i] = VarPhi_coh / n_coh + VarPhi_cal;

					i++;
				}
			}
		}
	}
	return (VarPhi);
}
