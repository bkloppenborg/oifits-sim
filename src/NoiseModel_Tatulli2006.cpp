/*
 * NoiseModel_Tatulli2006.cpp
 *
 *  Created on: Apr 12, 2011
 *      Author: fbaron (via. MROI simulator, rewritten for class by bkloppenborg)
 */

#include "NoiseModel_Tatulli2006.h"

NoiseModel_Tatulli2006::NoiseModel_Tatulli2006()
{
	// TODO Auto-generated constructor stub

}

NoiseModel_Tatulli2006::~NoiseModel_Tatulli2006()
{
	// TODO Auto-generated destructor stub
}

// Computes the squared unbiased error for the visibility squared
double NoiseModel_Tatulli2006::GetVis2Var(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, Baseline * baseline, UVPoint UV, int wavelength_num)
{
	// Here we assume that the noise is composed of two primary terms.  (1) variance in the coherent flux, and (2) variance in the photometric flux
	// We start computing these values by finding how many photons make it to the combiner:
	double n_photons = this->NumPhotons(array, target, combiner, spec_mode, wavelength_num);

	// Now, compute the photometric flux and it's variance.
	double p_flux = n_photons * combiner->GetPhotometryFluxFrac() * combiner->GetPhotometryThroughput();
	double p_flux_var = VarPhotometricFlux(combiner, p_flux);

	UV.Scale(spec_mode->mean_wavenumber[wavelength_num]);

	// Secondly we compute the coherent flux and it's variance.  Start by looking up the squared visibility:
	double vis2 = baseline->GetVis2(*target, UV);
	// and the number of photons making it to the fringe section of the combiner
	double c_photons = n_photons * combiner->GetFringeFluxFrac() * combiner->GetFringeThroughput() / combiner->GetNumSplits();
	// now compute the coherent flux and it's variation:
	double c_flux = CoherentFlux(array, combiner, c_photons, vis2);
	double c_flux_var = VarCoherentFlux(array, combiner, c_photons, vis2);

	// Now compute the variance of the vis2.  Note, we have two telescopes whose variances are assumed equal
	// so that's why we double the second term:
	double error2 = c_flux_var / (c_flux * c_flux) + 2 * p_flux_var / (p_flux * p_flux);

	return sqrt(error2);
}

// Compute the error for the closure phase using an unbiased estimator.  This function includes noise terms from
// photons and the detector only (i.e. NO atmospheric noise!)
double NoiseModel_Tatulli2006::GetT3PhaseVar(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, Triplet * triplet, UVPoint uv1, UVPoint uv2, int wavelength_num)
{
	// Here we assume that the noise may be computed from the detector and photon noise terms added in quadrature.  The paper
	// also assumes that the source is centro-symmetric.

	UVPoint uv3 = uv1 + uv2;

	// First look up the baselines and ask them for the visibility
	Baseline * bl1 = triplet->GetBaseline(0);
	Baseline * bl2 = triplet->GetBaseline(1);
	Baseline * bl3 = triplet->GetBaseline(2);

	// Scale the UV points
	uv1.Scale(spec_mode->mean_wavenumber[wavelength_num]);
	uv2.Scale(spec_mode->mean_wavenumber[wavelength_num]);
	uv3.Scale(spec_mode->mean_wavenumber[wavelength_num]);

	// Now look up the visibilities and bispectra
	double v_ij = bl1->GetVis2(*target, uv1);
	double v_jk = bl2->GetVis2(*target, uv2);
	double v_ik = bl3->GetVis2(*target, uv3);
	double t3amp = triplet->GetT3Amp(*target, uv1, uv2, uv3);

	double n_photons = this->NumPhotons(array, target, combiner, spec_mode, wavelength_num);
	double c_photons = n_photons * combiner->GetFringeFluxFrac() * combiner->GetFringeThroughput() / combiner->GetNumSplits();

	// Pre-compute a few values
	int nt[7] = {0}; 		// Number of telescopes, to the ith power
	double nph[7] = {0.0}; 	// Number of photons to the ith power
	double rn[7] = {0.0};	// Read noise to the ith power
	int npix[7] = {0};

	nt[1] = array->GetNumStations();
	nph[1] = c_photons;
	rn[1] = combiner->read_noise;
	npix[1] = combiner->GetNumFringePixels();

	for(int i = 2; i < 7; i++)
	{
		nt[i] = nt[i-1] * nt[1];
		nph[i] = nph[i-1] * nph[1];
		rn[i] = rn[i-1] * rn[1];
		npix[i] = npix[i-1] * npix[1];
	}

	// Now compute the uncertainty from photons and the detector:
	double var_photon = ComputeT3PhVar(nt, nph, v_ij, v_jk, v_ik, t3amp);
	double var_detect = ComputeT3DetVar(nt, nph, v_ij, v_jk, v_ik, t3amp, rn, npix);

	// Return the variance
	double error2 = var_photon + var_detect;
	return sqrt(error2);
}

// Compute the error for the quad closure -- TBD implement something good here !
double NoiseModel_Tatulli2006::GetT4PhaseVar(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, Quadruplet * quadruplet, 
UVPoint uv1, UVPoint uv2, UVPoint uv3, int wavelength_num)
{
	return 5.;
}

// Computes the photometric flux
double NoiseModel_Tatulli2006::VarPhotometricFlux(Combiner * combiner, double n_photons)
{
	double read_noise_sq = combiner->read_noise * combiner->read_noise;
	return n_photons + combiner->GetNumPhotometricPixels() * read_noise_sq;
}

double NoiseModel_Tatulli2006::CoherentFlux(Array * array, Combiner * combiner, double n_photons, double vis2)
{
	return n_photons * n_photons * vis2 / array->GetNumStations() + n_photons + combiner->GetNumFringePixels() * combiner->read_noise;
}

double NoiseModel_Tatulli2006::VarCoherentFlux(Array * array, Combiner * combiner, double n_photons, double vis2)
{
	// First, pre-compute some values that are used often

	// Read noise squared and to the fourth
	double rn2 = combiner->read_noise * combiner->read_noise;
	double rn4 = rn2 * rn2;
	// Number of photons
	double nph2 = n_photons * n_photons;
	double nph3 = nph2 * n_photons;
	// Number of pixels
	double npi = combiner->GetNumFringePixels();
	double npi2 = npi * npi;
	// Telescopes
	double nt = array->GetNumStations();
	double nt2 = nt * nt;

	// The variance in coherent flux comes from three terms (see paper, equations 6 and 8)
	// (1) photon noise
	double temp = 2 * nph3 * vis2 / nt2 + 4 * nph2 * vis2 / nt2 + nph2 + n_photons;
	// (2) detector noise
	temp += npi2 * rn4 + 3. * npi * rn4;
	// (3) coupled terms
	temp += 2.* npi * rn2 * n_photons + 2.* npi * rn2 * nph2 * vis2 / nt2;

	return temp;
}

double NoiseModel_Tatulli2006::ComputeT3PhVar(int nt[7], double nph[7], double v_ij, double v_jk, double v_ik, double t3amp)
{
	// Implementing equation 16.  Note, we factored out 1 / (2 t3amp^2)
	double temp = (nt[3] - 2 * t3amp) * nt[3] / nph[3]
	            + (nt[2] * (v_ij + v_jk + v_ik)) * nt[2] / nph[2]
	            - (v_ij*v_ij + v_jk*v_jk + v_ik*v_ik + 2*(v_ij*v_jk + v_ij*v_ik + v_jk*v_ik)) * nt[2] / nph[2]
	            + (nt[1] * (v_ij*v_jk + v_ij*v_ik + v_jk*v_ik)) * nt[1] / nph[1]
	            - (2 * t3amp * (v_ij + v_jk + v_ik)) * nt[1] / nph[1];

	// Now multiply by the factored out term:
	temp *= 1. / (2. * t3amp * t3amp);

	if(temp < 0)
		printf("whoops.\n");

	return temp;
}

double NoiseModel_Tatulli2006::ComputeT3DetVar(int nt[7], double nph[7], double v_ij, double v_jk, double v_ik, double t3amp, double rn[7], int npix[7])
{
	//TODO: this could be slighly improved by pre-computing nt/nph values before this function.

	// Implementing equation 17.  Note, we factored out 1 / (2 t3amp^2)
	double temp = (npix[3] * rn[6] + 3*npix[2]*rn[6]) * nt[6]/nph[6]
	            + (v_ij + v_jk + v_ik) * (3*npix[1]*rn[4] + npix[2]*rn[4]) * nt[4]/nph[4]
	            + (v_ij * v_jk + v_ij * v_ik + v_jk*v_ik)*npix[1]*rn[2] * nt[2]/nph[2];

	temp *= 1. / (2. * t3amp * t3amp);

	return temp;
}
