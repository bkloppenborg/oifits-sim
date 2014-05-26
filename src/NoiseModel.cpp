/*
 * NoiseModel.cpp
 *
 *  Created on: Apr 12, 2011
 *      Author: bkloppenborg
 */

#include "NoiseModel.h"


#include "Array.h"
#include "Target.h"
#include "Combiner.h"
#include "SpectralMode.h"
#include "Baseline.h"
#include "Triplet.h"
#include "Common.h"



NoiseModel::NoiseModel()
{
	// TODO Auto-generated constructor stub

}

NoiseModel::~NoiseModel()
{
	// TODO Auto-generated destructor stub
}

// Computes the Strehl ratio according to the wavelength (in meters), seeing
// at 500 nm, and the telescope diameter (in meters).
double NoiseModel::Strehl(double wavelength, double r0, double mirror_diameter)
{
	// First compute the seeing at the current wavelength.
	double r0lambda = r0 * pow( wavelength / 5e-7, 6.0 / 5.0);
	return exp(-0.134 * pow(mirror_diameter / r0lambda, 5.0 / 3.0));
}

// Compute the correct integration time according to the
// mean wavelength (in meters), seeing (in meters at 500 nm),
// and wind speed (m/s)
double NoiseModel::IntegrationTime(double median_wl, double r0, double wind_speed)
{
	// First compute the seeing at the current wavelength.
	double r0lambda = r0 * pow( median_wl / 5e-7 , 6.0 / 5.0);

	// formula taken from Bridget O'Donovan's Phd Thesis
	// see equation (1.23)
	return 0.314 * r0lambda / wind_speed;
}

// Computes the number of photons sent through the system per coherent integration period
// including throughput losses of the array only (combiner losses should be calculated elsewhere).
// Note: this function assumes all telescopes are of the same diameter and have the same gain
// also, this function doesn't split the photons between photometric and data channels, that should be done
// elsewhere.
double NoiseModel::NumPhotons(Array * array, Target * target, Combiner * combiner, SpectralMode * spec_mode, int wavelength_num)
{
	double median_wavelenth = spec_mode->median_wavelength;
	double wavelength = spec_mode->mean_wavelength[wavelength_num];
	double bandwidth = spec_mode->delta_wavelength[wavelength_num];

	// TODO: It would be good to write this data out in a memo, rather than
	// re-computing it each time.  These are all fast operations, at worst O(n_telescopes)
	// == linear.

	// TODO: Modify this function to permit different gains / sizes of telescopes.

	double integration = IntegrationTime(median_wavelenth, array->Get_r0(), array->GetWindSpeed());

	// TODO: right now we assume all stations are the same diameter. We will need to
	// change this if we permit unequal telescopes.
	double tele_diameter = array->GetStation(0)->diameter;
	double total_area = array->GetNumStations() * PI * pow(tele_diameter, 2.0);
	double strehl = Strehl(wavelength, array->Get_r0(), tele_diameter);

	// Essentially integrate the number of photons sent to the combiner:
	double n_photons = target->GetTargNPhotons(wavelength, bandwidth, total_area, integration);
	n_photons *= array->GetThroughput() * strehl;

	return n_photons;
}

