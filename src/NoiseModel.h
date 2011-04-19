/*
 * NoiseModel.h
 *
 *  Created on: Apr 12, 2011
 *      Author: bkloppenborg
 *
 * This is the base class for various noise models.
 */

#ifndef NOISEMODEL_H_
#define NOISEMODEL_H_

#include "Array.h"
#include "Baseline.h"
#include "Combiner.h"
#include "Target.h"
#include "Triplet.h"
#include "SpectralMode.h"
#include "UVPoint.h"

// A list of headers from other noise models.
//#include "NoiseModel_Tatulli2006.h"

class NoiseModel
{
public:
	NoiseModel();
	virtual ~NoiseModel();

	double Strehl(double wavelength, double r0, double mirror_diameter);
	double IntegrationTime(double med_wavenumber, double r0, double wind_speed);
	double NumPhotons(Array * array, Target * target, Combiner * combiner, SpectralMode * spec_mode, int wavelength_num);

	virtual double GetVis2Var(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, Baseline * baseline, UVPoint uv, int wavelength_num) = 0;
	virtual double GetT3PhaseVar(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, Triplet * triplet, UVPoint uv1, UVPoint uv2, int wavelength_num) = 0;
};

#endif /* NOISEMODEL_H_ */
