/*
 * NoiseModel_Tatulli2006.h
 *
 *  Created on: Apr 12, 2011
 *      Author: fbaron (via. MROI simulator, rewritten for class by bkloppenborg)
 *
 * This is a noise model implemented per Tatulli (2006) which consists of
 * three sources of noise:
 *  (1) Photon Noise
 *  (2) Detector Noise
 *  (3) Coupling Noise
 * these equations are suitable as first order approximations for the noise.
 *
 * Unless otherwise specified, all equations come from the aformentioned paper.
 *
 * Note: This differs from the MROI simulator in that it uses wavenumbers instead of
 * wavelengths.
 */

#ifndef NOISEMODEL_TATULLI2006_H_
#define NOISEMODEL_TATULLI2006_H_

#include "NoiseModel.h"

class NoiseModel_Tatulli2006 : public NoiseModel
{
public:
	NoiseModel_Tatulli2006();
	virtual ~NoiseModel_Tatulli2006();

public:
	double GetVis2Var(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, Baseline * baseline, UVPoint uv, int wavelength_num);
	double GetT3PhaseVar(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, Triplet * triplet, UVPoint uv1, UVPoint uv2, int wavelength_num);

private:
	double CoherentFlux(Array * array, Combiner * combiner, double n_photons, double vis2);
	double ComputeT3PhVar(int nt[7], double nph[7], double v_ij, double v_jk, double v_ik, double t3amp);
	double ComputeT3DetVar(int nt[7], double nph[7], double v_ij, double v_jk, double v_ik, double t3amp, double rn[7], int npix[7]);
	double VarCoherentFlux(Array * array, Combiner * combiner, double n_photons, double vis2);
	double VarPhotometricFlux(Combiner * combiner, double n_photons);
};

#endif /* NOISEMODEL_TATULLI2006_H_ */
