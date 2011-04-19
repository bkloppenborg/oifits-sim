/*
 * oifits-sim.h
 *
 *  Created on: Apr 8, 2011
 *      Author: bkloppenborg
 *
 *  Note:
 *  Everything is internally stored in radians, MKS units and must be converted
 *  during the import process.
 */

#include <string>
#include <vector>

class Target;
class Array;
class Combiner;
class SpectralMode;
class Observation;
class NoiseModel;

using namespace std;

void run_sim(Target * target, Array * array, Combiner * combiner, SpectralMode * spec, NoiseModel * noisemodel, vector<Observation*> observations, string output_filename);
