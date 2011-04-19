// Class defining the instrument
// Stores parameters needed for the VSI noise model

#ifndef COMBINER_H
#define COMBINER_H

#include <string>
using namespace std;

class Combiner
{
public:
	string name;
	double int_trans;	// Internal transmission
	double visibility;
	double read_noise;
	double quantum_efficiency;
	double v2_frac_cal_err;
	double phase_cal_err;
	double incoh_int_time;

	int n_pix_fringe;
	int n_pix_photometry;

	double flux_frac_photometry;
	double flux_frac_fringes;

	double throughput_photometry;
	double throughput_fringes;

	int n_splits;

public:
	Combiner();
	~Combiner();

	void ImportFile(string filename, string comment_chars);

	double	GetThroughput(void);
	string GetName(void);

	int 	GetNumPhotometricPixels();
	int		GetNumFringePixels();

	double	GetPhotometryFluxFrac();
	double	GetFringeFluxFrac();

	double 	GetPhotometryThroughput();
	double 	GetFringeThroughput();

	double 	GetNumSplits();
};

#endif							// #ifndef COMBINER_H
