
#include "random.h"

class VisSimParams;

/*
 * WGS84 Earth equatorial radius /m 
 */
static const double WGS_A0 = 6378137;

/*
 * WGS84 reference spheroid flattening factor 
 */
static const double WGS_F = 1.0 / 298.257223563;

// defining seed for the random number generating function
static Rand_t random_seed;

// Code to run the simulator.
void run_sim(const VisSimParams * p);
