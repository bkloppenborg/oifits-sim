
#include "random.h"

// Header files for other libraries
extern "C" {
    #include "exchange.h"
    #include "oifile.h"
}

class VisSimParams;


// defining seed for the random number generating function
static Rand_t random_seed;

// Code to run the simulator.
void run_sim(const VisSimParams * p);

// Code for building OIFITS data tables.
oi_array GetOIArray(Array & array);
oi_target GetOITarget(Source & target);
oi_wavelength GetOIWavelength(SpectralMode spec);
