
#include "random.h"

// Header files for other libraries
extern "C" {
    #include "exchange.h"
    #include "oifile.h"
}

class VisSimParams;
class Array;
class Source;
class SpectralMode;


// defining seed for the random number generating function
static Rand_t random_seed;

// Code to run the simulator.
void run_sim(const VisSimParams * p);
