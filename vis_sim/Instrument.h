// Class defining the instrument
// Stores parameters needed for the VSI noise model

#ifndef INSTRUMENT_H
#define INSTRUMENT_H

#include <string>
using namespace std;

class Instrument
{
  public:
    // AS 2010-06-22
    // adapted constructor to work with new type of instrument text file 
    // that contains more paramters
    Instrument(double throughput, double visibility,
        int Npix, double read_noise, double quantum_efficiency,
        double vsqFracCalErr, double cloPhaseCalErrDeg,
        double wind_speed, double ast_seeing, double incoh_time);
    
    Instrument(string filename, string comment_chars);
    
    double GlobalThroughput();

    // quantities needed for noise model
    double throughput;			//< global throughput, excluding detector QE
    double visibility;			//< instrumental visibility
    int Npix;					//< number of pixels to sample fringes
    double read_noise;			//< detector readout noise in units of e
    double quantum_efficiency;	//< detector QE
    double vsq_frac_cal_err;	//< fractional calibration error on squared visibility
    double clo_cal_err_deg;		//< calibration error on closure phase /degrees
    // AS 2010-06-22
    // added three new paramters to the instrument class
    double wind_speed;			//< seeing wind speed v in m/s
    double ast_seeing;			//< seeing Fried Parameter r_0 at 500nm wavelength in m
    double incoh_time;			//< incoherent integration time in s
    // Note that coherent integration time is 2t_0(lambda_obs),
    // where t_0 = 0.314 r_0/v
};

#endif							// #ifndef INSTRUMENT_H
