
/// \file Observation.cpp

#include <cmath>

#include "Observation.h"
#include <cstdio>

/// The minimum timespan
const double DAYINSECONDS = 24 * 60 * 60;
const double SAMEOBS = 5 * 60 / DAYINSECONDS; // 5 minutes in seconds.

bool SameObservation(Observation & A, Observation & B)
{
    double aJD = A.GetJD();
    double bJD = B.GetJD();
    
    //printf("time difference: %f", fabs(aJD - bJD));
    
    if(fabs(aJD - bJD) < SAMEOBS)
        return true;
        
    return false;
}

Observation::Observation()
{
    this->mMJD = 0;
    this->mTime = 0;
}

Observation::Observation(double mjd, double time)
{
    this->mMJD = mjd;
    this->mTime = time;
}

Observation::Observation(double mjd)
{
    this->mMJD = mjd;
    this->mTime = 0;
}

double Observation::GetJD()
{
    double jd = this->mMJD + 2400000.5; // + this->mTime / DAYINSECONDS;
    //printf("%f %f %f \n", this->mMJD, this->mTime/DAYINSECONDS, jd);
    return jd;
}

double Observation::GetHA(double targ_ra, double targ_dec, double array_lat, double array_long)
{
    double jd = this->GetJD();
    double lst = this->GetLocalSiderealTime(jd, 0, 0, array_long);
    
    return lst - targ_ra * 24.0/360;
}

double Observation::GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long)
{
    double lst = this->GetSiderealTime(jd_high, jd_low, ee) + array_long * 24 / 360;
    
	if(lst < 0)
        lst += 24;
        
    return lst;
}

// Get the Sidereal time as a double given the High and Low Julian Date
double  Observation::GetSiderealTime(double jd_high, double jd_low, double ee)
{
	// Code slightly modified from NOVAS to use the current SystemTime
	// Naval Observatory Vector Astrometry Subroutines (C Language Version 2.0)
	const double T0 = 2451545.00000000;
	double t_hi = 0;
	double t_lo = 0;
	double t = 0;
	double t2 = 0;
	double t3 = 0;
	double st = 0;
	double gst = 0;

	t_hi = (jd_high -  T0) / 36525.0;
	t_lo = jd_low / 36525.0;
	t = t_hi + t_lo;
	t2 = t * t;
	t3 = t2 * t;

	st =  ee - 6.2e-6 * t3 + 0.093104 * t2 + 67310.54841 + 8640184.812866 * t_lo  + 3155760000.0 * t_lo + 8640184.812866 * t_hi + 3155760000.0 * t_hi;

	gst = fmod ((st / 3600.0), 24.0);

	if (gst < 0.0)
		gst += 24.0;

	return gst;
}
