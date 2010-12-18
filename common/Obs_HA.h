#include "Observation.h"

#ifndef OBS_HA_H
#define OBS_HA_H

class Obs_HA : public Observation
{
  private:
    double  mHA;
    bool    mComputeHA;

  private:
    double  GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long);
    double  GetSiderealTime(double jd_high, double jd_low, double ee);    
    
  public:
    Obs_HA(Array * array, vector<Station*> stations, string exclude_baselines);
    Obs_HA(Array * array, double hour_angle, string telescopes, string exclude_baselines);
    Obs_HA(Array * array, double mjd, double time, string telescopes, string exclude_baselines);
  
    static vector <Observation*> ReadObservation_HA(Array * array, string filename, string comment_chars);
    static vector <Observation*> ReadObservation_Descriptive(Array * array, string filename, string comment_chars);
    
    oi_vis2 GetVis2(string ins_name, Source & source, vector<double> & wavenumbers);
    oi_t3   GetT3(string ins_name, Source & source, vector<double> & wavenumbers);
    
    double  GetHA(double targ_ra);
};

#endif // OBS_HA_H
