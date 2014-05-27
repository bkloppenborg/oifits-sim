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
  
    static vector <Observation*> MakeObservations(Array * array, double start, double stop, double every, string telescopes);
    
    static vector <Observation*> ReadObservation_HA(Array * array, vector < string > lines, int i);
    static vector <Observation*> ReadObservation_Descriptive(Array * array, vector < string > lines, int i);

    oi_vis2 GetVis2(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    oi_t3   GetT3(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    oi_t4   GetT4(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    
    double  GetHA(double targ_ra);
};

#endif // OBS_HA_H
