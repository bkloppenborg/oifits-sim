#include "Observation.h"

#ifndef OBS_OIFITS_H
#define OBS_OIFITS_H


class Obs_OIFITS : public Observation 
{
  private:

    
  public:
    static vector <Observation*> ReadObservation_OIFITS(Array * array, string filename);

    oi_vis2 GetVis2(string ins_name, Source & source, vector<double> & wavenumbers);
    oi_t3   GetT3(string ins_name, Source & source, vector<double> & wavenumbers);
};

#endif // OBS_OIFITS_H
