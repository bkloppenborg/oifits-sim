#include "Obs_OIFITS.h"

/// Reads an OIFITS data file and creates a series of observation based upon the data.
/// Note, presently this function only supports ONE array, combiner, spectral mode per OIFITS file.
vector <Observation*> Obs_OIFITS::ReadObservation_OIFITS(Array * array, string filename)
{
    // init some local variables
    vector<Observation*> observations;
    oi_fits data;
    int status = 0;
    
    // First read the file into memory
    read_oi_fits(filename.c_str(), &data, & status);
    if(status)
        throw std::runtime_error("Could not read OIFITS file.");
    
    // With any luck the array definition file with this program has IDs that match
    // what is listed in the OIFITS file so we simply read in UV coordinates and telescope IDs.
    
    return observations;
}


oi_vis2 Obs_OIFITS::GetVis2(string ins_name, Source & source, vector<double> & wavenumbers)
{
    // init local vars
    oi_vis2 vis2;
    
    return vis2;
}

oi_t3   Obs_OIFITS::GetT3(string ins_name, Source & source, vector<double> & wavenumbers)
{
    oi_t3 t3;
    
    return t3;
}
