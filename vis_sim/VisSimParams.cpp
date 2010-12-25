
#include "VisSimParams.h"

VisSimParams::VisSimParams()
{
	this->target_filename = "";
	this->array_filename = "";
	this->instrument_filename = "";
	this->SpectralMode_filename = "";
	this->observation_filename = "";
	this->oifits_filename = "";
	
	this->input_oifits_filename = "";
	
	this->bFromOIFITSFile = false;
}

bool VisSimParams::is_oifits_obs(void)
{
    return this->bFromOIFITSFile;
}

/// Determines whether or not all parameters required to run the simulator are present
///
/// Two cases to see if all parameters are specified.  If we are using an OIFITS file, we only need:
/// target_filename, SpectralMode_filename, input_oifits_filename
/// everything else is determined from the OIFITS file.  Otherwise we need
/// target_filename, array_filename, instrument_filename, SpectralMode_filename, HourAngles_filename
bool VisSimParams::have_all_params()
{
    
    // Using an OIFITS file as input
    if(this->bFromOIFITSFile)
    {
        if((this->instrument_filename != "") && (this->target_filename != "") 
            && (this->SpectralMode_filename != "") && (this->input_oifits_filename != "")
            && (this->array_filename != ""))
        {
            return true;
        }
    }
    else    // Using text files as input.
    {
        if ((this->target_filename != "") && (this->array_filename != "") 
            && (this->instrument_filename != "") && (this->SpectralMode_filename != "") 
            && (this->observation_filename != ""))
        {
            return true;
        }
    }

    return false;
}

string VisSimParams::GetFilename()
{
    if(bFromOIFITSFile)
        return this->input_oifits_filename;
        
    return observation_filename;
}
