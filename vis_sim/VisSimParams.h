#include <string>

using namespace std;

#ifndef VISSIMPARAMS_H
#define VISSIMPARAMS_H

/// \todo It would be nice to have this class print out the missing parameters with some information.

class VisSimParams
{
  public:
	string target_filename;
	string array_filename;
	string instrument_filename;
	string SpectralMode_filename;
	string observation_filename;
	string oifits_filename; // The output OIFITS filename
	
	string input_oifits_filename;
	
	bool bFromOIFITSFile;
	
    // Member Functions
  public:
    VisSimParams(void);
  
	bool have_all_params();
	
	string GetFilename(void);
	bool is_oifits_obs(void);
};

#endif //VISSIMPARAMS_H
