#include "Obs_OIFITS.h"

#include "fitsio.h"

/// Reads an OIFITS data file and creates a series of observation based upon the data.
/// Note, presently this function only supports ONE array, combiner, spectral mode per OIFITS file.
vector <Observation*> Obs_OIFITS::ReadObservation_OIFITS(Array * array, string filename)
{
    vector<Observation*> observations;
    
    // We don't attempt to read anything from the OIFITS file, just set the filename.
    observations.push_back(new Obs_OIFITS(array, filename) );
    
    return observations;

}

Obs_OIFITS::Obs_OIFITS(Array * array, string filename)
{
    this->mArray = array;
    this->mstrFilename = filename;
}

//
/// \todo This is a really hacky solution.  Come up with a better method.  Perhaps use
/// the Baseline::GetOI_Vis2_record routine?
oi_vis2 Obs_OIFITS::GetVis2(string ins_name, Source & source, vector<double> & wavenumbers)
{
    // Make a place to store the data as we are simulating it.
    vector<oi_vis2_record> vis2_data;

    // init some local vars:
    int nwave = int(wavenumbers.size());
    fitsfile * fptr;
    oi_fits data;
    int status = 0;
    
    // Open the input file as read-only data.
    fits_open_file(&fptr, this->mstrFilename.c_str(), READONLY, &status);
    if(status)
        throw std::runtime_error("Could not read OIFITS file.");
        
    // Now iterate through the OI_VIS2 tables, mirroring the sampling on the Source object 
    // at the specified wavenumbers;
    oi_vis2 vis2;
    UVPoint uv;
    oi_vis2_record input_record;
    Baseline * baseline;
    
    do
    {
        // Read in the next vis2 table
        read_next_oi_vis2(fptr, &vis2, &status);
        if(status)
            break;    
            
        for(int record_id = 0; record_id < vis2.numrec; record_id++)
        {
            // Use a local var to store some information
            input_record = vis2.record[record_id];
            oi_vis2_record output;
            
            baseline = mArray->GetBaseline(input_record.sta_index[0], input_record.sta_index[1]);
            
            // Copy some information over to the output record:
            output.target_id = 0;
            output.time = input_record.time;
            output.mjd = input_record.mjd;
            output.int_time = input_record.int_time;
            output.ucoord = input_record.ucoord;
            output.vcoord = input_record.vcoord;
            output.sta_index[0] = input_record.sta_index[0];
            output.sta_index[1] = input_record.sta_index[1];
            
            // Allocate memory for the vis2data, vis2error, and flag:
            output.vis2data = (double *) malloc(nwave * sizeof(double));
            output.vis2err = (double *) malloc(nwave * sizeof(double));
            output.flag = (char *) malloc(nwave * sizeof(char));
            
            // Now iterate over the wavenumbers
            for(int j = 0; j < nwave; j++)
            {            
                // Reset the UV coordinates
                uv.u = input_record.ucoord;
                uv.v = input_record.vcoord;
                uv.Scale(wavenumbers[j]);
                
                // Simulate the visibility based on the source.
                output.vis2data[j] = baseline->GetVis2(source, uv);
                // Copy the error from the input file.
                output.vis2err[j] = input_record.vis2err[j];
    			output.flag[j] = FALSE;
            }
            
            // Now append the data to the output vector
            vis2_data.push_back(output);
        }    

    }while(status == 0);
    
    // close the file
    fits_close_file(fptr, &status);
    oi_vis2_record vis2_record;
    
    // Now convert the vis2_data vector into a properly formatted OI_VIS2 table.
    oi_vis2 outvis2;
    int npow = int(vis2_data.size());
    string arrname = this->mArray->GetArrayName();
    
	outvis2.revision = 1;
	/// \bug The observation date is set to all zeros by default.  
	/// This is to ensure the user knows this is simulated data, but may not be compliant
	/// with the OIFITS format, or good "note taking"
	strncpy(outvis2.date_obs, "0000-00-00", 11);
	strncpy(outvis2.arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(outvis2.insname, ins_name.c_str(), FLEN_VALUE);
	outvis2.numrec = npow;
	outvis2.nwave = nwave;
	
	outvis2.record = (oi_vis2_record *) malloc(npow * sizeof(oi_vis2_record));	
	for(int i = 0; i < npow; i++)
	{
	    outvis2.record[i] = vis2_data.back();
	    vis2_data.pop_back();
	}
    
    
    return outvis2;
}

oi_t3   Obs_OIFITS::GetT3(string ins_name, Source & source, vector<double> & wavenumbers)
{
    oi_t3 t3;
    
    return t3;
}
