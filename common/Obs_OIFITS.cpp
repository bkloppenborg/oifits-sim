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
    
    /// \bug By default the OIFITS file is designated to have triplets.  Need to check this.
    this->mbHasTriplets = true;
    
    // Now, check the OIFITS file to see if it has T3 data.
}

//
/// \todo This is a really hacky solution.  Come up with a better method.  Perhaps use
/// the Baseline::GetOI_Vis2_record routine?
oi_vis2 Obs_OIFITS::GetVis2(string ins_name, Source & source, vector<double> & wavenumbers)
{
    // Make a place to store the data as we are simulating it.
    // Note, the vis2_data vector doesn't know how to allocate oi_vis2_record* entries so we have 
    // to do this manually below.
    vector<oi_vis2_record*> vis2_data;

    // init some local vars:
    int nwave = int(wavenumbers.size());
    fitsfile * fptr;
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
            oi_vis2_record * output = (oi_vis2_record*) malloc(sizeof(oi_vis2_record));
            
            baseline = mArray->GetBaseline(input_record.sta_index[0], input_record.sta_index[1]);
            
            if(baseline == NULL)
                printf("Indicies: %i %i \n", input_record.sta_index[0], input_record.sta_index[1]);
            
            // Copy some information over to the output record:
            /// \bug Target is set to 1 by default
            output->target_id = 1;
            output->time = input_record.time;
            output->mjd = input_record.mjd;
            output->int_time = input_record.int_time;
            output->ucoord = input_record.ucoord;
            output->vcoord = input_record.vcoord;
            output->sta_index[0] = input_record.sta_index[0];
            output->sta_index[1] = input_record.sta_index[1];
            
            // Allocate memory for the vis2data, vis2error, and flag:
            output->vis2data = (double *) malloc(nwave * sizeof(double));
            output->vis2err = (double *) malloc(nwave * sizeof(double));
            output->flag = (char *) malloc(nwave * sizeof(char));
            
            // Now iterate over the wavenumbers
            for(int j = 0; j < nwave; j++)
            {            
                // Reset the UV coordinates
                uv.u = input_record.ucoord;
                uv.v = input_record.vcoord;
                uv.Scale(wavenumbers[j]);
                
                // Simulate the visibility based on the source.
                output->vis2data[j] = baseline->GetVis2(source, uv);
                // Copy the error from the input file.
                output->vis2err[j] = input_record.vis2err[j];
    			output->flag[j] = FALSE;
            }
            
            // Now append the data to the output vector
            vis2_data.push_back(output);
        }    

    }while(status == 0);
    
    // close the file
    fits_close_file(fptr, &status);
    
    // Now convert the vis2_data vector into a properly formatted OI_VIS2 table.
    oi_vis2 * outvis2 = (oi_vis2*) malloc(sizeof(oi_vis2));
    int npow = int(vis2_data.size());
    string arrname = this->mArray->GetArrayName();
    
	outvis2->revision = 1;
	/// \bug The observation date is set to all zeros by default.  
	/// This is to ensure the user knows this is simulated data, but may not be compliant
	/// with the OIFITS format, or good "note taking"
	strncpy(outvis2->date_obs, "0000-00-00", 11);
	strncpy(outvis2->arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(outvis2->insname, ins_name.c_str(), FLEN_VALUE);
	outvis2->numrec = npow;
	outvis2->nwave = nwave;
	
	outvis2->record = (oi_vis2_record *) malloc(npow * sizeof(oi_vis2_record));	
	for(int i = 0; i < npow; i++)
	{
	    outvis2->record[i] = *vis2_data.back();
	    
	    // free memory and pop
	    delete vis2_data.back();
	    vis2_data.pop_back();
	}
    
    return *outvis2;
}

/// \todo This is a really hacky solution.  Try to find a better solution.
oi_t3   Obs_OIFITS::GetT3(string ins_name, Source & source, vector<double> & wavenumbers)
{
    // a place to store the data as we are simulating it.
    // Notice, the t3_data vector doesn't know how to deallocate the oi_t3_records
    // so we will need to manually deallocate memory below.
    vector<oi_t3_record*> t3_data;

    // init some local vars:
    int nwave = int(wavenumbers.size());
    fitsfile * fptr;
    int status = 0;
    
    // Open the input file as read-only data.
    fits_open_file(&fptr, this->mstrFilename.c_str(), READONLY, &status);
    if(status)
        throw std::runtime_error("Could not read OIFITS file.");
        
    // Now iterate through the OI_VIS2 tables, mirroring the sampling on the Source object 
    // at the specified wavenumbers;
    oi_t3 t3;
    UVPoint uv1;
    UVPoint uv2;
    UVPoint uv3;
    oi_t3_record input_record;
    Triplet * triplet;
    
    do
    {
        // Read in the next vis2 table
        read_next_oi_t3(fptr, &t3, &status);
        if(status)
            break;    
            
        for(int record_id = 0; record_id < t3.numrec; record_id++)
        {
            // Use a local var to store some information
            input_record = t3.record[record_id];
            oi_t3_record * output = (oi_t3_record*) malloc(sizeof(oi_t3_record));
            
            // Get the triplet
            triplet = mArray->GetTriplet(input_record.sta_index[0], input_record.sta_index[1], input_record.sta_index[2]);
            
            // Copy some information over to the output record:
            /// \bug Target is set to 1 by default
            output->target_id = 1;
            output->time = input_record.time;
            output->mjd = input_record.mjd;
            output->int_time = input_record.int_time;
            output->u1coord = input_record.u1coord;
            output->v1coord = input_record.v1coord;
            output->u2coord = input_record.u2coord;
            output->v2coord = input_record.v2coord;
            output->sta_index[0] = input_record.sta_index[0];
            output->sta_index[1] = input_record.sta_index[1];
            output->sta_index[2] = input_record.sta_index[2];
            
            // Allocate memory for the vis2data, vis2error, and flag:
            output->t3amp = (double *) malloc(nwave * sizeof(double));
            output->t3amperr = (double *) malloc(nwave * sizeof(double));
            output->t3phi = (double *) malloc(nwave * sizeof(double));
            output->t3phierr = (double *) malloc(nwave * sizeof(double));
            output->flag = (char *) malloc(nwave * sizeof(char));
                        	
            // Now iterate over the wavenumbers
            for(int j = 0; j < nwave; j++)
            {            
                // Reset the UV coordinates
                uv1.u = input_record.u1coord;
                uv1.v = input_record.v1coord;
                uv2.u = input_record.u2coord;
                uv2.v = input_record.v2coord;
                uv3.u = uv1.u + uv2.u;
                uv3.v = uv1.v + uv2.v;
                
                // Scale them
                uv1.Scale(wavenumbers[j]);
                uv2.Scale(wavenumbers[j]);
                uv3.Scale(wavenumbers[j]);

                // Simulate the bispectrum's amplitude based on the source image
                output->t3amp[j] = triplet->GetBisAmp(source, uv1, uv2, uv3);
                // Copy the error from the input file.
                output->t3amperr[j] = input_record.t3amperr[j];
                
                // Simulate the bispectrum's phase based on the source image.  Convert to Degrees!
                output->t3phi[j] = triplet->GetBisPhi(source, uv1, uv2, uv3) * 180 / PI;
                // Copy the error from the input file
                output->t3phierr[j] = input_record.t3phierr[j];
            
                // Lastly set the "ignore this data" flag to false
    			output->flag[j] = FALSE;
            }
            
            // Now append the data to the output vector
            t3_data.push_back(output);

        }    

    }while(status == 0);
    
    // close the file
    fits_close_file(fptr, &status);
    
    // Now convert the t3_data vector into a properly formatted OI_T3 table.
    oi_t3 * outt3 = (oi_t3*) malloc(sizeof(oi_t3));
    int ndata = int(t3_data.size());
    string arrname = this->mArray->GetArrayName();
    
	outt3->revision = 1;
	/// \bug The observation date is set to all zeros by default.  
	/// This is to ensure the user knows this is simulated data, but may not be compliant
	/// with the OIFITS format, or good "note taking"
	strncpy(outt3->date_obs, "0000-00-00", 11);
	strncpy(outt3->arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(outt3->insname, ins_name.c_str(), FLEN_VALUE);
	outt3->numrec = ndata;
	outt3->nwave = nwave;
	
	outt3->record = (oi_t3_record *) malloc(ndata * sizeof(oi_t3_record));		
	for(int i = 0; i < ndata; i++)
	{
	    outt3->record[i] = *t3_data.back();
	    
	    // Free memory and pop
	    //delete t3_data.back();
	    t3_data.pop_back();
	}
    
    // Note, this memory object contains pointers and should be freed the OIFITSLIB free_oi_t3.
    return *outt3;
}
