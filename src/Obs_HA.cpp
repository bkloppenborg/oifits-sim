#include "Obs_HA.h"
#include "Common.h"

#include "Array.h"
#include "Combiner.h"
#include "SpectralMode.h"
#include "Target.h"
#include "NoiseModel.h"

Obs_HA::Obs_HA(Array * array, vector<Station*> stations, string exclude_baselines)
{
    this->mArray = array;
    this->mStations = stations;
    this->mBaselines = this->FindBaselines(mStations, exclude_baselines);
}

// Construct an Observation object from the hour angle, and included/excluded telescopes.
Obs_HA::Obs_HA(Array * array, double hour_angle, string telescopes, string exclude_baselines)
{
    this->mArray = array;
    this->mHA = hour_angle;
    this->mComputeHA = false;

    // Now Make the baselines
    this->mStations = this->FindStations(telescopes);
    this->mBaselines = this->FindBaselines(mStations, exclude_baselines);
    this->mTriplets = this->FindTriplets(mStations, exclude_baselines);
    this->mQuadruplets = this->FindQuadruplets(mStations, exclude_baselines);

    if(mTriplets.size() > 0)
        this->mbHasTriplets = true;

    if(mQuadruplets.size() > 0)
        this->mbHasQuadruplets = true;
}

/// Construct an Observation object from the MJD, time, and included/excluded telescopes.
Obs_HA::Obs_HA(Array * array, double MJD, double time, string telescopes, string exclude_baselines)
{
    this->mArray = array;
    this->mHA = 0;
    // Compute the (full) Julian date from the Modified Julian Date (MJD)
    this->mJD = MJD + time / 24 + 2400000.5;
    this->mComputeHA = true;

    // Now Make the baselines
    this->mStations = this->FindStations(telescopes);
    this->mBaselines = this->FindBaselines(mStations, exclude_baselines);
    this->mTriplets = this->FindTriplets(mStations, exclude_baselines);
    this->mQuadruplets = this->FindQuadruplets(mStations, exclude_baselines);

    if(mTriplets.size() > 0)
        this->mbHasTriplets = true;

    if(mQuadruplets.size() > 0)
        this->mbHasQuadruplets = true;
}


vector <Observation*> Obs_HA::MakeObservations(Array * array, double start, double stop, double every, string telescopes)
{
    vector <Observation*> observations;
	double ha;
	double temp;

	// Enforce start < stop:
	if(start > stop)
	{
		temp = start;
		start = stop;
		stop = temp;
	}

	int intervals = int((stop - start) / every);
	cout << "There will be " << intervals << " observations." << endl;

	for(int i = 0; i < intervals; i++)
	{
		ha = start + i * every;
		observations.push_back(new Obs_HA(array, ha, array->GetAllStationNames(), "") );
	}

	return observations;
}

/// Reads in a file that consists of lines of hour angles with or without comments.
vector <Observation*> Obs_HA::ReadObservation_HA(Array * array, vector < string > lines, int i)
{
    vector <Observation*> observations;
	vector <string> results;

    // Now parse the file, make the observations
    double ha;
    for(unsigned int j = i+1; j < lines.size(); j++)
    {
        results.clear();
        results = SplitString(lines[i], '=');
        StripWhitespace(results);

        if(results[0] == "hour_angle")
        {
        	try
        	{

                ha = atof(results[1].c_str());

                // Make a new observation with all of the stations included.
                observations.push_back(new Obs_HA(array, ha, array->GetAllStationNames(), "") );
        	}
        	catch(...)
        	{
        		throw std::runtime_error("Invalid observation type field in observation file.");
        	}
        }
        else
        {
        	printf("Warning, detected non-keyword, %s, in observation definition file.  Ignoring\n", results[0].c_str());
        }
    }

    return observations;
}

/// Reads in a series of formatted lines that consist of keywords:
///     hour_angle = DECIMAL NUMBER
///     stations = STATION_NAMES
///     exclude = EXCLUDED_BASELINES
/// DECIMAL NUMBER is a decimal number denoting the hour angle of the observation.
/// STATION_NAMES is a CSV list of stations included in this observation in the following format:
///     S1, S2, ... , E1
/// white space is automatically stripped between commas.
/// EXCLUDED_BASELINES is a CSV list of excluded baselines consisting of two station names separated
/// by a hyphen, e.g.:
///     S1-S2, S2-E2
/// with the lower station index (as defined in the array definition file) always appearing first.
/// white space is automatically stripped between commas.  EXCLUDED_BASELINES may be blank or omitted.
vector <Observation*> Obs_HA::ReadObservation_Descriptive(Array * array, vector < string > lines, int i)
{
    // init local vars
    vector <Observation*> observations;
    vector<string> results;

    // local vars for an observation
    double hour_angle = 0;
    string stations = "";
    string excluded_baselines = "";

    // Flags to keep track of whether or not we have the corresponding data
    bool bStations = false;
    bool bHourAngle = false;
    bool bExclude = false;

    unsigned int j = 1;  // A counter for the number of observations.  Only used in error messages below.
    // Iterate over the lines in the file
    for(unsigned int k = i+1; k < lines.size(); k++)
    {
        // First split the line and strip out white space
        results.clear();
        results = SplitString(lines[k], '=');
        StripWhitespace(results);

        if(results[0] == "hour_angle")
        {
            // We allow for the "exclude" parameter to be omitted.  If we have already found
            // a "hour_angle" flag and a "stations" flag we need to write out an observation block
            if(bHourAngle && bStations)
            {
                // Set the excluded stations to be blank, mark it as found so we will write out an observation
                excluded_baselines = "";
                bExclude = true;

                // Step back one in "k" so at we may revisit this line
                k--;
            }
            else if(bHourAngle && !(bStations))
            {
                /// \todo Perhaps issue a warning, but keep going?
                // We've found another an additional hour angle, but no stations.
                // The file is not formatted correctly
                throw std::runtime_error("Found an 'hour_angle' specifier without stations in observation block " + j);
            }
            else
            {
                try
                {
                    hour_angle = atof(results[1].c_str());
                    bHourAngle = true;
                }
                catch (const std::exception&)
                {
                    throw std::runtime_error("Hour angle specification is not numeric in observation block " + j);
                }
            }
        }
        else if(results[0] == "telescopes")
        {

            if(results[1].size() < 2)
                throw std::runtime_error("Two few telescopes specified in observation entry " + j);

            stations = results[1];
            bStations = true;
        }
        else if(results[0] == "exclude")
        {
            excluded_baselines = results[1];
            bExclude = true;
        }
        else
        {
            // This is an unknown type, thrown an exception.
            throw std::runtime_error("Unknown parameter specified in observation block " + j);
        }

        if(bHourAngle && bStations && (bExclude || k == lines.size() - 1))
        {
            // Push the observation on to the list
            observations.push_back(new Obs_HA(array, hour_angle, stations, excluded_baselines) );

            // Increment the observation block counter, j
            j++;

            // Reset the station, hour angle, and exclude flags.
            bStations = false;
            bHourAngle = false;
            bExclude = false;
        }
    }

    return observations;
}



/// Creates an OIFITS oi_vis2 compliant entry for this observation.
oi_vis2 Obs_HA::GetVis2(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
    // init local vars
    oi_vis2 vis2;
    int npow = this->mBaselines.size();
    int nwave = spec_mode->mean_wavenumber.size();
    int iwave = 0;
    string arrname = array->GetArrayName();
    string ins_name = spec_mode->spec_mode;
    double wavenumber = 0;

    double ra = target->right_ascension;
    double dec = target->declination;
    UVPoint uv;
    double v2;
    double v2_err;

    // Allocate room for the vis2 records and their data.
	vis2.record = (oi_vis2_record *) malloc(npow * sizeof(oi_vis2_record));
	for (int i = 0; i < npow; i++)
	{
		vis2.record[i].vis2data = (double *) malloc(nwave * sizeof(double));
		vis2.record[i].vis2err = (double *) malloc(nwave * sizeof(double));
		vis2.record[i].flag = (char *) malloc(nwave * sizeof(char));
	}

	vis2.revision = 1;
	/// \bug The observation date is set to all zeros by default.
	/// This is to ensure the user knows this is simulated data, but may not be compliant
	/// with the OIFITS format, or good "note taking"
	strncpy(vis2.date_obs, "2014-01-01", 11);
	strncpy(vis2.arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(vis2.insname, ins_name.c_str(), FLEN_VALUE);
	vis2.numrec = npow;
	vis2.nwave = nwave;
	for (int i = 0; i < npow; i++)
	{
		vis2.record[i].target_id = target->GetTargetID();
		/// \bug The time is set to zero by default
		vis2.record[i].time = 0.0;
		vis2.record[i].mjd = this->mJD;
		/// \bug The integration time is set to 10 seconds by default.
		vis2.record[i].int_time = 10;

		// Compute the UV coordinates and record the station positions:
		uv = this->mBaselines[i]->UVcoords(this->GetHA(ra), dec);
		vis2.record[i].ucoord = uv.u;
		vis2.record[i].vcoord = uv.v;
		vis2.record[i].sta_index[0] = this->mBaselines[i]->GetStationID(0);
		vis2.record[i].sta_index[1] = this->mBaselines[i]->GetStationID(1);

		// Now compute the individual visibilities and uncertainties
		for (iwave = 0; iwave < nwave; iwave++)
		{
		    // look up the present wavenumber, and then find the data
		    wavenumber = spec_mode->mean_wavenumber[iwave];

		    // Get the squared visibility, and it's error.
		    v2 = this->mBaselines[i]->GetVis2(*target, mHA, wavenumber);
		    v2_err = noisemodel->GetVis2Var(array, combiner, spec_mode, target, this->mBaselines[i], uv, iwave);

		    // Put out the visibility
			vis2.record[i].vis2data[iwave] = v2 + v2_err * Rangauss(random_seed);
			vis2.record[i].vis2err[iwave] = v2_err;
			vis2.record[i].flag[iwave] = FALSE;
		}
	}

	return vis2;
}

/// Create an OIFITS-compliant t3 table for this observation.
oi_t3   Obs_HA::GetT3(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
    oi_t3 t3;
    int nTriplets = this->mTriplets.size();
    int nwave = spec_mode->mean_wavenumber.size();
    string arrname = array->GetArrayName();
    string ins_name = spec_mode->spec_mode;
    complex<double> bis;
    double phi_err;
    double wavenumber;

    UVPoint uv_AB;
    UVPoint uv_BC;

	t3.record = (oi_t3_record *) malloc(nTriplets * sizeof(oi_t3_record));
	for (int i = 0; i < nTriplets; i++)
	{
		t3.record[i].t3amp = (double *) malloc(nwave * sizeof(double));
		t3.record[i].t3amperr = (double *) malloc(nwave * sizeof(double));
		t3.record[i].t3phi = (double *) malloc(nwave * sizeof(double));
		t3.record[i].t3phierr = (double *) malloc(nwave * sizeof(double));
		t3.record[i].flag = (char *) malloc(nwave * sizeof(char));
	}
	t3.revision = 1;
	strncpy(t3.date_obs, "2014-01-01", FLEN_VALUE);
	strncpy(t3.arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(t3.insname, ins_name.c_str(), FLEN_VALUE);
	t3.numrec = nTriplets;
	t3.nwave = nwave;

	// Now copy the data into t3 records:
	for (int i = 0; i < nTriplets; i++)
	{

		t3.record[i].target_id = target->GetTargetID();
		/// \bug Time is set to zero by default.
		t3.record[i].time = 0.0;
		t3.record[i].mjd = this->mJD;
		/// \bug Integration time set to 10 seconds by default.
		t3.record[i].int_time = 10;

		// Get the UV coordinates for the AB and BC baselines
		uv_AB = mTriplets[i]->GetBaseline(0)->UVcoords(this->mHA, target->declination);
		uv_BC = mTriplets[i]->GetBaseline(1)->UVcoords(this->mHA, target->declination);

		t3.record[i].u1coord = uv_AB.u;
		t3.record[i].v1coord = uv_AB.v;
		t3.record[i].u2coord = uv_BC.u;
		t3.record[i].v2coord = uv_BC.v;
		t3.record[i].sta_index[0] = mTriplets[i]->GetStationID(0);
		t3.record[i].sta_index[1] = mTriplets[i]->GetStationID(1);
		t3.record[i].sta_index[2] = mTriplets[i]->GetStationID(2);

		for(int j = 0; j < nwave; j++)
		{
			wavenumber = spec_mode->mean_wavenumber[j];
		    bis = mTriplets[i]->GetT3(*target, this->mHA, wavenumber);
		    phi_err = noisemodel->GetT3PhaseVar(array, combiner, spec_mode, target, mTriplets[i], uv_AB, uv_BC, j);

			// assume circular noise cloud

		    // First save the amplitudes
			t3.record[i].t3amp[j] = abs(bis);
			t3.record[i].t3amperr[j] = sqrt(abs(bis) * abs(bis) * phi_err * phi_err);
			// Now save the phases.  Remember, the phase is in degrees rather than radians.
			t3.record[i].t3phi[j] = (arg(bis) + phi_err * Rangauss(random_seed)) * 180 / PI;
			t3.record[i].t3phierr[j] = phi_err * 180 / PI;
			t3.record[i].flag[j] = FALSE;
		}

	}

	return t3;
}

/// Create a t4 table for this observation.
oi_t4   Obs_HA::GetT4(Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
    oi_t4 t4;
    int nQuadruplets = this->mQuadruplets.size();
    int nwave = spec_mode->mean_wavenumber.size();
    string arrname = array->GetArrayName();
    string ins_name = spec_mode->spec_mode;
    complex<double> quad_clos;
    double phi_err;
    double wavenumber;

    UVPoint uv_AB;
    UVPoint uv_CD;
    UVPoint uv_AD;

	t4.record = (oi_t4_record *) malloc(nQuadruplets * sizeof(oi_t4_record));
	for (int i = 0; i < nQuadruplets; i++)
	{
		t4.record[i].t4amp = (double *) malloc(nwave * sizeof(double));
		t4.record[i].t4amperr = (double *) malloc(nwave * sizeof(double));
		t4.record[i].t4phi = (double *) malloc(nwave * sizeof(double));
		t4.record[i].t4phierr = (double *) malloc(nwave * sizeof(double));
		t4.record[i].flag = (char *) malloc(nwave * sizeof(char));
	}
	t4.revision = 1;
	/// \bug Observation date is set to 0000-00-00 by default
	strncpy(t4.date_obs, "2014-01-01", FLEN_VALUE);
	strncpy(t4.arrname, arrname.c_str(), FLEN_VALUE);
	strncpy(t4.insname, ins_name.c_str(), FLEN_VALUE);
	t4.numrec = nQuadruplets;
	t4.nwave = nwave;

	// Now copy the data into t4 records:
	for (int i = 0; i < nQuadruplets; i++)
	{
		t4.record[i].target_id = target->GetTargetID();
		/// \bug Time is set to zero by default.
		t4.record[i].time = 0.0;
		t4.record[i].mjd = this->mJD;
		/// \bug Integration time set to 10 seconds by default.
		t4.record[i].int_time = 10;

		// Get the UV coordinates for the AB and BC baselines
		uv_AB = mQuadruplets[i]->GetBaseline(0)->UVcoords(this->mHA, target->declination);
		uv_CD = mQuadruplets[i]->GetBaseline(1)->UVcoords(this->mHA, target->declination);
		uv_AD = mQuadruplets[i]->GetBaseline(2)->UVcoords(this->mHA, target->declination);

		t4.record[i].u1coord = uv_AB.u;
		t4.record[i].v1coord = uv_AB.v;
		t4.record[i].u2coord = uv_CD.u;
		t4.record[i].v2coord = uv_CD.v;
		t4.record[i].u3coord = uv_AD.u;
		t4.record[i].v3coord = uv_AD.v;

		t4.record[i].sta_index[0] = mQuadruplets[i]->GetStationID(0);
		t4.record[i].sta_index[1] = mQuadruplets[i]->GetStationID(1);
		t4.record[i].sta_index[2] = mQuadruplets[i]->GetStationID(2);
		t4.record[i].sta_index[3] = mQuadruplets[i]->GetStationID(3);

		for(int j = 0; j < nwave; j++)
		{
		  wavenumber = spec_mode->mean_wavenumber[j];
		  quad_clos = mQuadruplets[i]->GetT4(*target, this->mHA, wavenumber);
		  phi_err = noisemodel->GetT4PhaseVar(array, combiner, spec_mode, target, mQuadruplets[i], uv_AB, uv_CD, uv_AD, j);

		  // assume circular noise cloud

		  // First save the amplitudes
		  t4.record[i].t4amp[j] = abs(quad_clos);
		  t4.record[i].t4amperr[j] = sqrt(abs(quad_clos) * abs(quad_clos) * phi_err * phi_err);
		  // Now save the phases.  Remember, the phase is in degrees rather than radians.
		  t4.record[i].t4phi[j] = (arg(quad_clos) + phi_err * Rangauss(random_seed)) * 180 / PI;
		  t4.record[i].t4phierr[j] = phi_err * 180 / PI;
		  t4.record[i].flag[j] = FALSE;
		}

	}

	return t4;
}

double Obs_HA::GetHA(double targ_ra)
{
    // If the hour angle is already specified, just reutrn it directly.
    if(!this->mComputeHA)
        return this->mHA;

    double lst = this->GetLocalSiderealTime(this->mJD, 0, 0, this->mArray->GetLongitude());
    return lst - targ_ra * 12.0/PI;
}

double Obs_HA::GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long)
{
    double lst = this->GetSiderealTime(jd_high, jd_low, ee) + array_long * 12 / PI;

	if(lst < 0)
        lst += 24;

    return lst;
}

// Get the Sidereal time as a double given the High and Low Julian Date
double  Obs_HA::GetSiderealTime(double jd_high, double jd_low, double ee)
{
	// Code slightly modified from Naval Observatory Vector Astrometry Subroutines
	// (C Language Version 2.0)
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
