#include "Triplet.h"

#include "Array.h"
#include "Baseline.h"
#include "Station.h"
#include "Target.h"

#include <cstdio>
using namespace std;

Triplet::Triplet()
{
    // Just initialize one variable.
    this->name = "INVALID_TRIPLET";
}

Triplet::Triplet(Array * array, Station * station1, Station * station2, Station * station3)
{
    // Query the array for the baselines involved
    mStations[0] = station1;
    mStations[1] = station2;
    mStations[2] = station3;
    
    string baseline;
    baseline = station1->GetName() + '-' + station2->GetName();    // ab
    mBaselines[0] = array->GetBaseline(baseline);
    baseline = station2->GetName() + '-' + station3->GetName();    // bc
    mBaselines[1] = array->GetBaseline(baseline);
    baseline = station1->GetName() + '-' + station3->GetName();    // ac, note we do this backward to keep names ordered
    mBaselines[2] = array->GetBaseline(baseline);
    
    this->name = station1->GetName() + "-" + station2->GetName() + "-" + station3->GetName();
    
    // A useful comment if you wish to see all triplets created.
//    printf("Constructed %s from %s %s %s\n", this->name.c_str(), 
//        mBaselines[0]->GetName().c_str(),  mBaselines[1]->GetName().c_str(),  mBaselines[2]->GetName().c_str());
}

/// Returns the name of this triplet as the names of the stations separated by hyphens, e.g.:
///     S1-E1-E2
string Triplet::GetName(void)
{
    return this->name;
}

int     Triplet::GetStationID(int station_num)
{
    /// \todo Assert that station_num is less than 3.
    return this->mStations[station_num]->GetIndex();
}

//complex<double> Triplet::ComputeT3Error(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
//{
//    /// \bug T3pectrum error is set to 0.001 by default
//    return complex<double>(0.001, 0.001);
//}

complex<double> Triplet::ComputeT3(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
{
    // Get the visibilities on the baselines AB, BC, and CA.
    complex<double> AB = mBaselines[0]->GetVisibility(target, uv_ab);
    complex<double> BC = mBaselines[1]->GetVisibility(target, uv_bc);
    complex<double> AC = mBaselines[2]->GetVisibility(target, uv_ac);   // Note, this is AC NOT CA
    
    // Now compute the bispectrum.  We take the conjugate of AC to form CA.
    return AB * BC * conj(AC);
}

/// Returns a boolean to indicate if the baseline specified by bl_name is involved in this triplet.
bool    Triplet::ContainsBaseline(string bl_name)
{
    /// \todo We could probably implment this as a hash lookup instead of a string operation
    /// this should make the lookup slightly faster.
    for(int i = 0; i < 3; i++)
    {
        if(this->mBaselines[i]->GetName() == bl_name)
            return true;
    }

    return false;
}

/// Returns the norm (amplitude) of the bispectra
double  Triplet::GetT3Amp(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
{
    complex<double> bis = GetT3(target, uv_ab, uv_bc, uv_ac);
    return abs(bis);
}

// Wrapper function for the UV-coordinate taking versions.
double  Triplet::GetT3Amp(Target & target, double hour_angle, double wavenumber)
{
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);
    UVPoint uv_ac = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ac.Scale(wavenumber);

    return GetT3Amp(target, uv_ab, uv_bc, uv_ac);
}

// Returns the norm (amplitude) of the bispectrum
double  Triplet::GetT3AmpErr(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
{
    /// \bug Returns 0.001 by default
    return 0.001;
}
   

// Wrapper function for the UV-coordinate taking versions.
double  Triplet::GetT3AmpErr(Target & target, double hour_angle, double wavenumber)
{
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);
    UVPoint uv_ac = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ac.Scale(wavenumber);
    
    
     return GetT3AmpErr(target, uv_ab, uv_bc, uv_ac);
}


// Wrapper function for the UV-coordinate taking versions.
double  Triplet::GetT3Phi(Target & target, double hour_angle, double wavenumber)
{
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);
    UVPoint uv_ac = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ac.Scale(wavenumber);
    
    return GetT3Phi(target, uv_ab, uv_bc, uv_ac);
}

/// Returns the argument (phase) of the bispectra
double  Triplet::GetT3Phi(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
{
    complex<double> bis = GetT3(target, uv_ab, uv_bc, uv_ac);
    return arg(bis);
}


// Wrapper function for the UV-coordinate taking versions.
double  Triplet::GetT3PhiErr(Target & target, double hour_angle, double wavenumber)
{
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);
    UVPoint uv_ac = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ac.Scale(wavenumber);
    
    return GetT3PhiErr(target, uv_ab, uv_bc, uv_ac);
}

double  Triplet::GetT3PhiErr(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
{
    /// \bug Returns 0.001 by default.
    return 0.001;
}   
    
// Computes the bispectra from the three baselines in this triplet.
complex<double> Triplet::GetT3(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
{
    string hash_key = GetHashKey(target, uv_ab, uv_bc, uv_ac);
    complex <double> bis(0.0, 0.0);
    
    // First try looking up the value in the hash table
    if(mT3Values.find(hash_key) != mT3Values.end())
    {
        bis = mT3Values[hash_key];
    }
    else
    {
        // The value did not exist in the hash table, we need to compute and store it.
        bis = ComputeT3(target, uv_ab, uv_bc, uv_ac);
        mT3Values[hash_key] = bis;
    }
    
    return bis;
}

complex<double> Triplet::GetT3(Target & target, double hour_angle, double wavenumber)
{
	UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
	uv_ab.Scale(wavenumber);
	UVPoint uv_bc = mBaselines[1]->UVcoords(hour_angle, target.declination);
	uv_bc.Scale(wavenumber);
	UVPoint uv_ac = mBaselines[2]->UVcoords(hour_angle, target.declination);
	uv_ac.Scale(wavenumber);

	return GetT3(target, uv_ab, uv_bc, uv_ac);
}

// Computes a hash key from the source, hour angle, and wavenumber.
string  Triplet::GetHashKey(Target & target, UVPoint uv_ab, UVPoint uv_bc, UVPoint uv_ac)
{
    /// \todo It may be necessary for the doubles coming into this function to be cast into some 
    /// finite floating point format.
    
    /// \todo This function is in common with the Baseline class, need to factor this code.
    
    std::ostringstream sstream;
    sstream << target.GetName() << "-" << uv_ab.HashString() << '-' << uv_bc.HashString() << '-' << uv_ac.HashString();
    std::string str = sstream.str();
    return str;
}

Baseline * Triplet::GetBaseline(int baseline_num)
{
    /// \todo Assert that baseline_num is less than 3
    return mBaselines[baseline_num];
}

////////////////////////////////////////////////////////////////////
// Non Class Functions Below
////////////////////////////////////////////////////////////////////

/// Computes all possible triplets for the stations provided.
vector<Triplet> ComputeTriplets(Array * array, vector<Station> & stations)
{
    // init some local variables
    vector<Triplet> triplets;
    int num_stations = stations.size();

    // Compute all possible triplets.  This should total to N choose 3.
    for(int i = 0; i < num_stations - 2; i++)
    {
        for(int j = i + 1; j < num_stations - 1; j++)
        {
            for(int k = j + 1; k < num_stations; k++)
            {
                triplets.push_back( Triplet(array, &stations[i], &stations[j], &stations[k]) );
            }
        }
    }
    
    return triplets;
}

TripletHash ComputeTripletHash(vector<Triplet> & triplets)
{
    TripletHash hash;
    
    for(unsigned int i = 0; i < triplets.size(); i++)
    {
        hash.insert(TripletHash::value_type(triplets[i].GetName(), &triplets[i]) );
    }
    
    return hash;
}
