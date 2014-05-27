#include "Quadruplet.h"

#include "Array.h"
#include "Baseline.h"
#include "Station.h"
#include "Target.h"

#include <cstdio>
using namespace std;

Quadruplet::Quadruplet()
{
    // Just initialize one variable.
    this->name = "INVALID_QUADRUPLET";
}

Quadruplet::Quadruplet(Array * array, Station * station1, Station * station2, Station * station3, Station * station4)
{
    // Query the array for the baselines involved
    mStations[0] = station1;
    mStations[1] = station2;
    mStations[2] = station3;
    mStations[3] = station4;
    
    string baseline;
    baseline = station1->GetName() + '-' + station2->GetName();    // ab
    mBaselines[0] = array->GetBaseline(baseline);
    baseline = station2->GetName() + '-' + station3->GetName();    // cd
    mBaselines[1] = array->GetBaseline(baseline);
    baseline = station1->GetName() + '-' + station4->GetName();    // ad
    mBaselines[2] = array->GetBaseline(baseline);
    baseline = station2->GetName() + '-' + station3->GetName();    // bc instead of cb, to keep alphabetical order
    mBaselines[3] = array->GetBaseline(baseline);
    
    this->name = station1->GetName() + "-" + station2->GetName() + "-" + station3->GetName() + "-" + station4->GetName();
    
    // A useful comment if you wish to see all quadruplets created.
    printf("Constructed %s from %s %s %s %s\n", this->name.c_str(), mBaselines[0]->GetName().c_str(),  mBaselines[1]->GetName().c_str(),  mBaselines[2]->GetName().c_str(),  mBaselines[3]->GetName().c_str());
}

/// Returns the name of this quadruplet as the names of the stations separated by hyphens, e.g.:
///     S1-E1-E2
string Quadruplet::GetName(void)
{
    return this->name;
}

int     Quadruplet::GetStationID(int station_num)
{
    /// \todo Assert that station_num is less than 3.
    return this->mStations[station_num]->GetIndex();
}

complex<double> Quadruplet::ComputeT4(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc)
{
    // Get the visibilities on the baselines AB, CD, AD, BC
    complex<double> AB = mBaselines[0]->GetVisibility(target, uv_ab);
    complex<double> CD = mBaselines[1]->GetVisibility(target, uv_cd);
    complex<double> AD = mBaselines[2]->GetVisibility(target, uv_ad); 
    complex<double> BC = mBaselines[3]->GetVisibility(target, uv_bc); 
    
    // Now compute the complex quad closure 
    return ( AB * BC ) / ( AD * conj(BC) );
}

/// Returns a boolean to indicate if the baseline specified by bl_name is involved in this quadruplet.
bool    Quadruplet::ContainsBaseline(string bl_name)
{
    /// \todo We could probably implment this as a hash lookup instead of a string operation
    /// this should make the lookup slightly faster.
    for(int i = 0; i < 4; i++)
    {
        if(this->mBaselines[i]->GetName() == bl_name)
            return true;
    }

    return false;
}

/// Returns the norm (amplitude) of the quad closures
double  Quadruplet::GetT4Amp(Target & target,  UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc)
{
    complex<double> t4 = GetT4(target, uv_ab, uv_cd, uv_ad, uv_bc );
    return abs(t4);
}

// Wrapper function for the UV-coordinate taking versions.
double  Quadruplet::GetT4Amp(Target & target, double hour_angle, double wavenumber)
{
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_cd = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_cd.Scale(wavenumber);
    UVPoint uv_ad = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ad.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[3]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);

    return GetT4Amp(target, uv_ab, uv_cd, uv_ad, uv_bc );
}

// Returns the norm (amplitude) of the quad closure
double  Quadruplet::GetT4AmpErr(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc)
{
    /// \bug Returns 0.001 by default
    return 0.001;
}
   

// Wrapper function for the UV-coordinate taking versions.
double  Quadruplet::GetT4AmpErr(Target & target, double hour_angle, double wavenumber)
{

    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_cd = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_cd.Scale(wavenumber);
    UVPoint uv_ad = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ad.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[3]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);
    
     return GetT4AmpErr(target, uv_ab, uv_cd, uv_ad, uv_bc );
}


// Wrapper function for the UV-coordinate taking versions.
double  Quadruplet::GetT4Phi(Target & target, double hour_angle, double wavenumber)
{
  
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_cd = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_cd.Scale(wavenumber);
    UVPoint uv_ad = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ad.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[3]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);
    
    return GetT4Phi(target,  uv_ab, uv_cd, uv_ad, uv_bc );
}

/// Returns the argument (phase) of the quad closures
double  Quadruplet::GetT4Phi(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc)
{
    complex<double> t4 = GetT4(target,  uv_ab, uv_cd, uv_ad, uv_bc );
    return arg(t4);
}


// Wrapper function for the UV-coordinate taking versions.
double  Quadruplet::GetT4PhiErr(Target & target, double hour_angle, double wavenumber)
{
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_cd = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_cd.Scale(wavenumber);
    UVPoint uv_ad = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ad.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[3]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);
    
    return GetT4PhiErr(target,  uv_ab, uv_cd, uv_ad, uv_bc );
}

double  Quadruplet::GetT4PhiErr(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc)
{
    /// \bug Returns 0.001 by default.
    return 0.001;
}   
    
// Computes the quad closures from the three baselines in this quadruplet.
complex<double> Quadruplet::GetT4(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc)
{
  string hash_key = GetHashKey(target,  uv_ab, uv_cd, uv_ad, uv_bc);
  complex <double> t4(0.0, 0.0);
  
  // First try looking up the value in the hash table
  if(mT4Values.find(hash_key) != mT4Values.end())
    {
      t4 = mT4Values[hash_key];
    }
  else
    {
        // The value did not exist in the hash table, we need to compute and store it.
      t4 = ComputeT4(target, uv_ab, uv_cd, uv_ad, uv_bc);
      mT4Values[hash_key] = t4;
    }
  
    return t4;
}

complex<double> Quadruplet::GetT4(Target & target, double hour_angle, double wavenumber)
{
    UVPoint uv_ab = mBaselines[0]->UVcoords(hour_angle, target.declination);
    uv_ab.Scale(wavenumber);
    UVPoint uv_cd = mBaselines[1]->UVcoords(hour_angle, target.declination);
    uv_cd.Scale(wavenumber);
    UVPoint uv_ad = mBaselines[2]->UVcoords(hour_angle, target.declination);
    uv_ad.Scale(wavenumber);
    UVPoint uv_bc = mBaselines[3]->UVcoords(hour_angle, target.declination);
    uv_bc.Scale(wavenumber);

    return GetT4(target, uv_ab, uv_cd, uv_ad, uv_bc); 
}

// Computes a hash key from the source, hour angle, and wavenumber.
string  Quadruplet::GetHashKey(Target & target, UVPoint uv_ab, UVPoint uv_cd, UVPoint uv_ad, UVPoint uv_bc)
{
    /// \todo It may be necessary for the doubles coming into this function to be cast into some 
    /// finite floating point format.
    
    /// \todo This function is in common with the Baseline class, need to factor this code.
    
    std::ostringstream sstream;
    sstream << target.GetName() << "-" << uv_ab.HashString() << '-' << uv_cd.HashString() << '-' << uv_ad.HashString() <<  '-' << uv_bc.HashString();
    std::string str = sstream.str();
    return str;
}

Baseline * Quadruplet::GetBaseline(int baseline_num)
{
    /// \todo Assert that baseline_num is less than 4
    return mBaselines[baseline_num];
}

////////////////////////////////////////////////////////////////////
// Non Class Functions Below
////////////////////////////////////////////////////////////////////

/// Computes all possible quadruplets for the stations provided.
vector<Quadruplet> ComputeQuadruplets(Array * array, vector<Station> & stations)
{
    // init some local variables
    vector<Quadruplet> quadruplets;
    int num_stations = stations.size();

    // Compute all possible quadruplets.  This should total to N choose 4
    for(int i = 0; i < num_stations - 3; i++)
    {
        for(int j = i + 1; j < num_stations - 2; j++)
        {
            for(int k = j + 1; k < num_stations - 1; k++)
            {	      
				for(int l = k + 1; l < num_stations ; l++)
				{
					quadruplets.push_back( Quadruplet(array, &stations[i], &stations[j], &stations[k], &stations[l] ) );
				}
            }
        }
    }
    
    return quadruplets;
}

QuadrupletHash ComputeQuadrupletHash(vector<Quadruplet> & quadruplets)
{
    QuadrupletHash hash;
    
    for(unsigned int i = 0; i < quadruplets.size(); i++)
    {
        hash.insert(QuadrupletHash::value_type(quadruplets[i].GetName(), &quadruplets[i]) );
    }
    
    return hash;
}
