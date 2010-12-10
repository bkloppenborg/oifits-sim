#include "Simulator.h"

#include <iostream>
#include <fstream>
#include <stdexcept>

template<typename T, typename P>
T remove_if(T beg, T end, P pred)
{
    T dest = beg;
    for (T itr = beg;itr != end; ++itr)
        if (!pred(*itr))
            *(dest++) = *itr;
    return dest;
}

double sinc( double number )
{
  if (abs(number) < 1.0e-10)
    return 1.0;
  return sin(number) / number;
}

// A simple string splitting function for use with delimiters
void StringSplit(string str, string delim, vector<string> & results)
{
    int cutAt;
    while( (cutAt = str.find_first_of(delim)) != str.npos )
    {
        if(cutAt > 0)
            results.push_back(str.substr(0,cutAt));
        
        str = str.substr(cutAt+1);
    }
    
    // Push any remaining characters on to the back of the vector.
    if(str.length() > 0)
        results.push_back(str);
}

string StripWhitespace(string str)
{
    str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
    return str;
}


void StripWhitespace(vector<string> & strings)
{
    // Remove any remaining white space.
    for(unsigned int i = 0; i < strings.size(); i++)
    {
        strings[i] = StripWhitespace(strings[i]);
    }
}

/// Reads in a file, returns non-comment lines as a vector of strings
/// If the file cannot be opened, an exception is thrown with the message contained in error_message
vector<string> ReadFile(string filename, string comment_chars, string error_message)
{
    ifstream infile;
    infile.open(filename.c_str());
    vector < string > lines;   
    
    if (infile.is_open())
    {
        string line;

        while (!infile.eof())
        {
            getline(infile, line);
            while ((line.size() == 0 || comment_chars.find(line[0]) != string::npos)
                   && !infile.eof())
            {
                getline(infile, line);
            }
            if (!infile.eof())
                lines.push_back(line);
        }
        infile.close();
    }
    else
    {
        /// \exception runtime_error Error opening array file
        /// The array file could not be located.  It is likely that the user just specified an invalid path.
        throw std::runtime_error(error_message);
    }
    
    return lines;
}

// AS 2010-06-22
// added function from JSY's sidereal.c routine to compute the true
// geocentric coordinates 
// of the array location
/**
 * Convert position on Earth's surface from geodetic to geocentric coordinates.
 *
 * Uses WGS84 geoid parameters.
 *
 * @param lat     Geodetic latitude /radians.
 * @param height  Height above geoid /metres.
 * @param GeocLat     Return location for geocentric latitude /radians.
 * @param GeocRadius  Return location for geocentric radius /metres.
 */
void wgs84_to_geoc(double lat, double height, double *GeocLat, double *GeocRadius)
{
	double cosLat, sinLat, c, c0, s0, one_f_sq;

	cosLat = cos(lat);
	sinLat = sin(lat);
	one_f_sq = (1. - WGS_F) * (1. - WGS_F);
	c = pow(cosLat * cosLat + one_f_sq * sinLat * sinLat, -0.5);
	c0 = WGS_A0 * c + height;
	s0 = WGS_A0 * one_f_sq * c + height;

	*GeocLat = atan2(s0 * sinLat, c0 * cosLat);
	*GeocRadius = pow(c0 * c0 * cosLat * cosLat + s0 * s0 * sinLat * sinLat, 0.5);
}

