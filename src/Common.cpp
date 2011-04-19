
#include "Common.h"

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

void Swap(int *a, int *b) 
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

void Sort(int &a, int &b)
{
    if(a > b)
        Swap(&a, &b);
}

void Sort(int &a, int &b, int &c) 
{
    if (a > b) 
        Swap(&a,&b);
    if (a > c) 
        Swap(&a,&c);
    if (b > c) 
        Swap(&b,&c);
}

