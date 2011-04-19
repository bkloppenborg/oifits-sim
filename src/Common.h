/*
 * Constants.h
 *
 *  Created on: Apr 11, 2011
 *      Author: bkloppenborg
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// ============================================================================
// Name : simulator.cpp
// Version : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
// ============================================================================

/// \mainpage
/// The OIFITS simulator is a C++ application for simulating the observables
/// from interferometers using a realistic noise model.
///
/// This program uses the following libraries:
/// \ref OIFITSlib
///

// include file for the simulator
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <complex>
#include <vector>
using namespace std;

// First define some constants:
/// WGS84 Earth equatorial radius /m
static const double WGS_A0 = 6378137;

/// WGS84 reference spheroid flattening factor
static const double WGS_F = 1.0 / 298.257223563;

/// Global throughput.  Does not include the quantum efficiencey.
/// This value is used in FastOneRead function.
const double instrument_throughput = 0.1324 / 0.65;

const double instrument_visibility = 1.;        // visibility losses

const double PI = 3.141592653589793238462643;

const double milliarcsec = (PI / 180.0) / 3600000.0;
const complex<double> I(0.0, 1.0);

const complex<double> ZEROCOMP(0., 0.);

// Now define a few functions
void wgs84_to_geoc(double lat, double height, double *GeocLat, double *GeocRadius);
void testzern();

double phasemagrun(double magnitude);

double delaymagrun(double magnitude);

double simplyfringes(double magnitude);

double sinc(double number);

void    StringSplit(string str, string delim, vector<string> & results);
string  StripWhitespace(string str);
void    StripWhitespace(vector<string> & strings);
vector<string> ReadFile(string filename, string comment_chars, string error_message);

void Swap(int *a, int *b);
void Sort(int &a, int &b);
void Sort(int &a, int &b, int &c);

#endif // #ifndef SIMULATOR_H

#endif /* CONSTANTS_H_ */
