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

/// Global throughput.  Does not include the quantum efficiencey.
/// This value is used in FastOneRead function.
const double instrument_throughput = 0.1324 / 0.65; 

const double instrument_visibility = 1.;        // visibility losses

const double PI = 3.141592653589793238462643;

const double milliarcsec = (PI / 180.0) / 3600000.0;
const complex<double> I(0.0, 1.0);

const complex<double> ZEROCOMP(0., 0.);

void testzern();

double phasemagrun(double magnitude);

double delaymagrun(double magnitude);

double simplyfringes(double magnitude);

double sinc(double number);

void    StringSplit(string str, string delim, vector<string> results);
string  StripWhitespace(string str);
void    StripWhitespace(vector<string> strings);
vector<string> ReadFile(string filename, string comment_chars, string error_message);

#endif // #ifndef SIMULATOR_H
