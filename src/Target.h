/*
 * Target.h
 *
 *  Created on: Apr 8, 2011
 *      Author: bkloppenborg
 */

#ifndef TARGET_H_
#define TARGET_H_

#include <string>
#include "textio.hpp"
#include "Matrix.h"

// Header files for other libraries
extern "C" {
    #include "exchange.h"
}

using namespace std;

class Target
{
public:
	string name;
	double pixellation;
	double mag;
	char band;
	double temperature;
	double background;
	double background_ap;
	double declination;
	double right_ascension;

	Matrix < double > image;

	// A few computed values
	double flux;
	double flux_bg;

public:
	Target();
	virtual ~Target();

	void ImportFile(string filename, string comment_chars);

	void ParseFileOptions(char *argv[], int i, int argc);
	void ParseImageOptions(char *argv[], int i, int argc);
	void ParseOptions(char *argv[], int i, int argc);

	void SetImage(string filename);

	string GetName();
	int GetTargetID(void);

	double GetTargNPhotons(double wavelength, double bandwidth, double area, double time);
	double GetBackNPhotons(double wavelength, double bandwidth, double area, double time);

    oi_target  GetOITarget(void);

private:
    void InitFlux();
};

#endif /* TARGET_H_ */
