#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include "Instrument.h"

using std::vector;
using std::string;
using std::ifstream;
using std::cout;
using std::endl;

Instrument::Instrument(double throughput, double visibility, int Npix,
   double read_noise, double quantum_efficiency,
   double vsq_frac_cal_err, double clo_cal_err,
   double wind_speed, double ast_seeing, double incoh_time)
{
	this->throughput = throughput;
	this->visibility = visibility;
	this->Npix = Npix;
	this->read_noise = read_noise;
	this->quantum_efficiency = quantum_efficiency;
	this->vsq_frac_cal_err = vsq_frac_cal_err;
	this->clo_cal_err_deg = clo_cal_err;
	this->wind_speed = wind_speed;
	this->ast_seeing = ast_seeing;
	this->incoh_time = incoh_time;
}

Instrument::Instrument(const char *Inst_file)
{
	const string comments("\\/#~$&£%");
	vector < string > lines;	// stores non-blank, non-comment lines

	ifstream fil(Inst_file);
	if (fil.is_open())
	{
		string line;
		while (!fil.eof())
		{
			getline(fil, line);
			while ((line.size() == 0
				  || comments.find(line[0]) != string::npos) && !fil.eof())
				getline(fil, line);
			if (!fil.eof())
				lines.push_back(line);
		}
		fil.close();
	}
	else
	{
	    /// \except runtime_error "Error Opening Instrument File"
		throw std::runtime_error("Error opening instrument file");
	}

	if (lines.size() != 10)
		throw
		   std::runtime_error
		   ("Invalid number of parameters in the instrument file");

	throughput = atof(lines[0].c_str());
	visibility = atof(lines[1].c_str());
	Npix = atoi(lines[2].c_str());
	read_noise = atof(lines[3].c_str());
	quantum_efficiency = atof(lines[4].c_str());
	vsq_frac_cal_err = atof(lines[5].c_str());
	clo_cal_err_deg = atof(lines[6].c_str());
	wind_speed = atof(lines[7].c_str());
	ast_seeing = atof(lines[8].c_str());
	incoh_time = atof(lines[9].c_str());
	cout << "throughput = " << throughput << endl;	// 
	cout << "inst visibility = " << visibility << endl;	// 
	cout << "Npix = " << Npix << endl;	// 
	cout << "clo phase cal err = " << clo_cal_err_deg << " deg" << endl;	// 
	cout << "Wind speed = " << wind_speed << " m/s" << endl;
	cout << "Fried parameter@500nm = " << ast_seeing << " m" << endl;
	cout << "Incoherent integration time = " << incoh_time << " s" << endl;
}

/** Return total interferometer throughput, including detector QE. */
double Instrument::GlobalThroughput()
{
	return throughput * quantum_efficiency;
}
