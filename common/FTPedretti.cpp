//============================================================================
// Name        : simulator.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

// This module implements the group-delay tracking algorithm described in:
// Robust determination of optical path difference: fringe tracking
// at the Infrared Optical Telescope Array Interferometer, E. Pedretti et al.
// Applied Optics Vol. 44 No. 25

#include <iostream>

#include "FTPedretti.h"
#include "Simulator.h"
#include "SpectralMode.h"
#include "Modulator.h"

using std::cout;

FTPedretti::FTPedretti( int n_coh, int ndelta_coh, int n_incoh, SpectralMode* spectr , Modulator* modulator)  
{
	// Init classic variables
	this->nchannels = spectr->nchannels;
	this->spectr =  spectr;
	this->modulator = modulator;
	this->n_coh = n_coh; 
	this->ndelta_coh = ndelta_coh;
	this->n_incoh = n_incoh;
	
	// Init acquisition variables
	framecounter = 0;
	previous_time = 0.0; // assumes simulation starts a t(0)=0
	frametime.setsize(n_coh);
	frameduration.setsize(n_coh);
	modulator_position.setsize(n_coh);
	deltamod.setsize(n_coh);
	frames.setsize( nchannels , n_coh );
	
	
	ntilde.setsize( nchannels );
	X.setsize( nchannels - 1 );
	
	framecounter = 0;
	for( int ii=0 ; ii < n_coh ; ii++ ) frametime[ ii ] = 0. ;
	// Specific to Pedretti's algo
	delta_m12 = spectr->delta_wavenumber[ nchannels / 2 ];
	
}

void FTPedretti::StoreOneFrame( Row<int>& detectorframe, double time )
{
	// acquire fringes, detect potential problems (bad frames), convert to double precision
	
	for(int channel = 0 ; channel < nchannels ; channel++) frames[ channel ][ framecounter ] = (double)( detectorframe[ channel ] ) ; 
	
	// Storing data about the current frame
	frametime[ framecounter ] = (time + previous_time) / 2. ;
	frameduration[ framecounter ] = time - previous_time ;
	modulator_position[ framecounter ] = modulator->GetDelay( frametime[ framecounter ] );
	deltamod[ framecounter ] = modulator_position[ framecounter ] - modulator->GetDelay( previous_time );
	previous_time = time;
	//cout << framecounter << endl;
	if ( framecounter == (n_coh - 1) )
	{
		framecounter = 0 ; // reset frame counter for the next batch of coherent integration 	
		
		GDEstimate(); // launch group delay estimation on n_coh frames
		
	} else framecounter++;
	
}

void FTPedretti::DFT()
{
	
	for( int channel = 0 ; channel < nchannels ; channel++ )		
	{
		ntilde[ channel ]  = ZEROCOMP;
		
		for( int framenumber = 0 ; framenumber < n_coh ; framenumber++ )
		{
			ntilde[ channel ] += frames[ channel ][ framenumber ] * exp( - 2. * PI * I * spectr->mean_wavenumber[ channel ] * modulator_position[ framenumber ] ) ; 
		}
	}
	
	// The complex visibility is A + i B, with: 
	// A = C cos(D) + sin( C + D )
	// B = C sin(D) + cos( C + D ) 
	// C = 4*pi*wavenumber*scan_ampl*T
	// D = 2*pi*wavenumber*OPD
	// Computed with Maple:
	// f:=t->int(cos(2*Pi*s*(a*t+OPD))*exp(-2*I*Pi*s*a*t), t);
	// simplify(evalc(Im(f(T))/Re(f(T))));
	
	// Using arctan( B / A ) to get the phase only works under certain conditions
	// C >>> 1, preferably > 50, ensures this works adequately,
	// That means that the integration time should be T > 50/(4*pi*scan_ampl*wavenumber)
	// Then B/A = C cos(D)/C sin(D) = tan(D)
	// which gives OPD = atan(B/A) / ( 2*pi*wavenumber )
	// The Cross-spectra is also generally more precise for this reason (amongst others)
	// Note that this condition is independent of the value of the OPD
	
	
}


void FTPedretti::ComputeCrossSpectra()
{
	for(int channel = 0 ; channel < nchannels - 1 ; channel++)		
	{
		X[ channel ] =  ntilde[ channel + 1 ] * conj( ntilde[ channel ]) ;
	}
	
}

void FTPedretti::AverageSpectrum()
{
	X_avg = ZEROCOMP ;
	
	for(int channel = 0 ; channel < nchannels - 1 ; channel++)		
	{
		X_avg += X[ channel ] ;
	}
	
	X_avg /= nchannels;
}

void FTPedretti::GDEstimate()
{
	DFT();
	ComputeCrossSpectra();
	AverageSpectrum();
	
	double OPD_estimate = arg( X_avg ) / ( 2. * PI * delta_m12 ) ;
	cout << "Pedretti OPD: " << OPD_estimate * 1e6 << endl;
}
