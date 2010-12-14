//
//  ZernikeLib Class
//  - Generates a basis of n Zernike modes
//  - Orthonormalizes them using Gram-Schmidt method
//  - Stores the resulting discretised Zernike modes in a single tridimensional array
//
//
//  24/10/2008 F. Baron, version 1.0
//
//

// Next two are required only for the tests in main()

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "ZernikeLib.h"

using std::cout;

int factorial(int n) 
{
	return ( n < 2 ) ? 1 : n * factorial( n - 1 );
}

ZernikeLib::ZernikeLib(int number, int size) {
	this->number = number;
	this->size = size;
	Row< Matrix<double> > basis(number);
	double x, y, centre, rsq, r, phi, sum;
	centre = ( (double)size - 1.0 ) / 2.0;
	int maxnradterm = int( 1.5 + number / 2 );
	int n, l, m, nradterm;
	Row<int> radial_powers( maxnradterm );
	Row<double> radial_coefficients( maxnradterm );
	Matrix<double> temp( size , size );
	
	// Generate #number consecutive Zernike modes
	for(int ii = 0; ii < number; ii++) {
		n=int( ceil( -1.5 + .5 * sqrt( 9 + 8 * ii ) ) );
		l=2 * ii - n * ( n + 2 );
		m=abs( l );
		if( ( ( n - m ) & 1 ) == 0 ) // (n -m) is even
		{
			// Reset the radial power coefficients - we only have to compute them up to nradterm
			radial_powers = 0;
			radial_coefficients = 0. ;
			nradterm = 1 + ( ( n - m ) / 2 );
			// Compute the radial coefficients only once per ZernikeLib mode
			// DEBUG cout << "Evaluation of polynom: " << ii << " with n=" << n << " and m=" << m << endl;
			
			int sign = -1;
			for ( int s = 0; s < nradterm; s++ ) {
				sign *= -1;
				radial_powers[ s ] = n - 2 * s;
				radial_coefficients[ s ]  = sign * factorial( n - s ) / ( factorial( s ) * factorial( ( n + m ) / 2 - s) * factorial( ( n - m ) / 2 - s ) );
			}
			
			// Now loop on the x and y of the Zernike
			temp = 0.;
			for (int iy = 0; iy < size; iy++ )
			{
				for (int ix = 0; ix < size; ix++ )
				{
					x = ( ix - centre ) / ((double)size / 2.);
					y = ( iy - centre ) / ((double)size / 2.);
					rsq = x * x + y * y ;
					if (rsq >= 1.0) 
					{
						temp[ ix ][ iy ] = 0.;
					} 
					else
					{
						// Initialize the sum to storing the pixel value
						sum = 0.0;
						// Computation of the radial function
						r = sqrt( rsq );
						for ( int s = 0; s < nradterm; s++ ) sum += radial_coefficients[ s ] * pow( r , radial_powers[ s ] );
						
						// Computation of the trigonometric function
						// No need to normalise this, will be done during the Gram-Schmidt orthonormalization
						// check if ok with 0
						if(m !=0 )
						{
							// Then there are multiple modes sharing the same radial part
							// even ranks = cos, odd ranks = sin
							phi = atan2( y , x );
							if ( l > 0 )
								sum *= cos( m * phi );
							else
								sum *= sin( m * phi );
							
						}
						temp[ ix ][ iy ] = sum;
					}
					
				}
				
			}
			basis[ii]=temp;
			
		}
		
	}
	
	// Gram-Schmidt loop
	
	// Redimension modes vector to a size equal to "number"
	modes.setsize(number);	
	
	Matrix<double> orthog( size, size );
	double crossproduct;
	for(int ii=0 ; ii < number ; ii++ ) {
		// Starting from the previously computed Noll Basis
		orthog = basis[ii];
		// Orthogonalisation for modes greater than piston
		if(ii !=0 ) {
			for(int jj=0; jj<ii; jj++ ) {
				temp = modes[jj];
				crossproduct = scalprod( orthog , temp );
				orthog -= temp * crossproduct;
			}
		}
		
		// Normalisation
		orthog = orthog * ( 1./ sqrt( norm2( orthog ) ) );
		if(ii == 0) phase_to_a0 = orthog[ size/2 ][ size/2 ];
		// Affectation to the ZernikeLib class
		modes[ ii ] = orthog;

	}
	
}

ZernikeLib::~ZernikeLib() {
	// Nothing to do
}

void testzern() {
	int size = 256;
	int n = 20;
	ZernikeLib zernlib = ZernikeLib( n, size );
	Matrix<double> test1( size , size );
	Matrix<double> test2( size , size );
	double prod;
	
	cout << "Generating " << n << " Zernike modes of " << size << " x " << size << " pixels" << endl;
	
	// Test for the orthonormalisation
	for(int ii = 0; ii < n ; ii++) 
	{
		test1= zernlib.modes[ii];
		cout << "Mode "<< ii << endl;
		for(int jj = 0 ; jj <= ii ; jj++ ) 
		{
			test2 = zernlib.modes[jj];
			prod = abs( scalprod( test1 , test2 ) );
			if( prod < 1e-14) prod = 0.; // improves readability
			cout << " Scalar product of " << ii << " by "<< jj << " = " << prod  << endl;
		}
	}
	
	// Output modes to files 
	ofstream file;
	for(int ii = 0; ii < n ; ii++) {
		test1 = zernlib.modes[ii];
		file.open( "zernike.txt" );
		file << test1;
		file.close();
	//	mat2fits(test1, "zernike.txt");
		cout << "Mode " << ii << " is ready to be read from file zernike.txt" << endl;
		getchar();
	}
	
}
