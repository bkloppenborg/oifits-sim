// ============================================================================
// Name : matrix.cpp
// Version : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron and James Gordon
// ============================================================================

// Doxygen Items:
// / \file matrix.cpp
// / \todo Some work left to do. In particular to unifiy the double/complex
// routines into matrix.h

// This module defines some function used on 2D matrices or on vectors
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <fftw3.h>

#include "Matrix.h"
#include "Simulator.h"
#include "fitsio.h"


// ROUTINES FOR DOUBLE MATRICES
// / Calculates the norm (i.e. magnitude squared) for a row.
double norm2(Row < double >&r)
{
    double sum = 0.;

    for (int c = 0; c < r.size(); c++)
        sum += r[c] * r[c];
    return sum;

}

// / Calculates the norm (i.e. magnitude squared) for a matrix.
double norm2(Matrix < double >&m)
{
    double sum = 0.;

    for (int r = 0; r < m.GetRows(); r++)
        for (int c = 0; c < m.GetCols(); c++)
            sum += m[r][c] * m[r][c];
    return sum;

}

// / Computes the sum of all of the elements in a matrix.
double total(Matrix < double >&m)
{
    double sum = 0.;

    for (int r = 0; r < m.GetRows(); r++)
        for (int c = 0; c < m.GetCols(); c++)
            sum += m[r][c];
    return sum;

}

// / Calculates the scalar product of two matricies.
double scalprod(Matrix < double >&m1, Matrix < double >&m2)
{
    assert((m1.GetRows() == m2.GetRows()) && (m1.GetCols() == m2.GetCols()));
    double sum = 0;

    for (int r = 0; r < m1.GetRows(); r++)
        for (int c = 0; c < m1.GetCols(); c++)
            sum += m1[r][c] * m2[r][c];
    return sum;

}

// / Calculates the mean value of the elements in a matrix.
double mean(Matrix < double >&m)
{
    double sum = 0.;

    for (int r = 0; r < m.GetRows(); r++)
        for (int c = 0; c < m.GetCols(); c++)
            sum += m[r][c];
    return sum / double (m.GetRows() * m.GetCols());
}

// / Calculates the mean value of a row
double mean(Row < double >&row)
{
    double sum = 0.;

    for (int r = 0; r < row.size(); r++)
        sum += row[r];
    return sum / double (row.size());
}

// / Calculates the varaince of the elements in a matrix.
double variance(Matrix < double >&m)
{
    double sum = 0.;

    double average = mean(m);

    for (int r = 0; r < m.GetRows(); r++)
        for (int c = 0; c < m.GetCols(); c++)
            sum += (m[r][c] - average) * (m[r][c] - average);
    return sum / double (m.GetRows() * m.GetCols() - 1);
}

// / Calculates the variance for a row in a matrix.
double variance(Row < double >&row)
{
    double sum = 0.;

    double average = mean(row);

    for (int r = 0; r < row.size(); r++)
        sum += (row[r] - average) * (row[r] - average);
    return sum / double (row.size() - 1);
}

// ROUTINES FOR COMPLEX MATRICES

double norm2(Matrix < Complex > &m)
{
    double sum = 0.;

    for (int r = 0; r < m.GetRows(); r++)
        for (int c = 0; c < m.GetCols(); c++)
            sum += norm(m[r][c]);
    return sum;

}

Complex total(Matrix < Complex > &m)
{
    Complex sum(0, 0);

    for (int r = 0; r < m.GetRows(); r++)
        for (int c = 0; c < m.GetCols(); c++)
            sum += m[r][c];
    return sum;
}

Complex scalprod(Matrix < Complex > &m1, Matrix < Complex > &m2)
{
    assert((m1.GetRows() == m2.GetRows()) && (m1.GetCols() == m2.GetCols()));
    Complex sum(0, 0);

    for (int r = 0; r < m1.GetRows(); r++)
        for (int c = 0; c < m1.GetCols(); c++)
            sum += m1[r][c] * m2[r][c];
    return sum;
}

Complex scalprod(Matrix < Complex > &m1, Matrix < double >&m2)
{
    assert((m1.GetRows() == m2.GetRows()) && (m1.GetCols() == m2.GetCols()));
    Complex sum(0, 0);

    for (int r = 0; r < m1.GetRows(); r++)
        for (int c = 0; c < m1.GetCols(); c++)
            sum += m1[r][c] * m2[r][c];
    return sum;
}

Complex scalconj(Matrix < Complex > &m1, Matrix < Complex > &m2)        // total(A.B*)
{
    assert((m1.GetRows() == m2.GetRows()) && (m1.GetCols() == m2.GetCols()));
    Complex sum(0, 0);

    for (int r = 0; r < m1.GetRows(); r++)
        for (int c = 0; c < m1.GetCols(); c++)
            sum += m1[r][c] * conj(m2[r][c]);
    return sum;

}

void ZeroPad(Matrix < Complex > &unpadded, Matrix < Complex > &padded)  // for 
                                                                        // square 
                                                                        // matrices
{
    // Pad a matrix to padding x padding
    int unpadded_size = unpadded.GetCols();

    int padded_size = padded.GetCols();

    int low = (padded_size - unpadded_size) / 2;

    int high = (padded_size + unpadded_size) / 2;

    for (int ii = low; ii < high; ii++)
        for (int jj = low; jj < high; jj++)
            padded[ii][jj] = unpadded[ii - low][jj - low];

}

void ZeroUnPad(Matrix < Complex > &unpadded, Matrix < Complex > &padded)        // for 
                                                                                // square 
                                                                                // matrices
{
    // remove zero padding
    int unpadded_size = unpadded.GetCols();

    int padded_size = padded.GetCols();

    int low = (padded_size - unpadded_size) / 2;

    int high = (padded_size + unpadded_size) / 2;

    for (int ii = low; ii < high; ii++)
        for (int jj = low; jj < high; jj++)
            unpadded[ii - high][jj - high] = padded[ii][jj];
}

// Function to interpolate a 2d rectangular grid onto another 2d
// rectangular grid with arbitrary scaling, translation and
// rotation. Uses bilinear interpolation and uses periodic boundary
// conditions i.e. wraps the boundaries in the same way data from an
// FFT is wrapped.
// Taken from David Buscher and converted for use with this Matrix class by
// Fabien Baron
void Interp2d(Matrix < double >&dataIn, Matrix < double >&dataOut,
              double offsetX, double offsetY, double scale, double rotation)
{
    // Matrix<double> &dataIn = Input grid of data points from which to
    // extract a subgrid. The grid pitch is assumed to be the same
    // in both dimensions (this restriction only holds if a non-zero rotation
    // is required)
    // Matrix<double> &dataOut = Output gridded values.
    // double offsetX = X offset of (0,0) point in output grid in input grid
    // coordinate units.
    // If greater than nxIn, the coordinates are wrapped in the X dimension.
    // double offsetY = Y offset of (0,0) point in output grid in input grid
    // coordinate units.
    // If greater than nyIn, the coordinates are wrapped in the Y dimension
    // double scale = output grid pitch measured in units of the input grid
    // pitch
    // double rotation = rotation between input and output grids in radians.
    // Positive values rotate the output axes *clockwise* wrt the input axes

    int ix0, ix1, iy0, iy1;

    double x, y, intX0, intX1, intY0, intY1, fracX, fracY;

    int nxIn = dataIn.GetCols();

    int nyIn = dataIn.GetRows();

    int nxOut = dataOut.GetCols();

    int nyOut = dataOut.GetRows();

    double c = scale * cos(rotation);

    double s = scale * sin(rotation);

    for (int iy = 0; iy < nyOut; iy++)
    {
        for (int ix = 0; ix < nxOut; ix++)
        {
            x = c * ix + s * iy + offsetX;
            y = -s * ix + c * iy + offsetY;

            /* Get integer and fractional parts of x position */
            intX0 = floor(x);
            intX1 = intX0 + 1.0;
            fracX = x - intX0;

            /* Do modulo arithmetic correctly for negative numbers */
            ix0 = (int)(intX0 - nxIn * floor(intX0 / nxIn));
            ix1 = (int)(intX1 - nxIn * floor(intX1 / nxIn));

            /* Repeat for y */
            intY0 = floor(y);
            intY1 = intY0 + 1.0;
            fracY = y - intY0;
            iy0 = (int)(intY0 - nyIn * floor(intY0 / nyIn));
            iy1 = (int)(intY1 - nyIn * floor(intY1 / nyIn));

            dataOut[ix][iy] = fracX * (fracY * dataIn[ix1][iy1] + (1.0 - fracY)
                                       * dataIn[ix1][iy0]) + (1.0 -
                                                              fracX) * (fracY *
                                                                        dataIn
                                                                        [ix0]
                                                                        [iy1] +
                                                                        (1.0 -
                                                                         fracY)
                                                                        *
                                                                        dataIn
                                                                        [ix0]
                                                                        [iy0]);
        }
    }

}


// / Converts a matrix composed of doubles into a FITS file
void mat2fits(Matrix < double >&m, const char *filename)
{
    int status = 0;

    fitsfile *fptr;

    long fpixel = 1, naxis = 2, nelements;

    long naxes[2];

    // Initialise storage
    naxes[0] = (long)m.GetRows();
    naxes[1] = (long)m.GetCols();
    nelements = naxes[0] * naxes[1];

    double *ptrimg;

    // Create pointer image
    ptrimg = (double *)malloc(nelements * sizeof(double));

    for (int ii = 0; ii < naxes[0]; ii++)
        for (int jj = 0; jj < naxes[1]; jj++)
            ptrimg[ ii + jj * naxes[0] ] = m[ii][jj];

    // Create new file, write image, then close file
    if (status == 0)
        fits_create_file(&fptr, filename, &status);
    if (status == 0)
        fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    if (status == 0)
        fits_write_img(fptr, TDOUBLE, fpixel, nelements, &ptrimg[0], &status);
    if (status == 0)
        fits_close_file(fptr, &status);
    free(ptrimg);

    fits_report_error(stderr, status);

}

// / Converts a matrix composed of complex numbers into a FITS file
void mat2fits(Matrix < Complex > &m, char *filename)  // TBD: fix row major vs colum major
{
    int status = 0;

    fitsfile *fptr;

    long fpixel = 1, naxis = 2, nelements;

    long naxes[2];

    // Initialise storage
    naxes[0] = (long)m.GetRows();
    naxes[1] = (long)m.GetCols();
    nelements = naxes[0] * naxes[1];

    // Filenames
    char fitsreal[100];

    char fitsimag[100];

    for (int i = 0; i < 100; i++)
    {
        fitsreal[i] = '\0';
        fitsimag[i] = '\0';
    }

    strcpy(fitsreal, filename);
    strcpy(fitsimag, filename);
    strcat(fitsreal, "_r");
    strcat(fitsimag, "_i");

    double *ptrimgreal = (double *)malloc(nelements * sizeof(double));

    for (int ii = 0; ii < naxes[0]; ii++)
        for (int jj = 0; jj < naxes[1]; jj++)
            ptrimgreal[ii + jj * naxes[0]] = m[ii][jj].real();

    // Create new file, write image, then close file
    if (status == 0)
        fits_create_file(&fptr, fitsreal, &status);
    if (status == 0)
        fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    if (status == 0)
        fits_write_img(fptr, TDOUBLE, fpixel, nelements, &ptrimgreal[0],
                       &status);
    if (status == 0)
        fits_close_file(fptr, &status);

    free(ptrimgreal);

    double *ptrimgimag = (double *)malloc(nelements * sizeof(double));

    for (int ii = 0; ii < naxes[0]; ii++)
        for (int jj = 0; jj < naxes[1]; jj++)
            ptrimgimag[ii + jj * naxes[0]] = m[ii][jj].imag();

    // Create new file, write image, then close file
    if (status == 0)
        fits_create_file(&fptr, fitsimag, &status);
    if (status == 0)
        fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    if (status == 0)
        fits_write_img(fptr, TDOUBLE, fpixel, nelements, &ptrimgimag[0],
                       &status);
    if (status == 0)
        fits_close_file(fptr, &status);

    free(ptrimgimag);

    fits_report_error(stderr, status);

}

// / Converts a FITS file into a matrix of double values.
void fits2mat(const char *fname, Matrix < double >&m)
{
    int status = 0;

    char comment[200];          // should be enough for most comments

    fitsfile *fptr;             // pointer to the FITS file, defined in
                                // fitsio.h
    int nfound, anynull;

    long naxes[2], fpixel, npixels;

    float nullval;

    if (status == 0)
        fits_open_file(&fptr, fname, READONLY, &status);
    if (status == 0)
        fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);

    // The pixellation parameter has been moved to the SOURCE/Target definition file.
//    if (status == 0)            // Get the pixellation parameter
//    {
//        float xyint;

//        fits_read_key_flt(fptr, "PIXEL", &xyint, comment, &status);
//        if (status != 0)        // in the case where the keyword is missing
//        {
//            status = 0;
//            pixellation = 0.;

//            // The file should have a pixelation value.  Throw an exception.
//            /// \exception runtime_error The Target FITS file is missing PIXEL keyword!
//            /// This means that the input FITS file is missing the PIXEL(lation) parameter
//            /// which should be defined somewhere in the FITS header.
//            throw std::runtime_error("Target FITS File missing PIXEL keyword!");
//        }
//        else
//        {
//            pixellation = double (xyint);
//        }
//    }

    npixels = naxes[0] * naxes[1];      // number of pixels in the image
    fpixel = 1;
    nullval = 0;                // do not check for null values

    double *ptrimg = (double *)malloc(npixels * sizeof(double));

    if (status == 0)
        fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval, ptrimg,
                      &anynull, &status);
    if (status == 0)
        fits_close_file(fptr, &status);

    m.setsize(naxes[0], naxes[1]);
     for (int ii = 0; ii < naxes[0]; ii++)
        for (int jj = 0; jj < naxes[1]; jj++)
            m[ii][jj] = ptrimg[ii + jj * naxes[0]];
    free(ptrimg);

    // Report any errors
    fits_report_error(stderr, status);

}

// / \brief Compute the FFT of a matrix
// /
// / Note: this FFT is normalized.

void FFT(Matrix < Complex > &min, Matrix < Complex > &mout, int direction, int dccenter)
{
    assert(min.GetCols() == min.GetRows());     // square matrix
    assert(mout.GetCols() == mout.GetRows());
    assert(mout.GetCols() == min.GetRows());
    int nfft = min.GetCols();

    fftw_complex *in =
        (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nfft * nfft);
    fftw_complex *out =
        (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nfft * nfft);

    for (int iy = 0; iy < nfft; iy++)
    {
        for (int ix = 0; ix < nfft; ix++)
        {
            in[ix + nfft * iy][0] = min[ix][iy].real();
            in[ix + nfft * iy][1] = min[ix][iy].imag();
        }
    }

    fftw_plan p;

    if (direction >= 0)
        p = fftw_plan_dft_2d(nfft, nfft, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    else
        p = fftw_plan_dft_2d(nfft, nfft, in, out, FFTW_BACKWARD,
                             FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_free(in);

    if (dccenter == 0)
    {
        for (int iy = 0; iy < nfft; iy++)
            for (int ix = 0; ix < nfft; ix++)
                mout[ix][iy]
                    = Complex(out[ix + nfft * iy][0], out[ix + nfft * iy][1])
                    / double (nfft);
    }
    else
    {
        for (int iy = 0; iy < nfft / 2; iy++)
        {
            for (int ix = 0; ix < nfft / 2; ix++)
            {
                mout[ix][iy] =
                    Complex(out[nfft / 2 + ix + nfft * (iy + nfft / 2)][0],
                            out[nfft / 2 + ix +
                                nfft * (iy + nfft / 2)][1]) / double (nfft);
            }

            for (int ix = nfft / 2; ix < nfft; ix++)
            {
                mout[ix][iy] =
                    Complex(out[ix - nfft / 2 + nfft * (iy + nfft / 2)][0],
                            out[ix - nfft / 2 +
                                nfft * (iy + nfft / 2)][1]) / double (nfft);
            }
        }

        for (int iy = nfft / 2; iy < nfft; iy++)
        {
            for (int ix = 0; ix < nfft / 2; ix++)
            {
                mout[ix][iy] =
                    Complex(out[ix + nfft / 2 + nfft * (iy - nfft / 2)][0],
                            out[ix + nfft / 2 +
                                nfft * (iy - nfft / 2)][1]) / double (nfft);
            }

            for (int ix = nfft / 2; ix < nfft; ix++)
            {
                mout[ix][iy] =
                    Complex(out[ix - nfft / 2 + nfft * (iy - nfft / 2)][0],
                            out[ix - nfft / 2 +
                                nfft * (iy - nfft / 2)][1]) / double (nfft);
            }
        }
    }
    fftw_free(out);
}

