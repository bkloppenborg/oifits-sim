// ============================================================================
// Name : matrix.cpp
// Version : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron and James Gordon
// ============================================================================

// / \file matrix.h

// This file defines the matrix and vector classes used in the simulator

#ifndef MATRIX_H
#define MATRIX_H

#include <complex>
#include <cstring>
#include <assert.h>

using namespace std;

typedef std::complex < double >Complex;

template < class T > class Row
{
  public:
  Row(int n = 0):
    row(NULL)
    {
        setsize(n);
    }

    ~Row()
    {
        if (row)
            delete[]row;
    }

    Row(const Row & r):row(NULL)
    {
        setsize(r.ncols);
        for (int i = 0; i < ncols; i++)
            row[i] = r.row[i];
    }

    void setsize(int n)
    {
        if (row)
        {
            delete[]row;
        }
        if (n > 0)
        {
            row = new T[n];
            memset(row, 0, sizeof(T) * n / sizeof(char));       // init the
            // row at 0
        }
        else
            row = NULL;
        ncols = n;
    }

    int size()
    {
        return ncols;
    }

    T & operator[](int column)
    {
        assert(column < ncols);
        assert(column >= 0);
        return row[column];
    }

    Row & operator=(const Row & r)
    {
        if (this != &r)         // check for self-assignment
        {
            setsize(r.ncols);
            for (int i = 0; i < ncols; i++)
                row[i] = r.row[i];
        }
        return *this;
    }

    Row & operator=(const T s)
    {
        for (int i = 0; i < ncols; i++)
            row[i] = s;
        return *this;
    }

    const Row operator-(const T s)
    {
        Row theRow(ncols);

        for (int i = 0; i < ncols; i++)
            theRow[i] = row[i] - s;
        return theRow;
    }

    Row & operator+=(const T s)
    {
        for (int c = 0; c < ncols; c++)
            row[c] += s;
        return *this;
    }

    friend ostream & operator<<(ostream & os, const Row & r)
    {
        // os << "Ncols: " << r.ncols << " " << endl;
        for (int i = 0; i < r.ncols; i++)
        {
            os << r.row[i] << " ";
        }
        os << endl;
        return os;
    }

  private:
    int ncols;

    T *row;
};

template < class T > class Matrix
{
  public:
  Matrix(int rows = 0, int cols = 0):
    matrix(NULL)
    {
        setsize(rows, cols);
    }
    Matrix(const Matrix & m):matrix(NULL)
    {
        setsize(m.nrows, m.ncols);
        for (int r = 0; r < nrows; r++)
            matrix[r] = Row < T > (m.matrix[r]);
    }
    ~Matrix()
    {
        if (matrix)
            delete[]matrix;
    }

    void setsize(int rows, int cols)
    {
        if (matrix)
            delete[]matrix;
        if (cols > 0 && rows > 0)
        {
            matrix = new Row < T >[rows];
            for (int i = 0; i < rows; i++)
                matrix[i].setsize(cols);
        }
        else
            matrix = NULL;
        ncols = cols;
        nrows = rows;
    }

    int GetCols()
    {
        return ncols;
    }
    int GetRows()
    {
        return nrows;
    }
    int type()
    {
        return sizeof(T);
    }

    Row < T > &operator[](int index)
    {
        assert(index < nrows);
        assert(index >= 0);
        return matrix[index];
    }

    Matrix & operator=(const Matrix & m)
    {
        if (this != &m)         // check for self-assignement
        {
            setsize(m.nrows, m.ncols);
            for (int i = 0; i < nrows; i++)
                matrix[i] = Row < T > (m.matrix[i]);

        }
        return *this;
    }

    const Matrix operator+(const Matrix & m)
    {
        assert(ncols == m.ncols && nrows == m.nrows);
        Matrix theMatrix(nrows, ncols);

        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                theMatrix[r][c] = matrix[r][c] + m.matrix[r][c];
        return theMatrix;
    }

    const Matrix operator-(const Matrix & m)
    {
        assert(ncols == m.ncols && nrows == m.nrows);
        Matrix theMatrix(nrows, ncols);

        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                theMatrix[r][c] = matrix[r][c] - m.matrix[r][c];
        return theMatrix;
    }

    Matrix & operator+=(const Matrix & m)
    {
        assert(ncols == m.ncols && nrows == m.nrows);
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] += m.matrix[r][c];
        return *this;
    }

    Matrix & operator-=(const Matrix & m)
    {
        assert(ncols == m.ncols && nrows == m.nrows);
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] -= m.matrix[r][c];
        return *this;
    }

    const Matrix operator*(const double s)
    {
        Matrix theMatrix(nrows, ncols);

        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                theMatrix[r][c] = matrix[r][c] * s;
        return theMatrix;
    }

    const Matrix operator/(const double s)
    {
        Matrix theMatrix(nrows, ncols);

        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                theMatrix[r][c] = matrix[r][c] / s;
        return theMatrix;
    }

    const Matrix operator+(const double s)
    {
        Matrix theMatrix(nrows, ncols);

        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                theMatrix[r][c] = matrix[r][c] + s;
        return theMatrix;
    }

    const Matrix operator-(const double s)
    {
        Matrix theMatrix(nrows, ncols);

        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                theMatrix[r][c] = matrix[r][c] - s;
        return theMatrix;
    }

    Matrix & operator+=(const double s)
    {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] += s;
        return *this;
    }

    Matrix & operator-=(const double s)
    {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] -= s;
        return *this;
    }

    Matrix & operator*=(const double s)
    {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] *= s;
        return *this;
    }

    Matrix & operator/=(const double s)
    {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] /= s;
        return *this;
    }

    Matrix & operator=(const double s)
    {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] = s;
        return *this;
    }

    Matrix & operator*=(const Matrix & m)
    {
        assert(ncols == m.ncols && nrows == m.nrows);
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                matrix[r][c] *= m.matrix[r][c];
        return *this;
    }

    // Defining * as the element by element multiplication
    Matrix operator*(const Matrix & m)
    {
        assert(ncols == m.ncols && nrows == m.nrows);
        Matrix theMatrix(nrows, ncols);

        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < ncols; c++)
                theMatrix[r][c] = matrix[r][c] * m.matrix[r][c];
        return theMatrix;
    }

    friend ostream & operator<<(ostream & os, const Matrix & m)
    {
        // os << "Nrow: " << m.nrows<<" Ncols: "<<m.ncols<<" "<<endl;
        for (int r = 0; r < m.nrows; r++)
        {
            for (int c = 0; c < m.ncols; c++)
            {
                os << m.matrix[r][c] << " ";

            }
            os << endl;
        }
        return os;
    }

    friend istream & operator>>(istream & is, Matrix & m)
    {
        int rows, cols;

        is >> rows >> cols;
        m.setsize(rows, cols);
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
                is >> m[r][c];
        return is;
    }

  private:
    int ncols, nrows;

    Row < T > *matrix;

};

void Interp2d(Matrix < double >&dataIn, Matrix < double >&dataOut, double offsetX,
         double offsetY, double scale, double rotation);
double norm2(Matrix < double >&m);

double total(Matrix < double >&m);

double scalprod(Matrix < double >&m1, Matrix < double >&m2);

double norm2(Matrix < Complex > &m);

double mean(Matrix < double >&m);

double variance(Matrix < double >&m);

double variance(Row < double >&row);

double mean(Row < double >&m);

double norm2(Row < double >&r);

Complex total(Matrix < Complex > &m);

Complex scalprod(Matrix < Complex > &m1, Matrix < Complex > &m2);       // total(A.B)

Complex scalconj(Matrix < Complex > &m1, Matrix < Complex > &m2);       // total(A.B*)

Complex scalprod(Matrix < Complex > &m1, Matrix < double >&m2);

void FFT(Matrix < Complex > &min, Matrix < Complex > &mout, int direction, int dccenter);

void mat2fits(Matrix < double >&m, const char *filename);

void mat2fits(Matrix < Complex > &m, const char *filename);

void fits2mat(const char *fname, Matrix < double >&m);

void ZeroPad(Matrix < Complex > &unpadded, Matrix < Complex > &padded);

void ZeroUnPad(Matrix < Complex > &unpadded, Matrix < Complex > &padded);

#endif // #ifndef MATRIX_H
