/// \file ZernikeLib.h

#include "Matrix.h"

int factorial(int n);

class ZernikeLib
{
  public:
    int number;

    int size;

    Row < Matrix < double > > modes;

    double phase_to_a0;

    ZernikeLib(int number, int size);

    virtual ~ ZernikeLib();
};
