#include "Simulator.h"

double sinc( double number )
{
  if (abs(number) < 1.0e-10)
    return 1.0;
  return sin(number) / number;
}
