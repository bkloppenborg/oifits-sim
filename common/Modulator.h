/// \file Modulator.h
/// Header file for the modulator class

#ifndef MODULATOR_H
#define MODULATOR_H

/// \class Modulator
class Modulator
{
  public:
    Modulator(int modulationtype, double period, double amplitude,
              double dephasing);
    int modulationtype;

    double amplitude;

    double period;

    double dephasing;

    double GetDelay(double time);

    double GetDifferential(double time);

  private:
    double Square(double period, double amplitude, double dephasing, double time);
    
    double NoModulation(double time);

    double Triangle(double period, double amplitude, double dephasing,
                 double time);
    double Sawtooth(double period, double amplitude, double dephasing,
                    double time);
    double Ramp(double amplitude, double dephasing, double time);

};

#endif // MODULATOR_H
