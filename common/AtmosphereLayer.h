/// \file AtmosphereLayer.h
/// Header file for the AtmosphereLayer class.

#include <fftw3.h>

#include "Matrix.h"
#include "Simulator.h"

/// \class AtmosphereLayer simulator.h
/// \brief A class to represent an atmospheric layer.
class AtmosphereLayer
{
  public:
    AtmosphereLayer(int type, int generatorsize, int screensize,
                    double lambda0, double r0, double pixellation,
                    double windspeed, double winddirection, double outerscale,
                    double scintillation_diameter);
                    
    virtual ~ AtmosphereLayer();
                    
  public:
                    
    double Cn2(double);

    void InitScreen();

    void Update(double time);

    Matrix < double >phase_generator;   // frozen screen turbulence

    Matrix < double >phase_screen;      // current square phase screen
    // representing the turbulence just
    // above the telescope
    Matrix < double >amplitude_generator;       // frozen amplitude screen
    // (scintillation)
    Matrix < double >amplitude_screen;  // current square scintillation screen 
                                        // 
    // 
    // just above the telescope
    int type;                   // 0: Kolmogorov 1: Von Karman 2:
    // Greenwood-Tarazano
    int generatorsize;

    int screensize;             // Actual size of screen as determined by FFT
    // constraints
    double GetDelay(double time);

    double lambda0;             // wavelength for which the screen is computed 
                                // 
    // 
    // (for reference)
    double r0;                  // r0(lambda0) in meters

    double outerscale;          // von Karman outer scale

    double pixellation;         // number of pixels per meter

    double wind_direction;      // direction of the Wind in degrees

    double wind_speed;

    double previous_time;       // store the previous call time

    double scintillation_diameter;

};

