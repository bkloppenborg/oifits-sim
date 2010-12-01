
/// \class Fiber simulator.h
/// \brief A class to represent a single mono-mode optical fiber.
/// \details
/// Implements a single mono-mode optical fiber coupling light from free
/// space.  Optimum Gaussian mode has 1/e radius 0.9 times pupil radius
class Fiber
{
  public:
    Matrix < double >mode;

    Matrix < double >opdmode;

    Matrix < double >mask;

    Fiber(int screenSize, double diameter, double modeSize);

    Complex Couple(Beam & Beam1);
    double Opd(Beam & Beam1);
};
