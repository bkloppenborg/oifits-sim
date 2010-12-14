/// \class PinholdCombiner simulator.h
/// \brief A class to implement a pinhole combiner
class PinholeCombiner
{
  public:
    double normalisation;

    double amplitudes;

    void Couple(Beam & Beam1);

    PinholeCombiner(int screenSize, double diameter);
};
