/// \file DelayLine.h
/// Header file for the delay line class.


/// \class DelayLine DelayLine.h
/// \brief A class representing a delay line.
class DelayLine
{
  public:
    DelayLine(int type, double valueinit[3], double filtercoeffs[3], double time_resolution);
    int filter_type;            // 1st, 2nd order filter or instantaneous
    // delay line
    double coeffs[3];           // filter coefficients

    double init[3];             // initial conditions

    double t_current;

    double x_current;

    double v_current;

    double dt;

    double GetDelay(double time);

    double GetCurrentPosition(double time);
    
    void ComputePinhole(int size, int padded_size);

    void SetTargetPosition(double target_position);

  private:
    double target_position;

    double previous_calltime;

    double previous_position;

    double previous_speed;

    double previous_targetposition;

    double delay_equation_1st(double t, double x);

    double delay_equation_2nd(double t, double x, double v);
};
