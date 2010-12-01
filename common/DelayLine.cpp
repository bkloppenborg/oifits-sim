//============================================================================
// Name        : delaylines.cpp
// Version     : 1.0
// Copyright ( C ) 2008-2009 Fabien Baron
//============================================================================

#include "DelayLine.h"

// Delay lines module
// First and second order differential equation are solved 
// using the classic fourth-order Runge-Kutta method

DelayLine::DelayLine(int type, double valueinit[3], double filtercoeffs[3],
    double time_resolution)
{
  // Initialization of the delay line module
  // int type : 1st order or second order filter

  for (int i = 0; i <= 2; i++)
    this->init[i] = valueinit[i];
  for (int i = 0; i <= 2; i++)
    this->coeffs[i] = filtercoeffs[i];

  this->dt = time_resolution;
  this->filter_type = type;
  // Dummy initialize
  this->target_position = 0.0;
  this->previous_calltime = 0.0;
  this->previous_position = 0.0;
  this->previous_speed = 0.0;
  this->previous_targetposition = 0.0;
}

double DelayLine::delay_equation_1st(double t, double x)
{
  return -(coeffs[1] * x + coeffs[0] * t);
}

double DelayLine::delay_equation_2nd(double t, double x, double v)
{
  return -(coeffs[2] * v + coeffs[1] * x + coeffs[0] * t);
}

void DelayLine::SetTargetPosition(double new_target_position)
{
  // Set the new target
  previous_targetposition = target_position;
  target_position += new_target_position;
  // cout << "New target position: " << target_position << endl;
  //Change the equations (initial conditions don't change as the position is continuous)
  coeffs[0] = coeffs[0] - target_position;
}

double DelayLine::GetCurrentPosition(double current_calltime)
{
  double k1, k2, k3, k4;
  int intervals = int((current_calltime - previous_calltime) / dt);
  // Note: the intrinsic time resolution (here = dt) of the delay lines should be small enough
  // compared to the integration time so that intervals is an integer
  // (needs to be improved...)

  switch (filter_type)
  {
    case 0: // instantaneous delayline (would be good to have one of those !!!)
      x_current = target_position;
      break;

    case 1: // Update the position using Runge-Kutta for 1st order filter
      for (int ii = 0; ii < intervals; ii++)
      {
        k1 = dt * delay_equation_1st(t_current, x_current);
        k2 = dt * delay_equation_1st(t_current + dt / 2.0, x_current + k1
            / 2.0);
        k3 = dt * delay_equation_1st(t_current + dt / 2.0, x_current + k2
            / 2.0);
        k4 = dt * delay_equation_1st(t_current + dt, x_current + k3);
        x_current += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        t_current += dt;
      }
      break;

    case 2: // Update the position using Runge-Kutta for 2nd order filter
      for (int ii = 0; ii < intervals; ii++)
      {
        k1 = dt * delay_equation_2nd(t_current, x_current, v_current);
        k2 = dt * delay_equation_2nd(t_current + dt / 2., x_current
            + v_current * dt / 2., v_current + 0.5 * k1);
        k3 = dt * delay_equation_2nd(t_current + dt / 2., x_current
            + v_current * dt / 2. + 0.25 * k1 * dt, v_current + 0.5 * k2);
        k4 = dt * delay_equation_2nd(t_current + dt, x_current + v_current
            * dt + dt / 2. * k2, v_current + k3);
        x_current += (v_current + (k1 + k2 + k3) / 6.) * dt;
        v_current += (k1 + k2 + k2 + k3 + k3 + k4) / 6.;
        t_current += dt;
      }

  }

  //if (filter_type !=0) cout << "Time: " << current_calltime << "\t DL time: " << t_current << "\t DL pos: " << x_current <<endl;
  return x_current;
}






