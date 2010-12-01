// The station module

#include <cstdlib>
#include <stdexcept>

#include "Station.h"
#include "Simulator.h"

class AtmosphereLayer;

Station::Station(double coordx, double coordy, double coordz, AtmosphereLayer * atm,
                 double gain, double diameter)
{
    // :TODO: this->staname
    this->x = coordx;
    this->y = coordy;
    this->z = coordz;
    this->layer = atm;
    this->gain = gain;
    this->diameter = diameter;
}

// void Station::Update( double current_time )
// {
// layer->Update(current_time); // updates the atmosphere
// }

Station::~Station()
{
    // free(layer);
}


