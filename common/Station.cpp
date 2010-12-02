// The station module

#include <cstdlib>
#include <stdexcept>

#include "Station.h"

class AtmosphereLayer;

Station::Station(string station_name, double North, double East, double Up, AtmosphereLayer * atm, double gain, double diameter)
{
    this->staname = station_name;
    this->north = North;
    this->east = East;
    this->up = Up;
    this->gain = gain;
    this->diameter = diameter;
}

Station::~Station()
{
    // free(layer);
}

string Station::GetName()
{
    return this->staname;
}


