// The station module

#include <cstdlib>
#include <stdexcept>

#include "Station.h"

class AtmosphereLayer;

Station::Station()
{
    this->staname = "None";
    this->north = 0;
    this->east = 0;
    this->up = 0;
    this->gain = 0;
    this->diameter = 0;  
}

Station::Station(string station_name, double North, double East, double Up, double gain, double diameter)
{
    this->staname = station_name;
    this->north = North;
    this->east = East;
    this->up = Up;
    this->layer = NULL;
    this->gain = gain;
    this->diameter = diameter;
}

Station::Station(string station_name, double North, double East, double Up, AtmosphereLayer * atm, double gain, double diameter)
{
    this->staname = station_name;
    this->north = North;
    this->east = East;
    this->up = Up;
    this->layer = atm;
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


