oifits-sim
==========

OIFITS Simulator

(c) Brian Kloppenborg, Fabien Baron (2015)

# Description: The OIFITS simulator is a tool to assist in observation
planning or in reconstructed image artifact analysis.  The program has
the ability to either (1) simulate the data from an interferometer
given a list of observational hour angles and a FITS image of what the
source looks like, or (2) copy the UV coverage of an existing
observation and simulate the data given a FITS image of the source.

OIFITS-SIM is developped by Brian Kloppenborg. Originally based on the
MROI simulator written by Fabien Baron, it does not include all its
capabilities yet (though it brings new ones).

This repository contains Fabien Baron's own version of OIFITS-SIM, developped for PFI.

# Installation instructions

To check out and build the current stable version of this software do:

```
git clone https://github.com/fabienbaron/oifits-sim.git
cd oifits-sim
git submodule init
git submodule update
cd build
cmake ..
make
```

# Using the simulator

The simulator operates in two primary modes:

## Using with existing data

In this mode the program will simulate observations from the interferometer and COPY the
uncertainty estimates from an existing data file.  This mode is useful if you are exploring
artifacts in existing images as you can sample a model image and reconstruct the resulting 
OIFITS file.  To use this mode the user needs to specify the target, target image, output 
file, input OIFITS file, array, combiner, and spectral mode of the combiner.  In future versions
the last three parameters will be determined (in as much as possible) from the OIFITS file
but in the present version this wasn't possible.

    oifits-sim -t target_file.txt -i image_file.fits -o output.oifits -d existing.oifits -a array -c combiner -m mode
    
The simulator attempts to find information about the array, combiner and spectral mode from
the OIFITS file.  If it fails, it will prompt you for more information.

For example (running from the bin directory), we can emulate CHARA-MIRC's observations of the 2010-10 model image of epsilon Aurigae via:

    oifits-sim -t ../samples/target_epsAur.txt -i ../samples/2010-10.fits -o 2010-10_simulated.oifits -d ../samples/2010-10-eps_Aur-all-avg5.oifits -a CHARA -c MIRC -m LowH

## Simulating potential observations (for observation planning / proposals)

In this mode the program will simulate observations from an interferometer and estimate the
noise (see documentation for further description of noise estimation).  This mode is useful
during observing planning as the user can choose various array configurations and see the
sampling that results from the configuration.  The user may, of course, reconstruct images
from the simulated data to see which configuration results in "better" images.

To use this mode the user must specify the target, target image, output file, array, 
combiner, spectral mode, and observations:

    oifits-sim -t target_file.txt -i image_file.fits -o output.oifits -a array -c combiner -m spectral_mode -obs observations
      
Note: see the shortcut section below as arrays, combiners, and spectral modes may be
abbreviated.

Example of Hour Angle usage:

    oifits-sim -t ../samples/target_epsAur.txt -i ../samples/2010-10.fits -o 2010-10_simulated.oifits -a CHARA -obs ../samples/obs_ha_example.txt -c MIRC -m LowH

# Shortcuts

For the definitions of arrays, combiners, and spectral modes the software will automatically
search the files in oifits-sim/etc for the presence of a similarly named file.  Parameters in
these files can be overwidden (although they may not necessairly be used) on the command line
via -parameter_name=value.

Currently implemented arrays:

    CHARA

Currently implemented combiners and spectral modes

    MIRC	: Low_H

# Overriding (some) parameters:

Some (but not all) of the parameters in default configuration files can be manually overridden 
on the command line.  To override one of these parameters specify the options IMMEDIATELY after
the inclusion of the configuration file.  For example the following command overrides
the default value of r0 at CHARA (7.6 cm) and sets it to 3.1 cm:

    oifits-sim -a CHARA --r0 0.031 -t foo.txt ...
    ...
    r0 overridden from default value in array config file.  Now: 0.031
      
Conversely the following will do nothing as the --r0 parameter will be ignored.

    oifits-sim -a CHARA -t foo.txt ... --r0 0.031

All parameters that can be overriden are mentioned below.  See documentation for more details.

    Array:
        --r0 double
        --wind_speed double

    Target:
        --res double
        --mag double
        
# Bug Reporting and Feature Requests

Please use the [issue tracker on GitHub](https://github.com/fabienbaron/oifits-sim/issues)


# Licensing and acknowledgements

`oifits-sim` is free software, distributed under the GNU Lesser General Public License (Version 3).

If you use this software as part of a scientific publication, please cite the following work:

Kloppenborg, B.; Baron, F. (2012) "oifits-sim: an OIFITS-simulator" (Version 2). Available from https://github.com/fabiebaron/oifits-sim.
