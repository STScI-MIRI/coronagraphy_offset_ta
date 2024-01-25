"""
Offset TA Tools
Author: Jonathan Aguilar
Latest update: January 25, 2024

USAGE
-----
This library is meant to help users compute offsets in the case that they need to perform TA on a target other than their science target. Users should edit the file `compute_offsets.py` with the coordinates of their TA and Science targets, as well as the V3 angle of the telescope. The instructions for editing each section are written in comments in the script. 

Once the user has entered the requested values, they should run `compute_offsets.py` as a a script. `compute_offsets.py` will call tools from `offset_tools.py` and print out the appropriate X and Y offsets to enter into APT. `compute_offsets.py` and `offset_tools.py` must be located in the same folder.

It is suggested that users make a new copy of `compute_offsets.py` for each variation of TA star and science target.

Requirements:
- numpy
- matplotlib
- astropy
- pySIAF
"""
from astropy.coordinates import SkyCoord

import offset_tools


###############################
########## USER INPUT #########
###############################

# Star positions - make sure to enter all values.
# The coordinates should be given *at the time of observation*
# By default, we expect the coordinates to be given in RA, Dec system,
# in units of degrees, and in the ICRS frame.
# Users who are comfortable with astropy.coordiantes.SkyCoord may set them
# as they like.

# The "slew_to" variable stores the position of the final target of the observations
slew_to = {
    'label': 'B component',
    'position': SkyCoord( 
        272.81221024, 69.24905315,
        unit='deg',
        frame='icrs',
    )
}
# The "slew_from" variable stores the position of the star that is used for TA.
slew_from = {
    'label': 'A component',
    'position': SkyCoord( 
        272.81369648, 69.25014279,
        unit='deg',
        frame='icrs',
    )
}
# Telescope V3PA
# enter the PA angle of the *telescope* V3 axis, at the time of the observation
v3 = 320.074

# Choose a coronagraph by uncommenting one of these choices
coron_id = [
    # '1065',
    # '1140',
    '1550',
    # 'LYOT',
]

# Plotting - set to False if you don't want to show plots
show_plots = True

###############################
####### END USER INPUT ########
###############################
# Script takes over from here #
###############################

if __name__ == "__main__":
    offset_tools.compute_offsets(slew_from, slew_to, v3, coron_id, show_plots)

