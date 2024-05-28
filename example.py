"""
Offset TA Tools
Author: Jonathan Aguilar
Latest update: April 17, 2024

USAGE
-----

This script is provided as an example for importing compute_offsets as a module
into other code. Users should edit the `slew_to`, `slew_from`, `v3`, and
`coron_id` variables using the appropriate values for their observation. The
instructions for editing each section are written in comments in the script.

Once the user has entered the requested values, they should run `example.py` as
a a script. `example.py` will call tools from `compute_offsets.py` and print
out the appropriate X and Y offsets to enter into APT. `example.py` and
`compute_offsets.py` must be located in the same folder.

If the users wish to use the offsets for their own code, instead of just
printing them to screen, they can pass the argument `return_offsets=True` to
compute_offsets.compute_offsets(). It is suggested that users make a new copy
of this script for each variation of acquisition target, science target, and V3
angle.

Requirements:
- numpy
- matplotlib
- astropy
- pySIAF
- compute_offsets.py
"""
from astropy.coordinates import SkyCoord

from miri_coro_offset_ta import compute_offsets


###############################
###### BEGIN USER INPUT #######
###############################

# Star positions - make sure to enter all values.
# The coordinates should be given *at the time of observation*
# By default, we expect the coordinates to be given in RA, Dec system,
# in units of degrees, and in the ICRS frame.
# Users who are comfortable with astropy.coordiantes.SkyCoord may set them
# as they like.

# The "slew_to" variable stores the position of the final target of the observations.
# The "slew_from" variable stores the position of the star that is used for TA.
slew_to = {
    'label': 'B component',
    'position': SkyCoord( 
        272.81221024, 69.24905315,
        unit='deg',
        frame='icrs',
    )
}

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
v3pa = 320.074

# Choose a coronagraph by assigning one of the following to `coron_id`:
# 1065, 1140, 1550, LYOT
coron_id = '1550'

# Plotting - set to False if you don't want to show plots, or True if you do
# If False, other plotting commands are ignored.
# plot_full is a switch for plotting the Imager footprint as well
show_plots = True
plot_full = True

###############################
######## END USER INPUT #######
###############################
# Script takes over from here #
###############################

if __name__ == "__main__":
    dx, dy = compute_offsets.compute_offsets(
        slew_from, slew_to, v3pa, coron_id,
        verbose = False,
        show_plots = show_plots,
        plot_full = plot_full,
        return_offsets = True,
    )

    print("Offsets: ")
    print(f"\tdx: {dx:+4.5f}")
    print(f"\tdy: {dy:+4.5f}")

