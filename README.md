# Offset TA Tools
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
