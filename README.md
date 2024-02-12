# Offset TA Tools
Author: Jonathan Aguilar
Latest update: January 25, 2024

USAGE
-----
This library is meant to help users compute offsets in the case that they need to perform TA on a target other than their science target. Users should edit the `if __name__ == '__main__` section at the botto of `compute_offsets.py` with the coordinates of their TA and Science targets, as well as the V3 angle of the telescope. The instructions for editing each section are written in comments in the script. This will print out appropriate X and Y offsets that can be entered into APT, along with some diagnostic information. The file `example.py` shows how to use `compute_offsets.py` as an imported module.

It is suggested that users make a new copy of `compute_offsets.py` for each variation of TA star and science target.

Requirements:
- numpy
- matplotlib
- astropy
- pySIAF
