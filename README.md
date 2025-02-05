# Offset TA Tools
Author: Jonathan Aguilar Last update: March 6, 2024


## Requirements ##
- numpy
- matplotlib
- astropy
- pySIAF


## USAGE ##

This library is meant to help users compute offsets in the detector frame of
reference to enter into APT as Offset Special Requirements, in the case that
they need to perform TA on a target other than their science target. Users
should edit the `if __name__ == '__main__:` section at the bottom of
`compute_offsets.py` with the coordinates of their TA and Science targets, as
well as the position angle of the V3 axis of the telescope. The instructions for
editing each section are written in comments in the script. This will print out
appropriate X and Y offsets that can be entered into APT, along with some
diagnostic information. Alternatively, the file `example.py` shows how to use
`compute_offsets.py` as an imported module.

It is suggested that users make a new copy of `compute_offsets.py` or
`example.py` for each variation of TA star and science target.

It can also be imported as a module after downloading the source code from github:
```
cd $directory
pip install .
```

Then from a python terminal:
```
> from miri_coro_offset_ta import compute_offsets
```

### User-specified variables ###

#### `slew_to`, `slew_from` ####

`slew_to` and `slew_from` are dictionaries used to specify the science (SCI) and
acquisition (ACQ) targets, respectively.
- The `label` keyword is used for annotating output text and figures.
- The `position` keyword is used to store an astropy `SkyCoord` object
  ([SkyCoord
  documentation](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html)).
  This gives considerable flexibility in specifying the coordinates; for
  example, by providing a distance and proper motions, the user can propagate
  the positions using `SkyCoord.apply_space_motion()` to compute offsets for
  multiple epochs.


#### `v3pa` ####

`v3pa` refers to the position angle of the telescope's V3 axis *at the reference
position of the aperture used for the observation*. If that sounds confusing,
the short version is it corresponds to the angle in APT's `Special Requirements
-> PA -> PA Range` menu if the `V3PA` radio button is selected. It also
corresponds to the `ROLL_REF` header keyword (see the [JWST Keyword
Dictionary](https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/) for
keyword definitions). If you are using this library to plan observations in APT,
the `V3PA` field should match your `v3pa` variable.

This is not to be confused with the `PA_APER` header keyword, which corresponds
to the `Aperture PA Range` radio button and refers to the amount by which the
detector-aligned coordinate system is rotated with respect to the `V3` axis. It
also is not to be confused with the `PA_V3` header keyword, which refers to the
V3 position angle at the position of the telescope boresight. Due to spherical
trigonometric effects, the PA of the V3 axis varies across the telescope's focal
plane and varies strongly at high and low latitudes.

Offset slews are specified along the detector axes, in units of arcsec (see
https://jwst-docs.stsci.edu/jppom/special-requirements/general-special-requirements).
In order to convert between the detector coordinate system and two positions on
the sky, pySIAF requires information about the orientation of the telescope.
Here, we provide this information using a combination of the coronagraph used
(see `coron_id`), and position angle of the v3 axis of the telescope, measured
at the chosen coronagraph's reference position.

More details about the different coordinate systems used in describing positions
in the telescope can be found here:
https://jwst-docs.stsci.edu/jwst-observatory-characteristics/jwst-observatory-coordinate-system-and-field-of-regard/
.

#### `coron_id` ####

`coron_id` is a string that tells the script which of the four MIRI coronagraphs
will be used for the observation. Options are:
- '1065' -> '4QPM_1065'/'MIRIM_MASK1065'
- '1140' -> '4QPM_1140'/'MIRIM_MASK1140'
- '1550' -> '4QPM_1550'/'MIRIM_MASK1550'
- 'LYOT' -> '4QPM_LYOT'/'MIRIM_MASKLYOT'

#### `show_plots` ####

This is a switch to turn on (True) or off (False) the display of handy plots
that show the TA process from the points of view of the sky and detector.

#### `other_stars` ###

This parameter allows you to define other targets to plot in the field of view.
It takes a list, each entry of which is a dictionary of the same format as
`slew_from`/`slew_to`.

#### `plot_full` ###

If set to True, this will also plot the MIRI Imager footprint, in addition to
the coronagraph aperture. This is usefull if you have selected `SUBARRAY ->
FULL` in APT.

#### Scripting parameters ####

The function `compute_offsets` takes two more arguments that are more useful if
it is being imported into another script:
- `verbose`: a switch to print (True) or suppress (False) diagnostic text to the
  terminal. This text includes the offsets that should be entered into APT. Set
  to False if you don't want the printed output.
- `return_offsets`: a switch to return (True) or not (False) the numerical value
  of the x and y offset commands. Leave it as False if you intend to copy the
  offsets from the verbose output into APT, or set it to True if you want to
  capture the offsets in your code.



## FAQs for creating your APT program ##

### How do I choose my acquisition target? ###

To choose an acquisition target, you should consider the brightness, separation,
and position angle:
- Brightness: it should be bright enough to achieve high SNR in the TA filter
  without saturating (see the [ETC](https://jwst.etc.stsci.edu/)).
- Separation : it should be closer than the [visit-splitting
  distance](https://jwst-docs.stsci.edu/jwst-astronomers-proposal-tool-overview/additional-jwst-apt-functionality/apt-visit-splitting),
  which ranges between 30"-80" depending on the availability of guide stars for
  a particular target.
- Position angle: The acquisition target must be clear of diffraction spikes
  from nearby sources. This is especially important if the science target is
  very bright. [WebbPSF](https://webbpsf.readthedocs.io/) can be used to
  determine if your TA target will be clear of the diffraction spikes. Determine
  the orientation of your system for a given PA angle using `compute_offsets`,
  which will print the IDL positions of the sources you provide.

### How do I choose dates and V3PA angles? ###

To see available dates and V3PA angles, go to the Visit Planner window in APT
and find the `Reports` menu at the bottom. Select a visit, and then select
`Total Roll Analysis for Visit`. This will give you a plot of available V3 PA
angles against dates, as well as a table that can be read into a script.

### Do I need to calculate a separate offset for each roll? ###

Yes, unless your roll angle is very small or your acquisition target is very
close.
