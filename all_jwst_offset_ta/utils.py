import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon

from astropy.coordinates import SkyCoord
from astropy import units

import pysiaf
from pysiaf import Siaf


#------------------------------------------------------#
#----------------- Offset Computation -----------------#
#------------------------------------------------------#
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon

import ipywidgets as widgets

from astropy.coordinates import SkyCoord
from astropy import units

import pysiaf
from pysiaf import Siaf

# define user interface

instruments = ['NIRCAM', 'NIRSPEC', 'NIRISS', 'MIRI', 'FGS']

class ComputeOffsets():
    def __init__(self, initial_values={}):
        """
        Generate the Offset Computation GUI.
        Can initialize from a dictionary containing the entries:
        initial_values={
          'pa': 180.,
          'acq_ra' = 90., acq_dec = 90., sci_ra = 91., sci_dec=89.}
        """
        self._initial_values = initial_values
        self._current_values = {k: v for k, v in initial_values.items()}
        self.ui = self._make_ui()
        # if an initial dictionary is provided, run the computations.
        if initial_values != {}:
            self._compute_offsets()

    def _update_apers(self, *args, initial_values={}):
        self._acq_apers = [i for i in Siaf(self._instr_picker.value).apernames if '_TA' in i]
        self._sci_apers = [i for i in Siaf(self._instr_picker.value).apernames if '_TA' not in i]
        self._acq_aper_picker.options = self._acq_apers
        self._sci_aper_picker.options = self._sci_apers

        try:
            self._acq_aper_picker.value = initial_values.get('acq_aper', self._acq_apers[0]).upper()
        except IndexError:
            pass
        try:
            self._sci_aper_picker.value = initial_values.get('sci_aper', self._sci_apers[0]).upper()
        except IndexError:
            pass

    def _update_config_dict(self):
        self._current_values['instr'] = self._instr_picker.value
        self._current_values['sci_aper'] = self._sci_aper_picker.value
        self._current_values['pa'] = self._PA_setter.value
        self._current_values['acq_ra'] = self._acq_pos_widget.children[1].children[0].value
        self._current_values['acq_dec'] = self._acq_pos_widget.children[1].children[1].value
        self._current_values['sci_ra'] = self._sci_pos_widget.children[1].children[0].value
        self._current_values['sci_dec'] = self._sci_pos_widget.children[1].children[1].value
        self._current_values['other_stars'] = self._other_stars_widget.value

    def get_aper(self):
        aper = Siaf(self._instr_picker.value)[self._sci_aper_picker.value]
        return aper

    def _make_starpos_widget(self, title, initial_ra=0., initial_dec=0.):
        """Make a widget for getting a star's RA and Dec"""
        star_widget = widgets.VBox([
            widgets.Label(value=title),
            widgets.VBox([
                widgets.FloatText(value=initial_ra, description='RA [deg]', disabled=False),
                widgets.FloatText(value=initial_dec, description='Dec [deg]', disabled=False)
            ])
        ])
        return star_widget

    def _make_final_idl_widget(self, title, initial_x=0., initial_y=0.):
        """Make a widget for getting a star's RA and Dec"""
        star_widget = widgets.VBox([
            widgets.Label(value=title),
            widgets.VBox([
                widgets.FloatText(value=initial_x, description='Final IDL X', disabled=False),
                widgets.FloatText(value=initial_y, description='Final IDL Y', disabled=False)
            ])
        ])
        return star_widget


    def _parse_other_stars(self) -> list:
        """Parse the entries in the other_stars_widget and return a list"""
        other_stars = []
        if self._other_stars_widget.value != '':
            stars = self._other_stars_widget.value.split("\n")
            stars = [dict(zip(['label', 'position'], i.split(":"))) for i in stars]
            for star in stars:
                position = SkyCoord(
                    *[float(i) for i in star['position'].strip()[1:-1].split(",")],
                    frame='icrs', unit='deg',
                )
                star['position'] = position
            other_stars = stars
        return other_stars

    def _make_widgets(self, initial_values={}):
        """
        Container method for making and initializing the widgets
        """
        self._instr_picker = widgets.Dropdown(
            options = instruments,
            value = initial_values.get('instr', instruments[0]).upper(),
            description='Instrument'
        )
        self._instr_picker.observe(self._update_apers)
        # set self.acq_apers and self.sci_apers
        self._acq_aper_picker = widgets.Dropdown(description='ACQ aperture')
        self._sci_aper_picker = widgets.Dropdown(description='SCI aperture')
        self._update_apers(initial_values=initial_values)
        # Position Angle
        self._PA_setter = widgets.BoundedFloatText(
            value=initial_values.get("pa", 0),
            min=0.,
            max=360.,
            step=0.1,
            description='PA (deg):',
            disabled=False
        )
        self._slew_to_this_idl = self._make_final_idl_widget(
            "Slew to this IDL",
            initial_x = initial_values.get("final_idl_x", 0.),
            initial_y = initial_values.get("final_idl_y", 0.),
        )
        # star positions
        self._acq_pos_widget = self._make_starpos_widget(
            "ACQ target position",
            initial_ra = initial_values.get("acq_ra", 0.),
            initial_dec = initial_values.get("acq_dec", 0.),
        )
        self._sci_pos_widget = self._make_starpos_widget(
            "SCI target position",
            initial_ra = initial_values.get("sci_ra", 0.),
            initial_dec = initial_values.get("sci_dec", 0.)
        )
        self._other_stars_widget = widgets.Textarea(
            value=initial_values.get("other_stars", ''),
            placeholder='label : (ra.deg, dec.deg)',
            description='Other stars: ',
            disabled=False,
        )

        self._compute_offsets_button = widgets.Button(
            description='Compute offset',
            disabled=False,
            button_style='success'
        )
        self._compute_offsets_button.on_click(self._compute_offsets)
        self._plot_scene_button = widgets.Button(
            description = 'Plot scenes',
            disabled=False,
            button_style = 'info'
        )
        self._plot_scene_button.on_click(self._plot_scene)

        self._output_offset = widgets.Output()
        self._output_before = widgets.Output()
        self._output_after = widgets.Output()

    def _make_ui(self):
        self._make_widgets(self._initial_values)
        grid = widgets.GridspecLayout(
            n_rows=9, n_columns=3,
            style=dict(background='white')
        )
        grid[0, :] = widgets.Label(
            value="IDL Coordinate and Offset TA Calculator".upper(),
            layout = widgets.Layout(display='flex', justify_content='center'),
        )
        grid[1, 0] = self._instr_picker
        grid[2, 0] = self._sci_aper_picker
        grid[4, 0] = self._PA_setter
        grid[5:8, 0] = self._slew_to_this_idl

        grid[1:4, 1] = self._acq_pos_widget
        grid[4:7, 1] = self._sci_pos_widget

        grid[1:-1, 2] = self._other_stars_widget

        # place the buttons at the bottom of the central column
        grid[-1, 1] = widgets.HBox([self._compute_offsets_button, self._plot_scene_button])
        output_grid = widgets.GridspecLayout(
            n_rows=1, n_columns=3,
        )
        output_grid[0, 0] = self._output_offset
        output_grid[0, 1] = self._output_before
        output_grid[0, 2] = self._output_after
        ui = widgets.VBox([
            grid, output_grid
            # widgets.HBox([, self.output_before, self.output_after]),
        ])
        return ui

    def _plot_scene(self, *args):
        fig,axes = plt.subplots(nrows=1, ncols=2)
        fig = plot_aper_idl(
            self.get_aper(),
            self.idl_coords_after_ta,
            ax = axes[0],
            title='Before slew',
        )

        fig = plot_aper_idl(
            self.get_aper(),
            self.idl_coords_after_slew,
            ax = axes[1],
            title='After slew',
        )
        return fig

    def show_ui(self):
        return self.ui

    def _compute_offsets(self, *args):
        acq_pos = {
            'label': 'ACQ',
            'position': SkyCoord(
                self._acq_pos_widget.children[1].children[0].value,
                self._acq_pos_widget.children[1].children[1].value,
                frame='icrs', unit='deg',
            ),
        }
        sci_pos = {
            'label': 'SCI',
            'position': SkyCoord(
                self._sci_pos_widget.children[1].children[0].value,
                self._sci_pos_widget.children[1].children[1].value,
                frame='icrs', unit='deg',
            )
        }

        v3pa = self._PA_setter.value
        aperture = Siaf(self._instr_picker.value)[self._sci_aper_picker.value]

        other_stars = self._parse_other_stars()
        slew_to_idl = np.array([
            self._slew_to_this_idl.children[1].children[0].value,
            self._slew_to_this_idl.children[1].children[1].value,
        ])

        idl_coords = compute_idl_after_ta(
            acq_pos, sci_pos, v3pa, aperture,
            other_stars = other_stars,
        )
        self.idl_coords_after_ta = {i['label']: i['position'] for i in idl_coords}
        self.offset_to_sci = -self.idl_coords_after_ta['SCI'] + slew_to_idl
        self.idl_coords_after_slew = {
            k: v + self.offset_to_sci
            for k, v in self.idl_coords_after_ta.items()
        }

        self._output_offset.clear_output()
        self._output_before.clear_output()
        self._output_after.clear_output()
        outputstr = "Special Requirement -> Offset values\n"
        outputstr += "-"*(len(outputstr)-1) + "\n"
        outputstr += f"Offset X [arcsec]: {self.offset_to_sci[0]:+0.4f}\nOffset Y [arcsec]: {self.offset_to_sci[1]:+0.4f}"
        self._output_offset.append_stdout(outputstr)
        outputstr = f"IDL positions of stars after TA:\n"
        outputstr += "-"*(len(outputstr)-1) + "\n"
        outputstr += "\n".join(f"{k}\t{v[0]:>+10.4f}\t{v[1]:>+10.4f}" for k, v in self.idl_coords_after_ta.items())
        self._output_before.append_stdout(outputstr)
        outputstr = f"IDL positions of stars after slew:\n"
        outputstr += "-"*(len(outputstr)-1) + "\n"
        outputstr += "\n".join(f"{k}\t{v[0]:>+10.4f}\t{v[1]:>+10.4f}" for k, v in self.idl_coords_after_slew.items())
        self._output_after.append_stdout(outputstr)


#------------------------------------------------------#
#----------------- Offset Computation -----------------#
#------------------------------------------------------#
def compute_idl_after_ta(
    slew_from: dict,
    slew_to: dict,
    v3pa: float,
    sci_aper : str,
    other_stars : list = [],
) -> list[dict]:
    """
    Compute the IDL positions of all the given stars whent he ACQ target is at the reference position.

    Parameters
    ----------
    slew_from: dict
      A dictionary containing the label and position of the TA target, set by
      the user in compute_offsets.py
    slew_to: dict
      A dictionary containing the label and position of the science target, set by
      the user in compute_offsets.py
    v3pa: float
      The PA_V3 angle of the telescope for this observation
    coron_id : str
      Identifier for the desired coronagraphic subarray. Must be one of '1065',
      '1140', '1550', and 'LYOT'
    other_stars : list
      A list of dicts of other stars in the field, in the same format as slew_from/slew_to
    verbose : int = 1
      print diagnostics and offsets to screen.
      0 : nothing is printed
      1 : all output is printed
      2 : only the final IDl coordinates of the stars are printed

    Output
    ------
    Returns a list of dicts of floats

    """
    # make sure coron_id is valid
    star_positions = [slew_from, slew_to] + other_stars
    idl_coords = sky_to_idl(star_positions,
                            sci_aper,
                            v3pa)
    return idl_coords

def apply_offset(
    idl_coords : dict,
) -> dict:
    sci_idl = idl_coords.pop('SCI')
    offset_to_sci = -1 * sci_idl
    final_idl = {l: c + offset_to_sci for l, c in idl_coords.items()}
    return final_idl

def compute_offsets(
        slew_from: dict,
        slew_to: dict,
        v3pa: float,
        sci_aper : str,
        other_stars : list = [],
        verbose : int = 1,
        return_offsets : bool = False,
) -> np.ndarray :
    """
    Compute the slews for the TA sequences, print the offsets, and show the plots if requested.
    How it works:
    - Point the coronagraphic aperture (MASK or CORON) at the TA star by
      setting an attitude matrix for the V3PA value
    - The offset is the negative of the IDL coordinates of the SCI star
    - The rest of the machinery is basically just making verification plots

    Parameters
    ----------
    slew_from: dict
      A dictionary containing the label and position of the TA target, set by
      the user in compute_offsets.py
    slew_to: dict
      A dictionary containing the label and position of the science target, set by
      the user in compute_offsets.py
    v3pa: float
      The PA_V3 angle of the telescope for this observation
    coron_id : str
      Identifier for the desired coronagraphic subarray. Must be one of '1065',
      '1140', '1550', and 'LYOT'
    other_stars : list
      A list of dicts of other stars in the field, in the same format as slew_from/slew_to
    verbose : int = 1
      print diagnostics and offsets to screen.
      0 : nothing is printed
      1 : all output is printed
      2 : only the final IDl coordinates of the stars are printed
    show_plots : bool = True
      If True, display the diagnostic plots. If False, only print the offsets.
    return_offsets : bool = False
      If True, return an array of dx and dy offsets

    Output
    ------
    Prints offsets and shows plots. Returns a dict of floats

    """
    # make sure coron_id is valid
    star_positions = [slew_from, slew_to] + other_stars
    labels = {'ACQ': slew_from['label'],
              'SCI': slew_to['label']}


    # Offsets to science target
    sep = star_positions[0]['position'].separation(star_positions[1]['position']).to(units.arcsec)
    pa = star_positions[0]['position'].position_angle(star_positions[1]['position']).to(units.deg)

    idl_coords = sky_to_idl(star_positions,
                            sci_aper,
                            v3pa)
    # The offset you apply is as if you were moving the science target - i.e.
    # the negative of its position
    offset = -1*np.array(idl_coords[1]['position'])

    if verbose == 1:
        print_offset_information(
            slew_from=slew_from,
            slew_to=slew_to,
            idl_coords=idl_coords,
            sep=sep,
            pa=pa,
            offset=offset
        )
    if verbose == 2:
        len_label = max(len(star['label']) for star in idl_coords)
        print("IDL coordinates of all stars after slew:")
        for star in idl_coords:
            label = star['label']
            idl = star['position']
            print(f"{label:{len_label}s}:\t{idl[0]+offset[0]:+0.10f}, {idl[1]+offset[1]:+0.10f}")
        print("")

    # if show_plots == True:
    #     make_plots(
    #         all_apers,
    #         star_positions,
    #         v3pa,
    #         idl_coords,
    #         offset,
    #     )

    if return_offsets == True:
        return offset

def create_attmat(
        position : SkyCoord,
        aper : pysiaf.aperture.JwstAperture,
        pa : float,
        idl_offset : tuple[float, float] = (0., 0.),
        set_matrix : bool = False
) -> np.ndarray:
    """
    Create an attitude matrix for JWST when the reference point of a particular
    aperture is pointed at a given position for a specified PA

    Parameters
    ----------
    position : SkyCoord
      skycoord position on the sky
    aper : pysiaf.aperture.JwstAperture
      pySIAF-defined aperture object
    pa : float
      PA angle with respect to the V3 axis measured at the aperture reference
      point. This corresponds to the V3PA field in the APT PA range special
      requirement, and the ROLL_REF keyword in the data. This is *not* the
      PA_V3 keyword value, which is the PA angle of the V3 axis measured at the
      telescope boresight.
    idl_offset : tuple[float, float] = (0.0, 0.0)
      allows you to specify an arbitrary position in IDL coordinates that
      corresponds to the position
    set_matrix : bool = True
      if True, also set the matrix on the current aperture in addition to returning it
    Output
    ------
    attmat : np.ndarray
      matrix that pySIAF can use to specify the attitude of the telescope
    """
    v2, v3 = aper.idl_to_tel(idl_offset[0], idl_offset[1])
    # v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3,
                                                    ra=position.ra.deg,
                                                    dec=position.dec.deg,
                                                    pa=pa)
    if set_matrix == True:
        aper.set_attitude_matrix(attmat)
    return attmat


def sky_to_idl(
        stars : list[dict],
        aper : pysiaf.aperture.JwstAperture,
        pa : float,
        idl_offset : tuple[float, float] = (0., 0.)
) -> list[dict]:
    """
    Convert RA and Dec positions of a TA star and its target with an offset
    into a detector position (measured from the reference point, in arcsec)
    for a given PA of the V3 axis w.r.t. North.
    Assume the TA star is centered on the aperture (idl = (0, 0))

    Parameters
    ----------
    stars : list of {"label": label, "position": pos} elements to plot
      the first element is the ACQ target. The aperture is centerd on this target.
    aper : SIAF object for the aperture being used (e.g. MIRIM_MASK1550)
    pa : the PA in degrees of the V3 axis of the telescope (measured eastward of North) at the observation epoch

    Output
    ------
    idl_coords : dict {label: idl_tuple} of IDL coordinates for each target
      a dictionary of IDL x, y coordinates for the TA star and the science target
      the TA star coordinates should be very close to 0
    """
    acq_pos = stars[0]['position']
    attmat = create_attmat(acq_pos, aper, pa, idl_offset)
    aper.set_attitude_matrix(attmat)
    idl_coords = []
    # ta star - should be close to 0
    for star in stars:
        label = star['label']
        pos = star['position']
        idl_coords.append({'label': label, 'position': np.array(aper.sky_to_idl(pos.ra.deg, pos.dec.deg))})
    return idl_coords


def work_backwards(
        slew_from : dict,
        slew_to : dict,
        coron_id : str,
        v3pa : float,
        offset : tuple[float, float] = (0., 0.),
        other_stars : list = [],
) -> list[dict] :
    """
    Work backwards to find out where the slew_to star ended up, given sky
    coordinates for the targets, the commanded offset, and the v3pa value of
    the observation.

    Parameters
    ----------
    slew_from : dict
      label and coordinate of the ACQ target
    slew_to : dict
      label and coordinate of the SCI target
    coron_id : str
      Identifier for the desired coronagraphic subarray. Must be one of '1065',
      '1140', '1550', and 'LYOT'
    v3pa : float
      the v3pa angle of the telescope at the aperture reference position, in
      degrees
    offset : np.ndarray[float]
      the x and y offset commanded
    other_stars : list
      A list of dicts of other stars in the field that you might want to keep
      track of, in the same format as slew_from/slew_to

    Output
    ------
    idl_positions : list[dict]
      A list of positions in subarray IDL coordinates provided targets. Each
      list entry has format {'label': label, 'position': position}.

    """
    coron_id = coron_id.upper()
    coro = miri[f'MIRIM_CORON{coron_id}']
    mask = miri[f'MIRIM_MASK{coron_id}']

    star_positions = [slew_from, slew_to] + other_stars

    idl_coords = sky_to_idl(star_positions,
                            coro,
                            v3pa,
                            idl_offset=offset)

    return idl_coords

#--------------------------------------------#
#----------------- Plotting -----------------#
#--------------------------------------------#

def plot_before_offset_slew(
        aper_dict : dict,
        idl_coords : list,
        star_positions : list =[],
        ax = None,
):
    """Plot the scene on the detector when you're pointed at the acquisition target"""
    # plot 1 : POV of the detector
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout='constrained')
    else:
        fig = ax.get_figure()
    if star_positions != []:
        title = f"{star_positions[0]['label']} --> {star_positions[1]['label']}"
        fig.suptitle(title)
    ax.set_title("""Positions *before* offset slew""")
    frame = 'idl' # options are: tel (telescope), det (detector), sci (aperture)
    aper_dict['mask'].plot(ax=ax, label=False, frame=frame, c='C0')
    aper_dict['coro'].plot(ax=ax, label=False, frame=frame, c='C1', mark_ref=True)

    ax.scatter(0, 0,
               c='k',
               label=f"ACQ/{idl_coords[0]['label']}",
               marker='x',
               s=100)
    ax.scatter(*idl_coords[1]['position'],
               label=f"SCI/{idl_coords[1]['label']}",
               marker="*",
               c='k')
    for star in idl_coords[2:]:
        ax.scatter(*star['position'],
                   # c='k',
                   label=star['label'],
                   marker='.',
                   s=50)
    ax.add_artist(quad_boundaries(aper_dict['coro'], kwargs={'fc': 'grey'}))
    ax.legend()
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    return fig


def plot_aper_idl(
        aper : pysiaf.aperture.JwstAperture,
        star_positions : dict[str, np.ndarray] = {},
        offset : list | np.ndarray = [0., 0.],
        ax = None,
        title = '',
):
    """Plot the scene on the detector when you're pointed at the science target"""
    # plot 1 : POV of the detector
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout='constrained')
    else:
        fig = ax.get_figure()
    frame = 'idl' # options are: tel (telescope), det (detector), sci (aperture)
    if title != '':
        ax.set_title(title)
    aper.plot(ax=ax, label=False, frame=frame, c='k', mark_ref=True)

    offset = np.array(offset)
    star_positions = star_positions.copy()
    acq_pos = star_positions.pop('ACQ')
    ax.scatter(acq_pos[0] + offset[0], acq_pos[1] + offset[1],
               c='k',
               label=f"ACQ",
               marker='x',
               s=100)
    sci_pos = star_positions.pop("SCI")
    ax.scatter(*(sci_pos+ offset),
               label=f"SCI",
               marker="*",
               c='k')
    for star, position in star_positions.items():
        ax.scatter(*(position + offset),
                   # c='k',
                   label=star,
                   marker='.',
                   s=50)
    if aper.AperName[-4:] in ['1065', '1140', '1550']:
        ax.add_artist(quad_boundaries(aper, kwargs={'fc': 'grey'}))
    ax.legend()
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    return fig



def quad_boundaries(aperture, kwargs={}):
    """
    Generate a polygon to plot the 4QPM quadrant boundaries. Stolen from the JWST
    coronagraphic visibility tool

    Parameters
    ----------
    aperture: a pysiaf.Siaf aperture for the 1065, 1140, or 1550 coronagraph
    kwargs : {} arguments to pass to Polygon

    Output
    ------
    mask : matplotlib.patches.Polygon object
    """

    y_angle = np.deg2rad(aperture.V3IdlYAngle)
    corners_x, corners_y = aperture.corners(to_frame='idl')
    min_x, min_y = np.min(corners_x), np.min(corners_y)
    max_x, max_y = np.max(corners_x), np.max(corners_y)

    width_arcsec = 0.33
    x_verts0 = np.array([
        min_x,
        -width_arcsec,
        -width_arcsec,
        width_arcsec,
        width_arcsec,
        max_x,
        max_x,
        width_arcsec,
        width_arcsec,
        -width_arcsec,
        -width_arcsec,
        min_x
    ])
    y_verts0 = np.array([
        width_arcsec,
        width_arcsec,
        max_y,
        max_y,
        width_arcsec,
        width_arcsec,
        -width_arcsec,
        -width_arcsec,
        min_y,
        min_y,
        -width_arcsec,
        -width_arcsec
    ])
    x_verts = np.cos(y_angle) * x_verts0 + np.sin(y_angle) * y_verts0
    y_verts = -np.sin(y_angle) * x_verts0 + np.cos(y_angle) * y_verts0

    verts = np.concatenate([x_verts[:, np.newaxis], y_verts[:, np.newaxis]], axis=1)
    mask = Polygon(verts, alpha=0.5, **kwargs)
    return mask

