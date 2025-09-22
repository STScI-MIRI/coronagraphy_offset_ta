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

    def _get_aper(self):
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
            self._get_aper(),
            self.idl_coords_after_ta,
            ax = axes[0],
            title='Before slew',
        )

        fig = plot_aper_idl(
            self._get_aper(),
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

def print_offset_information(
        slew_from : dict,
        slew_to: dict,
        idl_coords : list,
        sep : units.Quantity,
        pa : units.Quantity,
        offset : np.ndarray,
) -> None:
    """
    Print the information about offsets to the user.

    Parameters
    ----------
    define your parameters

    Output
    ------
    None - prints to screen

    """
    print_output = []
    print_output.append("Computing offset command values to slew from:")
    print_output.append(f"\t{slew_from['label']}")
    print_output.append(f"\t\t RA: \t {slew_from['position'].ra.degree}")
    print_output.append(f"\t\t Dec: \t {slew_from['position'].dec.degree}")
    print_output.append(f"to:")
    print_output.append(f"\t{slew_to['label']}")
    print_output.append(f"\t\t RA: \t {slew_to['position'].ra.degree}")
    print_output.append(f"\t\t Dec: \t {slew_to['position'].dec.degree}")
    print_output.append(f"\n")
    print_output.append(f"Separation and PA: {sep.mas:0.2f} mas, {pa.degree:0.2f} deg\n")

    print_output.append(f"\n")
    print_output.append(f"After TA but before slewing, the position of the ACQ star should be close to (0, 0):")
    print_output.append(f"\t" + ', '.join(f"{i:+0.3e}" for i in idl_coords[0]['position']) + " arcsec")
    print_output.append(f"... and the position of the SCI star is:")
    print_output.append(f"\t" + ', '.join(f"{i:+0.3e}" for i in idl_coords[1]['position']) + " arcsec")

    print_output.append(f"\n")
    print_output.append(f"When the ACQ star is centered, the SCI star is at:")
    print_output.append(f"\tdX: {idl_coords[1]['position'][0]:+2.6f} arcsec")
    print_output.append(f"\tdY: {idl_coords[1]['position'][1]:+2.6f} arcsec")

    if len(idl_coords) > 2:
        print_output.append(f"\n")
        print_output.append(f"Before the slew, here are the positions of the other targets provided:")
        for ic in idl_coords[2:]:
            print_output.append(f"{ic['label']}")
            print_output.append(f"\tdX: {ic['position'][0]:+2.6f} arcsec")
            print_output.append(f"\tdY: {ic['position'][1]:+2.6f} arcsec")

    print_output.append(f"\n")
    print_output.append(f"Sanity check: on-sky angular separation should be the same distance as that of the slew.")
    # print in nice columns
    sep_as_str = f"{sep:0.6f}"
    slew_mag_str = f"{np.linalg.norm(offset) * units.arcsec :0.6f}"
    check_str = [['Separation', sep_as_str], ['Slew magnitude', slew_mag_str]]
    for row in check_str:
        print_output.append("{: >20} {: >20}".format(*row))

    if len(idl_coords) > 2:
        print_output.append(f"\n")
        print_output.append(f"After the slew, here are the positions of other targets provided:")
        for ic in idl_coords[2:]:
            print_output.append(f"{ic['label']}")
            print_output.append(f"\tdX: {ic['position'][0]+offset[0]:+2.6f} arcsec")
            print_output.append(f"\tdY: {ic['position'][1]+offset[1]:+2.6f} arcsec")

    print_output.append("\n")
    print_output.append("Therefore, the commanded offsets that will move the coronagraph from the ACQ star to the SCI are:")
    print_output.append(f"\tdX: {offset[0]:+4.6f} arcsec")
    print_output.append(f"\tdY: {offset[1]:+4.6f} arcsec")
    print_output.append("\n")

    for line in print_output:
        print(line)



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


def plot_detector_ta_sequence(
        aper_dict,
        ta_sequence,
        idl_coords,
        star_positions={},
        axes=None
) -> mpl.figure.Figure:
    """Plot the TA sequence as seen by the detector"""
    if axes is None:
        nrows = 1
        ncols = 4
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*ncols, 4*nrows),
                                 sharex=True, sharey=True, layout='constrained')
    else:
        fig = axes[0].get_figure()

    if star_positions != {}:
        title= f"{star_positions[0]['label']} --> {star_positions[1]['label']}\n"
    else:
        title=''
    fig.suptitle(title + f"TA sequence, as seen by the detector")

    # plot the SIAF apertures on every plot
    for ax in axes:
        aper_dict['mask'].plot(ax=ax, label=False, frame='det', fill=False, c='C0')
        ax.plot([], [], c='C0', label='Readout')
        aper_dict['coro'].plot(ax=ax, label=False, frame='det', mark_ref=True, fill=False, c='C1')
        ax.plot([], [], c='C1', label='Illuminated')
        # TA aperturesL
        for aper_id in ['UR', 'CUR']:
            aper_dict[aper_id].plot(ax=ax, label=False, frame='det', mark_ref=True, fill=False, c='C2')
        ax.plot([], [], c='C2', label='TA regions')

    # plot the positions of the stars at each step in the TA sequence
    # Outer TA
    ax = axes[0]
    ta_aper_id = "UR"
    ax.set_title("Step 1\n" + f"{ta_aper_id} TA region")

    # use the TA aperture object to convert coordinates
    ta_aper = aper_dict[ta_aper_id]
    acq_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][0]['position'])
    sci_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][1]['position'])
    acq_label = ta_sequence[ta_aper_id][0]['label']
    sci_label = ta_sequence[ta_aper_id][1]['label']

    ax.scatter(*acq_pos, 
               c='k', label=acq_label, marker='x', s=100)
    ax.scatter(*sci_pos,
               c='k', label=sci_label, marker='*', s=100)
    other_pos = {i['label']: ta_aper.idl_to_det(*i['position']) for i in ta_sequence[ta_aper_id][2:]}
    for label, pos in other_pos.items():
        ax.scatter(*pos, marker='.', s=50, label=label)

    # put the legend on this plot
    ax.legend(loc='best', ncol=1, fontsize='small', markerscale=0.7)

    # Inner TA
    ax = axes[1]
    ta_aper_id = 'CUR'
    ax.set_title("Step 2\n" + f"{ta_aper_id} TA region")
    # use the TA aperture object to convert coordinates
    ta_aper = aper_dict[ta_aper_id]
    ta_aper.plot(ax=ax, label=False, frame='det', mark_ref=True, fill=False, c='C2')

    acq_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][0]['position'])
    sci_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][1]['position'])
    ax.scatter(*acq_pos,
               c='k', label='ACQ', marker='x', s=100)
    ax.scatter(*sci_pos,
               c='k', label='SCI', marker='*', s=100)
    other_pos = [ta_aper.idl_to_det(*i['position']) for i in ta_sequence[ta_aper_id][2:]]
    for pos in other_pos:
        ax.scatter(*pos, marker='.', s=50)

#     # TA star centered
    ax = axes[2]
    ax.set_title("Step 3\n" + "TA star centered")
    # plot the final TA before the offset is applied
    aper = aper_dict['coro']
    ax.scatter(*aper.idl_to_det(*idl_coords[0]['position']),
               c='k', label='ACQ', marker='x', s=100)
    ax.scatter(*aper.idl_to_det(*idl_coords[1]['position']),
               c='k', label='SCI', marker='*', s=100)
    other_pos = [ta_aper.idl_to_det(*i['position']) for i in ta_sequence[ta_aper_id][2:]]
    for pos in other_pos:
        ax.scatter(*pos, marker='.', s=50)

    # Offset applied
    ax = axes[3]
    ax.set_title("Step 4\n" + "Offset applied")
    # apply the offset to the position
    aper  = aper_dict['coro']
    acq_label = ta_sequence['slew'][0]['label']
    acq_pos = aper.idl_to_det(*ta_sequence['slew'][0]['position'])
    sci_label = ta_sequence['slew'][1]['label']
    sci_pos = aper.idl_to_det(*ta_sequence['slew'][1]['position'])
    ax.scatter(*acq_pos, 
               c='k', label=f'ACQ/{acq_label}', marker='x', s=100)
    ax.scatter(*sci_pos, 
               c='k', label=f'SCI/{sci_label}', marker='*', s=100)
    for pos in ta_sequence['slew'][2:]:
        ax.scatter(*aper.idl_to_det(*pos['position']),
                   marker='.', s=50)

    for ax in axes:
        # plot customizations
        # ax.label_outer()
        ax.set_aspect('equal')
        ax.grid(True, ls='--', c='grey', alpha=0.5)

    return fig


def plot_sky_ta_sequence(aper_dict, star_positions, v3pa, offset, colors, axes=None):

    if axes is None:
        nrows = 1
        ncols = 4
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                 layout='constrained',
                                 sharex=True, sharey=True)
    else:
        fig = axes[0].get_figure()

    targ_label= f"{star_positions[0]['label']} --> {star_positions[1]['label']}"
    fig.suptitle(targ_label)

    for ax in axes.ravel():
        acq_pos = (star_positions[0]['position'].ra.deg, star_positions[0]['position'].dec.deg)
        sci_pos = (star_positions[1]['position'].ra.deg, star_positions[1]['position'].dec.deg)
        other_pos = [(star['position'].ra.deg, star['position'].dec.deg) for star in star_positions[2:]]
        ax.scatter(*acq_pos,
                   c='k', label=f"ACQ/{star_positions[0]['label']}", marker='x', s=100)
        ax.scatter(*sci_pos,
                   c='k', label=f"SCI/{star_positions[1]['label']}", marker='*', s=100)
        for pos in other_pos:
            ax.scatter(*pos, marker='.', s=50)


    # We start TA in the outer TA region
    ax = axes[0]
    ax.set_title(f"Step 1\nUR TA region")

    # center the attitude matrix at the Outer TA ROI
    attmat = create_attmat(star_positions[0]['position'], aper_dict['UR'], v3pa)
    formatting = dict(c=colors[0], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)


    # Continue to step 2 of TA, in the inner TA region
    ax = axes[1]
    ax.set_title(f"Step 2\nCUR TA region")

    # center the attitude matrix at the Inner TA ROI
    attmat = create_attmat(star_positions[0]['position'], aper_dict['CUR'], v3pa)
    formatting = dict(c=colors[1], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    # plot the final TA before the offset is applied
    ax = axes[2]
    ax.set_title("Step 3\nCentered")

    # center the attitude matrix on the coronagraph reference position
    attmat = create_attmat(star_positions[0]['position'], aper_dict['coro'], v3pa)
    formatting = dict(c=colors[2], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    # Plot the apertures and sources after the offset slew
    ax = axes[3]
    # compute the ra, dec of the offset from the TA position
    attmat = create_attmat(star_positions[0]['position'], aper_dict['coro'], v3pa)
    aper_dict['coro'].set_attitude_matrix(attmat)
    ra, dec = aper_dict['coro'].idl_to_sky(*(-offset))
    tel_sky = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    # compute the new attitude matrix at the slew position
    attmat = create_attmat(tel_sky, aper_dict['coro'], v3pa)
    ax.set_title(f"Step 4\nOffset applied")#\nTel-Targ sep: {tel_sky.separation(star_positions['SCI']).to('mas'):0.2e}")
    formatting = dict(c=colors[3], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    for ax in axes:
        # plot customizations
        ax.set_ylabel("Dec [deg]")
        ax.set_xlabel("RA [deg]")
        # ax.label_outer()
        ax.set_aspect("equal") 
        ax.grid(True, ls='--', c='grey', alpha=0.5)    
        # fix x-axis labels
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    # RA increases to the left
    # they share axes, so you just need to invert one of them
    axes.flat[0].invert_xaxis()

    return fig



def make_plots(
        aper_dict : dict,
        star_positions : list,
        v3pa : float,
        idl_coords : list,
        offset : np.ndarray,
) -> None:
    """
    Generate the plots

    Parameters
    ----------
    aper_dict : dict,
      a dict of the SIAF apertures to plot
    star_positions : dict,
      dictionary of the sky coordinates of the ACQ and SCI targets
    idl_coords : list,
      list of dicts of the idl coordinates of each star in the list, at each phase of the TA sequence
    offset : np.ndarray,
      the amount of the slew if you want to draw it I guess
    Output
    ------
    Displays 4 plots

    """

    figures = []

    fig1, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), layout='constrained')
    fig1 = plot_before_offset_slew(aper_dict, idl_coords, star_positions, ax=axes[0])
    fig1 = plot_after_offset_slew(aper_dict, idl_coords, offset, star_positions, ax=axes[1])
    figures.append(fig1)

    # Plot 2: The TA sequence on the detector
    # compute the positions in IDL coordinates at each step of the sequence
    ta_sequence = {}
    for aper_id in ['UR', 'CUR', 'coro']:
        aper = aper_dict[aper_id]
        ta_sequence[aper_id] = sky_to_idl(star_positions, 
                                          aper, 
                                          v3pa)
    ta_sequence['slew'] = [{'label': i['label'], 'position': np.array(i['position']) + offset} for i in ta_sequence['coro']]

    # Plot 2: The TA sequence in RA and Dec on a single plot
    colors = mpl.cm.plasma(np.linspace(0.2, 0.9, 4))
    fig2 = plot_sky_ta_sequence_one_axis(aper_dict, star_positions, v3pa, offset, colors)
    figures.append(fig2)

    # Plot 3: plot detector and sky POV on same figure
    fig3 = plot_observing_sequence(
        aper_dict, ta_sequence, v3pa, idl_coords,
        star_positions, offset
    )
    figures.append(fig3)

    # # now, show the plots in correct order
    for fig in figures[::-1]:
        fig.show()
    plt.show()



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

