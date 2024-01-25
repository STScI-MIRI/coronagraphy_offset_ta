"""
Offset TA Tools
Author: Jonathan Aguilar
Latest update: January 25, 2024

USAGE
-----

Users should enter the coordinates the TA target and Science target in `compute_offsets.py`, and run it as a script.
`compute_offsets.py` will call tools from this file and print out the appriate X and Y offsets to enter into APT.
`compute_offsets.py` and `offset_tools.py` must be located in the same folder.

It is suggested that users make a new copy of `compute_offsets.py` for each variation of TA star and science target.

Requirements:
- numpy
- pandas
- matplotlib
- astropy
- pySIAF
"""
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

from astropy.coordinates import SkyCoord
from astropy import units

import pysiaf
from pysiaf import Siaf


def sky_to_idl(ta_pos, targ_pos, aper, pa):
    """
    Convert RA and Dec positions of a TA star and its target with an offset
    into a detector position (measured from the reference point, in arcsec)
    for a given PA of the V3 axis w.r.t. North.
    Assume the TA star is centered on the aperture (idl = (0, 0))

    Parameters
    ----------
    ta_pos : SkyCoord object of the TA star position
    targ_pos : SkyCoord object of the target star position
    aper : SIAF object for the aperture being used (e.g. MIRIM_CORON1550)
    pa : the PA in degrees of the V3 axis of the telescope (measured eastward of North) at the observation epoch

    Output
    ------
    idl_coords : dict {'ta': tuple, 'targ': tuple}
      a dictionary of IDL x, y coordinates for the TA star and the science target
      the TA star coordinates should be very close to 0
    """
    v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=ta_pos.ra.deg, 
                                                    dec=ta_pos.dec.deg, 
                                                    pa=pa)
    aper.set_attitude_matrix(attmat)
    idl_coords = {}
    # ta star - should be close to 0
    idl_coords['ta'] = np.array(aper.sky_to_idl(ta_pos.ra.deg, ta_pos.dec.deg))
    # eps Mus
    idl_coords['targ'] = np.array(aper.sky_to_idl(targ_pos.ra.deg, targ_pos.dec.deg))
    return idl_coords


def plot_apers(ax, attmat, aper_dict, format_dict):
    """
    Helper function to plot the apertures for a given part of the TA sequence
    ax : axis to plot on
    attmat : attitude matrix
    aper_dict : dictionary of apertures to plot
    format_dict : aperture plot formatting parameters
    """
    for k, aper in aper_dict.items():
        aper.set_attitude_matrix(attmat)
        formatting = format_dict.copy()
        if k == 'mask':
            formatting['mark_ref'] = True
        if k == 'coro':
            # skip the illuminated region aperture, it's too crowded
            pass
        else:
            aper.plot(ax=ax, label=False, frame='sky', fill=False, **formatting)


def plot_before_offset_slew(aper_dict, idl_coords):
    """Plot the scene on the detector when you're pointed at the TA target"""
    # plot 1 : POV of the detector
    fig, ax = plt.subplots(1, 1, )
    ax.set_title("""Positions *before* offset slew""")
    frame = 'idl' # options are: tel (telescope), det (detector), sci (aperture)
    aper_dict['mask'].plot(ax=ax, label=False, frame=frame, c='C0')
    aper_dict['coro'].plot(ax=ax, label=False, frame=frame, c='C1', mark_ref=True)

    ax.scatter(0, 0, c='k', label='TACQ', marker='x', s=100)
    ax.scatter(*idl_coords['targ'], label="Target", marker="*", c='k')

    ax.legend()
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    return fig

def plot_detector_ta_sequence(aper_dict, ta_sequence, idl_coords):
    """Plot the TA sequence as seen by the detector"""
    nrows = 1
    ncols = 4
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*ncols, 4*nrows), sharex=True, sharey=True)
    fig.suptitle(f"TA sequence, as seen by the detector")

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
    ta_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id]['ta'])
    targ_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id]['targ'])

    ax.scatter(*ta_pos, 
               c='k', label='TA star', marker='x', s=100)
    ax.scatter(*targ_pos,
               c='k', label='Target', marker='*', s=100)

    # put the legend on this plot
    ax.legend(loc='best', ncol=1, fontsize='small', markerscale=0.7)

    # Inner TA
    ax = axes[1]
    ta_aper_id = 'CUR'
    ax.set_title("Step 2\n" + f"{ta_aper_id} TA region")
    # use the TA aperture object to convert coordinates
    ta_aper = aper_dict[ta_aper_id]
    ta_aper.plot(ax=ax, label=False, frame='det', mark_ref=True, fill=False, c='C2')

    ta_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id]['ta'])
    targ_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id]['targ'])

    ax.scatter(*ta_pos,
               c='k', label='TA star', marker='x', s=100)
    ax.scatter(*targ_pos,
               c='k', label='Target', marker='*', s=100)

    # TA star centered
    ax = axes[2]
    ax.set_title("Step 3\n" + "TA star centered")
    # plot the final TA before the offset is applied
    aper = aper_dict['coro']
    ax.scatter(*aper.idl_to_det(*idl_coords['ta']),
               c='k', label='TA star', marker='x', s=100)
    ax.scatter(*aper.idl_to_det(*idl_coords['targ']),
               c='k', label='Target', marker='*', s=100)

    # Offset applied
    ax = axes[3]
    ax.set_title("Step 4\n" + "Offset applied")
    # apply the offset to the position
#     ta_pos  = aper_dict['coro'].idl_to_det(*np.array(ta_sequence['UR']['ta']) + offset)
#     targ_pos = aper_dict['coro'].idl_to_det(*np.array(ta_sequence['UR']['targ']) + offset)
    aper  = aper_dict['coro']
    ta_pos = aper.idl_to_det(*ta_sequence['slew']['ta'])
    targ_pos = aper.idl_to_det(*ta_sequence['slew']['targ'])
    ax.scatter(*ta_pos, 
               c='k', label='TA star', marker='x', s=100)
    ax.scatter(*targ_pos, 
               c='k', label='Target', marker='*', s=100)


    for ax in axes:
        # plot customizations
        ax.label_outer()
        ax.set_aspect('equal')
        ax.grid(True, ls='--', c='grey', alpha=0.5)

    fig.tight_layout()
    return fig

def plot_sky_ta_sequence(aper_dict, star_positions, offset, colors):

    nrows = 4
    ncols = 1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             #figsize=(6*ncols, 10*nrows),
                             sharex=True, sharey=True)

    for ax in axes.ravel():
        ta_pos = (star_positions['TACQ'].ra.deg, star_positions['TACQ'].dec.deg)
        targ_pos = (star_positions['Target'].ra.deg, star_positions['Target'].dec.deg)    
        ax.scatter(*ta_pos,
                   c='k', label='TACQ', marker='x', s=100)
        ax.scatter(*targ_pos,
                   c='k', label='Target', marker='*', s=100)    



    # We start TA in the outer TA region
    ax = axes[0]
    ax.set_title(f"Step 1: UR TA region")

    # center the attitude matrix at the Outer TA ROI
    aper = aper_dict['UR']
    # the telescope is now pointing the *outer* TA region at the TA star
    v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=star_positions['TACQ'].ra.deg, 
                                                    dec=star_positions['TACQ'].dec.deg, 
                                                    pa=star_positions['v3'])
    formatting = dict(c=colors[0], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)


    # Continue to step 2 of TA, in the inner TA region
    ax = axes[1]
    ax.set_title(f"Step 2: CUR TA region")

    # center the attitude matrix at the Inner TA ROI
    aper = aper_dict['CUR']
    # the telescope is now pointing the *inner* TA region at the TA star
    v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing the inner TA region at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=star_positions['TACQ'].ra.deg, 
                                                    dec=star_positions['TACQ'].dec.deg, 
                                                    pa=star_positions['v3'])
    formatting = dict(c=colors[1], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    # plot the final TA before the offset is applied
    ax = axes[2]
    ax.set_title("Step 3: Centered")

    # center the attitude matrix on the coronagraph reference position
    aper = aper_dict['coro']
    # the telescope is now pointing the center of the coronagraph at the TA star
    v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=star_positions['TACQ'].ra.deg, 
                                                    dec=star_positions['TACQ'].dec.deg, 
                                                    pa=star_positions['v3'])
    formatting = dict(c=colors[2], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    # Plot the apertures and sources after the offset slew
    ax = axes[3]
    aper = aper_dict['coro']
    # note that the two methods of computing the slew below are equivalent
    # v2, v3 = aper.idl_to_tel(*(offset))
    # ra = star_positions['TACQ'].ra.degree
    # dec = star_positions['TACQ'].dec.degree
    v2, v3 = aper.reference_point('tel')
    # compute the RA and Dec of the pointing using the offset values you found earlier
    # note that you must CHANGE THE SIGN OF THE SLEW with respect to the previous plot
    # if you're right, you should end up right on the position of the star
    ra, dec = aper_dict['coro'].idl_to_sky(*(-offset))
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=ra, 
                                                    dec=dec, 
                                                    pa=star_positions['v3'])
    tel_sky = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    ax.set_title(f"Step 4: Offset applied\nTel-Targ sep: {tel_sky.separation(star_positions['Target']).to('mas'):0.2e}")
    formatting = dict(c=colors[3], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    for ax in axes:
        # plot customizations
        ax.set_ylabel("Dec [deg]")
        ax.set_xlabel("RA [deg]")
        ax.label_outer()
        ax.set_aspect("equal") 
        ax.grid(True, ls='--', c='grey', alpha=0.5)    
        # fix x-axis labels
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    fig.tight_layout()
    return fig

def plot_sky_ta_sequence_one_axis(aper_dict, star_positions, offset, colors):
    """Plot the TA sequence on the sky, all on one axis"""
    nrows = 1
    ncols = 1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 5))
    fig.suptitle(f"TA sequence, in RA/Dec")
    colors = mpl.cm.plasma(np.linspace(0.2, 0.9, 4))

    ta_pos = (star_positions['TACQ'].ra.deg, star_positions['TACQ'].dec.deg)
    targ_pos = (star_positions['Target'].ra.deg, star_positions['Target'].dec.deg)    
    ax.scatter(*ta_pos,
               c='k', label='TACQ', marker='x', s=100)
    ax.scatter(*targ_pos,
               c='k', label='Target', marker='*', s=100)    

    # We start TA in the outer TA region
    aper = aper_dict["UR"]
    # the telescope is now pointing the *outer* TA region at the TA star
    v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=star_positions['TACQ'].ra.deg, 
                                                    dec=star_positions['TACQ'].dec.deg, 
                                                    pa=star_positions['v3'])
    formatting = dict(c=colors[0], alpha=1, ls='dotted')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [], 
            **formatting,
            label='Step 1: Outer TA step')


    # Continue to step 2 of TA, in the inner TA region
    aper = aper_dict["CUR"]
    # the telescope is now pointing the *outer* TA region at the TA star
    v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=star_positions['TACQ'].ra.deg, 
                                                    dec=star_positions['TACQ'].dec.deg, 
                                                    pa=star_positions['v3'])
    formatting = dict(c=colors[1], alpha=1, ls='dashdot')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [], 
            **formatting,
            label='Step 2: Inner TA step')    

    # plot the final TA before the offset is applied

    # the telescope is now pointing the center of the coronagraph at the TA star
    aper = aper_dict['coro']
    v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3, 
                                                    ra=star_positions['TACQ'].ra.deg, 
                                                    dec=star_positions['TACQ'].dec.deg, 
                                                    pa=star_positions['v3'])
    formatting = dict(c=colors[2], alpha=1, ls='dashed')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [],
            **formatting,
            label='Step 3: Before offset')


    # the telescope now places the TA star at the commanded offset
    aper = aper_dict['coro']
    v2, v3 = aper.reference_point('tel')
    # note that you must CHANGE THE SIGN OF THE OFFSET to get the position of the reference point
    ra, dec = aper.idl_to_sky(*(-offset))
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3,
                                                    ra=ra,
                                                    dec=dec,
                                                    pa=star_positions['v3'])
    formatting = dict(c=colors[3], alpha=1, ls='solid')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [],
            **formatting,
            label='Step 4: After offset')


    # plot formatting
    ax.set_ylabel("Dec [deg]")
    ax.set_xlabel("RA [deg]")
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    # fix x-axis labels
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    fig.tight_layout()
    ax.legend()

    return fig


def make_plots(
        aper_dict : dict,
        star_positions : dict,
        idl_coords : dict,
        offset : np.ndarray,
) -> None:
    """
    Generate the plots

    Parameters
    ----------
    idl_coords : dictionary of the idl coordinates of each component at each phase of the TA sequence

    Output
    ------
    Displays 4 plots

    """
    colors = mpl.cm.plasma(np.linspace(0.2, 0.9, 4))

    fig1 = plot_before_offset_slew(aper_dict, idl_coords)

    # Plot 2: The TA sequence on the detector
    # compute the positions at each step of the sequence
    ta_sequence = {}
    for aper_id in ['UR', 'CUR', 'coro']:
        aper = aper_dict[aper_id]
        ta_sequence[aper_id] = sky_to_idl(star_positions['TACQ'], 
                                          star_positions['Target'], 
                                          aper, 
                                          star_positions['v3'])
    ta_sequence['slew'] = {k: np.array(v) + offset for k, v in ta_sequence['coro'].items()}
    fig2 = plot_detector_ta_sequence(aper_dict, ta_sequence, idl_coords)


    # Plot 3: The TA sequence in RA and Dec, split into separate plots
    fig3 = plot_sky_ta_sequence(aper_dict, star_positions, offset, colors)

    # Plot 4: The TA sequence in RA and Dec on a single plot
    fig4 = plot_sky_ta_sequence_one_axis(aper_dict, star_positions, offset, colors)

    # now, actually show the plots
    plt.show()


def compute_offsets(
        slew_from: dict,
        slew_to: dict,
        v3: float,
        coron_ids : list[str],
        show_plots : bool = True,
) -> np.ndarray :
    """
    Compute the slews for the TA sequences, print the offsets, and show the plots if requested

    Parameters
    ----------
    slew_from: dict
    slew_to: dict
    v3: float
    coron_ids : list[str]
    show_plots : bool = True

    Output
    ------
    Prints offsets and shows plots

    """

    coron_id = coron_ids[0]
    star_positions = {
        # the TA star
        'TACQ': slew_from['position'],
        # The star you will eventually slew to
        'Target': slew_to['position'],
        'v3': v3
    }

    # Offsets
    sep = star_positions['TACQ'].separation(star_positions['Target']).to(units.arcsec)
    pa = star_positions['TACQ'].position_angle(star_positions['Target']).to(units.deg)
    print("Separation and PA: ", f"{sep.mas:0.2f} mas, {pa.degree:0.2f} deg\n")


    # Siaf
    miri = Siaf("MIRI")
    # now that we have the MIRI object, let's get the 1550 coronagraph apertures used in 1618.
    # There are two relevant apertures: MIRIM_MASK[XXXX], which is the entire subarray, and
    # MIRIM_CORON[XXXX], which is just the portion that gets illuminated
    # let's combine all the SIAF objects in a dict for convenience

    all_apers = {}
    all_apers['UR'] = miri[f'MIRIM_TA{coron_id}_UR']
    all_apers['CUR'] = miri[f'MIRIM_TA{coron_id}_CUR']
    all_apers['coro'] = miri[f'MIRIM_CORON{coron_id}']
    all_apers['mask'] = miri[f'MIRIM_MASK{coron_id}']


    idl_coords = sky_to_idl(star_positions['TACQ'], 
                            star_positions['Target'],
                            all_apers['coro'],
                            star_positions['v3'])
    # The offset you apply is as if you were moving the science target - i.e.
    # the negative of its position
    offset = -1*np.array(idl_coords['targ'])

    print("Computing offset command values from:")
    print(f"\t{slew_from['label']}")
    print("\t\t RA:\t ", slew_from['position'].ra.degree)
    print("\t\t Dec:\t", slew_from['position'].dec.degree)
    print("to")
    print(f"\t{slew_to['label']}")
    print("\t\t RA:\t ", slew_to['position'].ra.degree)
    print("\t\t Dec:\t", slew_to['position'].dec.degree)

    print("\n")
    print("After TA but before slewing, the position of the TACQ star should be close to (0, 0):")
    print(f"\t", ', '.join(f"{i:+0.3e}" for i in idl_coords['ta']), "arcsec")
    print("... and the position of the Target star is:")
    print(f"\t", ', '.join(f"{i:+0.3e}" for i in idl_coords['targ']), "arcsec")
    print("")

    print("When the TACQ star is centered, the Target star is at:")
    print(f"\tdX: {idl_coords['targ'][0]:+2.6f} arcsec")
    print(f"\tdY: {idl_coords['targ'][1]:+2.6f} arcsec")
    print("\n")

    print("Therefore, the commanded offsets that will move the coronagraph from the TACQ star to the Target are:")
    print(f"\tdX: {offset[0]:+2.6f} arcsec")
    print(f"\tdY: {offset[1]:+2.6f} arcsec")

    print("\n")
    print("Sanity check: on-sky angular separation should be the same distance as the slew.")
    # print in nice columns
    sep_as_str = f"{sep:0.6f}"
    slew_mag_str = f"{np.linalg.norm(offset) * units.arcsec :0.6f}"
    check_str = [['Separation', sep_as_str], ['Slew magnitude', slew_mag_str]]
    for row in check_str:
        print("{: >20} {: >20}".format(*row))


    if show_plots == True:
        make_plots(
            all_apers,
            star_positions,
            idl_coords,
            offset,
        )
        
# if __name__ == "__main__":
#     compute_offsets(slew_from)
#     if show_plots == True:
#         make_plots()