"""
Miscellaneous plot methods
"""
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

from all_jwst_offset_ta import utils

def plot_wfss_traces(
    idl_coords : dict,
    aper : utils.pysiaf.JwstAperture,
    ax : mpl.axes.Axes = None,
    title : str = '',
    show_mirim_illum = True,
    plot_full : bool = False
):
    """
    Compute the idl coordinates overlaid on the aperture of interest. Add the WFSS traces as well.
    idl_coords : dict[label, tuple[ra, dec]]
      usually, one of utils.ComputeOffsets().idl_coords_after_ta or
      utils.ComputeOffsets().idl_coords_after_slew
    aper : usually the return value of utils.ComputeOffsets().get_aper()
    ax : axis to plot on
    title : title for the axis
    show_aper_title : if True, display the text MIRIM_ILLUM
    plot_full : if True, also plot the FULL detector array
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    ax.set_title(title)
    trace_up = 100 * aper.YSciScale
    trace_down = 300 * aper.YSciScale
    # trace_lims = {k: (v[1]-trace_down, v[1]+trace_up) for k, v in idl_coords.items()}
    
    for i, (k, coord) in enumerate(idl_coords.items()):
        width = 1
        height = trace_up + trace_down
        ll = (coord[0]-width/2, coord[1]-trace_down)
        trace = mpl.patches.Rectangle(ll, width, height, facecolor=f'C0', alpha=0.5)
        ax.add_patch(trace)
        ax.scatter(*coord, c=f'C{i}', marker='x',label=k)
    aper.plot(ax=ax, fill=False, mark_ref=False, frame='idl', label=True, zorder=-1)

    if plot_full:
        full_aper = utils.Siaf("MIRI")['MIRIM_FULL']
        # conver the corners to ILLUM IDL
        corners = aper.convert(*full_aper.corners("tel"), from_frame="tel", to_frame='idl')
        full_rect = mpl.patches.Rectangle(
            xy = np.min(corners, axis=1),
            width = np.diff(corners[0]).max(),
            height = np.diff(corners[1]).max(),
            fill=False, ec='gray'
        )
        ax.add_patch(full_rect)
    ax.legend(loc=(1.05, 0.3), title='Sources')
    return fig

def corners_to_other_frame(from_aper, to_aper, frame):
    new_corners = to_aper.convert(*from_aper.corners("tel"), from_frame="tel", to_frame='idl')
    return new_corners
