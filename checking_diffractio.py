from sim_setup import VPI_BPM
from masktopolygon import Masktopolygon
import nazca as nd
from diffractio import degrees, plt, um, np
import matplotlib.pyplot as plt
import math
from diffractio.scalar_masks_XZ import Scalar_mask_XZ
from diffractio.scalar_sources_X import Scalar_source_X
def norm(field, x):
    """Normalization of a field.

    Normalizing the field
    Args:
        field (list or array): Input field.
        x (int): Grid difference (delta) x0 parameter.
    Returns:
        normalized field (float)
    """

    return np.sqrt((sum(field * np.conjugate(field)))) * (x[1] - x[0])


def normalize(vec):
    """Normalizes an electric field from its maximum, use for visualization of the field inside the index profile
    Args:
        vec (vector): Field to be normalized
    Returns:
        norm (Array): Normalized field
    """

    vec = np.asarray(vec)
    norm = np.squeeze(vec / vec.max())

    return norm

def overlap(field_1, field_2, x):
    """Overlap between two fields.

    Args:
        field_1 (list or array): Input field 1.
        field_2 (list or array): Input field 2.
        x (list or array): Grid difference (delta) x0 parameter.
    Returns:
        overlap between two fields (float)
    """
    return sum(field_1 * np.conjugate(field_2) * (x[1] - x[0]))



if __name__ == '__main__':
    width = 0.45
    diffractio = Masktopolygon()
    angle = 0
    theta = 0
    wavelength = 1.55 * um
    substrate_index = 1.44
    core_index = 3
    filename = 'mask.png'
    polarization = 'TE'
    ic = diffractio.get_default_ic(width=width)
    wg_length = 50
    T = []

    with nd.Cell(name='cell ') as _cell:
        el = ic.strt(length=wg_length, width=width).put(0, 0)
        el.raise_pins(['a0', 'b0'], ['a0', 'b0'])

    sim = VPI_BPM(
        grid_size = 7000,
        beam_waist=1,
        core_index=core_index,
        substrate_index=substrate_index,
        image=None,
        wavelength=1.55,
        polarization='TE',
        filename='./mask.png',
        cell=_cell,
        pin='b0',
        mode=0,
        plotting=True,
        sim_mode='BPM',
        angle=0
    )


    amp_prof_input, amp_prof_output, x0, transmission = sim.run_bpm_mode(sim_mode='BPM',
        output_monitor_location=10, angle=0, plotting=True)



    sim.visualize(sim_mode = 'BPM',monitor_location=0)
    x, z = sim.cell2image(save=True)
    width = sim.get_ic_width()