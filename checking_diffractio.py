from sim_setup import VPI_BPM
from masktopolygon import Masktopolygon
import nazca as nd
from diffractio import degrees, plt, um, np
import matplotlib.pyplot as plt
import math
from diffractio.scalar_masks_XZ import Scalar_mask_XZ
from diffractio.scalar_sources_X import Scalar_source_X


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
    wg_length = 20
    T = []

    with nd.Cell(name='cell ') as _cell:
        el = ic.strt(length=wg_length, width=width).put(0, 0)
        el.raise_pins(['a0', 'b0'], ['a0', 'b0'])

    sim = VPI_BPM(
        grid_size = 5000,
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


    amp_prof_input, amp_prof_output, x0, transmission = sim.run_bpm_mode(
        output_monitor_location=wg_length, angle=0, plotting=True)
