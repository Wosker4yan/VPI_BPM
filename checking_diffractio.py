from sim_setup import VPI_BPM
from masktopolygon import Masktopolygon
import nazca as nd
from diffractio import degrees, plt, um, np

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
    width = 0.77
    diffractio = Masktopolygon()
    angle = 0
    theta = 0
    wavelength = 1.55 * um
    substrate_index = 1.44
    core_index = 2.5
    filename = 'mask.png'
    polarization = 'TM'
    ic = diffractio.get_default_ic(width=width)
    wg_length = 50
    T = []

    with nd.Cell(name='cell ') as _cell:
        el = ic.strt(length=wg_length, width=width).put(0, 0)
        el.raise_pins(['a0', 'b0'], ['a0', 'b0'])

    sim = VPI_BPM(
        grid_size = 4200,
        beam_waist=1,
        core_index=core_index,
        substrate_index=substrate_index,
        image=None,
        wavelength=1.55,
        thickness=0.28,
        polarization='TE',
        filename='./mask.png',
        cell=_cell,
        pin='b0',
        mode=0,
        plotting=True,
        sim_mode='BPM',
        angle=0
    )

    amp_prof_input, amp_prof_output, x0, transmission, E = sim.run_bpm_mode(
        output_monitor_location=18, angle=0, plotting=True)



    sim.visualize(sim_mode = 'BPM',monitor_location=25)
    x, z = sim.cell2image(save=True)
    width = sim.get_ic_width()

    sim.index_contour_plot()
    x, eff, E = sim.mode_and_neff()
    effective_index = sim.get_slab_index()
    plt.plot(E, abs(x))
    plt.show()