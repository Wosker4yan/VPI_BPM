from sim_setup import VPI_BPM
from masktopolygon import Masktopolygon
import nazca as nd
from diffractio import um
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


if __name__ == '__main__':
    width2 = 2
    input_wg_length = 5
    pitch_between_taper = 1.8
    length_input_taper = 50
    output_width = 2
    bb_length =20 #16.172
    width_bb = 6.7
    width1 = 0.8
    width =0.5
    output_taper_length = 50
    output_wg_length = 10
    diffractio = Masktopolygon()
    angle =0
    wavelength = 1.55*um
    clad_index = 1.44
    core_index = 1.92
    filename = 'mask.png'
    polarization = 'TE'
    thickness = 0.8
    ic = diffractio.get_default_ic(width=width)
    grid_size = 600 #for wpm keep it small

    T= []
    length_output = input_wg_length + length_input_taper + bb_length
    monitor_loc = length_output + output_wg_length + output_taper_length
    #monitor_loc = 16.172+input_wg_length+length_input_taper
    with nd.Cell(name='mmi') as _cell:

        #Input waveguide
        el1 = ic.strt(length=input_wg_length, width=width1).put()
        el1.raise_pins(['a0'])
        #waveguide to taper
        el2 = ic.taper(length=length_input_taper, width1=width1, width2=width2).put()
        # #big box region
        el3 = ic.strt(length=bb_length, width=width_bb).put()
        #
        # #1 output taper

        el4= ic.taper(length=output_taper_length, width1=width2, width2=width1).put(length_output,1.32)
        el5= ic.strt(length=output_wg_length, width=width1).put()
        el5.raise_pins(['b0'])
        #
        # #2 output taper
        el5= ic.taper(length=output_taper_length, width1=output_width, width2=width1).put(length_output,-1.32)
        el6= ic.strt(length=output_wg_length, width=width1).put()
        # el6.raise_pins(['b0'])

    sim = VPI_BPM(
        grid_size = 8000,
        beam_waist=1,
        core_index=core_index,
        substrate_index=clad_index,
        image=None,
        wavelength=1.55,
        thickness=0.25,
        polarization='TE',
        filename='./mask.png',
        cell=_cell,
        pin='b0',
        mode=0,
        plotting=True,
        sim_mode='BPM',
        angle=0
    )

    sim.cell2image()

    amp_prof_input,  amp_prof_output,  x0, transmission,E  = sim.run_bpm_mode(
        input_pin='a0',
        output_pin = 'b0',
        output_monitor_location=monitor_loc)

    sim.visualize(plotting = True, monitor_location=monitor_loc)

    # prof_input, prof_output, x0 = sim.get_complex_field(
    #     input_pin='a0',
    #     output_pin='b0',
    #     sim_mode = sim_mode,monitor_location=monitor_loc)

import numpy as np
x = np.linspace(500, 3000,6)
ff = np.linspace(0.4,0.65,6)

plt.plot(x, ff)
plt.ylabel("Filling factor")
plt.xlabel("GC teeth")
plt.tight_layout()
plt.show()