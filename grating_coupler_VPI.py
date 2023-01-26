from sim_setup import NazcaBPM
from nazcatopolygon import Nazcatopolygons
import nazca as nd
import numpy as np
from diffractio import um
import matplotlib.pyplot as plt
from tqdm import tqdm
import nazca.geometries as geom

if __name__ == '__main__':
    width =0.45
    grid_size =1000
    diffractio = Nazcatopolygons()
    angle =0
    theta = 0
    wavelength = 1.55*um
    clad_index = 1.44
    core_index = 3
    filename = 'mask.png'
    polarization = 'TE'
    ic = diffractio.get_default_ic(width=width)
    wg_length = 20
    T= []
    number_of_gc = 20
    period = 1


    with nd.Cell(name='mmi') as _cell:
        el0 = ic.strt(length=wg_length+1, width=0.45).put(0,0,0)
        for i in range(number_of_gc):
            el1 = ic.strt(length=0.6, width=0.45).put(0+i*period+1,0.25,0)
            #print(i*period)

        el0.raise_pins(['a0', 'b0'])

    sim = NazcaBPM(
        cell=_cell,
        grid_size=grid_size,
        beam_waist=1,
        thickness=0.22,
        core_index=core_index,
        clad_index=clad_index,
        substrate_index=clad_index,
        left_side_index=clad_index,
        right_side_index=clad_index,
        wavelength=wavelength,
        polarization=polarization,
        mode=0,
        filename='mask.png',
        pin='b0',
        plotting=False,
        sim_mode = 'WPM'
    )

    sim.plot_mask()
    #write a function for transversal profile for taking the mode

    sim.run_bpm_mode()
    sim.visualize()
