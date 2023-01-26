from sim_setup import VPI_BPM
from masktopolygon import Masktopolygon
import nazca as nd
import numpy as np
from diffractio import um
import matplotlib.pyplot as plt
from tqdm import tqdm
import nazca.geometries as geom

if __name__ == '__main__':
    width1 =0.45
    grid_size =1000
    diffractio = Masktopolygon()
    angle =0
    theta = 0
    wavelength = 1.55*um
    clad_index = 1.44
    core_index = 3
    filename = 'mask.png'
    polarization = 'TE'
    ic = diffractio.get_default_ic(width=width1)


    T = []
    number_of_gc = 25
    period = 0.2
    ff = 0.5
    duty_cycle = ff * period
    wg_length = period*number_of_gc  + period
    with nd.Cell(name='grating') as _cell:

        for i in range(number_of_gc):
            el1 = ic.strt(length=duty_cycle, width=0.45).put(0+i*period,0.45,0)

            #print(i*period)
        el = ic.strt(length=wg_length , width=0.45).put(0, 0, 0)
        el.raise_pins(['a0', 'b0'])






    sim = VPI_BPM(
        grid_size=grid_size,
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

    sim.plot_mask()
    #write a function for transversal profile for taking the mode

    sim.run_bpm_mode(output_monitor_location= wg_length, input_monitor_location = 0)
    sim.visualize(monitor_location = wg_length)






