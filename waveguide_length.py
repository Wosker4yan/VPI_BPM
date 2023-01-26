from sim_setup import VPI_BPM
from masktopolygon import Masktopolygon
import nazca as nd
import numpy as np
from diffractio import um
import matplotlib.pyplot as plt
from tqdm import tqdm

def norm(field, x):

    return sum((field * np.conjugate(field)) * (x[1] - x[0]))

def norm_2(field_1, field_2, x):

    return np.sqrt((sum(field_1 * np.conjugate(field_2))) * (x[1] - x[0]))

def overlap(field_1, field_2, x):

    return sum(field_1 * np.conjugate(field_2) * (x[1] - x[0]))

if __name__ == '__main__':
    width =0.75
    diffractio = Masktopolygon()
    angle =0
    theta = 0
    wavelength = 1.55*um
    clad_index = 1.44
    core_index = 3
    filename = 'mask.png'
    polarization = 'TE'
    ic = diffractio.get_default_ic(width=width)
    wg_length = np.linspace(5,100,25)
    grid_sizes = [3000]
    T= []

    for i,grid in tqdm(enumerate(grid_sizes)):
        T.append([])
        for j,length in tqdm(enumerate(wg_length)):
            with nd.Cell(name=f'Cell{i}{j}') as _cell:
                el = ic.strt(length=length, width=width).put()
                el.raise_pins(['a0', 'b0'], ['a0', 'b0'])

            grid_size = grid

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

            amp_prof_input, amp_prof_output, x0, transmission, E = sim.run_bpm_mode(
                output_monitor_location=length, plotting=False)



            T[-1].append(transmission)

    fig1, ax1 = plt.subplots()
    ax1.set_title('Straight Waveguide Length Sweep')
    ax1.boxplot(T, labels=grid_sizes)
    ax1.set_xlabel('Grid size')
    ax1.set_ylabel('Transmission TE [norm]')
    fig1.tight_layout()
    fig1.show()

    fig2, ax2 = plt.subplots()
    ax2.set_title('Straight waveguide length sweep')
    for i,trans in enumerate(T):
        ax2.plot(wg_length, trans, label = grid_sizes[i])

    ax2.set_xlabel('Waveguide length [um]')
    ax2.legend()
    ax2.set_ylabel('TE Transmission [norm]')
    fig2.tight_layout()
    fig2.show()

    fig, ax = plt.subplots()
    z = ax.contourf(wg_length, grid_sizes, np.log10(T))
    cbar = fig.colorbar(mappable=z, ax=ax)
    cbar.set_label('log$_{10}$(T$_{TM}$) [norm]')
    ax.set_xlabel('Length [$\mu$m]')
    ax.set_ylabel('Grid Size')
    fig.suptitle('Straight Waveguide')
    fig.tight_layout()
    fig.show()

import numpy as np
import matplotlib.pyplot as plt

#c = 0.1205
#a = 0.7705
#b = 0.001

a = 0.0255
b = 0.001
c = 0.126

#a = 0.527
#b = 0.001
#c = 0.677
x = np.linspace(0, 3021, 57)
y = (a * np.exp(b*x)) + c

plt.plot(x, y, '-r')




axes = plt.gca()

plt.xlabel('gratings')
plt.ylabel('FF')
plt.title('FF vs grating number')

plt.show()


for i in y:
    print(i)



