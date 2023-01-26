from vpipda.material import *
from vpipda.layout import *
from vpipda import units
import numpy as np
import matplotlib.pyplot as plt
from diffractio import degrees, plt, um,mm
from diffractio.scalar_masks_XZ import Scalar_mask_XZ
from diffractio.scalar_sources_X import Scalar_source_X
import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks
from devices import *

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


width = 0.75

sub_mat = Dielectric(1.45, color='cyan')
core_mat = Dielectric(3.4, color='teal')
core = Layer(core_mat, h=width, rot='90 deg')
wg2 = CrossSection(sub_mat, core)
wg2.plot()
plt.show()

#TODO
#LATER
grid_size = 5000

# for m in wg2.mode.values():
#     m.E.real.plot(contour=False, aspect='equal')
# plt.show()
#now when I have te and tm it's time to pass them to the software and see what happens good job
wg2.mesh.box = wg2.box.expand(5)
wg2.mesh.update(dx=1 * units.nm)
wg2.calc('1.55 um', nmodes=1, pol='TM')
x= np.linspace(-2,2,grid_size)
Ey = wg2.mode[0].E.y.real(x,0)
plt.plot(x,abs(Ey))
plt.show()




wg2.mesh.box = wg2.box.expand(5)
wg2.mesh.update(dx=1 * units.nm)
wg2.calc('1.55 um', nmodes=1, pol='TE')
x= np.linspace(-1,1,grid_size)
Ey = wg2.mode[0].E.x.real(x,0)
plt.plot(x,abs(Ey))
plt.show()
# Ey = wg.mode[0].H.y
# xx = '150nm'
# y = np.linspace(-1*um, 1*um, 5000)
# Ey_imag = Ey.real(0, y)
# Ey_imag_abs = abs(Ey.real)(0, y)
# plt.plot(y, Ey_imag_abs.m)
# plt.show()
effective_index = wg2.mode[0].neff()
mode = Ey.m
plt.plot(x, mode)
plt.show()

sub_mat = Dielectric(3, color='cyan')
air = Dielectric(1, color='white')
core_mat = Dielectric(1.445, color='teal')
core = Layer(core_mat, h=width, rot='0 deg')

waveguide = CrossSection(sub_mat, core)
waveguide.plot()
plt.show()
waveguide.mesh.box = waveguide.box.expand(5)
waveguide.mesh.update(dy=1 * units.nm)
waveguide.calc('1.55 um', nmodes=1, pol='TM')
effective_index = waveguide.mode[0].neff()

#find pick postion and shift the mode towards that positon
x0 = np.linspace(-4* um, 4* um, grid_size)
z0 = np.linspace(0 * um, 30 * um, grid_size)
wavelength = 1.55 * um


plt.plot(x, mode)
plt.show()


t0 = Scalar_mask_XZ(x=x0, z=z0, wavelength=wavelength,n_background=1)

# waveguide(x = 0,
#           z = 0,
#           width =0.6,
#           length =10,
#           effective_index = effective_index,
#           rotation = 0,
#           scalar_mask = t0)

grating_coupler(x=0,
                z=0,
                period=1,
                width=width,
                teeth_number=10,
                length_base=2,
                ff=0.5,
                etch_depth=0.25,
                index_teeth=effective_index,
                effective_index= effective_index,
                substrate_index=1.55,
                substrate_width  = 0,
                scalar_mask= t0
               )
#MMIM
u0 = Scalar_source_X(x=x0, wavelength=wavelength)

u0.user_defined_mode(x,mode, angle=0, neff=effective_index, x_location=0)

t0.incident_field(u0)

t0.BPM(verbose=True)
t0.draw(kind='intensity',
        normalize=True,
        logarithm=False,
        draw_borders=True,
        colormap_kind = "viridis"
        )
plt.show()

amp_prof1 = t0.profile_longitudinal(kind='intensity',
                                  x0=3*um,
                                  logarithm=False,
                                  draw=False,
                                  filename='')

plt.plot(z0,amp_prof1)
plt.show()

amp_prof = t0.profile_longitudinal(kind='intensity',
                                  x0=-3*um,
                                  logarithm=False,
                                  draw=False,
                                  filename='')

plt.plot(z0,amp_prof)
plt.show()


amp_prof = t0.profile_transversal(kind='intensity',
                                  z0=25*um,
                                  logarithm=False,
                                  draw=False,
                                  filename='')

plt.plot(x0,amp_prof)
plt.title('Output Profile')
plt.ylabel('Intensity')
plt.tight_layout()
plt.xlabel('um')
plt.show()



amp_prof = t0.profile_transversal(kind='intensity',
                                  z0=1*um,
                                  logarithm=False,
                                  draw=False,
                                  filename='')

plt.plot(x0,amp_prof)
plt.title('Input Profile')
plt.ylabel('Intensity')
plt.xlabel('um')
plt.show()


#TODO
#first simulating a simple waveguide
#for gc creat a functon
#drawing shapes first and converting them to the image
#launching mode
#getting coordinates for grid
#calculating mode at a cross section
#extracting neff from device designer


from scipy.fft import fft, fftfreq
fft_amp = fft(amp_prof1)
plt.plot(z0,fft_amp)
plt.show()
#
# import scipy.signal
#
# # apply a 3-pole lowpass filter at 0.1x Nyquist frequency
# b, a = scipy.signal.butter(3, 0.05)
# filtered = scipy.signal.filtfilt(b, a, amp_prof)
# plt.plot(z0, filtered)
# plt.show()


