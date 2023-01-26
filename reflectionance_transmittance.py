from vpipda.material import *
from vpipda.layout import *
from scipy.constants import c
from vpipda import units as u
import matplotlib.pyplot as plt
from diffractio import degrees, plt, um,mm


import numpy as np


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



period = 265e-9
number_grating = 4
cladding = HalfPlane(SiO2)
etch = 0.3
wavelengths =np.linspace(1.54, 1.55, 50)
R = []
T = []
for wavelength in wavelengths:


    freq = c/wavelength*1e-6
    slab = Rectangle(Si, w=0.5, h=0.45, rt=(0.45 / 2, 0.45), rot='0')

    wg = CrossSection([Air, cladding, slab])
    wg.name = 'Waveguide Cross section'
    wg.plot()
    plt.show()


    wg.mesh.box = Box2D((-1, -0.3), (1, 0.8))
    dx = dy = 0.05
    wg.mesh.update(dx, dy)

    wg.calc(str(wavelength) + ' um', nmodes=1)
    effective_indexTE_unetched= wg.mode[0].neff(freq  = str(freq) + ' THz')

    Ex = wg.mode[0].E.x

    y = np.linspace(-1, 1, 5000)

    x = np.linspace(-1, 1, 5000)
    y = '150nm'
    Ex_imag = Ex.real(x, y)
    Ex_imag_abs = abs(Ex.real)(x, y)



    slab = Rectangle(Si, w=0.5, h=0.45-etch, rt=(0.45 / 2, 0.45-etch), rot='0')

    wg = CrossSection([Air, cladding, slab])
    wg.name = 'Waveguide Cross section'
    wg.plot()
    plt.show()


    wg.mesh.box = Box2D((-1, -0.3), (1, 0.8))
    dx = dy = 0.05
    wg.mesh.update(dx, dy)

    wg.calc(str(wavelength) + ' um', nmodes=1)
    effective_indexTE_etched= wg.mode[0].neff(freq  = str(freq) + ' THz')

    Ex = wg.mode[0].E.x

    y = np.linspace(-1, 1, 5000)

    x = np.linspace(-1, 1, 5000)
    y = '150nm'
    Ex_imag = Ex.real(x, y)
    Ex_imag_abs = abs(Ex.real)(x, y)


    delta_N = effective_indexTE_unetched -effective_indexTE_etched
    N_eff = (effective_indexTE_unetched +effective_indexTE_etched)/2



    #delta_N  = 3.4575 - 1
    r_fresnel_low_high = -delta_N/(2*N_eff)
    r_fresnel_high_low = delta_N/(2*N_eff)


    wavelength = wavelength*1e-6
    Bragg_wavelength = 2*N_eff*period
    #period=wavelength/(2*N_eff)
    #bragg_period = period

    k = 2*np.pi/wavelength
    beta = k*N_eff


    kappa = 2*abs(r_fresnel_high_low)/period
    kappa = kappa

    delta_beta = beta - 2*np.pi*N_eff/Bragg_wavelength


    yotta = kappa**2 - delta_beta**2
    yotta = yotta *1e-9

    L = number_grating*period
    #COSH SINH
    cosh_part = np.cosh(yotta*L)
    sinh_part = np.sinh(yotta*L)
    delta_beta_yotta = 1j*delta_beta/yotta
    kappa_imag_part = -1j*kappa/yotta
    division_reflectance = cosh_part + delta_beta_yotta*sinh_part

    transmittance = cosh_part - delta_beta_yotta *sinh_part


    reflectance =kappa_imag_part*sinh_part/division_reflectance

    R.append(reflectance)
    T.append(transmittance)



norm_T = normalize(np.real(T))
norm_R = normalize(abs(np.real(R)))



plt.plot(wavelengths, 1-abs(np.real(R)), label=' transmission')
plt.plot(wavelengths, np.real(R), label='reflectance')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
#todo
#1redo this and understand regorously what is what
