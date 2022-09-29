import numpy as np
import matplotlib.pyplot as plt
from diffractio import degrees, plt, um,mm
from diffractio.scalar_masks_XZ import Scalar_mask_XZ
from diffractio.scalar_sources_X import Scalar_source_X
import numpy as np

field = np.load('Exfield.npy')

x0 = np.linspace(-10 * um, 10 * um, 1000)
z0 = np.linspace(0 * um, 10 * um, 1000)
wavelength = 1.55 * um

t0 = Scalar_mask_XZ(x=x0, z=z0, wavelength=wavelength)

t0.rectangle(
    r0=(0 * um, 0 * um),
    size=(2 * um, 15 * um),
    angle=0 * degrees,
    refraction_index=1.5)

u0 = Scalar_source_X(x=x0, wavelength=wavelength)

u0.user_defined_mode(x0, field, 0, 3, x_location=0)

t0.incident_field(u0)
t0.BPM(verbose=False)
t0.draw(kind='intensity',
        normalize=True,
        logarithm=0,
        draw_borders=True,
        min_incr=0.05,
        )
plt.show()