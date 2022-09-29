from vpipda.material import *
from vpipda.layout import *
from vpipda import units
import matplotlib.pyplot as plt
import numpy as np




# function for semicircle
def semicircle(r, h, k, down = False):
    x0 = h - r  # determine x start
    x1 = h + r  # determine x finish
    x = np.linspace(x0, x1, 10000)  # many points to solve for y

    # use numpy for array solving of the semicircle equation
    if down:
        y = -(k + np.sqrt(r**2 - (x - h)**2))
    else:
        y = (k + np.sqrt(r**2 - (x - h)**2))

    return x, y


wavelength = 1.3e-6
period = 0.52e-6
n1 = 1
n2 = 2
cladding = HalfPlane(SiO2)
slab = Rectangle(Si, w=0.5, h=0.45, rt=(0.45 / 2, 0.45), rot='0')
wg = CrossSection([Air, cladding, slab])
wg.name = 'Waveguide Cross section'
wg.plot()
plt.show()

wg.mesh.box = Box2D((-1, -0.3), (1, 0.8))
dx = dy = 0.05
wg.mesh.update(dx, dy)

wg.calc('1.55 um', nmodes=1)
effective_indexTE= wg.mode[0].neff(freq  = '193.414 THz')
Ex = wg.mode[0].E.x
y = np.linspace(-1, 1, 5000)
x = np.linspace(-1, 1, 5000)
y = '150nm'
Ex_imag = Ex.real(x, y)
Ex_imag_abs = abs(Ex.real)(x, y)

G = 2*np.pi/period

k_air = 2*np.pi*n1/wavelength
beta = 2*np.pi*effective_indexTE/wavelength
degrees = 180/np.pi
theta = np.arcsin((effective_indexTE-wavelength/period)/n1)
theta2 = np.arcsin((effective_indexTE-wavelength/period)/n2)

angle_substrate = theta2*degrees
angle_air = theta*degrees
#theta = theta_rad*180/np.pi
k_substrate = 2*np.pi*n2/wavelength
x, y = semicircle(k_air, 0, 0, down=False)  # function call
plt.plot(x, y,  c='red')
x, y = semicircle(k_substrate, 0, 0, down=True)  # function call
plt.plot(x, y,  c='orange')

offset =1.1*beta
plt.arrow(0, 0,0, offset,  color = 'black',  head_width=900000, head_length=900000)
plt.arrow(0, 0,offset, 0,  color = 'black',  head_width=900000, head_length=900000)

plt.arrow(0, 0,0, -offset,  color = 'black',  head_width=900000, head_length=900000)
plt.arrow(0, 0,-offset, 0,  color = 'black',  head_width=900000, head_length=900000)

plt.arrow(0, 0,k_air*np.sin(theta), k_air*np.cos(theta),  color = 'red',  head_width=0, head_length=0)
plt.text(-k_air*np.sin(theta)-3*k_air, k_air*np.cos(theta),r'$k_a=2*\pi*n_1/\lambda$', color='red')


plt.arrow(0, 0,beta, 0,  color = 'violet',   head_width=800000, head_length=800000)
plt.text(beta/2-beta/8,  0+0.1*k_air,r'$\beta=2*\pi*n_eff/\lambda$', color='violet')

plt.arrow(0, 0,k_substrate*np.sin(theta2), -k_substrate*np.cos(theta2),  color = 'orange',  head_width=0.003, head_length=0.003)

plt.text(-k_substrate*np.sin(theta2)-3*k_air, -k_substrate*np.cos(theta2),r'$k_s=2*\pi*n_2/\lambda$', color='orange')


plt.arrow(beta-G, k_air*np.cos(theta),G, 0,  color = 'grey',  head_width=0, head_length=0, linestyle=(5, (1,3)))
plt.text(beta/2-beta/8,  k_air*np.cos(theta)+0.1*k_air,r'$\kappa=2*\pi/\Lambda$', color='grey')

plt.arrow(beta-G, 0-k_substrate*np.cos(theta2),G, 0,  color = 'grey',  head_width=0, head_length=0, linestyle=(5, (1,3)))
plt.text(beta,  -k_substrate*np.cos(theta2),r'$m=1$', color='grey')

plt.text(0,  k_air*np.cos(theta)+0.65*k_air,r'$\theta=$' +str(angle_air), color='k')
plt.text(0,  -k_substrate*np.cos(theta2) -0.65*k_air ,r'$\theta=$' +str(angle_substrate), color='k')


plt.arrow(beta-G, 0-k_substrate*np.cos(theta2),0,  k_air*np.cos(theta)+k_substrate*np.cos(theta2),  color = 'grey',  head_width=0, head_length=0, linestyle=(5, (1,3)))


plt.arrow(beta, 0-k_substrate*np.cos(theta2),0,  k_air*np.cos(theta)+k_substrate*np.cos(theta2),  color = 'grey',  head_width=0, head_length=0, linestyle=(5, (1,3)))


offset_x = offset + 1.5*beta

plt.axis('off')
plt.title('K vector diagram for grating')
plt.show()

