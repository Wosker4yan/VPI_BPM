import numpy as np
import matplotlib.pyplot as plt

Period=100e-9
NG=200
L=NG*Period
width0=0.45
dwidth=0.02
width1=width0 - dwidth
width2=width0
loss_dBcm=3
loss=np.log(10)*loss_dBcm/10*100
span=30e-9
Npoints = 10000
neff_wavelength = 2.4105054717128764
delta_N = 0.8799012809505122