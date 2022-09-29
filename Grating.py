# Generated with SMOP  0.41
# .\Grating.m
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve



def TMM_Grating_RT(wavelength=None, Period=None, NG=None, n1=None, n2=None, loss=None):


    M = TMM_Grating_Matrix(wavelength, Period, NG, n1, n2, loss)
    q = len(wavelength)
    T = abs(np.ones(q, 1) / np.squeeze(M(1, 1, np.arange()))) ** 2
    R = abs(np.squeeze(M(2, 1, np.arange())) / np.squeeze(M(1, 1, np.arange()))) ** 2


def TMM_Grating_Matrix(wavelength=None, Period=None, NG=None, n1=None, n2=None, loss=None):


    l = Period / 2
    T_hw1 = TMM_HomoWG_Matrix(wavelength, l, n1, loss)
    T_is12 = TMM_IndexStep_Matrix(n1, n2)
    T_hw2 = TMM_HomoWG_Matrix(wavelength, l, n2, loss)
    T_is21 = TMM_IndexStep_Matrix(n2, n1)
    q = len(wavelength)
    Tp = np.zeros(2, q)
    T = np.copy(Tp)
    for i in np.arange(1, len(wavelength)).reshape(-1):
        Tp[np.arange(), np.arange(), i] = np.dot(
            np.dot(np.dot(T_hw2(np.arange(), np.arange(), i), T_is21(np.arange(), np.arange(), i)),
                   T_hw1(np.arange(), np.arange(), i)), T_is12(np.arange(), np.arange(), i))
        T[np.arange(), np.arange(), i] = Tp(np.arange(), np.arange(), i) ** NG
        # for an FP cavity, 1st order cavity, insert a high index region, n2.
        T[np.arange(), np.arange(), i] = np.dot(
            np.dot(np.dot(Tp(np.arange(), np.arange(), i) ** NG, (T_hw2(np.arange(), np.arange(), i)) ** 1),
                   Tp(np.arange(), np.arange(), i) ** NG), T_hw2(np.arange(), np.arange(), i))


def TMM_HomoWG_Matrix(wavelength=None, l=None, neff=None, loss=None):


    beta = np.dot(np.dot(2, np.pi), neff) / wavelength - np.dot(1j, loss) / 2

    def TMM_Grating_RT(wavelength=None, Period=None, NG=None, n1=None, n2=None, loss=None, *args, **kwargs):
        varargin = TMM_Grating_RT.varargin
        nargin = TMM_Grating_RT.nargin

        M = TMM_Grating_Matrix(wavelength, Period, NG, n1, n2, loss)
        q = np.length(wavelength)
        T = abs(np.ones(q, 1) / np.squeeze(M(1, 1, np.arange()))) ** 2
        R = abs(np.squeeze(M(2, 1, np.arange())) / np.squeeze(M(1, 1, np.arange()))) ** 2

    def TMM_Grating_Matrix(wavelength=None, Period=None, NG=None, n1=None, n2=None, loss=None, *args, **kwargs):
        varargin = TMM_Grating_Matrix.varargin
        nargin = TMM_Grating_Matrix.nargin

        l = Period / 2
        T_hw1 = TMM_HomoWG_Matrix(wavelength, l, n1, loss)
        T_is12 = TMM_IndexStep_Matrix(n1, n2)
        T_hw2 = TMM_HomoWG_Matrix(wavelength, l, n2, loss)
        T_is21 = TMM_IndexStep_Matrix(n2, n1)
        q = np.length(wavelength)
        Tp = np.zeros(2, q)
        T = np.copy(Tp)
        for i in np.arange(1, np.length(wavelength)).reshape(-1):
            Tp[np.arange(), np.arange(), i] = np.dot(
                np.dot(np.dot(T_hw2(np.arange(), np.arange(), i), T_is21(np.arange(), np.arange(), i)),
                       T_hw1(np.arange(), np.arange(), i)), T_is12(np.arange(), np.arange(), i))
            T[np.arange(), np.arange(), i] = Tp(np.arange(), np.arange(), i) ** NG
            # for an FP cavity, 1st order cavity, insert a high index region, n2.
            T[np.arange(), np.arange(), i] = np.dot(
                np.dot(np.dot(Tp(np.arange(), np.arange(), i) ** NG, (T_hw2(np.arange(), np.arange(), i)) ** 1),
                       Tp(np.arange(), np.arange(), i) ** NG), T_hw2(np.arange(), np.arange(), i))

    def TMM_HomoWG_Matrix(wavelength=None, l=None, neff=None, loss=None, *args, **kwargs):
        varargin = TMM_HomoWG_Matrix.varargin
        nargin = TMM_HomoWG_Matrix.nargin

        beta = np.dot(np.dot(2, np.pi), neff) / wavelength - np.dot(1j, loss) / 2

        T_hw = np.zeros(2, np.length(neff))
        T_hw[1, 1, np.arange()] = np.exp(np.dot(np.dot(1j, beta), l))
        T_hw[2, 2, np.arange()] = np.exp(np.dot(np.dot(- 1j, beta), l))

    def TMM_IndexStep_Matrix(n1=None, n2=None, *args, **kwargs):
        varargin = TMM_IndexStep_Matrix.varargin
        nargin = TMM_IndexStep_Matrix.nargin

        T_is = np.zeros(2,  np.length(n1))
        a = (n1 + n2) / (np.dot(2, np.sqrt(np.multiply(n1, n2))))
        b = (n1 - n2) / (np.dot(2, np.sqrt(np.multiply(n1, n2))))
        T_is[1, 1, np.arange()] = a
        T_is[1, 2, np.arange()] = b
        T_is[2, 1, np.arange()] = b
        T_is[2, 2, np.arange()] = a
        T_hw[1, 1, np.arange()] = np.exp(np.dot(np.dot(1j, beta), l))
        T_hw[2, 2, np.arange()] = np.exp(np.dot(np.dot(- 1j, beta), l))


def TMM_IndexStep_Matrix(n1=None, n2=None):


    T_is = np.zeros(2, len(n1))
    a = (n1 + n2) / (np.dot(2, np.sqrt(np.multiply(n1, n2))))
    b = (n1 - n2) / (np.dot(2, np.sqrt(np.multiply(n1, n2))))
    T_is[1, 1, np.arange()] = a
    T_is[1, 2, np.arange()] = b
    T_is[2, 1, np.arange()] = b
    T_is[2, 2, np.arange()] = a



Period = 1e-07
NG = 200

L = np.dot(NG, Period)

width0 = 0.45

dwidth = 0.02

width1 = width0 - dwidth
width2 = np.copy(width0)
loss_dBcm = 3

loss = np.dot(np.dot(np.log(10), loss_dBcm) / 10, 100)
span = 3e-08
Npoints = 10000

if 1 == 1:
    neff_wavelength = lambda w=None: 2.4379 - np.dot(1.1193, (np.dot(w, 1000000.0) - 1.554)) - np.dot(0.035, (
                np.dot(w, 1000000.0) - 1.554) ** 2)

    dneff_width = lambda w=None: np.dot(10.4285, (w - 0.5) ** 3) - np.dot(5.2487, (w - 0.5) ** 2) + np.dot(1.6142, (w - 0.5))

    f = lambda lambda_=None: lambda_ - np.dot(np.dot(Period, 2),
                                           (neff_wavelength(lambda_) + (dneff_width(width2) + dneff_width(width1)) / 2))
    wavelength0 = fsolve(f, 1.55e-06)
    wavelengths = wavelength0 + np.linspace(- span / 2, span / 2, Npoints)
    n1 = neff_wavelength(wavelengths) + dneff_width(width1)

    n2 = neff_wavelength(wavelengths) + dneff_width(width2)

    R, T = TMM_Grating_RT(wavelengths, Period, NG, n1, n2, loss)
    plt.plot(np.dot(wavelengths, 1000000.0), np.concat([R, T]), 'LineWidth', 3)
    plt.hold('all')
    plt.lot(np.dot(np.concat([wavelength0, wavelength0]), 1000000.0), np.concat([0, 1]), '--')

    plt.xlabel('Wavelength [\mum]')
    plt.ylabel('Response')
    plt.axis('tight')


