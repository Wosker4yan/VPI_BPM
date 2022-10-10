# -*- coding: utf-8 -*-
"""
Created on T Sep 30 11:55 2022
Calculating the xs sections and polarizations for the future simulations
Functions to calculate overlap
Normalizing
Getting the complex fields
@author: Vahram
"""

from vpipda.material import *
from vpipda.layout import *
from vpipda import units
import numpy as np
from diffractio.scalar_sources_X import Scalar_source_X
from PIL import Image
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from diffractio import degrees, plt, um
from diffractio.scalar_masks_XZ import Scalar_mask_XZ
from masktopolygon import Masktopolygon
import math
import sys, os

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def norm(field, x):
    """Normalization of a field.
    Normalizing the field
    Args:
        field (list or array): Input field.
        x (int): Grid difference (delta) x0 parameter.
    Returns:
        normalized field (float)
    """

    return np.sqrt((sum(field * np.conjugate(field)))) * (x[1] - x[0])





def nearest(vector, number):
    """Function from diffractio module to be used alongside with the profile
    Args:
        vector (list): vector
        number (float): number
    Returns:
        indexes
        values
        distances
    """

    # noinspection PyTypeChecker
    indexes = np.abs(vector - number).argmin()
    values = vector.flat[indexes]
    distances = values - number

    return indexes, values, distances


def profile(field, z, z0=0 * um):
    """Extracts the profile of the field at a desired location
    Args:
        field (array): Complex field
        z (vector): Direction  of propagation
        z0 (float): Location to extract the slice
    Returns:
        u (complex array): field slice at the desired location
    """

    index, value, distance = nearest(vector=z, number=z0)
    u = np.squeeze(field[:, index])

    return u


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


def overlap(field_1, field_2, x):
    """Overlap between two fields.
    Args:
        field_1 (list or array): Input field 1.
        field_2 (list or array): Input field 2.
        x (list or array): Grid difference (delta) x0 parameter.
    Returns:
        overlap between two fields (float)
    """

    return sum(field_1 * np.conjugate(field_2) * (x[1] - x[0]))


class VPI_BPM:

    def __init__(
            self,
            grid_size=8000,
            beam_waist=1,
            core_index=3.4,
            substrate_index=1.45,
            image=None,
            wavelength=1.55,
            polarization='TE',
            filename='./mask.png',
            cell=None,
            pin='b0',
            mode=0,
            plotting=True,
            sim_mode='BPM',
            angle=0
    ):

        """Initializes the Nazca BPM class.
        Args:
            grid_size (int): number of points for running the BPM simulation
            beam_waist (float): size of the beam waist launched into the structure
            core_index: refractive index of the cladding
            substrate_index (float): refractive index of the substrate
            image (float): mask that is given
            wavelength (float): wavelength of the light
            polarization (string): polarization for the mode calculation
            filename (string): name of the image saved in the same directory that the code is run
            cell (netlist.Cell): structure generated with nazca
            pin (string): Pin name from nazca
            mode (int): select modes of different orders
            plotting (bool): Set true/false for plotting.
            sim_mode (string):select for Beam propagation method or Wave propagation method
        Returns:
            None
        """

        self.grid_size = grid_size
        self.beam_waist = beam_waist
        self.core_index = core_index
        self.image = image
        self.wavelength = wavelength
        self.polarization = polarization
        self.mode = mode
        self.substrate_index = substrate_index
        self.cell = cell
        self.filename = filename
        self.pin = pin
        self.plotting = plotting
        self.sim_mode = sim_mode
        self.angle = angle


    def __str__(self):
        message = f'Grid size: {self.grid_size}\n' \
                  f'Beam waist: {self.beam_waist}\n' \
                  f'Core Index: {self.core_index}\n' \
                  f'Converted image from the cell: {self.image}\n' \
                  f'Wavelength: {self.wavelength}\n' \
                  f'Polarization: {self.polarization}\n' \
                  f'Order of the Fundamental Mode: {self.mode}\n' \
                  f'Index of the Substrate. Input width: {self.substrate_index}\n' \
                  f'Path to the Saved Image Converted from Nazca Cell: {self.filename}\n' \
                  f'Pin is : {self.pin}\n' \
                  f'Plotting is set to : {self.plotting}\n' \
                  f'Simulation is set to is set to : {self.sim_mode}\n'

        return message



    def get_ic_width(
            self,
            cell=None,
            pin=None,
    ):

        """Get width from pin.
        Args:
            cell (netlist): Cell generated from Nazca
            pin (string): selected pin
        Returns:
           (float) width of the selected pin
        """

        if pin is None:
            pin = self.pin

        if cell is None:
            cell = self.cell

        return cell.pin[pin].width

    def cell2image(
            self,
            cell=None,
            filename=None,
            grid_size=None,
            save = True,
            offset=2
    ):
        """Converting the image to cell. Saves the image on the current directory.
        Args:
            cell (netlist.Cell): structure generated with nazca
            filename (string): name of the image saved in the same directory that the code is run
            grid_size (int): number of points for running the BPM simulation
            offset (float): offset value for the grid dimension
            save (boolean): saving the plot of the mask
        Returns (float):
            grid_offset_z,  grid_offset_x: grid sizes and dimensions for later simulations
        """

        if filename is None:
            filename = self.filename
        if cell is None:
            cell = self.cell
        if grid_size is None:
            grid_size = self.grid_size

        diffractio = Masktopolygon()
        polygons = diffractio.get_polygons(cell)
        polygons = polygons['1/0/None']

        for item in polygons:
            xs, ys = zip(*item)
            plt.fill(xs, ys, 'k')

        x, y = max(max(polygons))

        grid_offset_z = x
        grid_offset_x = y + offset

        plt.ylim(-grid_offset_x, grid_offset_x)
        plt.xlim(0, grid_offset_z)
        plt.gca().invert_yaxis()
        plt.axis('off')
        if save:
            plt.savefig(filename, bbox_inches='tight', pad_inches=0)
            im = Image.open(filename)

        size = grid_size, grid_size

        if save:
            im_resized = im.resize(size, Image.ANTIALIAS)
            im_resized.save(filename, 'png')

        plt.close()

        return grid_offset_x, grid_offset_z

    def mode_and_neff(
            self,
            cell=None,
            filename=None,
            substrate_index=None,
            core_index=None,
            polarization=None,
            mode_order=0,
            grid_size=None,
            number_of_modes=1,
            wavelength=None,
            input_pin='a0'
    ):

        """Calculating the slab index and mode profile from VPI devices designer
        Args:
            cell :cell of the image for simulation
            filename (string): name of the image saved in the same directory that the code is run
            core_index (float): refractive index of the core ,
            substrate_index (float): substrate refractive index,
            core_index (float): core refractive index,
            polarization (string): polarization of light (TE, TM),
            mode (int): order of the fundamental mode,
            input_pin (str): input pin from nazca
            grid_size (int): size of the grid for the simulation, number of points in multilayer
            number_of_modes (int): select the number of modes to calculate
        """
        if filename is None:
            filename = self.filename
        if cell is None:
            cell = self.cell
        if substrate_index is None:
            filename = self.substrate_index
        if core_index is None:
            core_index = self.core_index
        if polarization is None:
            polarization = self.polarization
        if grid_size is None:
            grid_size = self.grid_size
        if input_pin is None:
            input_pin = self.input_pin
        if wavelength is None:
            wavelength = self.wavelength

        width = self.get_ic_width(pin=input_pin)

        grid_offset_x, grid_offset_z = self.cell2image(
            cell=cell,
            grid_size=grid_size,
            filename=filename,
            save = False
        )


        sub_mat = Dielectric(substrate_index, color='cyan')
        core_mat = Dielectric(core_index, color='teal')
        core = Layer(core_mat, h=width, rot='90 deg')

        waveguide = CrossSection(sub_mat, core)
        x = np.linspace(-grid_offset_x, grid_offset_x, grid_size)
        wavelength = wavelength*1e-6
        if polarization == 'TM':
            waveguide.mesh.box = waveguide.box.expand(5)
            waveguide.mesh.update(dx=1 * units.nm)
            waveguide.calc(str(wavelength) + ' um', nmodes=number_of_modes, pol='TM')
            E_TM = waveguide.mode[mode_order].E.y.real(x, 0)
            effective_index = waveguide.mode[mode_order].neff()
            E = E_TM.m
            return E, effective_index, x


        elif polarization == 'TE':
            waveguide.mesh.box = waveguide.box.expand(5)
            waveguide.mesh.update(dx=1 * units.nm)
            waveguide.calc(str(wavelength) + ' um', nmodes=number_of_modes, pol='TE')
            E_TE = waveguide.mode[mode_order].E.x.real(x, 0)
            effective_index = waveguide.mode[mode_order].neff()
            E = E_TE.m
            return E, effective_index, x


        else:
            print('POLARIZATION TYPE IS NOT SELECTED')

    def plot_mask(
            self,
            cell=None
    ):

        """Function for plotting the mask that is going to be simulated
        Args:
            cell (netlist.Cell): structure generated with nazca
        Plotting the mask.
            use for visualizing the structure that is going to be simulated.
        Returns:
              Plotting
        """

        if cell is None:
            cell = self.cell

        diffractio = Masktopolygon()
        polygons = diffractio.get_polygons(cell)
        polygons = polygons['1/0/None']

        for item in polygons:
            xs, ys = zip(*item)
            plt.fill(xs, ys, 'k')
        plt.title('Mask')
        plt.xlabel('X [um]')
        plt.ylabel('Y [um]')

        plt.tight_layout()
        plt.show()

    def run_bpm_mode(self,
                     cell=None,
                     pin=None,
                     grid_size=None,
                     filename=None,
                     substrate_index=None,
                     wavelength=None,
                     sim_mode=None,
                     plotting=True,
                     input_monitor_location=1.0,
                     output_monitor_location=10.0,
                     log_scale = False,
                     angle=0,
                     x_location=0,
                     input_pin='a0',
                     output_pin='b0'
                     ):

        """
        Run BPM with gaussian beam
        Args:
             cell (netlist.Cell): structure generated with nazca
             pin (string): Pin name from nazca
             grid_size (int): number of points for running the BPM simulation
             filename (string): name of the image saved in the same directory that the code is run
             substrate_index: refractive index of the substrate
             wavelength (float): wavelength of the light
             sim_mode (str): 'BPM' or 'WPM' for simulation mode
             plotting (bool): Set true/false for plotting
             input_monitor_location (float): location of the input monitor
             output_monitor_location (float): location of the output monitor
             log_scale (boolean): True if log scale is needed
             angle (float): angle of the launched beam
             x_location (float): location of the beam in X direction
             input_pin (str): input pin
             output_pin (str): output pin
        Returns (float):
            transmission, coefficient, overlap coefficient with the mode profile
        """

        if filename is None:
            filename = self.filename
        if grid_size is None:
            grid_size = self.grid_size
        if cell is None:
            cell = self.cell
        if pin is None:
            pin = self.pin
        if substrate_index is None:
            substrate_index = self.substrate_index
        if wavelength is None:
            wavelength = self.wavelength
        if plotting is None:
            plotting = self.plotting
        if sim_mode is None:
            sim_mode = self.sim_mode

        width = self.get_ic_width(cell=cell, pin=pin)
        width_2 = self.get_ic_width(cell=cell, pin='b0')

        grid_offset_x, grid_offset_z = self.cell2image(
            cell=cell,
            grid_size=grid_size,
            filename=filename
        )

        plt.close()

        E, neff, x = self.mode_and_neff(mode_order=0,input_pin=input_pin)

        ref_background = substrate_index
        E = normalize(E)
        x0 = np.linspace(-grid_offset_x, grid_offset_x, grid_size)  # minus plus !!!
        z0 = np.linspace(0, grid_offset_z, grid_size)
        u0 = Scalar_source_X(x=x0, wavelength=wavelength)
        u0.user_defined_mode(x, E, angle, neff, x_location)


        t0 = Scalar_mask_XZ(x=x0, z=z0, wavelength=wavelength)

        t0.image(filename=filename,
                 n_max=neff,
                 n_min=ref_background,
                 angle=0 * degrees,
                 invert=False)

        t0.incident_field(u0)
        if sim_mode == 'BPM':
            t0.BPM(verbose=False)
        else:
            t0.WPM(verbose=False)

        if plotting:
            t0.draw(kind='intensity',
                    normalize=True,
                    logarithm=log_scale,
                    draw_borders=True,
                    colormap_kind='viridis',
                    min_incr=0.05,
                    )
            plt.title('Structure and refractive index contrast')
            beam_waist_vline = width / grid_offset_x

            wg_width_vline = width / (2 * grid_offset_x)
            wg_width_vline2 = width_2 / (2 * grid_offset_x)
            y_mid = 0.5

            cbar = plt.colorbar()
            #plt.clim(0, 1)

            cbar.set_label('Normalized Intensity', rotation=270, size=10, labelpad=10)



            plt.axvline(0,
                        ymin=y_mid - beam_waist_vline,
                        ymax=y_mid + beam_waist_vline,
                        color='#ff4a4a',
                        label='source',
                        linewidth=4.0)

            plt.axvline(input_monitor_location,
                        ymin=y_mid - wg_width_vline,
                        ymax=y_mid + wg_width_vline,

                        color='deepskyblue',
                        label='input monitor',
                        linewidth=4.0)

            plt.axvline(output_monitor_location,
                        ymin=y_mid - wg_width_vline2,
                        ymax=y_mid + wg_width_vline2,
                        color='#F2A800',
                        label='output monitor',
                        linewidth=4.0)

            plt.title(f'{cell.cell_name} at $\lambda$={wavelength} $\mu$m')

            plt.legend(loc='upper left', prop={'size': 5.5})
            plt.tight_layout()
            plt.show()

        field = t0.u
        z = t0.z

        amp_prof_input_complex = profile(field, z, z0=input_monitor_location)
        amp_prof_output_complex = profile(field, z, z0=output_monitor_location)

        amp_prof_input = t0.profile_transversal(kind='intensity',
                                          z0=input_monitor_location,
                                          logarithm=False,
                                          draw=False,
                                          filename='')

        amp_prof_output = t0.profile_transversal(kind='intensity',
                                          z0=output_monitor_location,
                                          logarithm=False,
                                          draw=False,
                                          filename='')


        E, neff, x = self.mode_and_neff(
            mode_order=0,
            input_pin=output_pin
        )

        #TODO
        #CHECK OVERLAP INTEGRALS
        norm_amp_input = norm(amp_prof_input, x0)
        norm_E = norm(E, x0)

        norm_amp_output = norm(amp_prof_output, x0)


        overlap1 = overlap(amp_prof_input, E, x0) / (norm_amp_input * norm_E)
        overlap2 = overlap(amp_prof_output, E, x0) / (norm_amp_output * norm_E)

        sab = overlap2 / overlap1
        transmission = abs(sab) ** 2


        # norm_E = norm(E, x0)
        #
        # transmission = np.abs(overlap(norm_amp_input, norm_E, x0)) / np.abs(overlap(norm_amp_output, norm_E, x0))


        return amp_prof_input,  amp_prof_output,  x0, transmission


    def visualize(self,
                  plotting=None,
                  cell=None,
                  grid_size=None,
                  filename=None,
                  substrate_index=None,
                  wavelength=None,
                  sim_mode=None,
                  input_pin =None,
                  monitor_location=0.0,
                  angle=0,
                  x_location=0
                  ):

        """
        Visualizing the field profile and the index profile
        Args:
            plotting (bool): Select True for plotting
            cell (netlist.Cell): structure generated with nazca
            grid_size (int): number of points for running the BPM simulation
            filename (string): name of the image saved in the same directory
            substrate_index: refractive index of the substrate
            wavelength (float): wavelength of the light
            sim_mode (str): 'BPM' or 'WPM' for simulation mode
            monitor_location (float): position of the monitor for taking the field
            angle (float): Angle for the mode profile
            x_location (float): Location of the mode
        Returns:
            Transmission sweep, grid size
        """

        if plotting is None:
            plotting = self.plotting
        if grid_size is None:
            grid_size = self.grid_size
        if filename is None:
            filename = self.filename
        if cell is None:
            cell = self.cell
        if substrate_index is None:
            substrate_index = self.substrate_index
        if sim_mode is None:
            sim_mode = self.sim_mode
        if wavelength is None:
            wavelength = self.wavelength
        if input_pin is None:
            input_pin = self.pin

        self.get_ic_width()

        grid_offset_x, grid_offset_z = self.cell2image(
            cell=cell,
            grid_size=grid_size,
            filename=filename
        )

        plt.close()

        E, neff, x = self.mode_and_neff(
            mode_order=0,
            input_pin=input_pin
        )

        ref_background = substrate_index

        x0 = np.linspace(-grid_offset_x, grid_offset_x, grid_size)
        z0 = np.linspace(0 * um, grid_offset_z, grid_size)

        u0 = Scalar_source_X(x=x0, wavelength=wavelength*um)
        u0.user_defined_mode(x, E, angle, neff, x_location)
        t0 = Scalar_mask_XZ(x=x0, z=z0, wavelength=wavelength*um)

        t0.image(filename=filename,
                 n_max=neff,
                 n_min=ref_background,
                 angle=0 * degrees,
                 invert=False)

        t0.incident_field(u0)

        if sim_mode == 'BPM':
            t0.BPM(verbose=False)
        else:
            t0.WPM(verbose=False)

        amp_prof = t0.profile_transversal(kind='amplitude',
                                          z0=monitor_location,
                                          logarithm=False,
                                          draw=False,
                                          filename='')

        index_profile = t0.profile_transversal(kind='refraction_index',
                                               z0=monitor_location,
                                               logarithm=False,
                                               draw=False,
                                               filename='')
        if plotting:
            fig, ax1 = plt.subplots()

            color1 = 'tab:grey'
            color2 = 'tab:green'

            ax1.set_xlabel('x [$\mu$m]')
            ax1.set_ylabel('Index Profile', color='grey')
            ax1.plot(x0, index_profile, color=color1)
            ax1.set_title(f'z={monitor_location} $\mu$m')
            ax1.tick_params(axis='y', labelcolor='#78f213')
            ax1.set_facecolor("#000000")

            ax2 = ax1.twinx()
            ax2.set_ylabel('Amplitude [norm]', color='#78f213')
            ax2.tick_params(axis='y', labelcolor='#78f213')
            ax2.plot(x0, amp_prof, color=color2)
            ax2.set_facecolor("#000000")

            ax1.grid('both', color='#292929', linestyle='--', linewidth=0.5)
            plt.tight_layout()
            plt.show()

        return amp_prof, index_profile

    def plot_waveguide(self,
                       substrate_index=None,
                       core_index=None,
                       pin=None
                       ):


        """Function for plotting the waveguide
        Args:
            core_index (float): refractive index of the cladding
            substrate_index (float): refractive index of the substrate
            pin (str): input pin of the waveguide.
        """
        if substrate_index is None:
            filename = self.substrate_index
        if core_index is None:
            core_index = self.core_index
        if pin is None:
            pin = self.pin

        width = self.get_ic_width(cell=None, pin=pin)

        sub_mat = Dielectric(substrate_index, color='cyan')
        core_mat = Dielectric(core_index, color='teal')
        core = Layer(core_mat, h=width, rot='90 deg')
        waveguide = CrossSection(sub_mat, core)

        waveguide.plot()
        plt.title('The widht of the waveguide is '+str(width) +'$\mu$m')
        plt.tight_layout()
        plt.show()