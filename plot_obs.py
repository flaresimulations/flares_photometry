import numpy as np
import pandas as pd
from FLARE.photom import M_to_lum

def plot_Fink15(z, axs):

    data = np.genfromtxt(f"./Obs_data/uv_lum_Fink15_z{z}.txt", delimiter=',', skip_header=1)
    M, phi, phi_up, phi_low, uplims = data[:,0], data[:,1]*1e-3, data[:,2]*1e-3, data[:,3]*1e-3, data[:,4]
    phi_up = np.log10(phi + phi_up) - np.log10(phi)
    phi_low = np.log10(phi) - np.log10(phi - phi_low)

    axs.errorbar(M, np.log10(phi), yerr = [phi_low, phi_up], uplims = uplims, label = 'Finkelstein 2015', ls = 'None', fmt = '.', color = 'green')


def plot_Bouw15(z, axs):

    data = np.genfromtxt(f"./Obs_data/uv_lum_Bouw15_z{z}.txt", delimiter=',', skip_header=1)
    M, phi, phi_err, uplims = data[:,0], data[:,1], data[:,2], data[:,3]
    phi_up = np.log10(phi + phi_err) - np.log10(phi)
    phi_low = np.log10(phi) - np.log10(phi - phi_err)

    axs.errorbar(M, np.log10(phi), yerr = [phi_low, phi_up], uplims = uplims, label = 'Bouwens 2015', ls = 'None', fmt = '.', color = 'blue')


def plot_obs(z, axs):


    if z == 5:

        plot_Fink15(z, axs)
        plot_Bouw15(z, axs)

    if z == 6:

        plot_Fink15(z, axs)
        plot_Bouw15(z, axs)

        data = np.genfromtxt("./Obs_data/uv_lum_Atek18_z6.txt", delimiter=',', skip_header=1)
        M, phi, phi_err = data[:,0], data[:,1], data[:,2]
        axs.errorbar(M, phi, yerr = phi_err, label = 'Atek 2018', ls = 'None', fmt = '.', color = 'brown')

    if z == 7:

        plot_Fink15(z, axs)
        plot_Bouw15(z, axs)

    if z == 8:

        plot_Fink15(z, axs)
        plot_Bouw15(z, axs)

    if z == 10:

        plot_Bouw15(z, axs)

        data = np.genfromtxt(f"./Obs_data/uv_lum_Oesch18_z{z}.txt", delimiter=',', skip_header=1)
        M, phi, phi_up, phi_low, uplims = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
        phi_up = np.log10(phi + phi_up) - np.log10(phi)
        phi_low = np.log10(phi) - np.log10(phi - phi_low)

        axs.errorbar(M, np.log10(phi), yerr = [phi_low, phi_up], uplims = uplims, label = 'Oesch 2018', ls = 'None', fmt = '.', color = 'black')
