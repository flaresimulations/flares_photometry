import time
import sys
import os
import numpy as np
import pandas as pd
import h5py
from functools import partial
import schwimmbad
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from FLARE.photom import lum_to_M, M_to_lum
sns.set_context("paper")
sns.set_style(style='white')

from plot_obs import plot_obs
from FLARE.LF import evo

obs = ['FLARES', 'Finkelstein_obs', 'Ma2019', 'Mason15', 'Yung2018']
obs_color = ['brown', 'indigo', 'violet', 'grey', 'black']


def get_hist(ii, tag, bins, inp='GEAGLE', filter = 'FUV', Luminosity = 'Dustcorr'):

    if inp == 'GEAGLE':

        num = str(ii)

        if len(num) == 1:
            num =  '0'+num

        filename = '../photometry/out/GEAGLE_{}_sp_info.hdf5'.format(num)


    with h5py.File(filename,'r') as hf:

        lum = np.array(hf[F"{tag}/Subhalo/BPASS/SalpeterIMF/ModelI/Luminosity/{Luminosity}/{filter}"])

        tmp, edges = np.histogram(lum_to_M(lum), bins = bins)

        return tmp


def get_all(tag, bins = np.arange(-25, -16, 0.5), inp = 'GEAGLE', filter = 'FUV', Luminosity = 'Dustcorr'):

    if inp == 'GEAGLE':

        sims = np.arange(0,38)

        df = pd.read_csv('weight_files/weights_grid.txt')
        weights = np.array(df['weights'])

        calc = partial(get_hist, tag = tag, bins = bins, inp = inp, filter = filter, Luminosity = Luminosity)

        pool = schwimmbad.MultiPool(processes=12)
        dat = np.array(list(pool.map(calc, sims)))
        pool.close()


        hist = np.sum(dat, axis = 0)
        out = np.zeros(len(bins)-1)
        for ii, sim in enumerate(sims):

            out+=dat[ii]*weights[ii]

        return out, hist

    else:

        out = get_hist(00, tag, bins, inp = 'REF')

        return out


def plot_tng(M, z, axs):

    redshifts = np.array([5, 6, 7, 8])
    logphi = np.array([[-3.244, -3.107, -3.398], [-3.079, -3.025, -3.608], [-3.846, -3.418, -4.209], [-4.445, -4.111, -4.714]])
    Mstar = np.array([[-21.17, -20.95, -21.21], [-20.61, -20.52, -21.31], [-21.18, -20.58, -21.47], [-21.38, -20.86, -21.44]])
    alpha = np.array([[-1.924, -1.884, -1.941], [-1.876, -1.833, -2.042], [-2.133, -1.967, -2.279], [-2.280, -2.216, -2.455]])

    labels = ['TNG19:Model-A', 'TNG19:Model-B', 'TNG19:Model-C']
    colors = ['yellow', 'magenta', 'pink']
    if z in redshifts:

        ok = np.where(redshifts == z)[0][0]
        for ii in range(3):
            theta = [Mstar[ok][ii], logphi[ok][ii], alpha[ok][ii]]
            axs.plot(M, model(M, theta), label = labels[ii], color = colors[ii])


def model(M, theta):

    Mstar, log10phistar, alpha = theta
    delta = M - Mstar

    return np.log10(0.4*np.log(10)) + log10phistar - 0.4*delta*(alpha+1.) - (np.log10(np.e)) * (10**(-0.4*delta))


fig,axs = plt.subplots(nrows = 2, ncols = 3, figsize=(18, 10), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

h = 0.6777
tbins = np.arange(-24, -16, 0.2)
tbincen = (tbins[1:]+tbins[:-1])/2.
tbinwidth = tbins[1:] - tbins[:-1]
filter = 'FUV'

tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']

for ii in range(6):

    z = ii+5

    df = pd.read_csv('Magnitude_limits.txt')
    low = lum_to_M(np.array(df[filter])[ii])
    bins = -np.arange(-low, 25, 0.3)[::-1]
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]
    vol = (4/3)*np.pi*(14/h)**3

    out, hist = get_all(tags[ii], bins, Luminosity = 'Dustcorr')
    phi = out/(vol*binwidth)
    xerr = np.ones(len(phi))*binwidth[0]/2.
    yerr = np.sqrt(out)/(vol*binwidth)

    axs[ii].errorbar(bincen, np.log10(phi), yerr = [np.log10(phi)-np.log10(phi-yerr), np.log10(phi+yerr)-np.log10(phi)], xerr = xerr, label = 'FLARES', color = 'brown', alpha = 0.5)

    for kk, jj in enumerate(obs):

        this = eval(F'evo.{jj}()')

        if z in this.redshifts:

            ok = np.where(np.array(this.redshifts) == z)[0][0]

            theta = [this.M_star[ok], this.phi_star[ok], this.alpha[ok]]

            if jj == 'FLARES':
                label = 'FLARES fit'
            else:
                label = this.name

            axs[ii].plot(tbincen, model(tbincen, theta), label = label, color = obs_color[kk])

    plot_tng(tbincen, z, axs[ii])
    plot_obs(z, axs[ii])

    axs[ii].set_ylim(-8, -0.5)
    axs[ii].set_xlim(-24.5, -15.5)
    axs[ii].grid()
    axs[ii].legend(loc=2, frameon=False)

fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
fig.savefig('fit_compare.pdf')
plt.show()
