"""

    Fit UV LF for all redshifts and all simulations. Use mpi
    for parallelising different redshifts

"""

import sys
import numpy as np
import pandas as pd
import scipy
import matplotlib
matplotlib.use('Agg')
from functools import partial
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

from emcee.autocorr import integrated_time

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

import flares
from modules import get_lum_all, Schechter, DPL, fit_function
# import fit_bootstrap
from mpi4py import MPI
import seaborn as sns
sns.set_context("paper")


model = str(sys.argv[1])
zs = [5., 6., 7., 8., 9., 10.]
fl = flares.flares('./data/flares.hdf5')
h = 0.6777
vol = (4/3)*np.pi*(14/h)**3
refvol = 100**3
AGNdT9vol = 50**3
filters = 'FUV'

if model == 'Schechter':
    folder = 'fit_Sch'
elif model == 'DPL':
    folder = 'fit_DPL'

tags = fl.tags[::-1]

def fitdf(ii, tag, N_up, N, V, bins, model):

    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    #obs = [{'bin_edges': bins, 'N': N_up, 'volume': V, 'Ntot': N}]
    obs = [{'bin_edges': bins, 'N': N_up, 'volume': V, 'sigma': N_up/np.sqrt(N)}]
    # print("upscaled N", N_up,
    #       "\nN", N,
    #       "\nsigma", obs[0]['sigma'])
    print (obs)
    priors = {}

    if model == 'Schechter':
        model = models.Schechter_Mags()
        folder = 'fit_Sch'
        # out = fit_bootstrap.fitter(tag)
        # print (out)
        priors['log10phi*'] = scipy.stats.uniform(loc=-5.5, scale=3.0)
        priors['alpha'] = scipy.stats.uniform(loc=-3.2, scale=2)
        priors['D*'] = scipy.stats.uniform(loc = -23, scale = 4.0)
        # if tag == '010_z005p000':
        #     priors['log10phi*'] = scipy.stats.uniform(loc=out[0]-0.5, scale=1.)
        #     priors['alpha'] = scipy.stats.uniform(loc=out[1]-0.5, scale=1.)
        #     priors['D*'] = scipy.stats.uniform(loc=out[2]-0.5, scale=1.)
        # if tag == '006_z009p000':
        #     priors['D*'] = scipy.stats.uniform(loc=-20.95, scale=0.6)
        # if tag == '005_z010p000':
        #     priors['D*'] = scipy.stats.uniform(loc=-20.5, scale=0.8)


    elif model == 'DPL':
        model = models.DPL_Mags()
        folder = 'fit_DPL'
        priors['log10phi*'] = scipy.stats.uniform(loc=-5.5, scale=3.0)
        priors['alpha_1'] = scipy.stats.uniform(loc=-3.5, scale=2.)
        priors['alpha_2'] = scipy.stats.uniform(loc=-5.3, scale=3.)
        priors['D*'] = scipy.stats.uniform(loc = -23, scale = 4.)
        # if tag == '006_z009p000':
        #     priors['D*'] = scipy.stats.uniform(loc = -21.0, scale = 1.3)
        # if tag == '005_z010p000':
        #     priors['D*'] = scipy.stats.uniform(loc = -20.8, scale = 1.3)
    else:
        raise ValueError("Model not recognized")


    fitter = fitDF.fitter(obs, model=model, priors=priors, output_directory=folder)
    fitter.lnlikelihood = fitter.gaussian_lnlikelihood
    samples = fitter.fit(nsamples=int(6e4), burn=int(5e4), sample_save_ID=f'{folder}_z{int(zs[ii])}', use_autocorr=False, verbose=True)

    return samples


def fit(ii, tag, bins, model):

    print('tag:',tag)
    V = (3200)**3
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    out, hist_all, err = get_lum_all(tags[ii], bins = bins)
    ok = np.where(hist_all<=5)[0]
    out[ok] = 0.
    hist_all[ok] = 0.
    err[ok] = 0.

    phi = out/(vol*binwidth)
    err = err/(vol*binwidth)

    N = models.phi_to_N(phi,V,bins)

    samples = fitdf(ii, tag, N, hist_all, V, bins, model)

    return None


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0: print (F"MPI size = {size}")


for ii, z in enumerate(zs):
    # ii = ii
    if rank == ii:
        df = pd.read_csv('Magnitude_limits.txt')
        low = np.array(df[filters])[ii]

        bins = -np.arange(-low, 26, 0.5)[::-1]

        fit(ii, tags[ii], bins, model)


if rank == 0:

    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(8, 7), sharex=True, sharey=True, facecolor='w', edgecolor='k')

    norm = matplotlib.colors.Normalize(vmin=0.5, vmax=len(zs)+0.5)

    # choose a colormap
    c_m = matplotlib.cm.viridis_r

    # create a ScalarMappable and initialize a data structure
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    for ii, z in enumerate(zs):
        df = pd.read_csv('Magnitude_limits.txt')
        low = np.array(df[filters])[ii]
        bins = -np.arange(-low, 26, 0.5)[::-1]
        bincen = (bins[1:]+bins[:-1])/2.
        binwidth = bins[1:] - bins[:-1]
        parent_volume = (3200)**3

        out, hist, err = get_lum_all(tags[ii], bins=bins)

        Msim = out/(binwidth*vol)
        xerr = np.ones(len(out))*binwidth[0]/2.
        yerr = err/(vol*binwidth)
        mask = np.where(hist==1)[0]
        uplims = np.zeros(len(bincen))
        uplims[mask] = True
        y_lo = np.log10(Msim)-np.log10(Msim-yerr)
        y_up =  np.log10(Msim+yerr)-np.log10(Msim)
        y_lo[mask] = 4.


        observed = Msim*(binwidth*parent_volume)
        sigma = observed/np.sqrt(hist)

        yy = fit_function(model, observed, sigma, bins, z)

        axs.errorbar(bincen, np.log10(Msim), yerr=[y_lo, y_up], uplims=uplims, ls='', marker='o', color=s_m.to_rgba(ii+0.5))
        axs.plot(bincen, yy, lw = 2, ls = 'solid', color=s_m.to_rgba(ii+0.5))

    axs.grid(True, alpha=0.6)
    axs.set_xlim(-16.7, -25.2)
    axs.set_ylim(-9, -1.9)
    axs.set_ylabel(r'$\mathrm{log}_{10}(\Phi/(\mathrm{cMpc}^{-3}\mathrm{Mag}^{-1}))$', fontsize=22)
    axs.set_xlabel(r'$\mathrm{M}_{1500}$', fontsize=22)
    axs.minorticks_on()
    axs.tick_params(axis='x', which='minor', direction='in')
    axs.tick_params(axis='y', which='minor', direction='in')

    cbaxes = fig.add_axes([0.15, 0.25, 0.03, 0.3])
    fig.colorbar(s_m, cax=cbaxes)
    cbaxes.set_ylabel(r'$z$', fontsize = 24)
    cbaxes.set_yticks(np.arange(len(zs)))
    cbaxes.set_yticklabels(zs)
    cbaxes.invert_yaxis()
    for label in (axs.get_xticklabels() + axs.get_yticklabels() + cbaxes.get_yticklabels()):
        label.set_fontsize(22)


    plt.savefig(F'{folder}/fit_gaussian.pdf', bbox_inches='tight')
    plt.close()
