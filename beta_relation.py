"""

Figure 8, with argument 0

"""


import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import h5py
from FLARE.photom import lum_to_M
from plot_obs import plot_beta
import flares as fl

filters = 'FUV'
xlims = [-16.9, -24.7]
ylims = [-2.7,-1.4]

df = pd.read_csv('weight_files/weights_grid.txt')
weights = np.array(df['weights'])

zs = [5, 6, 7, 8, 9, 10]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']


plt_options = ['lum', 'Mstar']

input = plt_options[int(sys.argv[1])]

fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(13, 5), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

for ii, tag in enumerate(tags):

    df = pd.read_csv('Magnitude_limits.txt')
    low = np.array(df[filters])[ii]
    bins = -np.arange(-low, 25, 0.4)[::-1]
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    L = {}
    for f in ['FUV','NUV']:
        L[f] = np.array([])
        mstar, sfr = np.array([]), np.array([])
        w = np.array([])
        with h5py.File(f'data/flares.hdf5', 'r') as hf:
            for sim in hf.keys():
                tmp = np.array(hf[f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/{f}'])
                L[f] = np.hstack((L[f], tmp))
                w = np.append(w, np.ones(len(tmp))*weights[int(sim)])


    beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0

    y = beta
    if input == plt_options[0]:
        x = lum_to_M(L['FUV'])
    else:
        x = np.log10(get_data_all(tag, inp = 'FLARES', DF = False)*1e10)

    axs[ii].hexbin(x, y, gridsize=(35,13), bins='log', cmap='Greys_r', linewidths=0., mincnt=5, extent=[*[low,-24.5], *ylims], alpha=0.6, zorder=2)

    quantiles = [0.84,0.50,0.16]
    out = fl.binned_weighted_quantile(x, y, w, bins, quantiles)
    hist, binedges = np.histogram(x, bins)
    ok = np.where(hist>0)[0]
    ok1 = np.where(hist[ok]>3)[0][0]

    axs[ii].fill_between(bincen[ok][ok1:], out[:,2][ok][ok1:], out[:,0][ok][ok1:], color='black', alpha=0.5)
    axs[ii].plot(bincen[ok], out[:,1][ok], ls='dashed', color='black', alpha=.5, lw=2)
    axs[ii].plot(bincen[ok][ok1:], out[:,1][ok][ok1:], ls='-', color='black', alpha=.5, lw=2)

    # add observations
    plot_beta(zs[ii], axs[ii])

    axs[ii].set_ylim(ylims)
    axs[ii].set_xlim(xlims)
    axs[ii].minorticks_on()
    axs[ii].tick_params(axis='x', which='minor', direction='in')
    axs[ii].tick_params(axis='y', which='minor', direction='in')
    axs[ii].grid(True, alpha = 0.5)
    axs[ii].text(-20.1, -2.6, r'$z = {}$'.format(zs[ii]), fontsize = 12)

axs[1].legend(frameon=False,fontsize=10,loc=4)

fig.subplots_adjust(bottom=0.09, left = 0.05, wspace=0, hspace=0)
fig.text(0.005, 0.5, r'$\beta$', va='center', rotation='vertical', fontsize=15)
fig.text(0.48, 0.02, r'M$_{1500}$', va='center', fontsize=15)
fig.savefig(F'beta_{input}_z5_10.pdf', dpi = 300, bbox_inches='tight')


plt.show()
