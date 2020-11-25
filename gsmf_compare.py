"""

Comparing GSMF for all redshifts and all simulations; Figure 1

"""

import sys
import numpy as np
import pandas as pd
import scipy
import matplotlib
from functools import partial
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt

import flares
from modules import get_data_all, bluetides_gsmf

import seaborn as sns
sns.set_context("paper")

h = 0.6777
vol = (4/3)*np.pi*(14/h)**3
refvol = 100**3
AGNdT9vol = 50**3

zs = [5., 6., 7., 8., 9., 10.]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']

fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(9, 4), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

norm = matplotlib.colors.Normalize(vmin=0.5, vmax=len(zs)+0.5)

# choose a colormap
c_m = matplotlib.cm.viridis_r

# create a ScalarMappable and initialize a data structure
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

for ii, tag in enumerate(tags):

    axs[ii].text(7.9, -6.5, r'$z = {}$'.format(zs[ii]), fontsize = 10)

    bins = np.arange(7.5, 12.5, 0.25)
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    out, hist, err = get_data_all(tag, bins = bins, inp = 'FLARES', DF = True)

    ok = np.where(hist > 0)[0]
    hist = hist[ok]
    Msim = out/(binwidth*vol)
    xerr = np.ones(len(out))*binwidth[0]/2.
    yerr = err/(vol*binwidth)

    ok1 = np.where(hist >= 5)[0][-1]

    Mref = get_data_all(tags_ref[ii], bins = bins, inp = 'REF', DF = True)/(binwidth*refvol)
    okref = np.where(Mref>0)
    # Magn = get_data_all(tags_ref[ii], bins = bins, inp = 'AGNdT9', DF = True)/(binwidth*AGNdT9vol)
    # okagn = np.where(Magn>0)
    # if zs[ii]>=8:
    #     binBt, MBt = bluetides_gsmf(zs[ii])


    if zs[ii] == 10:

        axs[ii].plot(bincen[ok][:ok1+1], np.log10(Msim[ok][:ok1+1]), color = 'black', lw = 2, label = r'\textsc{Flares}')
        axs[ii].plot(bincen[ok][:ok1+1], np.log10(Msim[ok][:ok1+1]), color = s_m.to_rgba(ii+0.5), lw = 2, )
        axs[ii].plot(bincen[ok], np.log10(Msim[ok]), color = s_m.to_rgba(ii+0.5), lw = 2, ls = 'dashed')
        axs[ii].fill_between(bincen[ok], np.log10(Msim[ok]-yerr[ok]), np.log10(Msim[ok]+yerr[ok]), color=s_m.to_rgba(ii+0.5), alpha=0.3)
        axs[ii].plot(bincen[okref], np.log10(Mref[okref]), color = 'crimson', lw = 2, ls = 'dashed', label = r'\textsc{Eagle} $\mathrm{Ref}$')
        # axs[ii].plot(bincen[okagn], np.log10(Magn[okagn]), color = 'brown', lw =2, ls = 'dashed', label = r'\textsc{Eagle} $\mathrm{AGNdT9}$')
        # axs[ii].plot(binBt, MBt, color = 'blue', lw =2, ls = 'dashed', label = r'\textsc{BlueTides}')

    else:

        axs[ii].plot(bincen[ok][:ok1+1], np.log10(Msim[:ok1+1]), lw = 2, color = s_m.to_rgba(ii+0.5))
        axs[ii].plot(bincen[ok], np.log10(Msim[ok]), lw = 2, color = s_m.to_rgba(ii+0.5), ls = 'dashed')
        axs[ii].fill_between(bincen[ok], np.log10(Msim[ok]-yerr[ok]), np.log10(Msim[ok]+yerr[ok]), color=s_m.to_rgba(ii+0.5), alpha=0.3)
        axs[ii].plot(bincen[okref], np.log10(Mref[okref]), color = 'crimson', lw = 2, ls = 'dashed')
        # axs[ii].plot(bincen[okagn], np.log10(Magn[okagn]), color = 'brown', lw = 2, ls = 'dashed')
        # if zs[ii]>=8:
        #     axs[ii].plot(binBt, MBt, color = 'blue', lw =2, ls = 'dashed')

    for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
        label.set_fontsize(11)

    axs[ii].set_xlim((7.6, 11.4))
    axs[ii].set_ylim((-8.2, -0.8))
    axs[ii].set_xticks(np.arange(8., 11.5, 1))
    axs[ii].grid(True, alpha = 0.4)
    axs[ii].legend(frameon=False, fontsize = 11, numpoints=1)
    axs[ii].minorticks_on()
    axs[ii].tick_params(axis='x', which='minor', direction='in')
    axs[ii].tick_params(axis='y', which='minor', direction='in')

fig.subplots_adjust(bottom=0.11, left = 0.09, wspace=0, hspace=0)
fig.text(0.03, 0.5, r'$\mathrm{log}_{10}(\Phi/(\mathrm{cMpc}^{-3}\mathrm{dex}^{-1}))$', va='center', rotation='vertical', fontsize=12)
fig.text(0.445, 0.02, r'$\mathrm{log}_{10}(\mathrm{M}_{\star}/M_{\odot})$', va='center', fontsize=12)
plt.savefig("GSMF_z5_10.pdf", bbox_inches='tight')
plt.show()
