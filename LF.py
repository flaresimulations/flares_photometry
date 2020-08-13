"""

    Plots involving luminosity functions (Give the numbers
    as arguments to the call):
        0 - Figure 4
        1 - Figure 5
        2 - Figure 7

"""


import sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from functools import partial

from FLARE.LF import evo
from plot_obs import plot_UVLF
from modules import get_lum_all, bluetides_uvlf, Schechter, DPL, plot_tng, fit_function

import seaborn as sns
sns.set_context("paper")

h = 0.6777
parent_volume = 3200**3
vol = (4/3)*np.pi*(14/h)**3
refvol = 100**3
AGNdT9vol = 50**3
filters = 'FUV'

plt_options = ['Flares-Eagle', 'Flares-fit', 'Observations+Models']
function = 'Schechter'
function = 'DPL'
observations = ['Atek+2018', 'Bouwens+2015', 'Bouwens+2016', 'Bouwens+2017', 'Bowler+2020', 'Finkelstein+2015', 'McLeod+2015', 'Oesch+2018', 'Stefanon+2019']
sims = ['\\textsc{BlueTides}', '\\textsc{Illustris} \\textsc{Tng}', 'Ma+2019', 'Mason+2015', 'Yung+2018']
sims_sub = ['Ma2019', 'Mason15', 'Yung2018']
sims_label = ['Ma+2019', 'Mason+2015', 'Yung+2018']
sims_sub_style = [(0, (5, 1)), (0, (3, 1, 1, 1)), 'solid']

input = plt_options[int(sys.argv[1])]

zs = [5., 6., 7., 8., 9., 10.]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']


if input==plt_options[1]:
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(4, 6), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
elif input==plt_options[2]:
    fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(14, 9), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
    function=''
else:
    fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(13, 6), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
    function=''


norm = matplotlib.colors.Normalize(vmin=0.5, vmax=len(zs)+0.5)

# choose a colormap
c_m = matplotlib.cm.viridis_r

# create a ScalarMappable and initialize a data structure
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])
# dat = pd.DataFrame({})
for ii, z in enumerate(zs):
    df = pd.read_csv('Magnitude_limits.txt')
    low = np.array(df[filters])[ii]

    bins = -np.arange(-low, 26.5, 0.5)[::-1]
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    print ("\n tag: ", tags[ii])
    out, hist, err = get_lum_all(tags[ii], bins=bins)

    Msim = out/(binwidth*vol)
    #out_low, out_up = out_low/(binwidth*vol), out_up/(binwidth*vol)
    xerr = np.ones(len(out))*binwidth[0]/2.
    yerr = err/(vol*binwidth)
    ok = np.where(hist > 0)[0]

    observed = Msim*(binwidth*parent_volume)
    sigma = observed/np.sqrt(hist)
    # xx, yy, zz = np.zeros(15), np.zeros(15), np.zeros(15)
    # num = len(ok)
    # xx[:num], yy[:num], zz[:num] = bincen[ok], Msim[ok], yerr[ok]
    # dat[f"M{z}"] = xx
    # dat[f"phi{z}"] = yy
    # dat[F"err{z}"] = zz

    if input != plt_options[1]:

        axs[ii].text(-17.6, -6.5, r'$z = {}$'.format(z), fontsize = 14)

        for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
            label.set_fontsize(14)

        axs[ii].set_xlim((-16.9, -24.7))
        axs[ii].set_ylim((-9.3, -1.7))
        axs[ii].grid(True, alpha=0.4)
        axs[ii].minorticks_on()
        axs[ii].tick_params(axis='x', which='minor', direction='in')
        axs[ii].tick_params(axis='y', which='minor', direction='in')

    if input == plt_options[0]:

        if ii==0: M5, phi5 = bincen, Msim
        ok1 = np.where(hist[ok] >= 5)[0][0]
        Mintr, tmp, tmp1 = get_lum_all(tags[ii], bins=bins, Luminosity='Intrinsic')/(binwidth*vol)
        Mref = get_lum_all(tags_ref[ii], bins=bins, inp='REF')/(binwidth*refvol)
        okref = np.where(Mref>0)[0]
        # Magn = get_lum_all(tags_ref[ii], bins=bins, inp='AGNdT9')/(binwidth*AGNdT9vol)
        # okagn = np.where(Magn>0)[0]
        if z!=5:
            axs[ii].plot(M5, np.log10(phi5), lw=2, alpha=0.7, ls='dashed', color=s_m.to_rgba(0.5))

        if z==10:
            axs[ii].plot(bincen[ok][ok1:], np.log10(Msim[ok][ok1:]), lw=2, label=r'\textsc{Flares}', color='black', alpha=0.7)
            axs[ii].plot(bincen, np.log10(Mintr), lw=2, ls='dotted', label=r'\textsc{Flares} intrinsic',  color='black', alpha=0.7)
            axs[ii].plot(bincen[okref], np.log10(Mref[okref]), color='crimson', lw=2, ls='dashed', label=r'\textsc{Eagle} Ref')
            # axs[ii].plot(bincen[okagn], np.log10(Magn[okagn]), color='brown', lw=2, ls='dashed', label=r'\textsc{Eagle} AGNdT9')
            axs[ii].legend(frameon=False, fontsize = 12, numpoints=1, ncol = 2)

        axs[ii].plot(bincen[ok][ok1:], np.log10(Msim[ok][ok1:]), lw = 2, color=s_m.to_rgba(ii+0.5))
        axs[ii].plot(bincen[ok], np.log10(Msim[ok]), lw = 2, ls = 'dashed', color=s_m.to_rgba(ii+0.5))
        axs[ii].fill_between(bincen[ok], np.log10(Msim[ok]-yerr[ok]), np.log10(Msim[ok]+yerr[ok]), alpha=0.5, color=s_m.to_rgba(ii+0.5))
        axs[ii].plot(bincen, np.log10(Mintr), lw=2, ls='dotted', color=s_m.to_rgba(ii+0.5))
        axs[ii].plot(bincen[okref], np.log10(Mref[okref]), color='crimson', lw=2, ls='dashed')
        # axs[ii].plot(bincen[okagn], np.log10(Magn[okagn]), color='brown', lw=2, ls='dashed')


    else:
        mask = np.where(hist==1)[0]
        lolims = np.zeros(len(bincen))
        lolims[mask] = True

        y_lo = np.log10(Msim)-np.log10(Msim-yerr)
        y_up =  np.log10(Msim+yerr)-np.log10(Msim)
        #y_lo[mask] = (np.log10(Msim)-np.log10(Msim-out_low))[mask]
        #y_up[mask] =  (np.log10(Msim+out_up)-np.log10(Msim))[mask]
        #y_lo[mask] = np.log10(Msim[mask])-np.log10(out_low[mask])
        #y_up[mask] = np.log10(out_up[mask])-np.log10(Msim[mask])
        y_lo[mask] = 0.5
        y_up[mask] = 0.5
        ok = np.where(hist<5)[0]
        observed[ok]=0

    if input == plt_options[1]:
        print ("\n tag:", tags[ii])
        axs[0].errorbar(bincen, np.log10(Msim), yerr=[y_lo, y_up], lolims=lolims, ls='', marker='o', color=s_m.to_rgba(ii+0.5))
        axs[0].plot(bincen, fit_function('Schechter', observed, sigma, bins, z), lw=2, ls='solid', color=s_m.to_rgba(ii+0.5))
        axs[1].errorbar(bincen, np.log10(Msim), yerr=[y_lo, y_up], lolims=lolims, ls='', marker='o', color=s_m.to_rgba(ii+0.5))
        axs[1].plot(bincen, fit_function('DPL', observed, sigma, bins, z), lw=2, ls='solid', color=s_m.to_rgba(ii+0.5))


    elif input == plt_options[2]:
        axs[ii].errorbar(bincen, np.log10(Msim), yerr=[y_lo, y_up], lolims=lolims, ls='', marker='o', color=s_m.to_rgba(ii+0.5), zorder=100, alpha=0.75, markersize=8)
        # axs[ii].plot(bincen, fit_function('Schechter', observed, sigma, bins, z), lw=2, ls='solid', color=s_m.to_rgba(ii+0.5), zorder=100, alpha=0.75)
        # axs[ii].plot(bincen, fit_function('DPL', observed, sigma, bins, z), lw=2, ls='dashed', color=s_m.to_rgba(ii+0.5), zorder=100, alpha=0.75)

        plot_UVLF(int(z), axs[ii])
        if z <= 8.:
            plot_tng(bincen, z, axs[ii])

        if z >= 8.:
            binBt, MBt = bluetides_uvlf(int(z))
            axs[ii].plot(binBt, MBt, color='blue', lw=2, ls='solid', label=r'\textsc{BlueTides}', alpha=0.4)

        for kk, ll in enumerate(sims_sub):

            this = eval(F'evo.{ll}()')
            if z in this.redshifts:

                ok = np.where(np.array(this.redshifts) == z)[0][0]

                theta = [this.phi_star[ok], this.alpha[ok], this.M_star[ok]]

                axs[ii].plot(bincen, Schechter(bincen, theta), label=sims_label[kk], color='grey', ls=sims_sub_style[kk], lw=2, alpha=0.5)

if input == plt_options[1]:
    for ll in [0,1]:
        axs[ll].grid(True, alpha=0.6)
        axs[ll].set_xlim(-16.9, -24.7)
        axs[ll].set_ylim(-9.4, -1.7)
        #axs[ll].set_ylabel(r'$\mathrm{log}_{10}(\Phi/(\mathrm{cMpc}^{-3}\mathrm{Mag}^{-1}))$', fontsize=14)
        axs[ll].minorticks_on()
        axs[ll].tick_params(axis='x', which='minor', direction='in')
        axs[ll].tick_params(axis='y', which='minor', direction='in')
        for label in (axs[ll].get_xticklabels() + axs[ll].get_yticklabels()):
            label.set_fontsize(10)

    axs[1].set_xlabel(r'$\mathrm{M}_{1500}$', fontsize=12)
    fig.subplots_adjust(left=0.14, wspace=0, hspace=0)
    fig.text(0.02, 0.5, r'$\mathrm{log}_{10}(\Phi/(\mathrm{cMpc}^{-3}\mathrm{Mag}^{-1}))$', va='center', rotation='vertical', fontsize=12)
    cbaxes = fig.add_axes([0.2, 0.13, 0.015, 0.25])
    fig.colorbar(s_m, cax=cbaxes)
    cbaxes.set_ylabel(r'$z$', fontsize = 16)
    cbaxes.set_yticks(np.arange(len(zs)))
    cbaxes.set_yticklabels(zs)
    cbaxes.invert_yaxis()
    for label in (cbaxes.get_yticklabels()):
        label.set_fontsize(10)

if input != plt_options[1]:

    fig.subplots_adjust(bottom=0.09, left=0.08, wspace=0, hspace=0)

    fig.text(0.03, 0.5, r'$\mathrm{log}_{10}(\Phi/(\mathrm{cMpc}^{-3}\mathrm{Mag}^{-1}))$', va='center', rotation='vertical', fontsize=18)
    fig.text(0.48, 0.01, r'$\mathrm{M}_{1500}$', va='center', fontsize=18)

if input == plt_options[2]:
    lines, labels = [], []

    for ax in fig.axes:
        axLine, axLabel = ax.get_legend_handles_labels()
        lines.extend(axLine)
        labels.extend(axLabel)

    unique = np.unique(np.array(labels))
    unique_lines = []
    for ii, jj in enumerate(unique):
        unique_lines.append(lines[np.where(np.array(labels)==jj)[0][0]])

    observation_lines, sim_lines=[], []
    for ii, jj in enumerate(observations):
        observation_lines.append(unique_lines[np.where(np.array(unique)==jj)[0][0]])

    for ii, jj in enumerate(sims):
        sim_lines.append(unique_lines[np.where(np.array(unique)==jj)[0][0]])


    axs[0].legend(observation_lines[:3], list(observations)[:3], frameon=False, fontsize=14, numpoints=1, ncol=1, loc='lower left')
    axs[1].legend(observation_lines[3:6], list(observations)[3:6], frameon=False, fontsize=14, numpoints=1, ncol=1, loc='lower left')
    axs[2].legend(observation_lines[6:], list(observations)[6:], frameon=False, fontsize=14, numpoints=1, ncol=1, loc='lower left')
    axs[3].legend(sim_lines[:2], list(sims)[:2], frameon=False, fontsize=14, numpoints=1, ncol=1, loc='lower left')
    axs[4].legend(sim_lines[2:], list(sims)[2:], frameon=False, fontsize=14, numpoints=1, ncol=1, loc='lower left')

plt.savefig(F"LF_FUV_z5_10_{input}.pdf", bbox_inches='tight', dpi=300)
plt.show()
