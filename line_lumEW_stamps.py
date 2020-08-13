"""

    Plots the line LF and EW relationship of the 6 prominent
    nebular lines:
        CIII]1907, 1909
        [OII]3726, 3729
        [NeIII]3869, 3967
        Hbeta
        [OIII]4959, 5007
        Halpha

    The following arguments along with the script for the following plots:
        0 - line LF (Top plot in Figure 12)
        1 - EW vs Mstar (Middle plot in Figure 12)
        2 - EW vs UV luminosity (Bottom plot in Figure 12)
    
"""


import numpy as np
import pandas as pd
import matplotlib, sys
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from FLARE.photom import lum_to_M, M_to_lum
import FLARE.filters
import flares as fl
from modules import get_data_all, get_line_all, get_lum_all
import seaborn as sns
sns.set_context("paper")



filters = 'FUV'
conversion_fac = 2E15  #converting from ergs/s/Hz to ergs/s at FUV
h = 0.6777
vol = (4/3)*np.pi*(14/h)**3
quantiles = [0.84,0.50,0.16]

def LF_create(dataset, weights, l_fuvs, low, axs, title, color):

    ok = np.where(l_fuvs>low)
    low = np.min(np.concatenate(dataset[:,0])[ok])

    bins = np.arange(39, 47, 0.4)
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]
    LF_tot = np.zeros(len(bincen))
    err = np.zeros(len(bincen))
    LF = np.zeros(len(bincen))

    for ii in range(len(weights)):
        tmp, edges = np.histogram(np.log10(dataset[:,0][ii]), bins = bins)
        err+=np.square(np.sqrt(tmp)*weights[ii])
        LF+=tmp*weights[ii]
        LF_tot+=tmp

    LF = LF/(binwidth*vol)
    err = np.sqrt(err)/(binwidth*vol)

    axs.plot(bincen, np.log10(LF), ls = 'solid', color=color)
    axs.fill_between(bincen, np.log10(LF-err), np.log10(LF+err), alpha=0.3, color=color)
    axs.set_title(rF"{title}", fontsize=13)


def create_EW_plot(x, y, bins, weights, axs, title, color, ii):

    quantiles = [0.84,0.50,0.16]

    ok = np.where(10**y>0)
    out = fl.binned_weighted_quantile(x[ok],y[ok],weights[ok],bins,quantiles)

    bincen = (bins[1:]+bins[:-1])/2.
    yy, yy16, yy84 = out[:,1], out[:,0],out[:,2]
    # axs.fill_between(bincen, yy16, yy84, color=color, alpha=0.3)
    axs.scatter(bincen, yy, color=color, alpha=1.0/(ii/10 +1), s=100)
    axs.locator_params(axis='y', nbins=5)
    axs.set_title(rF"{title}", fontsize=13)

zs = [5., 6., 7., 8., 9., 10.]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']

plt_options = ['LFlines', 'EWlinesMstar', 'EWlineslum']

norm = matplotlib.colors.Normalize(vmin=0.5, vmax=len(zs)+0.5)
# choose a colormap
c_m = matplotlib.cm.viridis_r
# create a ScalarMappable and initialize a data structure
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

titles = ["H$\\alpha$", "H$\\beta$", "CIII]1907, 1909", "[OII]3726, 3729", "[OIII]4959, 5007"]

input = plt_options[int(sys.argv[1])]

if input == plt_options[0]:
    fig, axs = plt.subplots(nrows = 1, ncols = 5, figsize=(15,3), sharex=False, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()

    xlims = [[41.7,44.2], [41.2, 43.7], [40.8,43.3], [40.4,42.9],[41.7,44.1]]
    # [40.7,43]
    df = pd.read_csv('weight_files/weights_grid.txt')
    weights = np.array(df['weights'])
    sims = np.arange(len(weights))


    for ii, tag in enumerate(tags):

        z = float(tag[5:].replace('p','.'))
        print (z)
        df = pd.read_csv('Magnitude_limits.txt')
        low = np.array(df[filters])[ii]

        Halpha = get_line_all(tag, 'HI6563', inp = 'FLARES', LF = False)
        Hbeta = get_line_all(tag, 'HI4861', inp = 'FLARES', LF = False)
        CIII = get_line_all(tag, 'CIII1907', inp = 'FLARES', LF = False) + get_line_all(tag, 'CIII1909', inp = 'FLARES', LF = False)
        OII = get_line_all(tag, 'OII3726', inp = 'FLARES', LF = False) + get_line_all(tag, 'OII3729', inp = 'FLARES', LF = False)
        # NeIII = get_line_all(tag, 'NeIII3869', inp = 'FLARES', LF = False) + get_line_all(tag, 'NeIII3967', inp = 'FLARES', LF = False)
        OIII = get_line_all(tag, 'OIII4959', inp = 'FLARES', LF = False) + get_line_all(tag, 'OIII5007', inp = 'FLARES', LF = False)
        l_fuvs = np.concatenate(get_lum_all(tag, LF=False))

        LF_create(Halpha, weights, l_fuvs, low, axs[0], titles[0], s_m.to_rgba(ii+0.5))
        LF_create(Hbeta, weights, l_fuvs, low, axs[1], titles[1], s_m.to_rgba(ii+0.5))
        LF_create(CIII, weights, l_fuvs, low, axs[2], titles[2], s_m.to_rgba(ii+0.5))
        LF_create(OII, weights, l_fuvs, low, axs[3], titles[3], s_m.to_rgba(ii+0.5))
        # LF_create(NeIII, weights, l_fuvs, low, axs[4], titles[4], s_m.to_rgba(ii+0.5))
        LF_create(OIII, weights, l_fuvs, low, axs[4], titles[4], s_m.to_rgba(ii+0.5))


    for ii in range(5):
        axs[ii].set_xlim(xlims[ii])
        axs[ii].set_ylim(-9.,-1.8)
        for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
            label.set_fontsize(13)


    axs[0].set_ylabel(r'log$_{10}$($\phi$/(cMpc$^{-3}$dex$^{-1}$))', fontsize=15)
    fig.subplots_adjust(bottom=0.14, left=0.08, right=0.92,  wspace=0, hspace=0)
    fig.text(0.45, 0.001, r'log$_{10}$(L/(erg s$^{-1}$))', va='center', fontsize=15)
    cbaxes = fig.add_axes([0.925, 0.17, 0.007, 0.65])
    fig.colorbar(s_m, cax=cbaxes)
    cbaxes.set_ylabel(r'$z$', fontsize = 20)
    cbaxes.set_yticks(np.arange(len(zs)))
    cbaxes.set_yticklabels(zs)
    cbaxes.invert_yaxis()
    for label in (cbaxes.get_yticklabels()):
        label.set_fontsize(13)


if input != plt_options[0] :
    fig, axs = plt.subplots(nrows = 1, ncols = 5, figsize=(15,2.5), sharex=False, sharey=False, facecolor='w', edgecolor='k')
    axs = axs.ravel()
    # xlims =
    df = pd.read_csv('weight_files/weights_grid.txt')
    weights = np.array(df['weights'])
    sims = np.arange(len(weights))

    for ii, tag in enumerate(tags):

        z = float(tag[5:].replace('p','.'))
        print (z)
        df = pd.read_csv('Magnitude_limits.txt')
        low = np.array(df[filters])[ii]

        Halpha = np.concatenate(get_line_all(tag, 'HI6563', inp = 'FLARES', LF = False)[:,1])
        Hbeta = np.concatenate(get_line_all(tag, 'HI4861', inp = 'FLARES', LF = False)[:,1])
        CIII = np.concatenate(get_line_all(tag, 'CIII1907', inp = 'FLARES', LF = False)[:,1] + get_line_all(tag, 'CIII1909', inp = 'FLARES', LF = False)[:,1])
        OII = np.concatenate(get_line_all(tag, 'OII3726', inp = 'FLARES', LF = False)[:,1] + get_line_all(tag, 'OII3729', inp = 'FLARES', LF = False)[:,1])
        # NeIII = get_line_all(tag, 'NeIII3869', inp = 'FLARES', LF = False) + get_line_all(tag, 'NeIII3967', inp = 'FLARES', LF = False)
        OIII = np.concatenate(get_line_all(tag, 'OIII4959', inp = 'FLARES', LF = False)[:,1] + get_line_all(tag, 'OIII5007', inp = 'FLARES', LF = False)[:,1])
        if input == plt_options[1]:
            x_inp = get_data_all(tag, dataset = 'Mstar_30', DF=False)*1e10
            for kk in range(len(x_inp)): x_inp[kk] = np.log10(x_inp[kk])
            bins = np.arange(7.5, 11.5, 0.5)
            xlabel=r'log$_{10}$(M$_{\star}$/M$_{\odot}$)'
            xpos = 0.455
            bottom = 0.20

        elif input == plt_options[2]:
            x_inp = get_lum_all(tag, LF=False)
            for kk in range(len(x_inp)): x_inp[kk] = lum_to_M(x_inp[kk])
            bins = -np.arange(17, 25, 1)[::-1]
            xlabel = r'$\mathrm{M}_{1500}$'
            xpos = 0.47
            bottom = 0.21

            for kk in range(5): axs[kk].set_xlim((-17,-23.9))


        ws = np.array([])
        for jj in sims:
            ws = np.append(ws, np.ones(np.shape(x_inp[jj]))*weights[jj])
        x_inp = np.concatenate(x_inp)

        create_EW_plot(x_inp, np.log10(Halpha), bins, ws, axs[0], titles[0], s_m.to_rgba(ii+0.5), ii)
        create_EW_plot(x_inp, np.log10(Hbeta), bins, ws, axs[1], titles[1], s_m.to_rgba(ii+0.5), ii)
        create_EW_plot(x_inp, np.log10(CIII), bins, ws, axs[2], titles[2], s_m.to_rgba(ii+0.5), ii)
        create_EW_plot(x_inp, np.log10(OII), bins, ws, axs[3], titles[3], s_m.to_rgba(ii+0.5), ii)
        create_EW_plot(x_inp, np.log10(OIII), bins, ws, axs[4], titles[4], s_m.to_rgba(ii+0.5), ii)

    axs[0].set_ylabel(r'log$_{10}$(EW/$\AA$)', fontsize=15)
    fig.subplots_adjust(bottom=bottom, left=0.08)
    fig.text(xpos, 0.005, xlabel, va='center', fontsize=17)
    # cbaxes = fig.add_axes([0.95, 0.14, 0.004, 0.45])
    # fig.colorbar(s_m, cax=cbaxes)
    # cbaxes.set_ylabel(r'$z$', fontsize = 16)
    # cbaxes.set_yticks(np.arange(len(zs)))
    # cbaxes.set_yticklabels(zs)
    # cbaxes.invert_yaxis()
    # for label in (cbaxes.get_yticklabels()):
    #     label.set_fontsize(14)


for ii in range(5):
    for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
        label.set_fontsize(13)

    axs[ii].grid(True, alpha=0.4)
    axs[ii].minorticks_on()
    axs[ii].tick_params(axis='x', which='minor', direction='in')
    axs[ii].tick_params(axis='y', which='minor', direction='in')

plt.savefig(F'{input}_stamps.pdf', dpi=300, bbox_inches='tight')
plt.show()
