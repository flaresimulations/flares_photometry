"""

Plots in the Appendices

"""


import numpy as np
import pandas as pd
from functools import partial
import matplotlib, sys, h5py, schwimmbad
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from FLARE.photom import lum_to_M, M_to_lum
import FLARE.filters
import flares as fl
from plot_obs import plot_beta, plot_OIIIEW_Mstar, plot_OIIIEW_UV, plot_OIIIlum_UV
from modules import get_data_all, get_lum_all, get_line_all

import seaborn as sns
sns.set_context("paper")

def get_data(ii, dataset, tag):

    num = str(ii)
    if ii/10 < 1: num = '0'+num
    sim = "./data1/flares.hdf5"
    num = num+'/'

    with h5py.File(sim, 'r') as hf:

        data = np.array(hf[num+tag+'/Galaxy'].get(dataset))

    return data

def get_all(dataset, tag):

    df = pd.read_csv('weight_files/weights_grid.txt')
    weights = np.array(df['weights'])
    sims = np.arange(0,len(weights))

    calc = partial(get_data, tag = tag, dataset = dataset)

    pool = schwimmbad.MultiPool(processes=12)
    dat = np.array(list(pool.map(calc, sims)))
    pool.close()

    return dat



plt_options = ['BC_plot_beta', 'BC_plot_EW', 'beta_Extcurves', 'EW_Extcurves', 'Att_Extcurves']
kappa_BCs =     [0.0, 0.001,  0.1,  0.25,   0.5,  0.75,  1.,  1.25,  1.5,  1.75,    2.]
kappa_ISMs =  [0.1925,0.0775,0.0475,0.0175,0.0075,0.0063,0.005,0.0075,0.0025,0.0025,0.0025,0.005,0.005]

input = int(sys.argv[1])

if input%2!=1:
    tag = '010_z005p000'

else:
    tag = '007_z008p000'


# choose a colormap
if input<2:
    c_m = matplotlib.cm.cividis
else:
    c_m = matplotlib.cm.coolwarm

norm = matplotlib.colors.Normalize(vmin=0, vmax=3.)
# create a ScalarMappable and initialize a data structure
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])


df = pd.read_csv('weight_files/weights_grid.txt')
weights = np.array(df['weights'])
quantiles = [0.84,0.50,0.16]
conversion_fac = 2E15  #converting from ergs/s/Hz to ergs/s at FUV


if input == 0:
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(6, 4), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    for ii, kappa_BC in enumerate(kappa_BCs):
        xlims = [-16.9, -24.7]
        ylims = [-2.6,-1.]

        bins = -np.arange(16.8, 25, 0.4)[::-1]
        bincen = (bins[1:]+bins[:-1])/2.
        binwidth = bins[1:] - bins[:-1]

        L = {}
        for f in ['FUV','NUV']:
            L[f] = np.array([])
            w = np.array([])
            with h5py.File(f'data1/flares.hdf5', 'r') as hf:
                for sim in hf.keys():
                    tmp = np.array(hf[f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI_{kappa_BC}/{f}'])
                    L[f] = np.hstack((L[f], tmp))
                    w = np.append(w, np.ones(len(tmp))*weights[int(sim)])

        beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0

        x, y = lum_to_M(L['FUV']), beta

        out = fl.binned_weighted_quantile(x, y, w, bins, quantiles)
        hist, binedges = np.histogram(x, bins)
        ok = np.where(hist>0)[0]
        ok1 = np.where(hist[ok]>3)[0][0]

        axs.fill_between(bincen[ok][ok1:], out[:,0][ok][ok1:], out[:,2][ok][ok1:], color=s_m.to_rgba(kappa_BC), alpha=0.25)
        axs.plot(bincen[ok], out[:,1][ok], ls='dashed', color=s_m.to_rgba(kappa_BC), alpha=1, lw=2)
        axs.plot(bincen[ok][ok1:], out[:,1][ok][ok1:], ls='-', color=s_m.to_rgba(kappa_BC), alpha=1, lw=2, label=r"$\kappa_{\mathrm{BC}}=%s$"%(str(kappa_BC)))

    plot_beta(5, axs)
    axs.set_ylim(ylims)
    axs.set_xlim(xlims)
    axs.minorticks_on()
    axs.tick_params(axis='x', which='minor', direction='in')
    axs.tick_params(axis='y', which='minor', direction='in')
    axs.grid(True, alpha = 0.5)
    for label in (axs.get_xticklabels()+axs.get_yticklabels()):
        label.set_fontsize(12)
    axs.set_ylabel(r'$\beta$', fontsize=16)
    axs.set_xlabel(r'M$_{1500}$', fontsize=16)
    axs.legend(frameon=False, fontsize=11, ncol=3)


elif input == 1:
    fig, axs = plt.subplots(nrows = 1, ncols = 3, figsize=(13, 3), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    mstar = get_data_all(tag, dataset = 'Mstar_30', DF=False)
    z = 8

    ws = np.array([])
    for jj in range(40):
        ws = np.append(ws, np.ones(np.shape(mstar[jj]))*weights[jj])
    mstar = np.concatenate(mstar)

    for ii, kappa_BC in enumerate(kappa_BCs):
        OIII4959_lum = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{kappa_BC}/OIII4959/Luminosity', tag))
        OIII4959_EW = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{kappa_BC}/OIII4959/EW', tag))

        OIII5007_lum = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{kappa_BC}/OIII5007/Luminosity', tag))
        OIII5007_EW = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{kappa_BC}/OIII5007/EW', tag))

        Hbeta_lum = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{kappa_BC}/HI4861/Luminosity', tag))
        Hbeta_EW = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{kappa_BC}/HI4861/EW', tag))

        l_fuvs = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Luminosity/DustModelI_{kappa_BC}/FUV', tag))

        x, y, w = np.log10(mstar*1e10), OIII4959_EW+OIII5007_EW+Hbeta_EW, ws
        thisok = np.where(y>0)
        x = x[thisok]
        y = np.log10(y[thisok])
        w = w[thisok]
        tbins = np.arange(min(x), max(x)+0.4, 0.4)
        tbincen = (tbins[1:]+tbins[:-1])/2.
        tbinwidth = tbins[1:] - tbins[:-1]
        out = fl.binned_weighted_quantile(x,y,w,tbins,quantiles)
        hist, binedges = np.histogram(x, tbins)
        xx, yy, yy16, yy84 = tbins, out[:,1], out[:,0],out[:,2]
        tok = np.where(hist>0)[0]
        tok1 = np.where(hist[tok]>3)[0][-1]
        axs[0].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], alpha=0.25, color=s_m.to_rgba(kappa_BC))
        axs[0].plot(xx[tok], yy[tok], ls='dashed', color=s_m.to_rgba(kappa_BC), alpha=0.8, lw=2)
        axs[0].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color=s_m.to_rgba(kappa_BC), alpha=0.8, lw=2, label=r"$\kappa_{\mathrm{BC}}=%s$"%(str(kappa_BC)))
        xlims = [7.5, 10.5]
        ylims = [1.8, 3.9]
        axs[0].set_ylim(ylims)
        axs[0].set_xlim(xlims)



        x, y, w = np.log10(l_fuvs*conversion_fac), OIII4959_EW+OIII5007_EW+Hbeta_EW, ws
        thisok = np.where(y>0)
        x = x[thisok]
        y = np.log10(y[thisok])
        w = w[thisok]
        tbins = np.arange(min(x), max(x)+0.3, 0.3)
        tbincen = (tbins[1:]+tbins[:-1])/2.
        tbinwidth = tbins[1:] - tbins[:-1]
        out = fl.binned_weighted_quantile(x,y,w,tbins,quantiles)
        hist, binedges = np.histogram(x, tbins)
        xx, yy, yy16, yy84 = tbins, out[:,1], out[:,0],out[:,2]
        tok = np.where(hist>0)[0]
        tok1 = np.where(hist[tok]>3)[0][-1]
        axs[1].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color=s_m.to_rgba(kappa_BC), alpha=0.25)
        axs[1].plot(xx[tok], yy[tok], ls='dashed', color=s_m.to_rgba(kappa_BC), alpha=1.0, lw=2)
        axs[1].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color=s_m.to_rgba(kappa_BC), alpha=.8, lw=2)
        xlims = [42.9, 45.]
        ylims = [1.8, 3.9]
        axs[1].set_ylim(ylims)
        axs[1].set_xlim(xlims)



        x, y, w = np.log10(l_fuvs*conversion_fac), np.log10(OIII4959_lum+OIII5007_lum+Hbeta_lum) - np.log10(l_fuvs*conversion_fac), ws
        thisok = np.where(y > -5)
        x = x[thisok]
        y = y[thisok]
        w = w[thisok]
        tbins = np.arange(min(x), max(x)+0.3, 0.3)
        tbincen = (tbins[1:]+tbins[:-1])/2.
        tbinwidth = tbins[1:] - tbins[:-1]
        hist, binedges = np.histogram(x, tbins)
        out = fl.binned_weighted_quantile(x,y,w,tbins,quantiles)
        xx, yy, yy16, yy84 = tbins, out[:,1], out[:,0],out[:,2]
        tok = np.where(hist>0)[0]
        tok1 = np.where(hist[tok]>3)[0][-1]
        axs[2].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color=s_m.to_rgba(kappa_BC), alpha=0.25)
        axs[2].plot(xx[tok], yy[tok], ls='dashed', color=s_m.to_rgba(kappa_BC), alpha=1.0, lw=2)
        axs[2].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color=s_m.to_rgba(kappa_BC), alpha=1.0, lw=2)
        xlims = [42.9, 45.]
        ylims = [-2.7, 0.]
        axs[2].set_ylim(ylims)
        axs[2].set_xlim(xlims)


    plot_OIIIEW_Mstar(z, axs[0], errorbar=False)
    plot_OIIIEW_UV(z, axs[1], errorbar=False)
    plot_OIIIlum_UV(z, axs[2], errorbar=False)

    for ax in axs:
        for label in (ax.get_xticklabels()+ax.get_yticklabels()):
            label.set_fontsize(14)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='in')
        ax.tick_params(axis='y', which='minor', direction='in')
        ax.grid(True, alpha = 0.5)

    for ii in [0,1]:
        axs[ii].set_yticks(np.arange(2,3.6,0.5))
    axs[0].set_xlabel(r'log$_{10}$(M$_{\star}$/M$_{\odot}$)', fontsize=14)
    axs[1].set_xlabel(r'log$_{10}$(L$_{\mathrm{FUV}}$/(erg s$^{-1}$))', fontsize=14)
    axs[2].set_xlabel(r'log$_{10}$(L$_{\mathrm{FUV}}$/(erg s$^{-1}$))', fontsize=14)
    axs[0].set_ylabel(r'log$_{10}$(EW([OIII]+H$\beta$)/$\AA$)',  fontsize=14)
    axs[1].set_ylabel(r'log$_{10}$(EW([OIII]+H$\beta$)/$\AA$)', fontsize=14)
    axs[2].set_ylabel(r'log$_{10}$(L([OIII]+H$\beta$)/L$_{\mathrm{FUV}})$', fontsize=14)


    lines, labels = [], []
    for ax in fig.axes:
        axLine, axLabel = ax.get_legend_handles_labels()
        lines.extend(axLine)
        labels.extend(axLabel)
    lines = lines[:-4]
    labels = labels[:-4]

    axs[0].legend(lines[0:6], list(labels)[0:6], frameon=False, fontsize=11, numpoints=1, ncol=2)
    axs[1].legend(lines[6:11], list(labels)[6:11], frameon=False, fontsize=11, numpoints=1, ncol=2)
    axs[2].legend(lines[-2:], list(labels)[-2:], frameon=False, fontsize=11, numpoints=1, ncol=1, loc='lower right')

    fig.subplots_adjust(bottom=0.11, left = 0.05, wspace = 0.3, hspace=0)


elif input == 2:
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(5, 3), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    xlims = [-16.9, -24.7]
    ylims = [-2.8,-1.5]
    bins = -np.arange(16.8, 25, 0.4)[::-1]
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    ext_curves = ['Default', 'Calzetti_1.0', 'SMC_1.0', 'N18_1.0']
    labels = ['Default', 'Calzetti', 'SMC', 'N18']
    colors = ['black', 'brown', 'grey', 'orange']

    for ii, jj in enumerate(ext_curves):
        if ii == 0:
            datafolder = 'data/flares.hdf5'
            dataset = 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI'
        else:
            datafolder = 'data1/flares.hdf5'
            dataset = F'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI_{jj}'

        L = {}
        for f in ['FUV','NUV']:
            L[f] = np.array([])
            mstar, sfr = np.array([]), np.array([])
            w = np.array([])
            with h5py.File(datafolder, 'r') as hf:
                for sim in hf.keys():
                    tmp = np.array(hf[f'{sim}/{tag}/{dataset}/{f}'])
                    L[f] = np.hstack((L[f], tmp))
                    w = np.append(w, np.ones(len(tmp))*weights[int(sim)])

        beta = np.log10(L['FUV']/L['NUV'])/np.log10(1500./2500.) - 2.0
        x, y = lum_to_M(L['FUV']), beta

        quantiles = [0.84,0.50,0.16]
        out = fl.binned_weighted_quantile(x, y, w, bins, quantiles)
        hist, binedges = np.histogram(x, bins)
        ok = np.where(hist>0)[0]
        ok1 = np.where(hist[ok]>3)[0][0]

        axs.fill_between(bincen[ok][ok1:], out[:,0][ok][ok1:], out[:,2][ok][ok1:], color=colors[ii], alpha=0.3)
        axs.plot(bincen[ok], out[:,1][ok], ls='dashed', color=colors[ii], alpha=1, lw=2)
        axs.plot(bincen[ok][ok1:], out[:,1][ok][ok1:], ls='-', color=colors[ii], alpha=1, lw=2, label=labels[ii])

    plot_beta(5, axs)

    axs.set_ylim(ylims)
    axs.set_xlim(xlims)
    axs.minorticks_on()
    axs.tick_params(axis='x', which='minor', direction='in')
    axs.tick_params(axis='y', which='minor', direction='in')
    axs.grid(True, alpha = 0.5)
    axs.legend(frameon=False,fontsize=11,loc=4, ncol=2)
    for label in (axs.get_xticklabels()+axs.get_yticklabels()):
        label.set_fontsize(12)
    axs.set_ylabel(r'$\beta$', fontsize=15)
    axs.set_xlabel(r'M$_{1500}$', fontsize=15)

elif input == 3:
    fig, axs = plt.subplots(nrows = 1, ncols = 3, figsize=(13, 3), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    mstar = get_data_all(tag, dataset = 'Mstar_30', DF=False)
    z = 8

    fig, axs = plt.subplots(nrows = 1, ncols = 3, figsize=(13, 2.5), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    ext_curves = ['Default', 'Calzetti_1.0', 'SMC_1.0', 'N18_1.0']
    labels = ['Default', 'Calzetti', 'SMC', 'N18']
    colors = ['black', 'brown', 'grey', 'orange']

    ws = np.array([])
    for jj in range(40):
        ws = np.append(ws, np.ones(np.shape(mstar[jj]))*weights[jj])
    mstar = np.concatenate(mstar)

    for ii, curves in enumerate(ext_curves):
        if ii==0:
            dat1 = get_line_all(tag, 'OIII4959', inp = 'FLARES', LF = False)
            dat2 = get_line_all(tag, 'OIII5007', inp = 'FLARES', LF = False)
            dat3 = get_line_all(tag, 'HI4861', inp = 'FLARES', LF = False)
            l_fuvs = np.concatenate(get_lum_all(tag, LF=False))

            OIII4959_lum = np.concatenate(dat1[:,0])
            OIII4959_EW = np.concatenate(dat1[:,1])

            OIII5007_lum = np.concatenate(dat2[:,0])
            OIII5007_EW = np.concatenate(dat2[:,1])

            Hbeta_lum = np.concatenate(dat3[:,0])
            Hbeta_EW = np.concatenate(dat3[:,1])
        else:
            OIII4959_lum = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{curves}/OIII4959/Luminosity', tag))
            OIII4959_EW = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{curves}/OIII4959/EW', tag))

            OIII5007_lum = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{curves}/OIII5007/Luminosity', tag))
            OIII5007_EW = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{curves}/OIII5007/EW', tag))

            Hbeta_lum = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{curves}/HI4861/Luminosity', tag))
            Hbeta_EW = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Lines/DustModelI_{curves}/HI4861/EW', tag))

            l_fuvs = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Luminosity/DustModelI_{curves}/FUV', tag))

        x, y, w = np.log10(mstar*1e10), OIII4959_EW+OIII5007_EW+Hbeta_EW, ws
        thisok = np.where(y>0)
        x = x[thisok]
        y = np.log10(y[thisok])
        w = w[thisok]
        tbins = np.arange(min(x), max(x)+0.4, 0.4)
        tbincen = (tbins[1:]+tbins[:-1])/2.
        tbinwidth = tbins[1:] - tbins[:-1]
        out = fl.binned_weighted_quantile(x,y,w,tbins,quantiles)
        hist, binedges = np.histogram(x, tbins)
        xx, yy, yy16, yy84 = tbins, out[:,1], out[:,0],out[:,2]
        tok = np.where(hist>0)[0]
        tok1 = np.where(hist[tok]>3)[0][-1]
        axs[0].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], alpha=0.25, color=colors[ii])
        axs[0].plot(xx[tok], yy[tok], ls='dashed', color=colors[ii], alpha=0.8, lw=2)
        axs[0].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color=colors[ii], alpha=0.8, lw=2, label=labels[ii])
        xlims = [7.5, 10.5]
        ylims = [1.8, 3.9]
        axs[0].set_ylim(ylims)
        axs[0].set_xlim(xlims)



        x, y, w = np.log10(l_fuvs*conversion_fac), OIII4959_EW+OIII5007_EW+Hbeta_EW, ws
        thisok = np.where(y>0)
        x = x[thisok]
        y = np.log10(y[thisok])
        w = w[thisok]
        tbins = np.arange(min(x), max(x)+0.3, 0.3)
        tbincen = (tbins[1:]+tbins[:-1])/2.
        tbinwidth = tbins[1:] - tbins[:-1]
        out = fl.binned_weighted_quantile(x,y,w,tbins,quantiles)
        hist, binedges = np.histogram(x, tbins)
        xx, yy, yy16, yy84 = tbins, out[:,1], out[:,0],out[:,2]
        tok = np.where(hist>0)[0]
        tok1 = np.where(hist[tok]>3)[0][-1]
        axs[1].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color=colors[ii], alpha=0.25)
        axs[1].plot(xx[tok], yy[tok], ls='dashed', color=colors[ii], alpha=1.0, lw=2)
        axs[1].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color=colors[ii], alpha=.8, lw=2)
        xlims = [42.9, 45.]
        ylims = [1.8, 3.9]
        axs[1].set_ylim(ylims)
        axs[1].set_xlim(xlims)



        x, y, w = np.log10(l_fuvs*conversion_fac), np.log10(OIII4959_lum+OIII5007_lum+Hbeta_lum) - np.log10(l_fuvs*conversion_fac), ws
        thisok = np.where(y > -5)
        x = x[thisok]
        y = y[thisok]
        w = w[thisok]
        tbins = np.arange(min(x), max(x)+0.3, 0.3)
        tbincen = (tbins[1:]+tbins[:-1])/2.
        tbinwidth = tbins[1:] - tbins[:-1]
        hist, binedges = np.histogram(x, tbins)
        out = fl.binned_weighted_quantile(x,y,w,tbins,quantiles)
        xx, yy, yy16, yy84 = tbins, out[:,1], out[:,0],out[:,2]
        tok = np.where(hist>0)[0]
        tok1 = np.where(hist[tok]>3)[0][-1]
        axs[2].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color=colors[ii], alpha=0.25)
        axs[2].plot(xx[tok], yy[tok], ls='dashed', color=colors[ii], alpha=1.0, lw=2)
        axs[2].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color=colors[ii], alpha=1.0, lw=2)
        xlims = [42.9, 45.]
        ylims = [-2.7, 0.]
        axs[2].set_ylim(ylims)
        axs[2].set_xlim(xlims)


    plot_OIIIEW_Mstar(z, axs[0], errorbar=False)
    plot_OIIIEW_UV(z, axs[1], errorbar=False)
    plot_OIIIlum_UV(z, axs[2], errorbar=False)

    for ax in axs:
        for label in (ax.get_xticklabels()+ax.get_yticklabels()):
            label.set_fontsize(14)
        ax.minorticks_on()
        ax.tick_params(axis='x', which='minor', direction='in')
        ax.tick_params(axis='y', which='minor', direction='in')
        ax.grid(True, alpha = 0.5)

    for ii in [0,1]:
        axs[ii].set_yticks(np.arange(2,3.6,0.5))
    axs[0].set_xlabel(r'log$_{10}$(M$_{\star}$/M$_{\odot}$)', fontsize=14)
    axs[1].set_xlabel(r'log$_{10}$(L$_{\mathrm{FUV}}$/(erg s$^{-1}$))', fontsize=14)
    axs[2].set_xlabel(r'log$_{10}$(L$_{\mathrm{FUV}}$/(erg s$^{-1}$))', fontsize=14)
    axs[0].set_ylabel(r'log$_{10}$(EW([OIII]+H$\beta$)/$\AA$)',  fontsize=14)
    axs[1].set_ylabel(r'log$_{10}$(EW([OIII]+H$\beta$)/$\AA$)', fontsize=14)
    axs[2].set_ylabel(r'log$_{10}$(L([OIII]+H$\beta$)/L$_{\mathrm{FUV}})$', fontsize=14)


    lines, labels = [], []
    for ax in fig.axes:
        axLine, axLabel = ax.get_legend_handles_labels()
        lines.extend(axLine)
        labels.extend(axLabel)
    lines = lines[:-4]
    labels = labels[:-4]

    axs[0].legend(lines[:-2], list(labels)[:-2], frameon=False, fontsize=11, numpoints=1, ncol=2)
    axs[1].legend(lines[-2:], list(labels)[-2:], frameon=False, fontsize=11, numpoints=1, ncol=1, loc=2)

    fig.subplots_adjust(bottom=0.11, left = 0.05, wspace = 0.3, hspace=0)

elif input == 4:
    fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(5, 2.5), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    xlims = [-16.9, -24.7]
    ylims = [0,3]
    bins = -np.arange(16.8, 25, 0.4)[::-1]
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    ext_curves = ['Default', 'Calzetti_1.0', 'SMC_1.0', 'N18_1.0']
    labels = ['Default', 'Calzetti', 'SMC', 'N18']
    colors = ['black', 'brown', 'grey', 'orange']

    L_FUV_int = get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='Intrinsic')
    ws = np.array([])
    for jj in range(40):
        ws = np.append(ws, np.ones(np.shape(L_FUV_int[jj]))*weights[jj])
    L_FUV_int = np.concatenate(L_FUV_int)

    for ii, jj in enumerate(ext_curves):
        if ii == 0:
            L_FUV = np.concatenate(get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='DustModelI'))
        else:
            L_FUV = np.concatenate(get_all(F'BPASS_2.2.1/Chabrier300/Luminosity/DustModelI_{jj}/FUV', tag))

        att = -2.5*np.log10(L_FUV/L_FUV_int)

        x, y = lum_to_M(L_FUV_int), att

        out = fl.binned_weighted_quantile(x, y, ws, bins, quantiles)
        hist, binedges = np.histogram(x, bins)
        ok = np.where(hist>0)[0]
        ok1 = np.where(hist[ok]>3)[0][0]
        axs.fill_between(bincen[ok][ok1:], out[:,0][ok][ok1:], out[:,2][ok][ok1:], color=colors[ii], alpha=0.3)
        axs.plot(bincen[ok][ok1:], out[:,1][ok][ok1:], ls='-', color=colors[ii], alpha=.8, lw=2, label=labels[ii])
        axs.plot(bincen[ok], out[:,1][ok], ls='dashed', color=colors[ii], alpha=.8, lw=2)
        axs.set_ylim(ylims)
        axs.set_xlim(xlims)
        axs.minorticks_on()
        axs.tick_params(axis='x', which='minor', direction='in')
        axs.tick_params(axis='y', which='minor', direction='in')
        axs.grid(True, alpha = 0.5)
        axs.legend(frameon=False,fontsize=10, ncol=2, loc=2)
        for label in (axs.get_xticklabels()+axs.get_yticklabels()):
            label.set_fontsize(10)
        axs.set_ylabel(r'A$_{\mathrm{FUV}}$=-2.5 log$_{10}$(L$_{\mathrm{FUV}}^{\mathrm{Observed}}$/L$_{\mathrm{FUV}}^{\mathrm{Intrinsic}}$)', fontsize=11)
        axs.set_xlabel(r'M$_{1500}$(Intrinsic)', fontsize=11)

fig.savefig(F'App_{plt_options[input]}.pdf', dpi = 300, bbox_inches='tight')
plt.show()
