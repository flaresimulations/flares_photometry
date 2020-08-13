"""

    Plots the line luminosity & EW relation of [OIII]
    doublet at z=7, 8: Figure 13

"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from FLARE.photom import lum_to_M, M_to_lum
import FLARE.filters
import flares
from modules import get_data_all, get_line_all, get_lum_all
from plot_obs import plot_OIIIEW_Mstar, plot_OIIIEW_UV, plot_OIIIlum_UV
import seaborn as sns
sns.set_context("paper")


tags = ['008_z007p000', '007_z008p000']
filters = 'FUV'
conversion_fac = 2E15  #converting from ergs/s/Hz to ergs/s at FUV
h = 0.6777
vol = (4/3)*np.pi*(14/h)**3
quantiles = [0.84,0.50,0.16]

bins = np.arange(40, 47, 0.4)
bincen = (bins[1:]+bins[:-1])/2.
binwidth = bins[1:] - bins[:-1]

df = pd.read_csv('weight_files/weights_grid.txt')
weights = np.array(df['weights'])
sims = np.arange(len(weights))

fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(13, 6), sharex=False, sharey=False, facecolor='w', edgecolor='k')
axs = axs.ravel()

for ii, tag in enumerate(tags):

    z = float(tag[5:].replace('p','.'))
    print (z)

    dat1 = get_line_all(tag, 'OIII4959', inp = 'FLARES', LF = False)
    dat2 = get_line_all(tag, 'OIII5007', inp = 'FLARES', LF = False)
    dat3 = get_line_all(tag, 'HI4861', inp = 'FLARES', LF = False)
    l_fuvs = np.concatenate(get_lum_all(tag, LF=False))
    mstar = get_data_all(tag, dataset = 'Mstar_30', DF=False)

    ws = np.array([])
    for jj in sims:
        ws = np.append(ws, np.ones(np.shape(mstar[jj]))*weights[jj])
    mstar = np.concatenate(mstar)*1e10

    OIII4959_lum = np.concatenate(dat1[:,0])
    OIII4959_EW = np.concatenate(dat1[:,1])

    OIII5007_lum = np.concatenate(dat2[:,0])
    OIII5007_EW = np.concatenate(dat2[:,1])

    Hbeta_lum = np.concatenate(dat3[:,0])
    Hbeta_EW = np.concatenate(dat3[:,1])

    # df = pd.read_csv('Magnitude_limits.txt')
    # low = np.array(df[filters])[ii+2]
    # ok = np.where(lum_to_M(l_fuvs) < low)
    # mstar, l_fuvs, OIII4959_EW, OIII5007_EW, Hbeta_EW, OIII4959_lum, OIII5007_lum, Hbeta_lum, ws = mstar[ok], l_fuvs[ok], OIII4959_EW[ok], OIII5007_EW[ok], Hbeta_EW[ok], OIII4959_lum[ok], OIII5007_lum[ok], Hbeta_lum[ok], ws[ok]

    x, y, w = np.log10(mstar), OIII4959_EW+OIII5007_EW+Hbeta_EW, ws
    thisok = np.where(y>0)
    x = x[thisok]
    y = np.log10(y[thisok])
    w = w[thisok]
    tbins = np.arange(min(x), max(x)+0.4, 0.4)
    tbincen = (tbins[1:]+tbins[:-1])/2.
    tbinwidth = tbins[1:] - tbins[:-1]
    out = flares.binned_weighted_quantile(x,y,w,tbins,quantiles)
    hist, binedges = np.histogram(x, tbins)
    xx, yy, yy16, yy84 = tbincen, out[:,1], out[:,2],out[:,0]
    tok = np.where(hist>0)[0]
    tok1 = np.where(hist[tok]>3)[0][-1]
    axs[3*ii+0].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color='black', alpha=0.5)
    axs[3*ii+0].plot(xx[tok], yy[tok], ls='dashed', color='black', alpha=1.0, lw=1)
    axs[3*ii+0].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color='black', alpha=1.0, lw=1)
    xlims = [7.5, 10.5]
    ylims = [1.8, 3.5]
    axs[3*ii+0].hexbin(x, y, gridsize=(30,15), bins='log', cmap='Greys_r', linewidths=0., mincnt=5, extent=[*xlims, *ylims], alpha=0.6, zorder=2)
    plot_OIIIEW_Mstar(z, axs[3*ii+0])


    x, y, w = np.log10(l_fuvs*conversion_fac), OIII4959_EW+OIII5007_EW+Hbeta_EW, ws
    thisok = np.where(y>0)
    x = x[thisok]
    y = np.log10(y[thisok])
    w = w[thisok]
    tbins = np.arange(min(x), max(x)+0.3, 0.3)
    tbincen = (tbins[1:]+tbins[:-1])/2.
    tbinwidth = tbins[1:] - tbins[:-1]
    out = flares.binned_weighted_quantile(x,y,w,tbins,quantiles)
    hist, binedges = np.histogram(x, tbins)
    xx, yy, yy16, yy84 = tbincen, out[:,1], out[:,2],out[:,0]
    tok = np.where(hist>0)[0]
    tok1 = np.where(hist[tok]>3)[0][-1]
    axs[3*ii+1].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color='black', alpha=0.5)
    axs[3*ii+1].plot(xx[tok], yy[tok], ls='dashed', color='black', alpha=1.0, lw=1)
    axs[3*ii+1].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color='black', alpha=1.0, lw=1)
    xlims = [42.9, 45.]
    ylims = [1.8, 3.5]
    axs[3*ii+1].hexbin(x, y, gridsize=(30,15), bins='log', cmap='Greys_r', linewidths=0., mincnt=5, extent=[*xlims, *ylims], alpha=0.6, zorder=2)
    plot_OIIIEW_UV(z, axs[3*ii+1])


    x, y, w = np.log10(l_fuvs*conversion_fac), np.log10(OIII4959_lum+OIII5007_lum+Hbeta_lum) - np.log10(l_fuvs*conversion_fac), ws
    thisok = np.where(y > -5)
    x = x[thisok]
    y = y[thisok]
    w = w[thisok]
    tbins = np.arange(min(x), max(x)+0.3, 0.3)
    tbincen = (tbins[1:]+tbins[:-1])/2.
    tbinwidth = tbins[1:] - tbins[:-1]
    hist, binedges = np.histogram(x, tbins)
    out = flares.binned_weighted_quantile(x,y,w,tbins,quantiles)
    xx, yy, yy16, yy84 = tbincen, out[:,1], out[:,2],out[:,0]
    tok = np.where(hist>0)[0]
    tok1 = np.where(hist[tok]>3)[0][-1]
    axs[3*ii+2].fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color='black', alpha=0.5)
    axs[3*ii+2].plot(xx[tok], yy[tok], ls='dashed', color='black', alpha=1.0, lw=1)
    axs[3*ii+2].plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color='black', alpha=1.0, lw=1)
    xlims = [42.9, 45.]
    ylims = [-2.7, -0.5]
    axs[3*ii+2].hexbin(x, y, gridsize=(35,15), bins='log', cmap='Greys_r', linewidths=0., mincnt=3, extent=[*xlims, *ylims], alpha=0.6, zorder=2)
    plot_OIIIlum_UV(z, axs[3*ii+2])

for ii in [0,1,2]:
        axs[ii].set_xticklabels('')

for ax in axs:
    for label in (ax.get_xticklabels()+ax.get_yticklabels()):
        label.set_fontsize(14)
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', direction='in')
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.grid(True, alpha = 0.5)
for ii in [0,1]:
    axs[ii*3+0].set_ylim(1.8,3.7)
    axs[ii*3+0].set_xlim(7.5, 11.4)
    axs[ii*3+1].set_xlim(42.9, 45.7)
    axs[ii*3+1].set_ylim(1.8,3.6)
    axs[ii*3+2].set_xlim(42.9, 45.7)
    axs[ii*3+2].set_ylim(-2.7,-0.25)
for ii in [0,1,3,4]:
    axs[ii].set_yticks(np.arange(2,3.6,0.5))
axs[-3].set_xlabel(r'log$_{10}$(M$_{\star}$/M$_{\odot}$)', fontsize=14)
axs[-2].set_xlabel(r'log$_{10}$(L$_{\mathrm{FUV}}$/(erg s$^{-1}$))', fontsize=14)
axs[-1].set_xlabel(r'log$_{10}$(L$_{\mathrm{FUV}}$/(erg s$^{-1}$))', fontsize=14)

axs[1].legend(frameon=False, fontsize=12, loc=4)
axs[-2].legend(frameon=False, fontsize=12, loc=4)

fig.subplots_adjust(bottom=0.11, left = 0.05, wspace = 0.3, hspace=0)
fig.text(0.005, 0.5, r'log$_{10}$(EW$_{[OIII]4959,5007+H\beta}$/$\AA$)', va='center', rotation='vertical', fontsize=14)
fig.text(0.31, 0.5, r'log$_{10}$(EW$_{[OIII]4959,5007+H\beta}$/$\AA$)', va='center', rotation='vertical', fontsize=14)
fig.text(0.61, 0.5, r'log$_{10}$(L$_{[OIII]4959,5007+H\beta}$/L$_{\mathrm{FUV}})$', va='center', rotation='vertical', fontsize=14)

fig.text(0.23, 0.84, r'$z=7$', va='center', fontsize=15)
fig.text(0.54, 0.84, r'$z=7$', va='center', fontsize=15)
fig.text(0.85, 0.84, r'$z=7$', va='center', fontsize=15)
fig.text(0.23, 0.43, r'$z=8$', va='center', fontsize=15)
fig.text(0.54, 0.43, r'$z=8$', va='center', fontsize=15)
fig.text(0.85, 0.43, r'$z=8$', va='center', fontsize=15)

plt.savefig('line_lumOIII_z7_8.pdf', dpi=300, bbox_inches='tight')
plt.show()
