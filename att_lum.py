import sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from modules import get_lum_all, get_data_all
from FLARE.photom import lum_to_M, M_to_lum
import flares as fl
import seaborn as sns
sns.set_context("paper")

filters = 'FUV'
zs = [5, 6, 7, 8, 9, 10]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']

plt_options = ['Observed', 'Intrinsic', 'Stellar']

input = plt_options[int(sys.argv[1])]


xlims = [-16.9, -24.7]
ylims = [-0.5, 4.4]

df = pd.read_csv('weight_files/weights_grid.txt')
weights = np.array(df['weights'])



fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(13, 5), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()

for ii, tag in enumerate(tags):

    df = pd.read_csv('Magnitude_limits.txt')
    low = np.array(df[filters])[ii]
    bins = -np.arange(-low, 25, 0.4)[::-1]

    L_FUV = get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='DustModelI')
    L_NUV = get_lum_all(tag, LF = False, filter = 'NUV', Luminosity='DustModelI')
    L_FUV_int = get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='Intrinsic')
    Mstar_30 = get_data_all(tag, inp = 'FLARES', DF = False)
    # sfr_30 = get_data_all(tag, dataset = 'SFR_inst_30', inp = 'FLARES', DF = False)

    ws = np.array([])
    for jj in range(len(weights)):
        ws = np.append(ws, np.ones(np.shape(L_FUV[jj]))*weights[jj])
    L_FUV = np.concatenate(L_FUV)
    L_FUV_int = np.concatenate(L_FUV_int)
    L_NUV = np.concatenate(L_NUV)
    Mstar_30 = np.concatenate(Mstar_30)
    # sfr_30 = np.concatenate(sfr_30)


    ok = np.where(lum_to_M(L_FUV)<low)[0]
    L_FUV, L_NUV, L_FUV_int, Mstar_30 = L_FUV[ok], L_NUV[ok], L_FUV_int[ok], Mstar_30[ok]
    # sfr_30 = sfr_30[ok]/Mstar_30

    beta = np.log10(L_FUV/L_NUV)/np.log10(1500./2500.) - 2.0
    att = -2.5*np.log10(L_FUV/L_FUV_int)
    zpos = -21.3

    # if tag=='010_z005p000':
    #     from astropy.cosmology import Planck13
    #     from astropy import units as u
    #     tmp = 1/Planck13.H(5).decompose()
    #     print ((1.0/(3*tmp.to(u.yr))).value)
    #     check = np.logical_and(Mstar_30>1e10, att<0.4)
    #     print (sfr_30[check])
    #     print (beta[check])

    if input == plt_options[0]:
        x, y, z, w = lum_to_M(L_FUV), att, beta, ws[ok]
        gridsize=(50,21)
        extent=[*[low,-24.5], *ylims]
        add=''
        savename=F'att_lfuv_beta_z5_10.pdf'
        xlabel = r'M$_{1500}$'
    elif input == plt_options[1]:
        x, y, z, w = lum_to_M(L_FUV_int), att, beta, ws[ok]
        gridsize=(55,21)
        xlims = [-16.9, -25.7]
        extent=[*xlims, *ylims]
        savename=F'att_lfuvintr_beta_z5_10.pdf'
        xlabel = r'M$_{1500}\mathrm{(Intrinsic)}$'
    elif input == plt_options[2]:
        x, y, z, w = np.log10(Mstar_30), att, beta, ws[ok]
        gridsize=(55,21)
        bins = np.arange(7.25,12,0.5)
        xlims = [7.5, 11.3]
        extent=[*xlims, *ylims]
        zpos = 8.
        savename=F'att_Mstar_beta_z5_10.pdf'
        xlabel = r'$\mathrm{log}_{10}(\mathrm{M}_{\star}/\mathrm{M}_{\odot})$'

    if ii == 0:
        hb = axs[ii].hexbin(x, y, C=z, gridsize=gridsize, cmap=plt.cm.get_cmap('coolwarm'), mincnt=1, extent=extent, alpha=1, reduce_C_function=np.median, vmin=-2.4, vmax=-1.2)
    else:
        axs[ii].hexbin(x, y, C=z, gridsize=gridsize, cmap=plt.cm.get_cmap('coolwarm'), mincnt=1, extent=extent, alpha=1, reduce_C_function=np.median, vmin=-2.4, vmax=-1.1)


    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]
    quantiles = [0.84,0.50,0.16]
    out = fl.binned_weighted_quantile(x, y, w, bins, quantiles)
    hist, binedges = np.histogram(x, bins)
    ok = np.where(hist>0)[0]
    ok1 = np.where(hist[ok]>3)[0][0]
    axs[ii].fill_between(bincen[ok][ok1:], out[:,0][ok][ok1:], out[:,2][ok][ok1:], color='black', alpha=0.2)
    axs[ii].plot(bincen[ok][ok1:], out[:,1][ok][ok1:], ls='-', color='black', alpha=.5, lw=2)
    axs[ii].plot(bincen[ok], out[:,1][ok], ls='dashed', color='black', alpha=.5, lw=2)


    axs[ii].set_xlim(xlims)
    axs[ii].set_ylim(0.,3.7)
    axs[ii].minorticks_on()
    axs[ii].tick_params(axis='x', which='minor', direction='in')
    axs[ii].tick_params(axis='y', which='minor', direction='in')
    axs[ii].grid(True, alpha = 0.5)
    axs[ii].text(zpos, 3.2, r'$z = {}$'.format(zs[ii]), fontsize = 13)

    for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
        label.set_fontsize(11)

cbaxes = fig.add_axes([0.83, 0.15, 0.008, 0.3])
fig.colorbar(hb, cax=cbaxes)
cbaxes.set_ylabel(r'$\beta$', fontsize = 14)
for label in cbaxes.get_yticklabels():
    label.set_fontsize(12)

fig.subplots_adjust(bottom=0.11, left = 0.05, wspace=0, hspace=0)
fig.text(0.01, 0.5, r'A$_{\mathrm{FUV}}$=-2.5 log$_{10}$(L$_{\mathrm{FUV}}^{\mathrm{Observed}}$/L$_{\mathrm{FUV}}^{\mathrm{Intrinsic}}$)', va='center', rotation='vertical', fontsize=15)
fig.text(0.43, 0.03, xlabel, va='center', fontsize=15)

plt.savefig(savename, bbox_inches='tight', dpi=300)

plt.show()
