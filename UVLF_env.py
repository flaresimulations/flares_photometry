import numpy as np
import pandas as pd
import h5py, sys
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from astropy.cosmology import Planck13 as cosmo
from astropy import units as u

from modules import get_lum, get_lum_all
from FLARE.photom import lum_to_M, M_to_lum
import flares as fl

import seaborn as sns
sns.set_context("paper")

##fac = 0.1225##
h = 0.6777

rho_crit = (cosmo.critical_density(5).to(u.Msun/u.Mpc**3)).value
rho_crit_ref = (cosmo.critical_density(5.037).to(u.Msun/u.Mpc**3)).value

vol = (4/3)*np.pi*(14/h)**3
refvol = 100**3
AGNdT9vol = 50**3
filter = 'FUV'


tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']

# tags = ['010_z005p000', '008_z007p000', '006_z009p000']
# tags_ref = ['008_z005p037', '005_z007p050', '003_z008p988']

arr = np.arange(0,40)
density = np.array([-0.3, 0.3])
zs = [5., 6., 7., 8., 9., 10.]

norm = matplotlib.colors.Normalize(vmin=min(density), vmax=max(density))
dbins = np.arange(-0.3, 0.4, 0.1)
# choose a colormap
c_m = matplotlib.cm.plasma

# create a ScalarMappable and initialize a data structure
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

df = pd.read_csv('./weight_files/weights_grid.txt')
delta = np.array(df['log(1+delta)'])

fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(14, 9), sharex=True, sharey=True, facecolor='w', edgecolor='k')
axs = axs.ravel()
plt_options = ['UVLF_env', 'att_env']
input = plt_options[int(sys.argv[1])]

#
# for ii in range(len(dbins)-1):
#
#     ok = np.where(np.logical_and(delta >= dbins[ii], delta < dbins[ii+1]))[0]
#     if len(ok)>0:
#         jj = np.median(delta[ok])
#         mhist = np.zeros(len(bincen))
#         menc = 0
#         for kk in ok:
#             mhist += get_hist(kk, tags[0], bins = bins, LF = True)
#             menc += (10**(delta[kk])) * rho_crit * vol
#
#         yy = mhist/(menc)
#
#         axs.plot(bincen, np.log10(yy), lw = 2, ls = 'solid', color=s_m.to_rgba(jj))

#
# y = get_hist(0, tags_ref[0], bins = bins, inp = 'REF', LF = True)/(rho_crit_ref * refvol)
# axs.plot(bincen, np.log10(y), lw = 2, ls = 'dashed', color='red', label = 'EAGLE Ref')
#
# y = get_hist(0, tags_ref[0], bins = bins, inp = 'AGNdT9', LF = True)/(rho_crit_ref * AGNdT9vol)
# axs.plot(bincen, np.log10(y), lw = 2, ls = 'dotted', color='red', label = 'EAGLE AGNdT9')

for ii, jj in enumerate(tags):

    df = pd.read_csv('Magnitude_limits.txt')
    low = np.array(df[filter])[ii]
    bins = -np.arange(-low, 25, 0.5)[::-1]
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    if input == plt_options[1]:
        LFUV = get_lum_all(jj, LF = False, filter = 'FUV', Luminosity='DustModelI')
        LFUV_int = get_lum_all(jj, LF = False, filter = 'FUV', Luminosity='Intrinsic')
        xlabel = r'M$_{1500}$'#\mathrm{(Intrinsic)}$'
        ylabel = r'A$_{\mathrm{FUV}}$=-2.5 log$_{10}$(L$_{\mathrm{FUV}}^{\mathrm{Observed}}$/L$_{\mathrm{FUV}}^{\mathrm{Intrinsic}}$)'
        ylim = (0.,3.7)
        savename='att_env_obs.pdf'
        axs[ii].text(-21.3, 3.2, r'$z = {}$'.format(zs[ii]), fontsize = 13)

    for kk in range(len(dbins)-1):

        ok = np.where(np.logical_and(delta >= dbins[kk], delta < dbins[kk+1]))[0]

        if input == plt_options[0]:
            mhist = np.zeros(len(bincen))
            for ll in ok:
                mhist += get_lum(ll, jj, bins = bins, LF = True)

            yy = mhist/(binwidth*(len(ok)*vol))
            yyerr = np.sqrt(mhist)/(binwidth*(len(ok)*vol))
            nonzero = np.where(mhist>0)[0]
            axs[ii].errorbar(bincen[nonzero], np.log10(yy[nonzero]), lw=2, ls='solid', marker='o', yerr=[np.log10(yy[nonzero]) - np.log10(yy[nonzero]-yyerr[nonzero]), np.log10(yy[nonzero]+yyerr[nonzero]) - np.log10(yy[nonzero])], color=s_m.to_rgba((dbins[kk]+dbins[kk+1])/2))

        else:
            LFUV_this = np.concatenate(LFUV[ok])
            LFUV_int_this = np.concatenate(LFUV_int[ok])
            y = -2.5*np.log10(LFUV_this/LFUV_int_this)
            x = lum_to_M(LFUV_this)
            quantiles = [0.84,0.50,0.16]
            out = fl.binned_weighted_quantile(x, y, np.ones(len(x)), bins, quantiles)
            hist, binedges = np.histogram(x, bins)
            tok = np.where(hist>0)[0]

            axs[ii].errorbar(bincen[tok], out[:,1][tok], lw=2, ls='solid', marker='o', yerr=[out[:,1][tok] - out[:,0][tok], out[:,2][tok]-out[:,0][tok]], color=s_m.to_rgba((dbins[kk]+dbins[kk+1])/2))

    if input == plt_options[0]:
        out, hist, err = get_lum_all(jj, bins = bins)
        ok = np.where(hist > 0)[0]
        hist = hist[ok]
        Msim = out/(binwidth*vol)
        xerr = np.ones(len(out))*binwidth[0]/2.
        yerr = err/(vol*binwidth)


        axs[ii].errorbar(bincen[ok], np.log10(Msim[ok]), lw=2, ls='solid', marker='o', color='black', yerr=[np.log10(Msim[ok]) - np.log10(Msim[ok]-yerr[ok]), np.log10(Msim[ok]+yerr[ok]) - np.log10(Msim[ok])], label=r'$\mathrm{Composite}$ $\mathrm{Function}$')
        xlabel = r'$\mathrm{M}_{1500}$'
        ylabel = r'$\mathrm{log}_{10}(\Phi/(\mathrm{cMpc}^{-3}\mathrm{Mag}^{-1}))$'
        ylim=(-6.6, -1.3)
        savename = 'UVLF_env.pdf'
        axs[ii].text(-17.4, -5.4, r'$z = {}$'.format(zs[ii]), fontsize = 14)


    axs[ii].grid(True, alpha = 0.4)
    axs[ii].minorticks_on()
    axs[ii].tick_params(axis='x', which='minor', direction='in')
    axs[ii].tick_params(axis='y', which='minor', direction='in')

    for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
        label.set_fontsize(14)

    axs[ii].set_xlim(-17, -24.9)
    axs[ii].set_ylim(ylim)

    #axs[ii].set_ylabel(r'# galaxies/($\delta$ x $\rho_{crit}$(z) x Vol)', fontsize=20)
    #axs[ii].set_xlabel(r'M$_{FUV}$', fontsize=20)

axs[4].legend(frameon=False, fontsize=14)
cbaxes = fig.add_axes([0.93, 0.12, 0.008, 0.3])
fig.colorbar(s_m, cax=cbaxes)
cbaxes.set_ylabel(r'$\mathrm{log}_{10}(1+\delta)$', fontsize = 14)
for label in cbaxes.get_yticklabels():
    label.set_fontsize(14)


fig.subplots_adjust(bottom=0.1, left = 0.08, right = 0.999, wspace=0, hspace=0)

fig.text(0.02, 0.5, ylabel, va='center', rotation='vertical', fontsize=16)
fig.text(0.51, 0.03, xlabel, va='center', fontsize=16)

plt.savefig(savename, bbox_inches='tight')
plt.show()
