"""

    Plots the EW relation of [CIII] doublet at z=7: Figure 15

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

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(6, 5), sharex=False, sharey=False, facecolor='w', edgecolor='k')


obs = {}

obs['EGS-zs8-1'] = {}
obs['EGS-zs8-1']['EW'] = [22., 2.2]
obs['EGS-zs8-1']['L1500'] = 29.46  # M = -22.06
obs['EGS-zs8-1']['z'] = 7.73
obs['EGS-zs8-1']['ref'] = r'Stark+2017'

# Stark et al. https://arxiv.org/pdf/1408.3649.pdf
obs['GN-108036'] = {}
obs['GN-108036']['EW'] = [7.6, 2.8]
obs['GN-108036']['L1500'] = 29.31  # J=25.2 unlensed M=-21.7
obs['GN-108036']['z'] = 7.2
obs['GN-108036']['ref'] = r'Stark+2015'

obs['A383-5.2'] = {}
obs['A383-5.2']['EW'] = [12.7, 3.5]
obs['A383-5.2']['L1500'] = 28.37  # J=25.2 magnified 7.4 M (lensed) = -21.5
obs['A383-5.2']['z'] = 6.0265
obs['A383-5.2']['M*'] = [9.5,[0.1,0.1]]
obs['A383-5.2']['ref'] = r'Stark+2015'

# Hutchison et al. 2019 https://arxiv.org/pdf/1905.08812.pdf
obs['z7-GND-42912'] = {}
obs['z7-GND-42912']['EW'] = [16.23, 2.32]
obs['z7-GND-42912']['z'] = 7.5032
obs['z7-GND-42912']['L1500'] = np.log10(M_to_lum(-21.58))
obs['z7-GND-42912']['ref'] = r'Hutchison+2019'


tag = '008_z007p000'
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

z = float(tag[5:].replace('p','.'))
print (z)


for id, marker in zip(obs.keys(), ['o','^','s','p','D','H','d','8','*']):

    o = obs[id]

    # c = 'k'
    c = 'forestgreen'

    ax.scatter(o['L1500']+np.log10(conversion_fac),np.log10(o['EW'][0]), marker=marker, s = 15, c = c, zorder = 5, label = r'$\rm '+id+'\ ('+o['ref']+')\ z='+str(o['z'])+'$')

    if o['EW'][1]>0:

        ax.plot([o['L1500']+np.log10(conversion_fac)]*2, [np.log10(o['EW'][0]-o['EW'][1]),np.log10(o['EW'][0]+o['EW'][1])], c = c, lw = 1)

    else:

        ax.arrow(o['L1500']+np.log10(conversion_fac), np.log10(o['EW'][0]), 0.0, -0.07, color=c, head_width=0.05, head_length=0.01)


dat1 = get_line_all(tag, 'CIII1907', inp = 'FLARES', LF = False, Type = 'Intrinsic')
dat2 = get_line_all(tag, 'CIII1909', inp = 'FLARES', LF = False, Type = 'Intrinsic')
l_fuvs = get_lum_all(tag, LF=False)

ws = np.array([])
for jj in sims:
    ws = np.append(ws, np.ones(np.shape(l_fuvs[jj]))*weights[jj])
l_fuvs = np.concatenate(l_fuvs)*conversion_fac

CIII1907_lum = np.concatenate(dat1[:,0])
CIII1907_EW = np.concatenate(dat1[:,1])

CIII1909_lum = np.concatenate(dat2[:,0])
CIII1909_EW = np.concatenate(dat2[:,1])


x, y, w = np.log10(l_fuvs), CIII1907_EW+CIII1909_EW, ws
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
ax.fill_between(xx[tok][:tok1+1], yy16[tok][:tok1+1], yy84[tok][:tok1+1], color='black', alpha=0.5)
ax.plot(xx[tok], yy[tok], ls='dashed', color='black', alpha=1.0, lw=1)
ax.plot(xx[tok][:tok1+1], yy[tok][:tok1+1], ls='-', color='black', alpha=1.0, lw=1, label='FLARES z={}'.format(z))
xlims = [42.9, 45.4]
ylims = [0.1, 3.5]
ax.hexbin(x, y, gridsize=(45,35), bins='log', cmap='Greys_r', linewidths=0., mincnt=3, extent=[*xlims, *ylims], alpha=0.6, zorder=2)

ax.set_xlim(43.2,45.5)
ax.set_ylim(0.25,2.2)
for label in (ax.get_xticklabels()+ax.get_yticklabels()):
    label.set_fontsize(14)
ax.minorticks_on()
ax.tick_params(axis='x', which='minor', direction='in')
ax.tick_params(axis='y', which='minor', direction='in')
ax.grid(True, alpha = 0.5)
ax.legend(frameon=False, fontsize=12, loc=1)
ax.set_xlabel(r'log$_{10}$(L$_{\mathrm{FUV}}$/(erg s$^{-1}$))', fontsize=14)
ax.set_ylabel(r'log$_{10}$(EW$_{[CIII],CIII]}$/$\AA$)', fontsize=14)
plt.savefig(F'line_lum_CIII_z{int(z)}.pdf', dpi=300, bbox_inches='tight')

plt.show()
