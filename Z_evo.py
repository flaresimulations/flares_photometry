"""

    Plots the redshift evolution of metallicity as a function
    of stellar mass, Figure 2

"""

import numpy as np
import pandas as pd
from functools import partial
import h5py, matplotlib, schwimmbad
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from modules import get_data_all
import flares
import seaborn as sns
sns.set_context("paper")

#Faisst et al.(2016) metallicity data
def Faisst_Z(ax):

    logmstar = np.array([8., 9., 8.75, 9.5, 10.])
    xerr = np.array([0.2, 0.2, 0.25, 0.5, 0.2])
    xuplim = np.array([1, 1, 0, 0, 1])

    OH = np.array([8.30, 8.26, 8.05, 8.12, 8.32])
    Z = 10**(OH - 12.) * (0.02/(10**(8.69-12)))
    yerr = np.array([[0.32, 0.51, 0.74, 0.52, 0.74], [0.26, 0.31, 0.56, 0.33, 0.40]])
    yerr[0]=OH-yerr[0]
    yerr[1]+=OH
    yerr = 10**(yerr-12.)
    yerr = yerr * (0.02/(10**(8.69-12)))
    yerr[0] = np.log10(Z) - np.log10(yerr[0])
    yerr[1] = np.log10(yerr[1]) - np.log10(Z)

    print (np.log10(Z), yerr)

    ax.errorbar(logmstar, np.log10(Z), xerr=xerr, yerr=yerr, xuplims=xuplim, label='Faisst+2016 [$z\sim 5$]', color='black', fmt='o', alpha=0.5, markersize=7)

def Troncosco(ax):

    log10M0, K0 = 11.59, 8.44 #average metallicity fit
    log10Mstar = np.arange(7., 12., 0.5)
    OH = 10**(-0.0864 * (log10Mstar - log10M0)**2 + K0 - 12.)
    Z = OH * (0.02/(10**(8.69-12)))
    ax.plot(log10Mstar, np.log10(Z), color='black', ls='dotted', label='Troncosco+2014 [$z\sim 3.4$]')


def get_data(ii, tag, inp = 'FLARES'):

    num = str(ii)
    if inp == 'FLARES':
        if len(num) == 1:
            num =  '0'+num

        sim = rF"../flares_pipeline/data/FLARES_{num}_sp_info.hdf5"
        num = F"{num}/"


    else:
        sim = rF"../flares_pipeline/data/EAGLE_{inp}_sp_info.hdf5"
        num=""


    with h5py.File(sim, 'r') as hf:
        S_len = np.array(hf[tag+'/Galaxy'].get('S_Length'), dtype = np.int64)
        G_len = np.array(hf[tag+'/Galaxy'].get('G_Length'), dtype = np.int64)
        S_mass = np.array(hf[tag+'/Particle'].get('S_MassInitial'), dtype = np.float64)*1e10
        G_mass = np.array(hf[tag+'/Particle'].get('G_Mass'), dtype = np.float64)*1e10
        S_Z = np.array(hf[tag+'/Particle'].get('S_Z_smooth'), dtype = np.float64)
        G_Z = np.array(hf[tag+'/Particle'].get('G_Z_smooth'), dtype = np.float64)

    begin = np.zeros(len(S_len), dtype = np.int64)
    end = np.zeros(len(S_len), dtype = np.int64)
    begin[1:] = np.cumsum(S_len)[:-1]
    end = np.cumsum(S_len)

    gbegin = np.zeros(len(G_len), dtype = np.int64)
    gend = np.zeros(len(G_len), dtype = np.int64)
    gbegin[1:] = np.cumsum(G_len)[:-1]
    gend = np.cumsum(G_len)

    G_met = np.zeros(len(begin))
    S_met = np.zeros(len(begin))
    for jj in range(len(begin)):
        S_met[jj] = np.sum(S_mass[begin[jj]:end[jj]]*S_Z[begin[jj]:end[jj]])/np.sum(S_mass[begin[jj]:end[jj]])
        G_met[jj] = np.sum(G_mass[gbegin[jj]:gend[jj]]*G_Z[gbegin[jj]:gend[jj]])/np.sum(G_mass[gbegin[jj]:gend[jj]])

    return S_met, G_met

def get_Z(tag):

    df = pd.read_csv('weight_files/weights_grid.txt')
    weights = np.array(df['weights'])
    sims = np.arange(0,len(weights))

    calc = partial(get_data, tag = tag)

    pool = schwimmbad.MultiPool(processes=12)
    dat = np.array(list(pool.map(calc, sims)))
    pool.close()

    return dat

tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
quantiles = [0.84,0.50,0.16]
df = pd.read_csv('weight_files/weights_grid.txt')
weights = np.array(df['weights'])
sims = np.arange(len(weights))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5), sharex=True, sharey=True, facecolor='w', edgecolor='k')

norm = matplotlib.colors.Normalize(vmin=0.5, vmax=len(tags)+0.5)
# choose a colormap
c_m = matplotlib.cm.viridis_r
# create a ScalarMappable and initialize a data structure
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])

bins = np.arange(7.5, 12, 0.5)
bincen = (bins[1:]+bins[:-1])/2.

yerr1_lo, yerr1_up, yerr2_lo, yerr2_up = np.array([]), np.array([]), np.array([]), np.array([])

for ii, tag in enumerate(tags):

    z = float(tag[5:].replace('p','.'))
    mstar = get_data_all(tag, inp = 'FLARES', DF = False)
    ws = np.array([])
    for jj in sims:
        ws = np.append(ws, np.ones(np.shape(mstar[jj]))*weights[jj])
    mstar = np.log10(np.concatenate(mstar)*1e10)

    data = get_Z(tag)

    S_met = np.log10(np.concatenate(data[:,0]))
    G_met = np.log10(np.concatenate(data[:,1]))

    out = flares.binned_weighted_quantile(mstar, G_met,ws,bins,quantiles)
    xx, yy, yy84, yy16 = bincen, out[:,1], out[:,0],out[:,2]
    hist, edges=np.histogram(mstar, bins)
    ok = np.where(hist>=5)[0]
    ax.errorbar(xx[ok], yy[ok], ls='-', color=s_m.to_rgba(ii+0.5), alpha=.9, lw=1, fmt='s')
    yerr1_lo = np.append(yerr1_lo, np.max(yy[ok]-yy16[ok]))
    yerr1_up = np.append(yerr1_up, np.max(yy84[ok]-yy[ok]))


    out = flares.binned_weighted_quantile(mstar, S_met,ws,bins,quantiles)
    xx, yy, yy84, yy16 = bincen, out[:,1], out[:,0],out[:,2]
    ax.errorbar(xx[ok], yy[ok], ls='-', color=s_m.to_rgba(ii+0.5), alpha=0.35, lw=2, fmt='D')
    yerr2_lo = np.append(yerr2_lo, np.max(yy[ok]-yy16[ok]))
    yerr2_up = np.append(yerr2_up, np.max(yy84[ok]-yy[ok]))

ax.errorbar(9.85, -3.5, yerr=[[np.max(yerr1_lo)], [np.max(yerr1_up)]], color='red', fmt='s', markersize=10, alpha=0.8)
ax.errorbar(10.15, -3.5, yerr=[[np.max(yerr2_lo)], [np.max(yerr2_up)]], color='red', fmt='D', markersize=10, alpha=0.35)

Troncosco(ax)
Faisst_Z(ax)


ax.set_xlim(7.5, 11.5)
for label in (ax.get_xticklabels()+ax.get_yticklabels()):
    label.set_fontsize(14)
ax.minorticks_on()
ax.tick_params(axis='x', which='minor', direction='in')
ax.tick_params(axis='y', which='minor', direction='in')
ax.grid(True, alpha = 0.5)
ax.legend(frameon=False, fontsize=12, loc=2)
ax.set_xlabel(r'log$_{10}$(M$_{\star}$/M$_{\odot}$)', fontsize=14)
ax.set_ylabel(r'log$_{10}$(Z)', fontsize=14)

cbaxes = fig.add_axes([0.75, 0.15, 0.02, 0.4])
fig.colorbar(s_m, cax=cbaxes)
cbaxes.set_ylabel(r'$z$', fontsize = 17)
cbaxes.set_yticks(np.arange(len(tags)))
cbaxes.set_yticklabels(np.arange(5,11))
cbaxes.invert_yaxis()
for label in (cbaxes.get_yticklabels()):
    label.set_fontsize(13)

plt.savefig('Z_star_gas.pdf', dpi=300, bbox_inches='tight')

plt.show()
