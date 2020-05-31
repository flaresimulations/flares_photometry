import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from FLARE.photom import lum_to_M, M_to_lum
import FLARE.filters
import flares
from modules import get_data_all, get_line_all, get_lum_all

#de Barros+19 data
def deBarros_sch(logL):

    log10phistar, alpha, log10Lstar = [-4.05, -2.22, 43.45]
    delta = logL - log10Lstar

    phi = np.log(10) * 10**log10phistar * np.exp(-10**delta) * 10**(delta*(alpha + 1))

    return np.log10(phi)


tag = '008_z007p000'
filters = 'FUV'
h = 0.6777
vol = (4/3)*np.pi*(14/h)**3
quantiles = [0.84,0.50,0.16]

bins = np.arange(40, 47, 0.4)
bincen = (bins[1:]+bins[:-1])/2.
binwidth = bins[1:] - bins[:-1]

df = pd.read_csv('weight_files/weights_grid.txt')
weights = np.array(df['weights'])
sims = np.arange(len(weights))

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(5, 3), sharex=False, sharey=False, facecolor='w', edgecolor='k')
c = 'indianred'

dat1 = get_line_all(tag, 'OIII4959', inp = 'FLARES', LF = False)
dat2 = get_line_all(tag, 'OIII5007', inp = 'FLARES', LF = False)
dat3 = get_line_all(tag, 'HI4861', inp = 'FLARES', LF = False)

LFOIII = np.zeros(len(bincen))
err = np.zeros(len(bincen))
LFOIII_tot = np.zeros(len(bincen))
for jj, sim in enumerate(sims):
    tmp, edges = np.histogram(np.log10(dat1[:,0][jj] + dat2[:,0][jj] + dat3[:,0][jj]), bins = bins)

    err+=np.square(np.sqrt(tmp)*weights[jj])
    LFOIII+=tmp*weights[jj]
    LFOIII_tot+=tmp

LFOIII = LFOIII/(binwidth*vol)
err = np.sqrt(err)/(binwidth*vol)
axs.plot(bincen, np.log10(LFOIII), color='black', ls = 'solid')
axs.fill_between(bincen, np.log10(LFOIII-err), np.log10(LFOIII+err), color='black', alpha=0.15)

L, phi16, phi84 = np.loadtxt('Obs_data/OIII_LF_z8_wSigmaInt_1SigmaContour.txt').T
axs.fill_between(L, phi16, phi84, color='indianred', alpha=0.5, label=r'de Barros 2019')
axs.plot(bincen, deBarros_sch(bincen), color='indianred', ls='solid')



for label in (axs.get_xticklabels()+axs.get_yticklabels()):
    label.set_fontsize(12)
axs.legend(frameon=False, fontsize=12, loc=1)
axs.minorticks_on()
axs.tick_params(axis='x', which='minor', direction='in')
axs.tick_params(axis='y', which='minor', direction='in')
axs.grid(True, alpha = 0.5)

axs.set_xlim(42, 44.2)
axs.set_ylim(-6.,-1.8)
axs.set_xticks(np.arange(42, 44.1, 0.5))

axs.set_xlabel(r'log$_{10}$(L([OIII]+H$\beta$)/(erg s$^{-1}$))', fontsize=13)
axs.set_ylabel(r'log$_{10}$($\phi$/(cMpc$^{-3}$dex$^{-1}$))', fontsize=13)

plt.savefig('LFOIII_z8.pdf', dpi=300, bbox_inches='tight')
plt.show()
