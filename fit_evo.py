import numpy as np
from collections import OrderedDict
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from modules import get_fit_params

import seaborn as sns
sns.set_context("paper")

zs = [5., 6., 7., 8., 9., 10.]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']

fig, axs = plt.subplots(nrows = 2, ncols = 2, figsize=(5,5), sharex=True, sharey=False, facecolor='w', edgecolor='k')
axs = axs.ravel()

Schechter = ['log10phi*', 'D*', 'alpha']
DPL = ['log10phi*', 'D*', 'alpha_1', 'alpha_2']


for ii, z in enumerate(zs):

    median_fit, fit84, fit16 = get_fit_params('Schechter', z)
    for jj, kk in enumerate(Schechter):
        axs[jj].errorbar(z, median_fit[kk], yerr=[[fit16[kk]],[fit84[kk]]], ls='', marker='o', color = 'black')

    median_fit, fit84, fit16 = get_fit_params('DPL', z)
    for jj, kk in enumerate(DPL):
        axs[jj].errorbar(z, median_fit[kk], yerr=[[fit16[kk]],[fit84[kk]]], ls='', marker='s', color = 'grey')

axs[1].yaxis.set_label_position("right")
axs[1].yaxis.tick_right()
axs[-1].yaxis.set_label_position("right")
axs[-1].yaxis.tick_right()
axs[2].set_xlabel(r'$z$', fontsize=14)
axs[3].set_xlabel(r'$z$', fontsize=14)
ylabels = [r'log$_{10}(\phi^*)$', r'M$^*$', r'$\alpha$', r'$\beta$']
for ii in range(4):
    axs[ii].set_ylabel(ylabels[ii], fontsize=14)
    axs[ii].grid(True, alpha=0.6)
    axs[ii].minorticks_on()
    axs[ii].tick_params(axis='x', which='minor', direction='in')
    axs[ii].tick_params(axis='y', which='minor', direction='in')
    for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
        label.set_fontsize(13)

fig.subplots_adjust(wspace=0, hspace=0)
plt.savefig('fit_params_evo.pdf', dpi=300, bbox_inches='tight')
plt.show()
