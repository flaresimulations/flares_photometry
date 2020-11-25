"""

Plots the evolution of fit parameters; Figure 6

"""

import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from modules import get_fit_params
from FLARE.LF import evo

import seaborn as sns
sns.set_context("paper")


def get_params(label):

    if label == 'Finkelstein 2015':
        redshifts = np.array([5, 6, 7, 8])
        phi = np.array([8.95e-4, 1.86e-4, 1.57e-4, 0.72e-4])
        logphi_err = np.array([[1.31e-4, 0.80e-4, 0.95e-4, 0.65e-4], [1.92e-4, 0.94e-4, 1.49e-4, 2.52e-4]])
        logphi_err[0], logphi_err[1] = np.log10(phi/(phi-logphi_err[0])), np.log10((phi+logphi_err[1])/phi)
        logphi = np.log10(phi)
        Mstar = np.array([-20.81, -21.13, -21.03, -20.89])
        Mstar_err = np.array([[0.12, 0.31, 0.50, 1.08], [0.13, 0.25, 0.37, 0.74]])
        alpha = np.array([-1.67, -2.02, -2.03, -2.36])
        alpha_err = np.array([[0.06, 0.10, 0.20, 0.40], [0.05, 0.10, 0.21, 0.54]])

        theta = np.array([redshifts, logphi, Mstar, alpha])
        theta_err = np.array([logphi_err, Mstar_err, alpha_err])

    elif label == 'Bouwens 2015':
        redshifts = np.array([4.9, 5.9, 6.8, 7.9, 10.4])
        phi = np.array([0.74, 0.50, 0.29, 0.21, 0.008])*1e-4
        logphi_err = np.array([[0.14, 0.16, 0.12, 0.11, 0.003], [0.18, 0.22, 0.21, 0.23, 0.004]])*1e-4
        logphi_err[0], logphi_err[1] = np.log10(phi/(phi-logphi_err[0])), np.log10((phi+logphi_err[1])/phi)
        logphi = np.log10(phi)
        Mstar = np.array([-21.17, -20.94, -20.87, -20.63, -20.92])
        Mstar_err = np.array([[0.12, 0.20, 0.26, 0.36, 0.0], [0.12, 0.20, 0.26, 0.36, 0.0]])
        alpha = np.array([-1.76, -1.87, -2.06, -2.02, -2.27])
        alpha_err = np.array([[0.05, 0.10, 0.13, 0.23, 0.], [0.05, 0.10, 0.13, 0.23, 0.]])

        theta = np.array([redshifts, logphi, Mstar, alpha])
        theta_err = np.array([logphi_err, Mstar_err, alpha_err])

    elif label == 'Illustris TNG':
        redshifts = np.array([5, 6, 7, 8])
        logphi = np.array([-3.398, -3.608, -4.209, -4.714])
        Mstar = np.array([-21.21, -21.31, -21.47, -21.44])
        alpha = np.array([-1.941, -2.042, -2.279, -2.455])

        theta = np.array([redshifts, logphi, Mstar, alpha])

    elif label in ['bluetides', 'Mason15', 'Yung2018', 'Ma2019']:
        this = eval(F'evo.{label}()')

        theta = np.array([this.redshifts, this.phi_star, this.M_star, this.alpha])

    elif label == 'Bowler':
        redshifts = np.array([5, 6, 7, 8, 9])
        # logphi = np.log10(np.array([2.5e-4, 1.9e-4, 2.2e-4, 4.83e-4, 2.85e-4]))
        # Mstar = np.array([-21.40, -21.20, -20.61, -19.80, -19.67])
        # alpha = np.array([-2.00, -2.10, -2.19, -1.96, -2.10])
        # beta = np.array([-4.8, -5.1, -4.6, -3.98, -3.75])

        phi = np.array([3.24e-4, 1.93e-4, 3.81e-4, 4.83e-4, 2.85e-4])
        phi_err = np.array([0.40, 0.45, 1.36, 2.25, 1.39])*1e-4
        logphi_err = np.array([np.log10(phi/(phi-phi_err)), np.log10((phi+phi_err)/phi)])
        logphi = np.log10(phi)
        Mstar = np.array([-21.48, -21.25, -20.27, -19.80, -19.67])
        Mstar_err = np.array([0.06, 0.10, 0.19, 0.26, 0.33])
        alpha = np.array([-1.88, -2.07, -2.13, -1.96, -2.10])
        alpha_err = np.array([0.04, 0.07, 0.08, 0.15, 99.0])
        beta = np.array([-5.3, -5.23, -4.23, -3.98, -3.75])
        beta_err = np.array([0.19, 0.18, 0.17, 0.14, 0.22])

        theta = np.array([redshifts, logphi, Mstar, alpha, beta])
        theta_err = np.array([logphi_err, Mstar_err, alpha_err, beta_err])


    return theta



zs = [5., 6., 7., 8., 9., 10.]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']

fig, axs = plt.subplots(nrows = 4, ncols = 1, figsize=(3,6), sharex=True, sharey=False, facecolor='w', edgecolor='k')
axs = axs.ravel()

Schechter = ['log10phi*', 'D*', 'alpha']
DPL = ['log10phi*', 'D*', 'alpha_1', 'alpha_2']

obs = ['Bouwens 2015', 'Finkelstein 2015', 'bluetides', 'Illustris TNG', 'Mason15', 'Yung2018', 'Ma2019', 'Bowler']
labels = ['Bouwens+2015', 'Finkelstein+2015', '\\textsc{BlueTides}', '\\textsc{Illustris} \\textsc{Tng}', 'Mason+2015', 'Yung+2018', 'Ma+2019', 'Bowler+2020']
colors = ['green', 'blue', 'red', 'violet', 'yellow', 'magenta', 'orange', 'brown']

# obs = ['Bouwens 2015', 'Finkelstein 2015', 'Bowler']
# labels = ['Bouwens+2015', 'Finkelstein+2015', 'Bowler+2020']

for ii, z in enumerate(zs):

    median_fit, fit84, fit16 = get_fit_params('Schechter', z)
    for jj, kk in enumerate(Schechter):
        axs[jj].errorbar(z, median_fit[kk], yerr=[[fit16[kk]],[fit84[kk]]], ls='', marker='D', color = 'black', markersize=6, alpha=0.6)

    median_fit, fit84, fit16 = get_fit_params('DPL', z)
    for jj, kk in enumerate(DPL):
        axs[jj].errorbar(z, median_fit[kk], yerr=[[fit16[kk]],[fit84[kk]]], ls='', marker='s', color = 'grey', markersize=6, alpha=0.6)

for ii, jj in enumerate(obs):
    theta = get_params(jj)
    z = theta[0]
    ok = np.logical_and(z>=4.9, z<=10.5)
    if jj != 'Bowler':
        for kk, ll in enumerate(Schechter):
            axs[kk].scatter(z[ok], theta[1+kk][ok], label=labels[ii], s=25, color=colors[ii], marker='o', alpha=0.6)
            # yerr=[theta_err[kk][0][ok], theta_err[kk][1][ok]]
            #
            # axs[kk].errorbar(z[ok], theta[1+kk][ok], yerr = yerr, label=labels[ii], markersize=6, color=colors[ii], marker='o', alpha=0.6)

    else:
        for kk, ll in enumerate(DPL):
            axs[kk].scatter(z[ok], theta[1+kk][ok], label=labels[ii], s=25, color=colors[ii], marker='o', alpha=0.6)
            # if ll == 'log10phi*':
            #     yerr = [theta_err[kk][0][ok], theta_err[kk][1][ok]]
            # else:
            #     yerr = theta_err[kk][ok]
            #
            #
            # axs[kk].errorbar(z[ok], theta[1+kk][ok], yerr = yerr, label=labels[ii], markersize=6, color=colors[ii], marker='o', alpha=0.6)


# axs[1].yaxis.set_label_position("right")
# axs[1].yaxis.tick_right()
# axs[-1].yaxis.set_label_position("right")
# axs[-1].yaxis.tick_right()
# axs[2].set_xlabel(r'$z$', fontsize=14)
axs[3].set_xlabel(r'$z$', fontsize=14)
ylabels = [r'log$_{10}(\phi^*)$', r'M$^*$', r'$\alpha$', r'$\beta$']
for ii in range(4):
    axs[ii].set_xlim(4.7,10.5)
    axs[1].set_ylim(-22.5,-18.8)
    axs[2].set_ylim(-3.2,-1.5)
    axs[ii].set_ylabel(ylabels[ii], fontsize=12)
    axs[ii].grid(True, alpha=0.6)
    axs[ii].minorticks_on()
    axs[ii].tick_params(axis='x', which='minor', direction='in')
    axs[ii].tick_params(axis='y', which='minor', direction='in')
    for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
        label.set_fontsize(10)

lines, labels = [], []
for ax in fig.axes:
    axLine, axLabel = ax.get_legend_handles_labels()
    lines.extend(axLine)
    labels.extend(axLabel)

unique = np.unique(np.array(labels))
unique_lines = []
for ii, jj in enumerate(unique):
    unique_lines.append(lines[np.where(np.array(labels)==jj)[0][0]])
fig.subplots_adjust(wspace=0, hspace=0, top=0.9)
axs[0].legend(unique_lines, unique, bbox_to_anchor=(-0.3, 0.97), loc = 'lower left', ncol=2, fontsize=10, frameon=False)

plt.savefig('fit_params_evo.pdf', dpi=300, bbox_inches='tight')
plt.show()
