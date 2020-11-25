"""

    Plots the obscured and unobscured sfr relations:
    0 - SFR distribution function (Figure 16)
    1 - SFR density function (Figure 17)

"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from modules import get_lum_all, get_data_all
from FLARE.photom import lum_to_M, M_to_lum

import seaborn as sns
sns.set_context("paper")


def get_obs_sfrf():

    obs_df = {}

    obs_df['katsianis_bouwens'] = {}
    obs_df['katsianis_bouwens']['z'] = np.array([8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5])
    obs_df['katsianis_bouwens']['log10SFR'] = np.log10([43.269 ,21.704 ,10.891 ,5.469  ,2.850  ,1.803  ,0.902  ,0.359  ,73.533 ,41.186 ,
                                          23.070 ,12.921 ,7.239  ,4.235  ,2.674  ,1.687  ,0.534  ,0.171  ,141.748 ,77.951 ,
                                          42.862 ,23.585 ,12.974 ,7.132  ,3.921  ,1.860  ,0.742  ,0.309  ,382.081 ,208.215 ,
                                          113.462 ,61.828 ,33.695 ,18.369 ,10.001 ,5.452  ,2.974  ,1.280  ,0.512  ,0.203])

    obs_df['katsianis_bouwens']['phi'] = 1e-2 * np.array([0.0010 ,0.0026 ,0.0116 ,0.0120 ,0.0662 ,0.1066 ,0.2120 ,0.5480 ,
                                                         0.0002 ,0.0062 ,0.0090 ,0.0362 ,0.0578 ,0.1224 ,0.1697 ,0.3212 ,
                                                         1.0925 ,1.5901 ,0.0004 ,0.0028 ,0.0100 ,0.0330 ,0.0598 ,0.1305 ,
                                                         0.2330 ,0.3554 ,1.2496 ,2.5517 ,0.0004 ,0.0012 ,0.0063 ,0.0189 ,
                                                         0.0495 ,0.1270 ,0.1925 ,0.2486 ,0.3900 ,0.8343 ,1.6080 ,4.5640])

    obs_df['katsianis_bouwens']['sigma'] = 1e-2 * np.array([0.0006 ,0.0010 ,0.0030 ,0.0050 ,0.0208 ,0.0452 ,0.0680 ,0.2080 ,
                                                           0.0004 ,0.0017 ,0.0028 ,0.0064 ,0.0114 ,0.0187 ,0.0331 ,0.0894 ,
                                                           0.2731 ,0.5499 ,0.0004 ,0.0012 ,0.0024 ,0.0047 ,0.0077 ,0.0015 ,
                                                           0.0026 ,0.0598 ,0.2581 ,0.7857 ,0.0004 ,0.0006 ,0.0015 ,0.0026 ,
                                                           0.0047 ,0.0086 ,0.0125 ,0.0175 ,0.0319 ,0.0101 ,0.0331 ,0.0133])


    # correct Salpeter->Chabrier IMF
    obs_df['katsianis_bouwens']['log10SFR'] = np.log10(10**obs_df['katsianis_bouwens']['log10SFR'] * 0.63)


    name = 'mashian'
    ## NOTE: SFRs require recalibrating by a factor of 0.63
    ## due to the updated Kennicutt & Evans+12 calibrations
    out = {'z': [4.9,5.9,6.8,7.9],
           'log10phi*': [-3.25,-3.45,-3.67,-3.79],
           'log10SFR*': [1.75,1.62,1.54,1.31],
           'alpha': [-1.59,-1.62,-1.76,-1.79]
           }

    # correct Salpeter->Chabrier IMF
    # out['log10SFR*'] = np.log10(10**np.array(out['log10SFR*']) * 0.63)

    obs_df[name] = pd.DataFrame(out)


    name = 'smit12'

    obs_df[name] = pd.read_csv('Obs_data/smit12.csv',
                              delim_whitespace=True,
                              skiprows=3,
                              header=None,
                              names = ['z','log10SFR','phi','sigma'])

    # correct Salpeter->Chabrier IMF
    obs_df[name]['log10SFR'] = np.log10(10**obs_df[name]['log10SFR'] * 0.63)


    return obs_df



h = 0.6777
parent_volume = 3200**3
vol = (4/3)*np.pi*(14/h)**3

plt_options = ['SFRF', 'sfrd']
inp = plt_options[int(sys.argv[1])]

zs = [5., 6., 7., 8., 9., 10.]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']

df = pd.read_csv('weight_files/weights_grid.txt')
weights = np.array(df['weights'])



if inp == plt_options[0]:
    fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize=(10, 4), sharex=True, sharey=True, facecolor='w', edgecolor='k')
    axs = axs.ravel()
    obs_df = get_obs_sfrf()

    norm = matplotlib.colors.Normalize(vmin=0.5, vmax=len(zs)+0.5)
    # choose a colormap
    c_m = matplotlib.cm.viridis_r
    # create a ScalarMappable and initialize a data structure
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])

    bins = np.arange(-1, 4, 0.4)
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    for ii, tag in enumerate(tags):

        z = float(tag[5:].replace('p','.'))
        axs[ii].text(-0.5, -5.7, r'$z = {}$'.format(z), fontsize = 12)

        sfr_30 = get_data_all(tag, dataset = 'SFR/SFR_100', inp = 'FLARES', DF = False)
        L_FUV = get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='DustModelI')
        L_FUV_int = get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='Intrinsic')

        sfr_tot = np.zeros(len(bincen))
        sfr_tot_err = np.zeros(len(bincen))

        sfr_obsc = np.zeros(len(bincen))
        sfr_unobsc = np.zeros(len(bincen))

        sfr_obsc_err = np.zeros(len(bincen))
        sfr_unobsc_err = np.zeros(len(bincen))

        for jj in range(len(weights)):

            tmp, binedges = np.histogram(np.log10(sfr_30[jj]), bins=bins)
            sfr_tot+=tmp*weights[jj]
            sfr_tot_err+=np.square(np.sqrt(tmp)*weights[jj])

            this_frac = L_FUV[jj]/L_FUV_int[jj]

            tmp, binedges = np.histogram(np.log10(this_frac*sfr_30[jj]), bins=bins)
            sfr_unobsc+=tmp*weights[jj]
            sfr_unobsc_err+=np.square(np.sqrt(tmp)*weights[jj])

            tmp, binedges = np.histogram(np.log10((1-this_frac)*sfr_30[jj]), bins=bins)
            sfr_obsc+=tmp*weights[jj]
            sfr_obsc_err+=np.square(np.sqrt(tmp)*weights[jj])

        sfr_tot_err = np.sqrt(sfr_tot_err)
        sfr_unobsc_err = np.sqrt(sfr_unobsc_err)
        sfr_obsc_err = np.sqrt(sfr_obsc_err)


        sfr_tot = sfr_tot/(binwidth*vol)
        sfr_tot_err = sfr_tot_err/(vol*binwidth)
        yerr_tot = np.log10(sfr_tot) - np.log10(sfr_tot-sfr_tot_err), np.log10(sfr_tot+sfr_tot_err) - np.log10(sfr_tot)

        sfr_unobsc = sfr_unobsc/(binwidth*vol)
        sfr_unobsc_err = sfr_unobsc_err/(vol*binwidth)
        yerr_unobsc = np.log10(sfr_unobsc) - np.log10(sfr_unobsc-sfr_unobsc_err), np.log10(sfr_unobsc+sfr_unobsc_err) - np.log10(sfr_unobsc)

        sfr_obsc = sfr_obsc/(binwidth*vol)
        sfr_obsc_err = sfr_obsc_err/(vol*binwidth)
        yerr_obsc = np.log10(sfr_obsc) - np.log10(sfr_obsc-sfr_obsc_err), np.log10(sfr_obsc+sfr_obsc_err) - np.log10(sfr_obsc)

        if z!=10:
            axs[ii].errorbar(bincen, np.log10(sfr_tot), ls='solid', color=s_m.to_rgba(ii+0.5))
            axs[ii].errorbar(bincen, np.log10(sfr_unobsc), yerr = yerr_unobsc, ls='dashed', color=s_m.to_rgba(ii+0.5))
            axs[ii].errorbar(bincen, np.log10(sfr_obsc), yerr = yerr_obsc, ls='dotted', color=s_m.to_rgba(ii+0.5))
        else:
            axs[ii].errorbar(bincen, np.log10(sfr_tot), ls='solid', color=s_m.to_rgba(ii+0.5), label='Total')
            axs[ii].errorbar(bincen, np.log10(sfr_unobsc), yerr = yerr_unobsc, ls='dashed', color=s_m.to_rgba(ii+0.5), label='Unobscured/UV')
            axs[ii].errorbar(bincen, np.log10(sfr_obsc), yerr = yerr_obsc, ls='dotted', color=s_m.to_rgba(ii+0.5), label='Obscured/IR')



        ## ---- Observations
        for author,label,ms in zip(['smit12','katsianis_bouwens'],['Smit+12','Katsianis+17'],['o','d','s','<']):

            mask = (obs_df[author]['z'] < (z + 0.5)) & (obs_df[author]['z'] > (z - 0.5))

            if np.sum(mask) > 0:
                phi = obs_df[author]['phi'][mask]
                lo = np.log10(phi) - np.log10(phi - obs_df[author]['sigma'][mask])
                hi = np.log10(phi + obs_df[author]['sigma'][mask]) - np.log10(phi)

                axs[ii].errorbar(obs_df[author]['log10SFR'][mask],
                            np.log10(obs_df[author]['phi'][mask]),
                            yerr=[lo,hi], ls='none',marker=ms,color='grey',label=label, alpha=0.3)


        axs[ii].grid(True, alpha=0.6)
        axs[ii].set_xlim((-1,3.5))
        axs[ii].minorticks_on()
        axs[ii].tick_params(axis='x', which='minor', direction='in')
        axs[ii].tick_params(axis='y', which='minor', direction='in')
        for label in (axs[ii].get_xticklabels() + axs[ii].get_yticklabels()):
            label.set_fontsize(12)

    axs[0].legend(frameon=False, fontsize=10)
    axs[-1].legend(frameon=False, fontsize=10)

    fig.subplots_adjust(bottom=0.1, left=0.08, wspace=0, hspace=0)
    fig.text(0.02, 0.5, r'$\mathrm{log}_{10}(\Phi/(\mathrm{cMpc}^{-3}\mathrm{dex}^{-1}))$', va='center', rotation='vertical', fontsize=14)
    fig.text(0.41, 0.001, r'$\mathrm{log}_{10}(\mathrm{SFR}/\mathrm{M}_{\odot}\mathrm{yr}^{-1})$', va='center', fontsize=14)
    plt.savefig(F"sfr_obsc_unobsc.pdf", bbox_inches='tight', dpi=300)
    plt.show()

elif inp == plt_options[1]:

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(4, 3), sharex=True, sharey=True, facecolor='w', edgecolor='k')

    sfrd_obsc = np.zeros(len(tags))
    sfrd_unobsc = np.zeros(len(tags))
    sfrd_tot = np.zeros(len(tags))
    zs = np.zeros(len(tags))

    for ii, tag in enumerate(tags):

        sfr_30 = get_data_all(tag, dataset = 'SFR/SFR_100', inp = 'FLARES', DF = False)
        L_FUV = get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='DustModelI')
        L_FUV_int = get_lum_all(tag, LF = False, filter = 'FUV', Luminosity='Intrinsic')


        zs[ii]=float(tag[5:].replace('p','.'))

        for jj in range(len(weights)):
            tmp = sfr_30[jj]
            ok = np.where(tmp>=0.1)[0]
            #tmp = tmp/(8e27)
            att = L_FUV[jj]/L_FUV_int[jj]
            sfrd_tot[ii]+=np.sum((tmp[ok]*weights[jj])/vol)
            sfrd_unobsc[ii]+=np.sum((tmp[ok]*att[ok]*weights[jj])/vol)
            sfrd_obsc[ii]+=np.sum((tmp[ok]*(1.-att[ok])*weights[jj])/vol)

    ax.plot(zs, np.log10(sfrd_tot), marker='o', color='black', label='Total')
    ax.plot(zs, np.log10(sfrd_unobsc), marker='o', color='green', label='Unobscured/UV')
    ax.plot(zs, np.log10(sfrd_obsc), marker='o', color='red', label='Obscured/IR')

    print ("Obscured fraction = ", sfrd_obsc/sfrd_tot)

    ax.errorbar([4.9, 5.9, 6.8, 7.9, 10.4], [-1.85, -2.05, -2.17, -2.48, -3.28], yerr=[[0.06,0.06,0.06,0.07,0.45], [0.06,0.06,0.06,0.07,0.36]], marker='s', ls = 'None', alpha=0.5, color='green', label='Bouwens+2020 (Unobscured)')

    #ax.errorbar([5.25], np.log10([7.46e-2]), xerr=[0.75], yerr=[[np.log10([7.46e-2])-np.log10([4.71e-2])], [np.log10([1.36e-1])-np.log10([7.46e-2])]], color='red', alpha=0.5, label='Gruppioni+2020 (sub-mm)')
    #ax.plot(zs, np.log10(0.015 * ((1+zs)**2.7)/(1+((1+zs)/2.9)**5.6)), ls='dashed', color='black', label=r'Madua $\&$ Dickinson 2014')

    ax.errorbar([4.5, 5.5], [-2.44, -2.67], yerr= [0.25,0.25], lolims=[1,1], marker='s', ls = 'None', color='red', alpha=0.5, label='Khusanova+2020')

    for label in (ax.get_xticklabels()+ax.get_yticklabels()):
        label.set_fontsize(11)
    ax.legend(frameon=False, fontsize=10, loc=3)
    ax.minorticks_on()
    ax.tick_params(axis='x', which='minor', direction='in')
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.grid(True, alpha = 0.5)
    ax.set_ylim(-4.5,-1)
    ax.set_xlim(4.3,10.5)

    ax.set_xlabel(r'$z$', fontsize=14)
    ax.set_ylabel(r'$\mathrm{log}_{10}(\mathrm{SFRD}/(M_{\odot}\mathrm{yr}^{-1}\mathrm{Mpc}^3))$', fontsize=12)
    plt.savefig('sfrd.pdf', bbox_inches='tight', dpi=300)
    plt.show()
