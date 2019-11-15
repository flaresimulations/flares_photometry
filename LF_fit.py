import numpy as np
import pandas as pd
import schwimmbad

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

from functools import partial
import h5py
import scipy
import matplotlib.pyplot as plt
from FLARE.photom import lum_to_M, M_to_lum
from plot_obs import plot_obs

def fit(M, theta):

    log10phistar, alpha, Mstar = theta
    delta = M - Mstar

    return np.log10(0.4*np.log(10)) + log10phistar - 0.4*delta*(alpha+1.) - (np.log10(np.e)) * (10**(-0.4*delta))


def get_hist(ii, tag, bins, inp='GEAGLE', filter = 'FUV'):

    if inp == 'GEAGLE':

        num = str(ii)

        if len(num) == 1:
            num =  '0'+num

        filename = '../photometry/out/GEAGLE_{}_sp_info.hdf5'.format(num)


    with h5py.File(filename,'r') as hf:

        lum = np.array(hf[F"{tag}/Subhalo/BPASS/SalpeterIMF/ModelI/Luminosity/Dustcorr/{filter}"])

        tmp, edges = np.histogram(lum_to_M(lum), bins = bins)

        return tmp


def get_all(tag, bins = np.arange(-25, -16, 0.5), inp = 'GEAGLE', filter = 'FUV'):

    if inp == 'GEAGLE':

        sims = np.arange(0,38)

        df = pd.read_csv('weight_files/weights_grid.txt')
        weights = np.array(df['weights'])

        calc = partial(get_hist, tag = tag, bins = bins, inp = inp, filter = filter)

        pool = schwimmbad.MultiPool(processes=12)
        dat = np.array(list(pool.map(calc, sims)))
        pool.close()


        hist = np.sum(dat, axis = 0)
        out = np.zeros(len(bins)-1)
        for ii, sim in enumerate(sims):

            out+=dat[ii]*weights[ii]

        return out, hist

    else:

        out = get_hist(00, tag, bins, inp = 'REF')

        return out


# plt.style.use('simple')
ID = 'fit_out'
h=0.6777

zs = [5, 6, 7, 8, 9, 10]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']

bins = np.arange(-25, -17., 0.4)#np.arange(np.log10(M_to_lum(-17.7)), np.log10(M_to_lum(-23.4)), 0.2)
bincen = (bins[1:]+bins[:-1])/2.
binwidth = bins[1:] - bins[:-1]
vol = (4/3)*np.pi*(14/h)**3

out, hist = get_all(tags[3], bins)

phi = out/(vol*binwidth)

N_sample = phi*(binwidth*(3200**3))

observations = [{'bin_edges': bins, 'N': N_sample, 'volume': 3200**3}]


print(observations)

# ----- Define Priors manually...
model = models.Schechter_Mags()

priors = {}
priors['log10phi*'] = scipy.stats.uniform(loc = -5., scale = 4.)
priors['alpha'] = scipy.stats.uniform(loc = -2.5, scale = 1)
priors['D*'] = scipy.stats.uniform(loc = -23, scale = 4.)



# model = models.DoubleSchechter_Mags()
#
# priors = {}
# priors['log10phi*_1'] = scipy.stats.uniform(loc = -4, scale = 7.0)
# priors['alpha_1'] = scipy.stats.uniform(loc = -2.0, scale = 3.0)
# priors['log10phi*_2'] = scipy.stats.uniform(loc = -5, scale = 7.0)
# priors['alpha_2'] = scipy.stats.uniform(loc = -3.0, scale = 3.0)
# priors['D*'] = scipy.stats.uniform(loc = -23., scale = 5.0)

# -------------------- fit sampled LF and plot median fit

fitter = fitDF.fitter(observations, model=model, priors=priors, output_directory = ID)
fitter.fit(nwalkers = 50, nsamples = 500, burn = 100)
# fitter.fit(nsamples = 2000, burn = 500, sample_save_ID = 'a_different_ID') # to save the samples as something other than samples.p


# -------------------- make simple analysis plots

a = analyse.analyse(ID = ID, model = model, observations=observations)
fig = a.triangle(hist2d = True, ccolor='0.5')
plt.savefig(f'{fit_out}/triangle.png')
plt.close()

fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(8, 8), sharex=True, sharey=True,
                            facecolor='w', edgecolor='k')

axs.plot(bincen, fit(bincen, model.sp.values()), label = 'fit')
axs.plot(bincen, np.log10(phi), label = 'data')
plot_obs(z, axs)
axs.legend(fontsize = 15, frameon = False)
plt.savefig(f'{fit_out}/fit_fitdf.png')
plt.show()
