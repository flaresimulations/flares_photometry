import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from functools import partial
import schwimmbad
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import FLARE.filters


def get_data(ii, tag, inp='GEAGLE', filter = 'FUV', Luminosity = 'Dustcorr'):

    if inp == 'GEAGLE':

        num = str(ii)

        if len(num) == 1:
            num =  '0'+num

        photfile = '../photometry/out/GEAGLE_{}_sp_info.hdf5'.format(num)
        origfile = '../los/data1/GEAGLE_{}_sp_info.hdf5'.format(num)

    with h5py.File(origfile,'r') as hf:
        glen = np.array(hf[F"{tag}/Subhalo/G_Length"])
        slen = np.array(hf[F"{tag}/Subhalo/S_Length"])

    with h5py.File(photfile,'r') as hf:
        if np.isscalar(filter):
            lum = np.array(hf[F"{tag}/Subhalo/BPASS/SalpeterIMF/ModelI/Luminosity/{Luminosity}/{filter[8:]}"])
        else:
            lum = np.zeros((len(glen), len(filter)))
            for ii, jj in enumerate(filter):
                lum[:,ii] = np.array(hf[F"{tag}/Subhalo/BPASS/SalpeterIMF/ModelI/Luminosity/{Luminosity}/{jj[8:]}"])


    return lum, glen+slen


def get_all(tag, inp = 'GEAGLE', filter = 'FUV', Luminosity = 'Dustcorr'):

    if inp == 'GEAGLE':

        sims = np.arange(0,38)

        calc = partial(get_data, tag = tag, inp = inp, filter = filter, Luminosity = Luminosity)

        pool = schwimmbad.MultiPool(processes=12)
        dat = np.array(list(pool.map(calc, sims)))
        pool.close()

        return dat


###Defining input values
fac = 0.0315
h = 0.6777
z = [5, 6, 7, 8, 9, 10]
tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075', '003_z008p988', '002_z009p993']
filters = FLARE.filters.TH[:-1]

df_filter_lims = pd.DataFrame({'z': z})
percentile = np.empty((len(z), len(filters)))


for ii, jj in enumerate(z):

    data = get_all(tags[ii], Luminosity = 'Dustcorr', filter = filters)

    for kk in range(38):

        if kk == 0:
            Mlum = data[kk][0]
            part = data[kk][1]
        else:
            Mlum = np.append(Mlum, data[kk][0], axis = 0)
            part = np.append(part, data[kk][1])

    ok = np.where(part == 100)

    percentile[ii] = np.percentile(Mlum[ok], 5, axis = 0)

for ii, jj in enumerate(filters):

    df_filter_lims[jj[8:]] = percentile[:, ii]

print (df_filter_lims)

df_filter_lims.to_csv('Magnitude_limits.txt')
