import numpy as np
import pandas as pd
import schwimmbad
from functools import partial
import h5py
from FLARE.photom import lum_to_M, M_to_lum
import fitDF.models as models
from collections import OrderedDict

def histogram(edges, counts):

    left,right = edges[:-1],edges[1:]
    X = np.array([left,right]).T.flatten()
    Y = np.array([counts,counts]).T.flatten()

    return X, Y

def plot_posteriors(data):

    import corner, json
    s = []
    with open(data) as f:
        samples = json.load(f)

    if 'DPL' in data:
        tmp1 = np.array(samples['alpha_1'])
        tmp2 = np.array(samples['alpha_2'])
        ok = np.where(tmp2>tmp1)
        tmp = tmp2[ok]
        tmp2[ok]= tmp1[ok]
        tmp1[ok]=tmp
        samples['alpha_2']=tmp2
        samples['alpha_1']=tmp1

    for ii in samples.keys():
        s = np.vstack(np.array([samples[ii] for ii in samples.keys()]).T)
    figure = corner.corner(s)

def Schechter(M, theta):

    log10phistar, alpha, Mstar = theta
    delta = M - Mstar

    return np.log10(0.4*np.log(10)) + log10phistar - 0.4*delta*(alpha+1.) - (np.log10(np.e)) * (10**(-0.4*delta))

def DPL(M, theta):

    log10phistar, alpha_1, alpha_2, Mstar = theta
    delta = M - Mstar

    return -np.log10(10**(0.4*(delta)*(alpha_1+1)) + 10**(0.4*(delta)*(alpha_2+1))) + log10phistar


def get_data(ii, tag, dataset = 'Mstar_30', bins = np.arange(7.5,12,0.5), DF=False, inp = 'FLARES'):

    if inp == 'FLARES':

        num = str(ii)
        if ii/10 < 1: num = '0'+num

        sim = rF"./data/flares.hdf5"
        num = num+'/'


    else:

        sim = F"./data/EAGLE_{inp}_sp_info.hdf5"
        num = ''

    with h5py.File(sim, 'r') as hf:

        data = np.array(hf[F'{num}{tag}/Galaxy'].get(dataset))
        if DF == True:

            tmp, edges = np.histogram(np.log10(data*1e10), bins = bins)

            return tmp

        else:return data


def get_lum(ii, tag, bins = np.arange(-26,-16,0.5), inp='FLARES', filter = 'FUV', LF = True, Luminosity='DustModelI'):

    if inp == 'FLARES':

        num = str(ii)

        if len(num) == 1:
            num =  '0'+num

        filename = rF"./data/flares.hdf5"
        num = num+'/'


    else:

        filename = F'./data/EAGLE_{inp}_sp_info.hdf5'
        num = ''

    with h5py.File(filename,'r') as hf:

        lum = np.array(hf[F"{num}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{Luminosity}/{filter}"])

    if LF == True:

        tmp, edges = np.histogram(lum_to_M(lum), bins = bins)

        return tmp

    else: return lum


def get_line(ii, tag, line = 'HI6563', inp = 'FLARES', LF = False, bins = np.arange(42,46,0.5), Type = 'DustModelI'):

    if inp == 'FLARES':

        num = str(ii)

        if len(num) == 1:
            num =  '0'+num

        filename = './data/flares.hdf5'
        num = num+'/'

    else:

        filename = F'./data/EAGLE_{inp}_sp_info.hdf5'
        num = ''

    with h5py.File(filename,'r') as hf:

        if LF:
            lum = np.array(hf[F"{num}{tag}/Galaxy/BPASS_2.2.1/Chabrier300/Lines/{Type}/{line}/Luminosity"])
            tmp, edges = np.histogram(lum_to_M(lum), bins = bins)
            return tmp

        else:
            lum = np.array(hf[F"{num}{tag}/Galaxy/BPASS_2.2.1/Chabrier300/Lines/{Type}/{line}/Luminosity"])
            EW = np.array(hf[F"{num}{tag}/Galaxy/BPASS_2.2.1/Chabrier300/Lines/{Type}/{line}/EW"])
            return lum, EW


def get_data_all(tag, dataset = 'Mstar_30', bins = np.arange(7.5, 12, 0.5), inp = 'FLARES', DF = False):

    if inp == 'FLARES':

        df = pd.read_csv('weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        sims = np.arange(0,len(weights))

        calc = partial(get_data, tag = tag, dataset = dataset, bins = bins, inp = inp, DF = DF)

        pool = schwimmbad.MultiPool(processes=12)
        dat = np.array(list(pool.map(calc, sims)))
        pool.close()

        if DF:
            hist = np.sum(dat, axis = 0)
            out = np.zeros(len(bins)-1)
            err = np.zeros(len(bins)-1)

            for ii, sim in enumerate(sims):
                out+=dat[ii]*weights[ii]
                err+=np.square(np.sqrt(dat[ii])*weights[ii])

            return out, hist, np.sqrt(err)

        else: return dat

    else:

        out = get_data(00, tag = tag, bins = bins, inp = inp, DF = DF)

        return out


def get_lum_all(tag, bins = np.arange(-25, -16, 0.5), inp = 'FLARES', LF = True, filter = 'FUV', Luminosity='DustModelI'):

    if inp == 'FLARES':

        df = pd.read_csv('weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        sims = np.arange(0,len(weights))

        calc = partial(get_lum, tag = tag, bins = bins, inp = inp, LF = LF, filter = filter, Luminosity = Luminosity)

        pool = schwimmbad.MultiPool(processes=12)
        dat = np.array(list(pool.map(calc, sims)))
        pool.close()

        if LF:
            hist = np.sum(dat, axis = 0)
            out = np.zeros(len(bins)-1)
            err = np.zeros(len(bins)-1)
            out_up = np.zeros(len(bins)-1)
            out_low = np.zeros(len(bins)-1)
            for ii, sim in enumerate(sims):
                out+=dat[ii]*weights[ii]
                err+=np.square(np.sqrt(dat[ii])*weights[ii])

            return out, hist, np.sqrt(err)

        else: return dat

    else:

        out = get_lum(00, tag = tag, bins = bins, inp = inp, filter = filter, Luminosity = Luminosity)

        return out


def get_line_all(tag, line, inp = 'FLARES', LF = True, bins = np.arange(40, 46, 0.5), Type='DustModelI'):

    if inp == 'FLARES':

        df = pd.read_csv('weight_files/weights_grid.txt')
        weights = np.array(df['weights'])
        sims = np.arange(0,len(weights))

        calc = partial(get_line, tag = tag, line = line, inp = inp, LF = LF, bins = bins, Type = Type)

        pool = schwimmbad.MultiPool(processes=12)
        dat = np.array(list(pool.map(calc, sims)))
        pool.close()

        if LF:
            hist = np.sum(dat, axis = 0)
            out = np.zeros(len(bins)-1)
            err = np.zeros(len(bins)-1)

            for ii, sim in enumerate(sims):
                out+=dat[ii]*weights[ii]
                err+=np.square(np.sqrt(dat[ii])*weights[ii])

            return out, hist, np.sqrt(err)

        else: return dat

    else:

        out = get_line(00, tag = tag, line = line, inp = inp, LF = LF, bins = bins, Type = Type)


        return out


# Function to output the number density for the functional form
def fit_function(model, observed, sigma, bins, z):
    bincen = (bins[1:]+bins[:-1])/2.
    binwidth = bins[1:] - bins[:-1]

    if model == 'Schechter':
        yy = partial(Schechter, M=bincen)
    elif model == 'DPL':
        yy = partial(DPL, M=bincen)

    median_fit, fit84, fit16 = get_fit_params(model, z)
    BIC_value(model, observed, sigma, bins, median_fit)

    return yy(theta=median_fit.values())

def get_fit_params(model, z):

    if model == 'Schechter':
        folder = 'fit_Sch'
    elif model == 'DPL':
        folder = 'fit_DPL'
    import json
    with open(F"{folder}/{folder}_z{int(z)}.json") as f:
        samples = OrderedDict(json.load(f))


    if model == 'DPL':
        tmp1 = np.array(samples['alpha_1'])
        tmp2 = np.array(samples['alpha_2'])
        ok = np.where(tmp2>tmp1)
        tmp = tmp2[ok]
        tmp2[ok]= tmp1[ok]
        tmp1[ok]=tmp
        samples['alpha_2']=tmp2
        samples['alpha_1']=tmp1

    median_fit, fit84, fit16 = {}, {}, {}
    for ip, p in enumerate(samples.keys()):

        median_fit[p] = np.percentile(samples[p], 50)
        fit84[p] = np.percentile(samples[p], 84) - np.percentile(samples[p], 50)
        fit16[p] = np.percentile(samples[p], 50) - np.percentile(samples[p], 16)

    print ("median fit: ", median_fit)
    print ('84 percentile:', fit84)
    print ('16 percentile:', fit16)

    # np.savez(F'fit_{model}_{z}', median_fit=median_fit, fit84=fit84, fit16=fit16)

    return median_fit, fit84, fit16

def BIC_value(model, observed, sigma, bins, median_fit):
    if model == 'Schechter':
        func = models.Schechter_Mags(sp=median_fit)
    elif model == 'DPL':
        func = models.DPL_Mags(sp=median_fit)

    parent_volume = (3200)**3

    expected = func.N(parent_volume, bins)
    s = np.logical_and(observed>0, expected>0)
    binwidth = bins[1:] - bins[:-1]
    sig = sigma[s]/(parent_volume*binwidth[s])
    lnlike = -0.5 * np.sum((observed[s]/(parent_volume*binwidth[s]) - expected[s]/(parent_volume*binwidth[s])) ** 2 / sig**2 + np.log(sig**2))
    k = len(median_fit)

    BIC = k*np.log(np.sum(s)) - 2*lnlike
    AIC = 2*k - 2*lnlike
    print (F"{model} BIC: ", BIC)
    print (F"{model} AIC: ", AIC)


def bluetides_gsmf(z):

    if z == 8:
        mbin = np.arange(8.1, 10.3, 0.2)
        phi = np.array([-2.76, -2.97, -3.19, -3.42, -3.66, -3.92, -4.2, -4.5, -4.87, -5.17, -5.59, -5.96])

    elif z == 9:
        mbin = np.arange(8.1, 9.9, 0.2)
        phi = np.array([-3.22, -3.45, -3.71, -3.98, -4.27, -4.58, -4.93, -5.28, -5.65, -6.07])

    elif z == 10:
        mbin = np.arange(8.1, 9.5, 0.2)
        phi = np.array([-3.72, -4.01, -4.3, -4.61, -4.96, -5.34, -5.65, -6.11])

    else: raise ValueError("z not in list")

    return mbin, phi

def bluetides_uvlf(z):

    if z == 8:
        Mbin = -np.arange(17.25, 23, 0.5)
        phi = -np.array([2.44, 2.64, 2.82, 3.03, 3.25, 3.5, 3.85, 4.27, 4.8, 5.35, 6.07, 6.89])

    elif z == 9:
        Mbin = -np.arange(17.25, 22.75, 0.5)
        phi = -np.array([2.77, 2.97, 3.18, 3.41, 3.65, 3.93, 4.29, 4.76, 5.37, 6.06, 6.61])

    elif z == 10:
        Mbin = -np.arange(17.25, 22.25, 0.5)
        phi = -np.array([3.02, 3.25, 3.49, 3.76, 4.04, 4.35, 4.73, 5.28, 5.91, 6.59])

    else: raise ValueError("z not in list")

    return Mbin, phi


def plot_tng(M, z, axs):

    redshifts = np.array([5, 6, 7, 8])
    logphi = np.array([[-3.244, -3.107, -3.398], [-3.079, -3.025, -3.608], [-3.846, -3.418, -4.209], [-4.445, -4.111, -4.714]])
    Mstar = np.array([[-21.17, -20.95, -21.21], [-20.61, -20.52, -21.31], [-21.18, -20.58, -21.47], [-21.38, -20.86, -21.44]])
    alpha = np.array([[-1.924, -1.884, -1.941], [-1.876, -1.833, -2.042], [-2.133, -1.967, -2.279], [-2.280, -2.216, -2.455]])

    labels = ['TNG19:Model-A', 'TNG19:Model-B', 'TNG19:Model-C']
    colors = ['brown', 'magenta', 'pink']
    if z in redshifts:

        ok = np.where(redshifts == z)[0][0]

        theta = [logphi[ok][2], alpha[ok][2], Mstar[ok][2]]
        axs.plot(M, Schechter(M, theta), label=r'\textsc{Illustris} \textsc{Tng}', color='red', lw=2, alpha=0.4)
