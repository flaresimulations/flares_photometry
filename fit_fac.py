import numpy as np
import pandas as pd
from scipy import interpolate

def chisquare(obs, exp, err):

    return ((obs - exp)**2)/(err**2)

def compare(Mobs, phiobs, phiobserr, uplims, xdata, ydata):

    ok = np.where(ydata > 0)

    f = interpolate.interp1d(xdata[ok], ydata[ok], fill_value="extrapolate")

    yy = f(Mobs)

    w = np.ones(len(uplims))
    w[np.where(uplims==1)] = 1e-4

    chi = chisquare(phiobs, yy, phiobserr)

    chi = np.sum(chi*w)/(len(chi)-1)

    return chi




h = 0.6777
z = 5



data = np.genfromtxt('Obs_data/uv_lum_Fink15_z{}.txt'.format(z), delimiter=',', skip_header=1)
M, phi, phi_up, phi_low, uplims = data[:,0], data[:,1]*1e-3, data[:,2]*1e-3, data[:,3]*1e-3, data[:,4]
phi_up = np.log10(phi + phi_up) - np.log10(phi)
phi_low = np.log10(phi) - np.log10(phi - phi_low)
phi_err = (phi_low+phi_up)/2


data = np.load(F'UVLF_z{z}.npz')

bincen = data['UVMag']
hists = data['UVLF']
facs = data['fac']

binwidth = bincen[1] - bincen[0]
vol = (4/3)*np.pi*(14/h)**3

phi_resim = hists


rchi = np.array([])

for ii in range(len(facs)):

    rchi = np.append(rchi, compare(M, phi, phi_err, uplims, bincen, phi_resim[ii]))
