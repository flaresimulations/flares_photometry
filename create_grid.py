

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

from mpi4py import MPI
from modules import get_lum_all

h=0.6777
sims = np.arange(0,38)

z = 5

tags = ['010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000']
tags_ref = ['008_z005p037', '006_z005p971', '005_z007p050', '004_z008p075']
vol = (4/3)*np.pi*(14/h)**3

bins = np.arange(-24, -16, 0.1)
bincen = (bins[1:]+bins[:-1])/2.
binwidth=bins[1:]-bins[:-1]
facs = np.arange(0.025, 0.07, 0.0005)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    hists = np.zeros((len(facs),len(bincen)))
    hist = np.zeros((len(facs),len(bincen)))

else:
    hists = np.empty((len(facs),len(bincen)))
    hist = np.zeros((len(facs),len(bincen)))

part = int(len(facs)/size)

if rank!=size-1:
    thisok = np.arange(rank*part, (rank+1)*part, 1).astype(int)

else:
    thisok = np.arange(rank*part, len(facs), 1).astype(int)

for ii, jj in enumerate(thisok):

    data = get_lum_all(facs[jj], tags[0], bins = bins, LF = True)
    hist[jj] = data/(binwidth*vol)

comm.Reduce(hist, hists, op=MPI.SUM, root=0)

if rank == 0:
    np.savez('UVLF_z{}.npz'.format(z), UVLF = hists, UVMag = bincen, fac = facs)
