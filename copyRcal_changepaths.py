import numpy, h5py, os
from PnSC_main import *
from PnSC_h5io import *

h5f=h5py.File(mm.h5path, mode='r')
exppaths=experimentgrppaths(h5f)

for expp in exppaths:
    exp=expp.rpartition('/')[2]
    if not exp.startswith('Rcal'):
        continue
    rcalpath='/'.join((expp, 'analysis', 'Res_TempCal'))
    if not rcalpath in h5f:
        continue
    print expp
    ds=h5f[rcalpath]
    for count, (l, v) in enumerate(zip(plotinfo, ds[:, 0])):
        if v>0:
            l+=[[(v, rcalpath.strip('/').partition('/')[2].partition('/')[0])]]

            if 'Rovaluesaveragedwith' in ds.attrs.keys():
                p1=ds.attrs['Ro']
                p2=ds.attrs['Rovaluesaveragedwith']
                l[-1]+=[(h5f[p1][count], p1.strip('/').partition('/')[2].partition('/')[0])]
                l[-1]+=[(h5f[p2][count], p2.strip('/').partition('/')[2].partition('/')[0])]

