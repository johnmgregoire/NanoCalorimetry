import numpy, h5py, os
from PnSC_main import *
from PnSC_h5io import *

savef='E:/CHESS2010PnSC/cycleplots'
h5f=h5py.File(mm.h5path, mode='r')
exppaths=experimentgrppaths(h5f)
for expp in exppaths:
    exp=expp.rpartition('/')[2]
    hppaths=experimenthppaths(h5f, exp)
    for hpp in hppaths:
        data=CreateHeatProgSegDictList(mm.h5path, exp, hpp.rpartition('/')[2])
        i=SegmentCyclePlot(None, data)
        p1=os.path.join(savef,exp)
        if not os.path.exists(p1):
            os.mkdir(p1)
        i.fig.savefig(os.path.join(p1, hpp.rpartition('/')[2])+'.png')
h5f.close()
print 'plotting done'
