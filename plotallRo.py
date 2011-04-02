import numpy, h5py, os
from PnSC_main import *
from PnSC_h5io import *

savef='E:/CHESS2010PnSC/Roplots'
h5f=h5py.File(mm.h5path, mode='r')
exppaths=experimentgrppaths(h5f)
plotinfo=[[] for i in range(25)]
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

#E:\CHESS2010PnSC\AuSiCu_pnsc_all.h5
h5f.close()
offset=0.
for count, celllist in enumerate(plotinfo):
    print 'plotting', count
    if len(celllist)==0:
        continue
    xlabs=[]
    for count2, heatlist in enumerate(celllist):
        temp=heatlist.pop(0)
        xlabs+=[temp[1]]
        pylab.plot(count2, temp[0], 'ro')
        if len(heatlist)>0:
            i=numpy.argmin([heatlist[0][0], heatlist[1][0]])
            j=-1*i+1
            if heatlist[i][0]>0:
                pylab.plot(count2, heatlist[i][0], 'b^')
                pylab.text(count2, heatlist[i][0]-offset, heatlist[i][1],rotation='vertical',va='top',ha='center',alpha=.7,color='b')
            
            pylab.plot(count2, heatlist[j][0], 'gv')
            pylab.text(count2, heatlist[j][0]+offset, heatlist[j][1],rotation='vertical',va='bottom',ha='center',alpha=.7,color='g')
    
    pylab.xticks(range(len(xlabs)))
    pylab.gca().set_xticklabels(xlabs)
    pylab.ylabel('Ro (Ohms)')
    pylab.title('cell %d' %(count+1))
    pylab.xlim(-.5, len(xlabs)-.5)
    pylab.savefig(os.path.join(savef,('Ro_cell%02d' %(count+1))+'.png'))
    pylab.clf()

print 'plotting done'
