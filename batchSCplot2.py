import time
import os
import sys
import numpy
import h5py
#from PnSC_ui import *
#from PnSC_dataimport import *
from PnSC_SCui import *
#from PnSC_math import *
from PnSC_h5io import *
from PnSC_main import *
from matplotlib.ticker import FuncFormatter

def myexpformat(x, pos):
    for ndigs in range(2):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

rootfolder='E:/CHESS2010PnSC'

nskip=100
#name, (cycles,daqHz)
exp_rec=[\
#('heat1a',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat1b',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat1d',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat1c',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat1_cell21_other',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat2', [(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -180*nskip))]),\
#('heat2_cell21_other',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat3',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -80*nskip))]),\
#('heat3_cell21_other',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -80*nskip))]),\
#('heat4a',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -100*nskip))]),\
#('heat4b',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -100*nskip))]),\
#('heat4_cell21_other',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
('heat4c',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (3, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat6a',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat6b',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat6c',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (3, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat6_cell21_other',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat7',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (3, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip))]),\
#('heat8',[(2, ('cycletime', 'sampletemperature', 'sampletemperature'),('samplecurrent', 'sampleheatrate', 'samplepowerperrate'), (nskip, -1*nskip)), (4, ('cycletime', 'sampletemperature', 'sampletemperature'),('sampletemperature', 'sampleheatrate', 'samplepowerperrate'), (nskip, -25*nskip))]),\
]
p='E:/CHESS2010PnSC/AuSiCu_pnsc_all.h5'
#p=mm.h5path
#f=h5py.File(p, mode='r+')
#f=h5py.File(p, mode='r')
savef='E:/CHESS2010PnSC/SCplots'

forcecopyrecipes=False
for exp, plotlist in exp_rec:
    m=len(plotlist)
    f, hppaths=experimenthppaths(p, exp)
    f.close()
    for hpp in hppaths:
        try:
        #if True:
            pylab.figure(figsize=(10, 8))
            h5hpname=hpp.rpartition('/')[2]
            f, g=gethpgroup(p, exp, h5hpname=h5hpname)
            cell=g.attrs['CELLNUMBER']
            f.close()
            hpsdl=CreateHeatProgSegDictList(p, exp, h5hpname)
            titles=['cell %02d' %cell, exp]
            segcols=['b', 'g', 'r', 'c', 'm']
            for count, (segind, xkeys, ykeys, indswithinseg) in enumerate(plotlist):
                n=len(xkeys)
                for count2, (xk, yk) in enumerate(zip(xkeys, ykeys)):
                    if count2==0:
                        pylab.subplot(n, m, count2*m+count+1)
                        for segcount, d in enumerate(hpsdl):
                            c='k'
                            for tempcount, (tempsegind, tempxkeys, tempykeys, tempindswithinseg) in enumerate(plotlist):
                                if segcount==tempsegind:
                                    c=segcols[tempcount]
                            if xk in d.keys() and yk in d.keys():
                                for x, y in zip(d[xk], d[yk]):
                                    pylab.plot(x, y, c+'.', markersize=1)
                        pylab.xlabel(xk)
                        pylab.ylabel(yk)
                        if count<len(titles):
                            pylab.title(titles[count])
                    else:
                        xcyc=hpsdl[segind][xk][:, indswithinseg[0]:indswithinseg[1]]
                        ycyc=hpsdl[segind][yk][:, indswithinseg[0]:indswithinseg[1]]
                        pylab.subplot(n, m, count2*m+count+1)
                        for x, y in zip(xcyc, ycyc):
                            pylab.plot(x, y, segcols[count]+'.', markersize=1)
                        pylab.xlabel(xk)
                        pylab.ylabel(yk+',seg.ind.%d' %segind)
                        yc=ycyc.mean(axis=0)
                        t1=numpy.median(yc)
                        t2=.01*yc.std()
                        fincl=((yc>t1-t2)&(yc<t1+t2)).sum()*1./len(yc)
                        delf=0.
                        while fincl<.7 or delf>.005:
                            t2*=1.2
                            temp=((yc>t1-t2)&(yc<t1+t2)).sum()*1./len(yc)
                            delf=temp-fincl
                            fincl=temp
                        if (t1-2*t2)<yc.min():
                            t3=yc.min()
                        else:
                            t3=t1-t2
                        if (t1+2*t2)>yc.max():
                            t4=yc.max()
                        else:
                            t4=t1+t2
                        pylab.ylim(t3, t4)
                        #pylab.ylim(yc[yc>(t1-t2)].min(), yc[yc<(t1+t2)].max())
                        if count2>0:
                            pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
            p1=os.path.join(savef,exp)
            if not os.path.exists(p1):
                os.mkdir(p1)
            pylab.subplots_adjust(left=.1, right=.97, top=.93, bottom=.08, wspace=.25, hspace=.25)
#            pylab.show()
#            break
            pylab.savefig(os.path.join(p1,'SCplot_%s_cell%02d' %(exp, cell)+'.png'))
            pylab.clf()
        except:
            print 'error in ', exp, h5hpname
