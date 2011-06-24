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


nskip=100
#name, (cycles,daqHz)
exp_rec=[\
('heat1a',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat1b',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat1d',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat1c',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat1_cell21_other',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat2', [(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -180*nskip))]),\
('heat2_cell21_other',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat3',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -80*nskip))]),\
('heat3_cell21_other',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -80*nskip))]),\
('heat4a',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -100*nskip))]),\
('heat4b',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -100*nskip))]),\
('heat4_cell21_other',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat4c',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (3, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat6a',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat6b',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat6c',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (3, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat6_cell21_other',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat7',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (3, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip))]),\
('heat8',[(2, 'sampletemperature', 'samplepowerperrate', (nskip, -1*nskip)), (4, 'sampletemperature', 'samplepowerperrate', (nskip, -25*nskip))]),\
]

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'
#p=mm.h5path
#f=h5py.File(p, mode='r+')
#f=h5py.File(p, mode='r')
savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_May17'
plotinfo0=[[] for i in range(25)]
plotinfo1=[[] for i in range(25)]
plotinfo=[plotinfo0, plotinfo1]
forcecopyrecipes=False
for exp, plotlist in exp_rec:
    m=len(plotlist)
    f, hppaths=experimenthppaths(p, exp)
    f.close()
    for hpp in hppaths:
        h5hpname=hpp.rpartition('/')[2]
        f, g=gethpgroup(p, exp, h5hpname=h5hpname)
        cell=g.attrs['CELLNUMBER']
        f.close()
        hpsdl=CreateHeatProgSegDictList(p, exp, h5hpname)
        titles=['cell %02d' %cell, exp]
        segcols=['b', 'g', 'r', 'c', 'm']
        for count, ((segind, xk, yk, indswithinseg), pi) in enumerate(zip(plotlist, plotinfo)):
            xc=hpsdl[segind][xk][0, indswithinseg[0]:indswithinseg[1]]
            yc=hpsdl[segind][yk][0, indswithinseg[0]:indswithinseg[1]]
            pi[cell-1]+=[(xc, yc, exp)]

hight=600
for cellcount, (pi0, pi1) in enumerate(zip(plotinfo0, plotinfo1)):
    print cellcount, len(pi0)
    if len(pi0)==0:
        continue

    titles=['heating cell%02d' %(cellcount+1), 'cooling cell%02d' %(cellcount+1)]
    pylab.figure(figsize=(10, 8))
    for count, pi in enumerate([pi0, pi1]):
        pylab.subplot(2, 1, count+1)
        yc=numpy.float32([])
        for x, y, l in pi:
            pylab.plot(x, y, '.', markersize=1, label=l)
            yc=numpy.append(yc, y[x<hight])
        pylab.legend(loc=2)    
        if count<len(titles):
            pylab.title(titles[count])
        pylab.xlim(0, hight)
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
        pylab.xlabel('T (C)')
        pylab.ylabel('P / dT/dt')
        pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
    pylab.subplots_adjust(left=.1, right=.97, top=.93, bottom=.08, wspace=.25, hspace=.25)
#    pylab.show()
#    idialog=messageDialog(title='continue')
#    if not idialog.exec_():
#        break
#    break
    pylab.savefig(os.path.join(savef,'SCcellplot_cell%02d' %(cellcount+1)+'.png'))
    pylab.clf()
