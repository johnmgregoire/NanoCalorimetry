import time, copy
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
import scipy.integrate

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'

def myexpformat(x, pos):
    for ndigs in range(2):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

def make_ticklabels_invisible(ax, x=True, y=True):
    if x:
        for tl in ax.get_xticklabels():
            tl.set_visible(False)
    if y:
        for tl in ax.get_yticklabels():
            tl.set_visible(False)
    
    
    
nskip=100

cycleindex=0

savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'

cols=['k', 'b', 'g', 'r', 'c', 'm', 'y']*3

expdicts=[\
dict([('name', 'heat1a'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat1b'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat1c'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat1d'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat2'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat3'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat4a'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat4b'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat4c'), ('heatseg', 2), ('coolseg', 3)]), \
dict([('name', 'heat6a'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat6b'), ('heatseg', 2), ('coolseg', 4)]), \
dict([('name', 'heat6c'), ('heatseg', 2), ('coolseg', 3)]), \
dict([('name', 'heat7'), ('heatseg', 2), ('coolseg', 3)]), \
dict([('name', 'heat8'), ('heatseg', 2), ('coolseg', 4)]), \
]


#axl=get25pylabaxes(horizslowaxis=True)
#taxl=[ax.twinx() for ax in axl]

tcrcells=numpy.array([1, 9, 11, 16, 17, 24])
f, g=getcalanalysis(p, 'R_TCR_RT')
tcr=g['Res_TempCal'][:, :][:, 2]
f.close()

fom=numpy.zeros((25, len(expdicts)), dtype='float32')
metadictlist=[]
allsegdict=[]
for count, ed in enumerate(expdicts):
    exp=ed['name']
    f, hppaths=experimenthppaths(p, exp)
    f.close()
    hpsdl=None
    saveh5hpname=None
    for hpp in hppaths:
        h5hpname=hpp.rpartition('/')[2]
        f, g=gethpgroup(p, exp, h5hpname=h5hpname)
        cell=g.attrs['CELLNUMBER']
        f.close()
        hpsdl=CreateHeatProgSegDictList(p, exp, h5hpname)
        saveh5hpname=h5hpname
        if not hpsdl is None:
            ed['cell']=cell
            T=hpsdl[ed['heatseg']]['sampletemperature'][cycleindex][nskip:-1*nskip]
            maxT=T.max()
            T=numpy.concatenate([T, hpsdl[ed['coolseg']]['sampletemperature'][cycleindex][nskip:-1*nskip]])
            i=0
            dt=[]
            while len(dt)<2:
                try:
                    dt=hpsdl[i]['cycletime'][cycleindex][0:2]
                    i+=1
                except:
                    dt=[]
                    break
            if len(dt)<0:
                continue
            dt=dt[1]-dt[0]
            for critt in [500, 600, 700, 900, 1100, 1300, 1420, 1500]:
                timeabove=len(numpy.where(T>(1.*critt))[0])*dt
                ed['secabove%dC' %critt]=timeabove
                
            ed['maxT']=maxT
            
            ed['edcount']=count
            metadictlist+=[copy.deepcopy(ed)]
            if count==0:
                fom[cell-1, count]=1.
            if count+1<fom.shape[1]:
                fom[cell-1, count+1]=max(maxT, numpy.max(fom[cell-1, :count+1]))
            #fom[cell-1]+=timeabove
    
if 1:
    f=h5py.File(p, mode='r+')
    for ed in metadictlist:
        if not 'maxT' in ed.keys():
            continue
        exp=ed['name']
        try:
            g=f['calbycellmetadata'][`ed['cell']`][exp]
        except:
            continue
        for k, v in ed.iteritems():
            if k=='maxT' or k.startswith('secabove'):
                print ed['cell'],  k, v
                g.attrs[k]=v
    f.close()
    
tcrcalc=lambda Tm:(Tm>1410. and (6.7e-7*(Tm-1410.), ) or (0, ))[0]+9.47e-4
tcrall=numpy.array([[tcrcalc(t) for t in arr] for arr in fom])

#for count, (ax, tcrl) in enumerate(zip(axl, tcrall)):
#    ax.plot(numpy.where(tcrl>0)[0], tcrl[tcrl>0], '.')
#    pylab.ylabel('cell%d' %(count+1))
#    pylab.xlabel('heat program')

savebool=False
if savebool:
    f=h5py.File(p, mode='r+')
    for count, (ed, tcr_cells) in enumerate(zip(expdicts, tcrall.T)):
        exp=ed['name']
        try:
            rcp=f['Calorimetry'][exp].attrs['Res_TempCalPath']
            f[rcp][:, 2]=tcr_cells[:]
        except:
            print 'saving skipped for ', exp
    f.close()
    
    
#pylab.show()
print 'done'
