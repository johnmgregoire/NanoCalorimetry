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

selectcell=11
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

def paramave_T(d, T, Twin=10., hpsdkey='sampleheatrate'):
    #i=numpy.argmin((T-d['sampletemperature'])**2)
    Ta=d['sampletemperature'][cycleindex]
    x=numpy.where((Ta>=T-Twin)&(Ta<=T+Twin))[0]
    prev=numpy.array([not (t-1 in x) for t in x])
    previ=numpy.where(prev)[0]
    if len(previ)==0:
        return 0.
    stopi=numpy.append(previ[1:],len(x))
    longestbunchind=numpy.argmax(stopi-previ)
    inds=x[previ[longestbunchind]:stopi[longestbunchind]]
    return d[hpsdkey][cycleindex][inds].mean()
    
    

def findenthalpyandpinacles(segdict, critenth=1.e-5, dTmin=.4, Tmeanmin=100.):
    T=segdict['sampletemperature'][cycleindex]
    C=segdict['sampleheatcapacity'][cycleindex]
    nci=numpy.where((C[:-1]>0.)&(C[1:]<=0.))[0]#neg crossings
    pci=numpy.where((C[1:]>0.)&(C[:-1]<=0.))[0]#pos crossings
    ci=numpy.sort(numpy.concatenate([nci, pci]))
    ans=[]
    for i, j in zip(ci[:-1], ci[1:]):
        enth=scipy.integrate.trapz(C[i:j], T[i:j])
        if numpy.abs(enth)>critenth and (T[j]-T[i])>dTmin:
            itemp=numpy.argmax(numpy.abs(C[i:j]))
            Tmean=scipy.integrate.trapz(C[i:j]*T[i:j], T[i:j])/scipy.integrate.trapz(C[i:j], T[i:j])
            if Tmean<Tmeanmin:
                continue
            ans+=[dict([('enthalpy', enth), ('T_Cmax', T[i:j][itemp]), ('Cmax', C[i:j][itemp]), ('Tweightedmean', Tmean), ('cycindstart', i), ('cycindstop', j)])]
    return ans
    
    
    
    
nskip=100

cycleindex=0

#p=mm.h5path
#f=h5py.File(p, mode='r+')
#f=h5py.File(p, mode='r')
savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'


plotTlim=(50., 700.)


metadictlist=[]
allsegdict=[]
f=h5py.File(p, mode='r')
cg=f['calbycellmetadata'][`selectcell`]
for mg in cg.itervalues():
    if isinstance(mg, h5py.Group) and 'name' in mg.attrs.keys():
        d={}
        for k, v in mg.attrs.iteritems():
            d[k]=v
#        if selectcell==1 and d['name'].startswith('heat1'):#heat1a was botched and heat1b we don't know cooling rate and the XRd for heat0 was questionable anyway
#            continue
        metadictlist+=[d]
        allsegdict+=[CreateHeatProgSegDictList(p, d['name'], d['h5hpname'])]


#xrddictlist=[]
#if 'xrdbycell' in f and `selectcell` in f['xrdbycell']:
#    cg=f['xrdbycell'][`selectcell`]
#    for mg in cg.itervalues():
#        if isinstance(mg, h5py.Group):
#            d={}
#            for k, v in mg.attrs.iteritems():
#                d[k]=v
#            xrddictlist+=[d]
f.close()

orderarray=numpy.abs(numpy.array([metadict['heatrate_170C500C'] for metadict in metadictlist]))

sortinds=numpy.argsort(orderarray)
cols=['b', (160./256.,160./256.,0), 'r', 'g', 'c', 'm', 'k']


## plotting series of heat ramps
mult=1.e6
nplots=len(orderarray)
pylab.figure(figsize=(8, 8))
axl=[pylab.subplot(nplots, 1, nplots)]
for i in range(1, nplots):
    #ax=pylab.subplot2grid((n, 3), (n-1-i, 0), colspan=2, sharex=axl[0], sharey=axl[0])
    ax=pylab.subplot(nplots, 1, nplots-i, sharex=axl[0], sharey=axl[0])
    pylab.setp(ax.get_xticklabels(), visible=False)
    axl+=[ax]

f=h5py.File(p, mode='r')
x=[]
y=[]
x2=[]
y2=[]
for count, i in  enumerate(sortinds):
    hpsdl=allsegdict[i]
    metadict=metadictlist[i]
    pave=paramave_T(hpsdl[metadict['heatseg']], 400., Twin=100., hpsdkey='samplepower')
    dtave=paramave_T(hpsdl[metadict['heatseg']], 400., Twin=100., hpsdkey='sampleheatrate')
    cdtave=numpy.abs(paramave_T(hpsdl[metadict['coolseg']], 400., Twin=100., hpsdkey='sampleheatrate'))
    
    if dtave>0:
        y+=[dtave/pave]
        s=f['Calorimetry'][metadict['name']]['measurement/HeatProgram'][metadict['h5hpname']].attrs['ambient_atmosphere']
        press=(s.endswith('mT') and (eval(s.strip().partition('mT')[0]), ) or (0., ))[0]
        x+=[press]
        print metadict['name'], pave, dtave, press
        x2+=[press]
        y2+=[cdtave]
    if metadict['name']=='heat4a':
        break
f.close()

pylab.clf()
pylab.plot(x, y, 'ro')
pylab.xlabel('He pressure (mT)', fontsize=14)
pylab.ylabel('heat rate per power, dT/dt / P (K/s/W), at 300-500C', fontsize=14)
pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
pylab.twinx()
pylab.plot(x2, y2, 'bo')
pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
pylab.ylabel('cool rate (K/s), at 300-500C', fontsize=14)
pylab.show()
