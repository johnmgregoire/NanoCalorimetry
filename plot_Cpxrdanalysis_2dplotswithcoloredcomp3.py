import os
os.chdir('C:/Users/JohnnyG/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
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
from scipy.interpolate import griddata
import matplotlib.cm as cm

selectcell=25
p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'

#CHOICE='amfrac'
CHOICE='Hgperamfrac'

ylabdict={\
'amfrac':'amorphous fraction after cooling', \
'amfracinB':'amorph / amorph+phaseB', \
'Hgperamfrac':r'volumetric glass transition enthalpy (J/cm$^3$)', \
}

mdboolfcndict={\
'amfrac':lambda md:'prevcoolrate_320C' in md.keys(), \
'amfracinB':lambda md:'prevcoolrate_320C' in md.keys(), \
'Hgperamfrac':lambda md:'prevcoolrate_320C' in md.keys() and 'Cpregions_enthalpy' in md.keys() and numpy.abs(md['prevcoolrate_320C'])>1.e3, \
}

xdboolfcndict={\
'amfrac':True, \
'amfracinB':True, \
'Hgperamfrac':lambda xd:'amfrac' in xd.keys() and xd['amfrac']>0.02, \
}

fomfcndict={\
'amfrac':lambda x, y, z:y, \
'amfracinB':lambda x, y, z:y/(y+z), \
'Hgperamfrac':lambda x, y, z:y/z*1.e-6/(3.6*.8*1.e-15), \
}


xfcndict={\
'amfrac':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
'amfracinB':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
'Hgperamfrac':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
}

yfcndict={\
'amfrac':lambda md, xd:xd['amfrac'], \
'amfracinB':lambda md, xd:xd['amfrac'], \
'Hgperamfrac':lambda md, xd:(md['Cpregions_glassind']<0 and (0., ) or (md['Cpregions_enthalpy'][md['Cpregions_glassind']], ))[0], \
}

zfcndict={\
'amfrac':lambda md, xd:xd['othfrac'], \
'amfracinB':lambda md, xd:xd['othfrac'], \
'Hgperamfrac':lambda md, xd:xd['amfrac']*md['nm'], \
}


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


cycleindex=0

#p=mm.h5path
#f=h5py.File(p, mode='r+')
#f=h5py.File(p, mode='r')
savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'


plotTlim=(50., 700.)

f=h5py.File(p, mode='r')
comp=f['CompThick/atfrac'][:, :]
nm=f['CompThick/nm'][:]

metadictlists=[]
for selectcell in range(1, 26):
    metadictlist=[]
    #allsegdict=[]
    if `selectcell` in f['calbycellmetadata']:
        cg=f['calbycellmetadata'][`selectcell`]
        for mg in cg.itervalues():
            if isinstance(mg, h5py.Group) and 'Cpregions_enthalpy' in mg.attrs.keys():
                d={}
                for k, v in mg.attrs.iteritems():
                    d[k]=v
        #        if selectcell==1 and d['name'].startswith('heat1'):#heat1a was botched and heat1b we don't know cooling rate and the XRd for heat0 was questionable anyway
        #            continue
                metadictlist+=[d]
                #allsegdict+=[CreateHeatProgSegDictList(p, d['name'], d['h5hpname'])]
    metadictlists+=[metadictlist]

xrddictlists=[]
for selectcell in range(1, 26):
    xrddictlist=[]
    if 'xrdbycell' in f and `selectcell` in f['xrdbycell']:
        cg=f['xrdbycell'][`selectcell`]
        for mg in cg.itervalues():
            if isinstance(mg, h5py.Group):
                d={}
                for k, v in mg.attrs.iteritems():
                    d[k]=v
                xrddictlist+=[d]
    xrddictlists+=[xrddictlist]
f.close()


x_cell=[]
y_cell=[]
z_cell=[]
c_cell=[]
cell_cell=[]
x=[]
y=[]
z=[]
d=[]
c=[]
cells=[]
for cell, cv, nmv, mdl, xdl in zip(range(1, 26), comp, nm, metadictlists, xrddictlists):
    xtemp=[]
    ytemp=[]
    ztemp=[]
    ctemp=[]
    celltemp=-1
    for md in mdl:
        if not mdboolfcndict[CHOICE](md):
            continue
        md['nm']=nmv
        md['comp']=cv
        pcal=md['prevname'][:5]
        cal=md['name'][:5]
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
            for xd in xdl:
                if not xdboolfcndict[CHOICE](xd):
                    continue
                if xd['name']==pcal:#use xrd that happened after the prev scan
                    
                    x+=[xfcndict[CHOICE](md, xd)]
                    y+=[yfcndict[CHOICE](md, xd)]
                    z+=[zfcndict[CHOICE](md, xd)]
                    c+=[cv]
                    cells+=[cell]
                    
                    
                    xtemp+=[xfcndict[CHOICE](md, xd)]
                    ytemp+=[yfcndict[CHOICE](md, xd)]
                    ztemp+=[zfcndict[CHOICE](md, xd)]
                    ctemp+=[cv]
                    celltemp=cell
                    break
    if len(ctemp)>0:
        x_cell+=[xtemp]
        y_cell+=[ytemp]
        z_cell+=[ztemp]
        c_cell+=[ctemp]
        cell_cell+=[celltemp]

c=numpy.array(c)
y=numpy.array(y)
x=numpy.array(x)
z=numpy.array(z)


pylab.figure(figsize=(10, 6))
ax=pylab.subplot(111)
stp = TernaryPlot(ax, ellabels=['Au', 'Si', 'Cu']) 
celllist=[1, 2]+range(4, 21)+[22]+[24, 25]
minlist=[c.min() for c in comp[numpy.array(celllist)-1, :].T]
rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
colors=stp.color_comp_calc(comp[numpy.array(celllist)-1, :], rangelist=rangelist)

pylab.clf()
for cell, xv, yv, zv in zip(cell_cell, x_cell, y_cell, z_cell):
    xv=numpy.array(xv)
    yv=numpy.array(yv)
    zv=numpy.array(zv)
    sortarr=numpy.argsort(xv)
    #pylab.plot(xv[sortarr], yv[sortarr], '.', color=colors[celllist.index(cell)], markersize=16)
    #pylab.plot(xv[sortarr], yv[sortarr], '-', color=colors[celllist.index(cell)], linewidth=1, alpha=.5)
    
    a=fomfcndict[CHOICE](xv[sortarr], yv[sortarr], zv[sortarr])
    if 0:
        pylab.plot(xv[sortarr], a, '.', color=colors[celllist.index(cell)], markersize=16)
    else:
        for xx, aa in zip(xv[sortarr], a):
            pylab.text(xx, aa, `cell`, ha='center', va='center', color=colors[celllist.index(cell)], fontsize=12)
    pylab.plot(xv[sortarr], a, '-', color=colors[celllist.index(cell)], linewidth=1, alpha=.5)
    
    print cell, xv[sortarr], a

pylab.gca().set_xscale('log')

xmin=x.min()-.05*(x.max()-x.min())
xmax=x.max()+.05*(x.max()-x.min())

ylabel=ylabdict[CHOICE]
a=fomfcndict[CHOICE](x, y, z)

ymin=a.min()-.05*(a.max()-a.min())
ymax=a.max()+.05*(a.max()-a.min())

pylab.ylabel(ylabel, fontsize=14)
pylab.xlabel('cooling rate at 320C (K/s)', fontsize=14)
pylab.title('colored by composition', fontsize=14)
#pylab.ylabel(, fontsize=14)
pylab.xlim(xmin, xmax)
pylab.ylim(ymin, ymax)


pylab.show()
