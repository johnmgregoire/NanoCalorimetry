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
c_cell=[]
cell_cell=[]
x=[]
y=[]
z=[]
d=[]
c=[]
cells=[]
for cell, cv, mdl, xdl in zip(range(1, 26), comp, metadictlists, xrddictlists):
    xtemp=[]
    ytemp=[]
    ctemp=[]
    celltemp=-1
    for md in mdl:
        if not 'prevcoolrate_400C' in md.keys():
            continue
        pcal=md['prevname'][:5]
        cal=md['name'][:5]
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
            for xd in xdl:
                if xd['name']==pcal:#use xrd that happened after the prev scan
                    
                    x+=[numpy.abs(md['prevcoolrate_400C'])]
                    y+=[xd['amfrac']]
                    c+=[cv]
                    cells+=[cell]
                    
                    
                    xtemp+=[numpy.abs(md['prevcoolrate_400C'])]
                    ytemp+=[xd['amfrac']]
                    ctemp+=[cv]
                    celltemp=cell
                    break
    if len(ctemp)>0:
        x_cell+=[xtemp]
        y_cell+=[ytemp]
        c_cell+=[ctemp]
        cell_cell+=[celltemp]

c=numpy.array(c)
y=numpy.array(y)
x=numpy.array(x)

pylab.figure(figsize=(10, 6))
ax=pylab.subplot(111)
stp = TernaryPlot(ax, ellabels=['Au', 'Si', 'Cu']) 
celllist=[1, 2]+range(4, 21)+[22]+[24, 25]
minlist=[c.min() for c in comp[numpy.array(celllist)-1, :].T]
rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
colors=stp.color_comp_calc(comp[numpy.array(celllist)-1, :], rangelist=rangelist)

pylab.clf()
for cell, xv, yv in zip(cell_cell, x_cell, y_cell):
    xv=numpy.array(xv)
    yv=numpy.array(yv)
    sortarr=numpy.argsort(xv)
    #pylab.plot(xv[sortarr], yv[sortarr], '.', color=colors[celllist.index(cell)], markersize=16)
    for xx, yy in zip(xv, yv):
        pylab.text(xx, yy, `cell`, ha='center', va='center', color=colors[celllist.index(cell)], fontsize=12)
        
    pylab.plot(xv[sortarr], yv[sortarr], '-', color=colors[celllist.index(cell)], linewidth=1, alpha=.5)

pylab.gca().set_xscale('log')

xmin=x.min()-.05*(x.max()-x.min())
xmax=x.max()+.05*(x.max()-x.min())
ymin=y.min()-.05*(y.max()-y.min())
ymax=y.max()+.05*(y.max()-y.min())

pylab.ylabel('amorphous fraction after cooling')
pylab.xlabel('cooling rate (K/s)', fontsize=14)
pylab.title('colored by composition', fontsize=14)
#pylab.ylabel(, fontsize=14)
pylab.xlim(xmin, xmax)
pylab.ylim(ymin, ymax)


#pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'HeatRatestack_cell%02d_%s.png' %(selectcell, metadict['name'])))

pylab.show()
