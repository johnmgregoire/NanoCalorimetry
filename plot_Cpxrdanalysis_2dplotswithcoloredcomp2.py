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


#following 8 lines are common to several belwo subroutines
z=[[] for cv in comp]
y=[[] for cv in comp]
x=[[] for cv in comp]
for i, (xv, yv, zv, cv, mdl, xdl) in enumerate(zip(x, y, z, comp, metadictlists, xrddictlists)):
    for md in mdl:
        pcal=md['prevname'][:5]
        cal=md['name'][:5]
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
        
        
        
##prev cool rate vs amfrac
#            if 'prevcoolrate_400C' in md.keys():
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[0.]
#                        yv+=[xd['amfrac']]
#                        xv+=[numpy.abs(md['prevcoolrate_400C'])]
#                        print 'cell ', i+1
#                        break
#data=[(cell, xv, yv, zv, cv) for cell, xva, yva, zva, cv in zip(range(1,26),x, y, z, comp) for xv, yv, zv in zip(xva, yva, zva)]
#cells=numpy.array(map(operator.itemgetter(0), data))
#xplot=numpy.array(map(operator.itemgetter(1), data))
#yplot=numpy.array(map(operator.itemgetter(2), data))
#garb=numpy.array(map(operator.itemgetter(3), data))
#ca=numpy.array(map(operator.itemgetter(4), data))
#xlab='cooling rate (K/s)'
#ylab='amorphous fraction after cooling'
##****


##melt enthalpy - xtal enthalpy / fcc, including when xtal enthalpy is zero
#            if 'Cpregions_xtalind' in md.keys() and 'Cpregions_meltind' in md.keys() and md['Cpregions_meltind']>=0:
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[md['Cpregions_enthalpy'][md['Cpregions_meltind']]]
#                        yv+=[xd['fccfrac']]
#                        xv+=[(md['Cpregions_xtalind']<0 and (0, ) or (md['Cpregions_enthalpy'][md['Cpregions_xtalind']], ))[0]]
#                        print 'cell ', i+1
#                        break
#data=[(cell, xv, yv, zv, (zv-xv), cv) for cell, xva, yva, zva, cv in zip(range(1,26),x, y, z, comp) for xv, yv, zv in zip(xva, yva, zva)]
#cells=numpy.array(map(operator.itemgetter(0), data))
#garb=map(operator.itemgetter(1), data)
#yplot=numpy.float32(map(operator.itemgetter(2), data))
#garb=map(operator.itemgetter(3), data)
#xplot=numpy.float32(map(operator.itemgetter(4), data))
#ca=numpy.array(map(operator.itemgetter(5), data))
#xlab='melt enthalpy - xtal enthalpy'
#ylab='fcc frac'
##****

#melt enthalpy - xtal enthalpy / fcc, including when xtal enthalpy is zero
            if 'Cpregions_xtalind' in md.keys() and 'Cpregions_meltind' in md.keys() and md['Cpregions_meltind']>=0:
                for xd in xdl:
                    if xd['name']==pcal:#use xrd that happened after the prev scan
                        zv+=[md['Cpregions_enthalpy'][md['Cpregions_meltind']]]
                        yv+=[xd['othfrac']]
                        xv+=[(md['Cpregions_xtalind']<0 and (0, ) or (md['Cpregions_enthalpy'][md['Cpregions_xtalind']], ))[0]]
                        print 'cell ', i+1
                        break
data=[(cell, xv, yv, zv, (zv-xv), cv) for cell, xva, yva, zva, cv in zip(range(1,26),x, y, z, comp) for xv, yv, zv in zip(xva, yva, zva)]
cells=numpy.array(map(operator.itemgetter(0), data))
garb=map(operator.itemgetter(1), data)
yplot=numpy.float32(map(operator.itemgetter(2), data))
garb=map(operator.itemgetter(3), data)
xplot=numpy.float32(map(operator.itemgetter(4), data))
ca=numpy.array(map(operator.itemgetter(5), data))
xlab='melt enthalpy - xtal enthalpy'
ylab='oth frac'
#****


##melt enthalpy - xtal enthalpy / fcc, including when xtal enthalpy is zero
#            if 'heatrate_170C500C' in md.keys() and 'Cpregions_glassind' in md.keys() and md['Cpregions_glassind']>=0:
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[md['Cpregions_enthalpy'][md['Cpregions_glassind']]]
#                        yv+=[md['Cpregions_Tweightedmean'][md['Cpregions_glassind']]]
#                        xv+=[md['heatrate_170C500C']]
#                        print 'cell ', i+1
#                        break
#data=[(cell, xv, yv, zv, 0., cv) for cell, xva, yva, zva, cv in zip(range(1,26),x, y, z, comp) for xv, yv, zv in zip(xva, yva, zva)]
#cells=numpy.array(map(operator.itemgetter(0), data))
#xplot=numpy.array(map(operator.itemgetter(1), data))
#yplot=numpy.array(map(operator.itemgetter(2), data))
#garb=numpy.array(map(operator.itemgetter(3), data))
#garb=numpy.array(map(operator.itemgetter(4), data))
#ca=numpy.array(map(operator.itemgetter(5), data))
#xlab='heat rate (K/s)'
#ylab='Tg weighted mean (K)'
##****

#start common plotting routine
pylab.figure(figsize=(10, 6))
ax=pylab.subplot(111)
stp = TernaryPlot(ax, ellabels=['Au', 'Si', 'Cu']) 
minlist=[c.min() for c in ca.T]
rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
colors=stp.color_comp_calc(ca, rangelist=rangelist)

pylab.clf()
for cell in set(cells):
    i=numpy.where(cells==cell)
    x_cell=numpy.array(xplot[i])
    y_cell=numpy.array(yplot[i])
    sortarr=numpy.argsort(x_cell)
    x_cell=x_cell[sortarr]
    y_cell=y_cell[sortarr]
    col=colors[i][0]#should all be the same color so take the 1st
    pylab.plot(x_cell, y_cell, '.', color=col, markersize=16)
#    for xx, yy in zip(x_cell, y_cell):
#        pylab.text(xx, yy, `cell`, ha='center', va='center', color=col, fontsize=12)

    pylab.plot(x_cell, y_cell, '-', color=col, linewidth=1, alpha=.5)

#pylab.gca().set_xscale('log')

xmin=xplot.min()-.05*(xplot.max()-xplot.min())
xmax=xplot.max()+.05*(xplot.max()-xplot.min())
ymin=yplot.min()-.05*(yplot.max()-yplot.min())
ymax=yplot.max()+.05*(yplot.max()-yplot.min())


pylab.ylabel(ylab, fontsize=14)
pylab.xlabel(xlab, fontsize=14)
pylab.title('colored by composition', fontsize=14)
#pylab.ylabel(, fontsize=14)
pylab.xlim(xmin, xmax)
pylab.ylim(ymin, ymax)


#pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'HeatRatestack_cell%02d_%s.png' %(selectcell, metadict['name'])))

pylab.show()
