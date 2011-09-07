import os
os.chdir('C:/Users/JohnnyG/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
import time, copy
import os
import sys
import numpy
import h5py
import pylab
#from PnSC_ui import *
#from PnSC_dataimport import *
#from PnSC_SCui import *
#from PnSC_math import *
#from PnSC_h5io import *
#from PnSC_main import *
from matplotlib.ticker import FuncFormatter
import scipy.integrate
from scipy.interpolate import griddata
import matplotlib.cm as cm

sp='C:/Users/JohnnyG/Documents/HarvardWork/MG/plots_Nov2010sample/plots_Aug2011/22cellarrays'
readarr=lambda l:numpy.load(os.path.join(sp, l+'.npy'))
Hg=readarr('Hg_kJpermol_maxamfrac')
comp=readarr('comp_AuSiCu')
cell=readarr('cells')
ccr=readarr('critcoolingrate_pureaminB_universalexp_maxamfracinBasintercept')
SiCu=readarr('aveSitoCu_compB')
amfrac=readarr('amfrac')

pylab.figure(figsize=(10, 6))
ax=pylab.subplot(111)
stp = TernaryPlot(ax, ellabels=['Au', 'Si', 'Cu']) 
minlist=[c.min() for c in comp.T]
rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
colors=stp.color_comp_calc(comp, rangelist=rangelist)

def colplot(ax, col, x, y, shape='o', **kwargs):
    inds=numpy.where(numpy.logical_not(numpy.isnan(x))&numpy.logical_not(numpy.isnan(y)))[0]
    if 1:
        inds2=numpy.argsort(x[inds])
        ax.plot(x[inds][inds2], y[inds][inds2], 'k-')
    for i in inds:
        ax.plot(x[i], y[i], shape, c=col[i], **kwargs)



pylab.clf()
ax1=pylab.subplot(111)
ax1.set_xscale('log')
ax2=ax1.twinx()
if 1:
    colplot(ax1, colors, ccr, Hg, 'o', ms=10)
    colplot(ax2, colors, ccr, SiCu, 'D', ms=10)
elif 0:
    colplot(ax1, colors, SiCu, ccr, 'o', ms=10)
    colplot(ax2, colors, SiCu, Hg, 'D', ms=10)
    ax1.set_yscale('log')

#ax1.set_ylabel('glass transition enthalpy (kJ/mol)', fontsize=14)
#ax1.set_xlabel('critical cooling rate (K/s)', fontsize=14)
#ax2.set_ylabel('Si:Cu ratio of glass-former', fontsize=14)


pylab.show()
