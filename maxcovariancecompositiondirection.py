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


cycleindex=0

#p=mm.h5path
#f=h5py.File(p, mode='r+')
#f=h5py.File(p, mode='r')
savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'


plotTlim=(50., 700.)

f=h5py.File(p, mode='r')
#metadictlist=[]
#allsegdict=[]
#cg=f['calbycellmetadata'][`selectcell`]
#for mg in cg.itervalues():
#    if isinstance(mg, h5py.Group) and 'Cpregions_enthalpy' in mg.attrs.keys():
#        d={}
#        for k, v in mg.attrs.iteritems():
#            d[k]=v
##        if selectcell==1 and d['name'].startswith('heat1'):#heat1a was botched and heat1b we don't know cooling rate and the XRd for heat0 was questionable anyway
##            continue
#        metadictlist+=[d]
#        allsegdict+=[CreateHeatProgSegDictList(p, d['name'], d['h5hpname'])]

comp=f['CompThick/atfrac'][:, :]
cell_comp_weight1_weight2_fomlist=[]
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
    if len(xrddictlist)>0:
        amfrac=numpy.array([d['amfrac'] for d in xrddictlist])
        cell_comp_weight1_weight2_fomlist+=[(selectcell, comp[selectcell-1, :], amfrac, [1./len(amfrac)]*len(amfrac), amfrac)]
f.close()

xl=[]
yl=[]
weightl1=[]
weightl2=[]
foml=[]
for cell, comp, weight1, weight2, fomlist in cell_comp_weight1_weight2_fomlist:
    xl+=[comp[0]]*len(fomlist)
    yl+=[comp[1]]*len(fomlist)
    weightl1+=list(weight1)
    weightl2+=list(weight2)
    foml+=list(fomlist)
xl=numpy.float64(xl)
yl=numpy.float64(yl)
weightl1=numpy.float64(weightl1)
weightl2=numpy.float64(weightl2)
foml=numpy.float64(foml)
fomdev=foml-foml.mean()

##find composition that makes maximum covariance of amfrac and disstance from composition
def minfcn(ab, x, y, wts, fomdevfrommean):
    a=ab[0]
    b=ab[1]
#    xr=(xp+a*yp-a*b)/(1.+a**2)
#    yr=a*xr+b
#    b-=yr.mean()
#    d=xp+a*yp-a*b
    d=(((x-a)**2+(y-b)**2+(x+y-a-b)**2)/2.)**.5
    return -1./(fomdevfrommean*wts*(d-d.mean())).sum()
xp=1.-xl-yl/2.
yp=3**.5*yl/2.

bestans=None
cov=None
for ag, bg in zip(xl, yl):
    ans, evals, code=scipy.optimize.fmin_tnc(minfcn, [ag, bg], args=(xl, yl, weightl1*weightl2, fomdev), approx_grad=1)#, bounds=[(0., 1.), (0., 1.)])
    if cov is None or minfcn(ans, xl, yl, weightl1*weightl2, fomdev)<cov:
        bestans=ans
        cov=minfcn(ans, xl, yl, weightl1*weightl2, fomdev)
        print ans, cov
a=ans[0]
b=ans[1]
#***

#show optimal composition and distances
print a, b
xr=1.-a-b/2.
yr=3**.5*b/2.
pylab.plot(xp, yp, 'b.')
pylab.plot([xr], [yr], 'r.')
for x1, y1 in zip(xp, yp):
    pylab.plot([x1, xr], [y1, yr], 'k-')

pylab.gca().set_aspect(1.)
pylab.show()

##find line that when points are collapsed onto it, the covariance between amfrac and distance along line is maximized - slope is maximized and for each slope the intercept that is centers line on data is chosen
#bminfcn=lambda b, xp, yp, a: ((yp-a*xp-b)**2).sum()
#def minfcn(a, xp, yp, wts, fomdevfrommean):
#    b=scipy.optimize.fmin(bminfcn, 1., args=(xp, yp, a))
#    xr=(xp+a*yp-a*b)/(1.+a**2)
#    yr=a*xr+b
#    d=(xp+a*yp-a*b)/(1.+a**2)**.5
#    return numpy.abs(1./(fomdevfrommean*wts*(d-d.mean())).sum())
#xp=1.-xl-yl/2.
#yp=3**.5*yl/2.
#ans=scipy.optimize.fmin(minfcn, -.4, args=(xp, yp, weightl1*weightl2, fomdev))#, maxfun=50)
#a=ans[0]
##***

##find line that spaces out the projected points
#bminfcn=lambda b, xp, yp, a: ((yp-a*xp-b)**2).sum()
#def minfcn(a, xp, yp, wts, fomdevfrommean):
#    b=scipy.optimize.fmin(bminfcn, 1., args=(xp, yp, a))
#    xr=(xp+a*yp-a*b)/(1.+a**2)
#    yr=a*xr+b
#    d=(xp+a*yp-a*b)/(1.+a**2)**.5
#    return (wts/numpy.array([numpy.min((d[d!=v]-v)**2) for i, v in enumerate(d)])).sum()
#
#xp=1.-xl-yl/2.
#yp=3**.5*yl/2.
#ans=scipy.optimize.fmin(minfcn, -.4, args=(xp, yp, weightl2, fomdev))#, maxfun=50)
#a=ans[0]
##****

#line collapse plot for above routines that have optimal lines
#b=scipy.optimize.fmin(bminfcn, 1., args=(xp, yp, a))[0]
#print a, b
#xr=(xp+a*yp-a*b)/(1.+a**2)
#yr=a*xr+b
#pylab.plot(xp, yp, 'b.')
#pylab.plot(xr, yr, 'r.')
#for x1, y1, x2, y2 in zip(xp, yp, xr, yr):
#    pylab.plot([x1, x2], [y1, y2], 'k-')


#def minfcn(ab, x, y, wts, fomdevfrommean):
#    a=ab[0]
#    b=ab[1]
#    d=(((x-a)**2+(y-b)**2+(x+y-a-b)**2)/2.)**.5
#    return 1./(fomdevfrommean*(d-d.mean())*wts).sum()
#ans=scipy.optimize.fmin_tnc(minfcn, [.4, .3], args=(xl, yl, weightl, fomdev), bounds=[(0., 1.), (0., 1.)], approx_grad=1)
#a, b=ans[0]
#axp=1.-a-b/2.
#ayp=3**.5*b/2.
#xp=1.-xl-yl/2.
#yp=3**.5*yl/2.
#pylab.plot(xp, yp, 'b.')
#pylab.plot(axp, ayp, 'r.')
#pylab.gca().set_aspect(1.)
#pylab.show()
