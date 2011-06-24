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

dlist=[]
for selectcell, cv, mdl, xdl in zip(range(1, 26), comp, metadictlists, xrddictlists):
    dlisttemp=[]
    for md in mdl:
        if not 'prevcoolrate_400C' in md.keys() or numpy.abs(md['prevcoolrate_400C'])<1.e3:
            continue
        pcal=md['prevname'][:5]
        cal=md['name'][:5]
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
            for xd in xdl:
                if xd['name']==pcal:#use xrd that happened after the prev scan
                    dlisttemp+=[dict([('cell', selectcell), ('comp', comp[selectcell-1, :]), ('amfrac', xd['amfrac']), ('prevcoolrate_400C', md['prevcoolrate_400C'])])]
                    break
    for d in dlisttemp:
        d['numdataforcell']=len(dlisttemp)
        dlist+=[d]

f.close()

xl=numpy.float64([d['comp'][0] for d in dlist])
yl=numpy.float64([d['comp'][1] for d in dlist])
weights1=numpy.float64([d['amfrac']/d['numdataforcell'] for d in dlist])
weights2=numpy.float64([d['amfrac'] for d in dlist])

foml1=[d['amfrac'] for d in dlist]
foml1=numpy.float64(foml1)
fomdev1=foml1-foml1.mean()


foml2=[d['amfrac']/numpy.abs(d['prevcoolrate_400C']) for d in dlist]
foml2=numpy.float64(foml2)
fomdev2=foml2-foml2.mean()

##find composition that makes maximum covariance of amfrac and disstance from composition
def minfcn(ab, x, y, wts, fomdevfrommean):
    a=ab[0]
    b=ab[1]
    d=(((x-a)**2+(y-b)**2+(x+y-a-b)**2)/2.)**.5
    return ((fomdevfrommean*wts*(d-d.mean())).sum())
xp=1.-xl-yl/2.
yp=3**.5*yl/2.

bestans=None
cov=None
for ag, bg in zip(xl, yl):
    fomd=fomdev2
    ans, evals, code=scipy.optimize.fmin_tnc(minfcn, [ag, bg], args=(xl, yl, weights1, fomd), approx_grad=1, disp=0, bounds=[(0., 1.), (0., 1.)])
    #ans, evals, code=scipy.optimize.leastsq(minfcn, numpy.float64([ag, bg]), args=(xl, yl, weights1, fomd))
    if cov is None or minfcn(ans, xl, yl, weights1, fomd)<cov:
        bestans=ans
        cov=minfcn(ans, xl, yl, weights1, fomd)
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
#ans=scipy.optimize.fmin(minfcn, -.4, args=(xp, yp, weights1, fomdev2))#, maxfun=50)
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
#ans=scipy.optimize.fmin(minfcn, -.4, args=(xp, yp, weights2, fomdev1))#, maxfun=50)
#a=ans[0]
##****

##line collapse plot for above routines that have optimal lines
#b=scipy.optimize.fmin(bminfcn, 1., args=(xp, yp, a))[0]
#print a, b
#xr=(xp+a*yp-a*b)/(1.+a**2)
#yr=a*xr+b
#pylab.plot(xp, yp, 'b.')
#pylab.plot(xr, yr, 'r.')
#for x1, y1, x2, y2 in zip(xp, yp, xr, yr):
#    pylab.plot([x1, x2], [y1, y2], 'k-')
#pylab.gca().set_aspect(1.)

pylab.show()
