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


x=[]
y=[]
z=[]
d=[]
c=[]

for cv, mdl, xdl in zip(comp, metadictlists, xrddictlists):
    for md in mdl:
        if not 'prevcoolrate_400C' in md.keys():
            continue
        pcal=md['prevname'][:5]
        cal=md['name'][:5]
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
            for xd in xdl:
                if xd['name']==pcal:#use xrd that happened after the prev scan
                    z+=[xd['amfrac']]
                    x+=[numpy.abs(md['prevcoolrate_400C'])]
                    c+=[cv]
                    break

c=numpy.array(c)
z=numpy.array(z)
x=numpy.array(x)

#results of optimal composition line
a,b=0.783984375,  -0.1214
xp=1.-c[:, 0]-c[:, 1]/2.
yp=3**.5*c[:, 1]/2.
xr=(xp+a*yp-a*b)/(1.+a**2)
yr=a*xr+b
d=(xp+a*yp-a*b)/(1.+a**2)**.5
y=d-d.min()
endpointcomp=[]
for i in [numpy.argmin(d), numpy.argmax(d)]:
    si=yp[i]*2./3**.5
    au=1.-xp[i]-si/2.
    cu=1.-si-au
    endpointcomp+=[[au, si, cu]]
ylab=r'composition, Au$_{%.2f}$Si$_{%.2f}$Cu$_{%.2f}$ to Au$_{%.2f}$Si$_{%.2f}$Cu$_{%.2f}$' %tuple(endpointcomp[0]+endpointcomp[1])
#****

##distance from select comp
#a, b=0.238083046206,  0.36854
#d=(((c[:, 0]-a)**2+(c[:, 1]-b)**2+(c[:, 0]+c[:, 1]-a-b)**2)/2.)**.5
#y=d
#ylab=r'composition, distance from Au$_{%.2f}$Si$_{%.2f}$Cu$_{%.2f}$' %(a, b, 1.-a-b)
##****

##si content
#y=c[:, 1]
#ylab='at. frac. Si'
##****


xmin=x.min()-.05*(x.max()-x.min())
xmax=x.max()+.05*(x.max()-x.min())
ymin=y.min()-.05*(y.max()-y.min())
ymax=y.max()+.05*(y.max()-y.min())


#npts = 200
#xi = numpy.linspace(xmin, xmax,npts)
#yi = numpy.linspace(ymin, ymax,npts)
#zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
#zi=numpy.ma.masked_where(numpy.isnan(zi), zi)
# contour the gridded data, plotting dots at the randomly spaced data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = pylab.contourf(xi,yi,zi,15,cmap=cm.jet)


# plot data points.
pylab.scatter(x,y,marker='o',c=z,s=25, cmap=cm.jet)
#pylab.scatter(x,y,marker='s',c='b',s=5)
pylab.gca().set_xscale('log')
pylab.colorbar() # draw colorbar

orderarray=numpy.abs(numpy.array([metadict['prevcoolrate_400C'] for metadict in metadictlist]))

sortinds=numpy.argsort(orderarray)
cols=['b', (160./256.,160./256.,0), 'r', 'g', 'c', 'm', 'k']

pylab.title('amorphous fraction after cooling')
pylab.xlabel('cooling rate (K/s)', fontsize=14)
pylab.ylabel(ylab, fontsize=14)
#pylab.ylabel(, fontsize=14)
pylab.xlim(xmin, xmax)
pylab.ylim(ymin, ymax)
#pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'HeatRatestack_cell%02d_%s.png' %(selectcell, metadict['name'])))

pylab.show()
