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
#allheatnames=['heat1a', 'heat1b', 'heat1c', 'heat1d', 'heat2', 'heat3', 'heat4a', 'heat4b','heat4c','heat6a', 'heat6b', 'heat6c','heat7','heat8']
allheatnames=['heat1b', 'heat1c', 'heat1d', 'heat2', 'heat3', 'heat4a', 'heat4b','heat4c','heat6a', 'heat6b', 'heat6c','heat8']
CHOICE='amfracinB'
#CHOICE='Hgperamfrac'
#CHOICE='timeabove'
#CHOICE='coolratevsfcc'
#CHOICE='heatnumvsfcc'
#CHOICE='coolratevsfcccts'
#CHOICE='totalxrdmass_fcc'

xlabdict={\
'amfrac':'cooling rate at 320C (K/s)', \
'amfracinB':'cooling rate at 320C (K/s)', \
'Hgperamfrac':'cooling rate at 320C (K/s)', \
'timeabove':'time above 900C (s)', \
'coolratevsfcc':'cooling rate at 320C (K/s)', \
'coolratevsfcccts':'cooling rate at 320C (K/s)', \
'heatnumvsfcc':'heat number', \
'totalxrdmass_fcc':'fcc counts', \
'totalxrdmass_oth':'oth counts', \
'totalxrdmass_am':'am counts', \
}

ylabdict={\
'amfrac':'amorphous fraction after cooling', \
'amfracinB':'amorphous fractionof\nglass-forming component', \
'Hgperamfrac':r'volumetric glass transition enthalpy (J/cm$^3$)', \
'timeabove':'fcc frac/max', \
'coolratevsfcc':'fcc frac/max', \
'coolratevsfcccts':'fcc counts/max', \
'heatnumvsfcc':'fcc frac/max', \
'totalxrdmass_fcc':'sum_phases: xrdcts/aveZ2  /nm', \
'totalxrdmass_oth':'sum_phases: xrdcts/aveZ2  /nm', \
'totalxrdmass_am':'sum_phases: xrdcts/aveZ2  /nm', \
}


xlogbool={\
'amfrac':True, \
'amfracinB':True, \
'Hgperamfrac':True, \
'timeabove':False, \
'coolratevsfcc':True, \
'coolratevsfcccts':True, \
'heatnumvsfcc':False, \
'totalxrdmass_fcc':False, \
'totalxrdmass_oth':False, \
'totalxrdmass_am':False, \
}

mdboolfcndict={\
'amfrac':lambda md:'prevcoolrate_320C' in md.keys(), \
'amfracinB':lambda md:'prevcoolrate_320C' in md.keys(), \
'Hgperamfrac':lambda md:'prevcoolrate_320C' in md.keys() and 'Cpregions_enthalpy' in md.keys() and numpy.abs(md['prevcoolrate_320C'])>1.e3, \
'timeabove':lambda md:'secabove900C' in md.keys(), \
'coolratevsfcc':lambda md:'prevcoolrate_320C' in md.keys(), \
'coolratevsfcccts':lambda md:'prevcoolrate_320C' in md.keys(), \
'heatnumvsfcc':lambda md:True, \
'totalxrdmass_fcc':lambda md:True, \
'totalxrdmass_oth':lambda md:True, \
'totalxrdmass_am':lambda md:True, \
}

xdboolfcndict={\
'amfrac':lambda xd:True, \
'amfracinB':lambda xd:True, \
'Hgperamfrac':lambda xd:'amfrac' in xd.keys() and xd['amfrac']>0.02, \
'timeabove':lambda xd:'fccfrac' in xd.keys(), \
'coolratevsfcc':lambda xd:'fccfrac' in xd.keys(), \
'coolratevsfcccts':lambda xd:'fccctsinrange' in xd.keys(), \
'heatnumvsfcc':lambda xd:'fccfrac' in xd.keys(), \
'totalxrdmass_fcc':lambda xd:'fccctsinrange' in xd.keys(), \
'totalxrdmass_oth':lambda xd:'fccctsinrange' in xd.keys(), \
'totalxrdmass_am':lambda xd:'fccctsinrange' in xd.keys(), \
}

fomfcndict={\
'amfrac':lambda x, y, z:y, \
'amfracinB':lambda x, y, z:y/(y+z), \
'Hgperamfrac':lambda x, y, z:y/z*1.e-6/(3.6*.8*1.e-15), \
'timeabove':lambda x, y, z:y/y.max(), \
'coolratevsfcc':lambda x, y, z:y/y.max(), \
'coolratevsfcccts':lambda x, y, z:y/y.max(), \
'heatnumvsfcc':lambda x, y, z:y/y.max(), \
'totalxrdmass_fcc':lambda x, y, z:y/z, \
'totalxrdmass_oth':lambda x, y, z:y/z, \
'totalxrdmass_am':lambda x, y, z:y/z, \
}


xfcndict={\
'amfrac':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
'amfracinB':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
'Hgperamfrac':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
'timeabove':lambda md, xd:md['secabove900C'], \
'coolratevsfcc':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
'coolratevsfcccts':lambda md, xd:numpy.abs(md['prevcoolrate_320C']), \
'heatnumvsfcc':lambda md, xd:allheatnames.index(md['name']), \
'totalxrdmass_fcc':lambda md, xd:xd['fccctsinrange'], \
'totalxrdmass_oth':lambda md, xd:xd['othctsinrange'], \
'totalxrdmass_am':lambda md, xd:xd['amctsinrange'], \
}


totxrdfcn=lambda md, xd: xd['fccctsinrange']/(md['compfcc']*numpy.array([79., 14., 29.])**2).sum()+(xd['othctsinrange']+xd['amctsinrange'])/(md['compoth']*numpy.array([79., 14., 29.])**2).sum()

yfcndict={\
'amfrac':lambda md, xd:xd['amfrac'], \
'amfracinB':lambda md, xd:xd['amfrac'], \
'Hgperamfrac':lambda md, xd:(md['Cpregions_glassind']<0 and (0., ) or (md['Cpregions_enthalpy'][md['Cpregions_glassind']], ))[0], \
'timeabove':lambda md, xd:xd['fccfrac'], \
'coolratevsfcc':lambda md, xd:xd['fccfrac'], \
'coolratevsfcccts':lambda md, xd:xd['fccctsinrange'], \
'heatnumvsfcc':lambda md, xd:xd['fccfrac'], \
'totalxrdmass_fcc':totxrdfcn, \
'totalxrdmass_oth':totxrdfcn, \
'totalxrdmass_am':totxrdfcn, \
}

zfcndict={\
'amfrac':lambda md, xd:xd['othfrac'], \
'amfracinB':lambda md, xd:xd['othfrac'], \
'Hgperamfrac':lambda md, xd:xd['amfrac']*md['nm'], \
'timeabove':lambda md, xd:1., \
'coolratevsfcc':lambda md, xd:1., \
'coolratevsfcccts':lambda md, xd:1., \
'heatnumvsfcc':lambda md, xd:1., \
'totalxrdmass_fcc':lambda md, xd:md['nm'], \
'totalxrdmass_oth':lambda md, xd:md['nm'], \
'totalxrdmass_am':lambda md, xd:md['nm'], \
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
compfcc=f['CompThick/atfrac_fcc'][:, :]
compoth=f['CompThick/atfrac_oth'][:, :]
nm=f['CompThick/nm'][:]

metadictlists=[]
for selectcell in range(1, 26):
    metadictlist=[]
    #allsegdict=[]
    if `selectcell` in f['calbycellmetadata']:
        cg=f['calbycellmetadata'][`selectcell`]
        for mg in cg.itervalues():
            if isinstance(mg, h5py.Group) and mg.attrs['name'] in allheatnames:
                d={}
                for k, v in mg.attrs.iteritems():
                    d[k]=v
                
        #        if selectcell==1 and d['name'].startswith('heat1'):#heat1a was botched and heat1b we don't know cooling rate and the XRd for heat0 was questionable anyway
        #            continue
                metadictlist+=[d]
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

for cell, cv, cfv, cov, nmv, mdl, xdl in zip(range(1, 26), comp, compfcc, compoth, nm, metadictlists, xrddictlists):
    print cell
    xtemp=[]
    ytemp=[]
    ztemp=[]
    ctemp=[]
    celltemp=-1
    print len(mdl)
    for md in mdl:
        print '$', mdboolfcndict[CHOICE](md)
        if not mdboolfcndict[CHOICE](md):
            continue
        md['nm']=nmv
        md['comp']=cv
        md['compfcc']=cfv
        md['compoth']=cov
        pcal=md['prevname'][:5]
        cal=md['name'][:5]
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
            for xd in xdl:
                print '*', xdboolfcndict[CHOICE](xd)
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


pylab.figure(figsize=(6, 4))
ax=pylab.subplot(111)
stp = TernaryPlot(ax, ellabels=['Au', 'Si', 'Cu']) 
celllist=[1, 2]+range(4, 21)+[22]+[24, 25]
minlist=[c.min() for c in comp[numpy.array(celllist)-1, :].T]
rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
colors=stp.color_comp_calc(comp[numpy.array(celllist)-1, :], rangelist=rangelist)

pylab.clf()
fommin=None
fommax=None
for cell, xv, yv, zv in zip(cell_cell, x_cell, y_cell, z_cell):
    xv=numpy.array(xv)
    yv=numpy.array(yv)
    zv=numpy.array(zv)
    sortarr=numpy.argsort(xv)
    #pylab.plot(xv[sortarr], yv[sortarr], '.', color=colors[celllist.index(cell)], markersize=16)
    #pylab.plot(xv[sortarr], yv[sortarr], '-', color=colors[celllist.index(cell)], linewidth=1, alpha=.5)
    
    a=fomfcndict[CHOICE](xv[sortarr], yv[sortarr], zv[sortarr])
    fommin=(fommin is None and (a.min(), ) or (min(fommin, a.min()),))[0]
    fommax=(fommax is None and (a.max(), ) or (max(fommax, a.max()),))[0]
    if 1:
        pylab.plot(xv[sortarr], a, '.', color=colors[celllist.index(cell)], markersize=10)
    else:
        for xx, aa in zip(xv[sortarr], a):
            pylab.text(xx, aa, `cell`, ha='center', va='center', color=colors[celllist.index(cell)], fontsize=9)
    pylab.plot(xv[sortarr], a, '-', color=colors[celllist.index(cell)], linewidth=1, alpha=.5)
    
    print cell, xv[sortarr], a

if xlogbool[CHOICE]:
    pylab.gca().set_xscale('log')

xmin=x.min()-.05*(x.max()-x.min())
xmax=x.max()+.05*(x.max()-x.min())

xlabel=xlabdict[CHOICE]
ylabel=ylabdict[CHOICE]

ymin=fommin-.05*(fommax-fommin)
ymax=fommax+.05*(fommax-fommin)

pylab.ylabel(ylabel, fontsize=14)
pylab.xlabel(xlabel, fontsize=14)
#pylab.title('colored by composition', fontsize=14)
#pylab.ylabel(, fontsize=14)
pylab.xlim(xmin, xmax)
pylab.ylim(ymin, ymax)

pylab.plot([3.e2, 3.e3], [.6, 1.], 'k-', lw=1.6)
pylab.text(3.*10**2.5, .75, '0.4/decade', color='k', va='bottom', ha='center', rotation=27, fontsize=14)
pylab.subplots_adjust(left=.17, bottom=.15)

pylab.show()
