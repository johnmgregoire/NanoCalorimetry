from scipy.special import erf as erf
from scipy.special import erfinv as erfinv
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




def area_gauss(p, x=numpy.inf):
    return (numpy.pi/2.)**.5*p[1]*p[2]*(1.+erf((x-p[0])/(2.**.5*p[1])))

def area_lorentz(p, x=numpy.inf):
    return p[1]*p[2]*(numpy.arctan((x-p[0])/p[1])+numpy.pi/2.)
    
def area_gausslorentz(p, x=numpy.inf):
    return p[3]*area_gauss(p, x)+(1.-p[3])*area_lorentz(p, x)

def x_gaussarea(area, p):
    return p[0]+p[1]*2.**.5*erfinv(area/((numpy.pi/2.)**.5*p[1]*p[2])-1.)

def x_lorentzarea(area, p):
    return p[0]+p[1]*numpy.tan(area/(p[1]*p[2])-numpy.pi/2.)

def x_gausslorentzarea(area, p):
    guess=p[3]*x_gaussarea(area*p[3], p)+(1.-p[3])*x_lorentzarea(area*(1.-p[3]), p)
    #print '**', area, area_gausslorentz(p, guess)
    fcn=lambda x, p:area-area_gausslorentz(p, x)
    ans=scipy.optimize.fsolve(fcn, guess, args=(p))
    #print ans
    return ans[0]
    
    
    
p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'



cycleindex=0

#p=mm.h5path
#f=h5py.File(p, mode='r+')
#f=h5py.File(p, mode='r')
savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2/fitglass'


plotTlim=(50., 700.)

f=h5py.File(p, mode='r')
comp=f['CompThick/atfrac'][:, :]

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


dictlists=[]
for selectcell in range(1, 26):
    dictlist=[]
    #allsegdict=[]
    if `selectcell` in f['calbycellmetadata']:
        cg=f['calbycellmetadata'][`selectcell`]
        for mg in cg.itervalues():
            if isinstance(mg, h5py.Group) and 'Cpregions_glassind' in mg.attrs.keys() and mg.attrs['Cpregions_glassind']>=0:
                d={}
                for k, v in mg.attrs.iteritems():
                    d[k]=v
                hpsdl=CreateHeatProgSegDictList(p, d['name'], d['h5hpname'])
                d['sampleheatcapacity']=hpsdl[d['heatseg']]['sampleheatcapacity'][cycleindex]
                d['sampletemperature']=hpsdl[d['heatseg']]['sampletemperature'][cycleindex]
                dictlist+=[d]
                
    dictlists+=[dictlist]

#xrddictlists=[]
#for selectcell in range(1, 26):
#    xrddictlist=[]
#    if 'xrdbycell' in f and `selectcell` in f['xrdbycell']:
#        cg=f['xrdbycell'][`selectcell`]
#        for mg in cg.itervalues():
#            if isinstance(mg, h5py.Group):
#                d={}
#                for k, v in mg.attrs.iteritems():
#                    d[k]=v
#                xrddictlist+=[d]
#    xrddictlists+=[xrddictlist]
#f.close()

#dlist=[]
#for selectcell, cv, mdl, xdl in zip(range(1, 26), comp, metadictlists, xrddictlists):
#    dlisttemp=[]
#    for md in mdl:
#        if not 'Cpregions_glassind' in md.keys() or md['Cpregions_glassind']<0:
#            continue
#        pcal=md['prevname'][:5]
#        cal=md['name'][:5]
#        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
#            for xd in xdl:
#                if xd['name']==pcal:#use xrd that happened after the prev scan
#                    dlisttemp+=[dict([('cell', selectcell), ('comp', comp[selectcell-1, :]), ('amfrac', xd['amfrac']), ('prevcoolrate_400C', md['prevcoolrate_400C'])])]
#                    break
#    for d in dlisttemp:
#        d['numdataforcell']=len(dlisttemp)
#        dlist+=[d]

f.close()

for cell, dlist in zip(range(1, 26), dictlists):
    for d in dlist:
        gi=d['Cpregions_glassind']
        C=d['sampleheatcapacity'][d['Cpregions_cycindstart'][gi]:d['Cpregions_cycindstop'][gi]]
        T=d['sampletemperature'][d['Cpregions_cycindstart'][gi]:d['Cpregions_cycindstop'][gi]]

        twm=d['Cpregions_Tweightedmean'][gi]
        itwm=numpy.argmin((T-twm)**2)
        critenth=d['Cpregions_enthalpy'][gi]*.05

        cmax=numpy.max(C[:itwm])
        imid=numpy.argmin((C[:itwm]-cmax/2.)**2)
        tsigest=numpy.abs(twm-T[imid])
        x, y=T[:itwm], C[:itwm]

        fitclass=fitfcns()
        fcn=fitclass.genfit(Gaussian, [twm, tsigest, cmax], (x, y))
        gl=True
        if gl:
            g=fitclass.finalparams
            fcn=fitclass.genfit(GaussLorentz, list(g)+[.8], (x, y))

        ans=fitclass.finalparams
        ans[1]*=numpy.sign(ans[1])
        if gl:
            ans[3]=min(max(ans[3], 0.), 1.)


        if gl:
            xs=x_gausslorentzarea(critenth, ans)
        else:
            xs=x_gaussarea(critenth, ans)
        yf=fcn(x)
        #print ans, totarea, totarea*.5, xs, area_gausslorentz(ans, x=xs)

        pylab.plot(x, y, 'b.')
        pylab.plot(x, yf, 'k-')
        pylab.title('crit enthalpy is %.2e, T=%.1f' %(critenth, xs))
        try:
            pylab.fill(numpy.concatenate([x[x<xs], x[x<xs][[-1, 0]]]), numpy.concatenate([yf[x<xs], [0, 0]]), 'r')
        except:
            pass

        pylab.savefig(os.path.join(savef, 'cell%02d_%s.png' %(cell, d['name'])))

        #pylab.show()
        pylab.clf()



    
    
    
    
###test with clean data
#x=numpy.linspace(0, 1, 100)
##y=Gaussian([.5, .1, 10.], x)
#y=GaussLorentz([.5, .1, 10., .6], x)
#fitclass=fitfcns()
#
#fcn=fitclass.genfit(Gaussian, [.4, .05, 12.], (x, y))
#g=fitclass.finalparams
#fcn=fitclass.genfit(GaussLorentz, list(g)+[.5], (x, y))
#
##fcn=fitclass.genfit(Gaussian, [.5, .05, 10.], (x, y))
#ans=fitclass.finalparams
#ans[1]*=numpy.sign(ans[1])
#ans[3]=min(max(ans[3], 0.), 1.)
#
#totarea=area_gausslorentz(ans)
#xs=x_gausslorentzarea(totarea*.03, ans)
##totarea=area_gauss(ans)
##xs=x_gaussarea(totarea*.5, ans)
#yf=fcn(x)
##print ans, totarea, totarea*.5, xs, area_gausslorentz(ans, x=xs)
#
#pylab.subplot(211)
#pylab.plot(x, y, 'b.')
#pylab.plot(x, yf, 'k-')
#pylab.fill(numpy.concatenate([x[x<xs], x[x<xs][[-1, 0]]]), numpy.concatenate([yf[x<xs], [0, 0]]), 'r')
#
#y2=y-1.
#fcn2=fitclass.genfit(Gaussian, [.4, .05, 12.], (x, y2))
#g2=fitclass.finalparams
##fcn2=fitclass.genfit(Lorentzian, [.4, .05, 12.], (x, y2))
##l2=fitclass.finalparams
#fcn2=fitclass.genfit(GaussLorentz, list(g2)+[.5], (x, y2))
#ans2=fitclass.finalparams
#ans2[1]*=numpy.sign(ans2[1])
#ans2[3]=min(max(ans2[3], 0.), 1.)
#
#totarea2=area_gausslorentz(ans2)
#xs2=x_gausslorentzarea(totarea2*.03, ans2)
#yf2=fcn2(x)
##print ans2, totarea2, totarea2*.5, xs2, area_gausslorentz(ans2, x=xs2)
#pylab.subplot(212)
#pylab.plot(x, y2, 'b.')
#pylab.plot(x, yf2, 'k-')
#pylab.fill(numpy.concatenate([x[x<xs2], x[x<xs2][[-1, 0]]]), numpy.concatenate([yf2[x<xs2], [0, 0]]), 'r')
#
#pylab.show()
