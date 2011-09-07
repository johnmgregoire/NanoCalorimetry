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

celllist=[1, 2]+range(4, 21)+[22]+[24, 25]
celllist=[11]
for selectcell in celllist:
    print selectcell
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

    def heatrate_T(d, T, Twin=10.):
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
        return d['sampleheatrate'][cycleindex][inds].mean()
        
        

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
    savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_Aug28/'
    


    plotTlim=(50., 700.)


    metadictlist=[]
    allsegdict=[]
    f=h5py.File(p, mode='r')
    if not `selectcell` in f['calbycellmetadata']:
        print 'no cal data for ', selectcell
        f.close()
        continue
    cg=f['calbycellmetadata'][`selectcell`]
    for mg in cg.itervalues():
        if isinstance(mg, h5py.Group) and 'Cpregions_enthalpy' in mg.attrs.keys():
            d={}
            for k, v in mg.attrs.iteritems():
                d[k]=v
    #        if selectcell==1 and d['name'].startswith('heat1'):#heat1a was botched and heat1b we don't know cooling rate and the XRd for heat0 was questionable anyway
    #            continue
            metadictlist+=[d]
            allsegdict+=[CreateHeatProgSegDictList(p, d['name'], d['h5hpname'])]

    xrddictlist=[]
    if 'xrdbycell' in f and `selectcell` in f['xrdbycell']:
        cg=f['xrdbycell'][`selectcell`]
        for mg in cg.itervalues():
            if isinstance(mg, h5py.Group):
                d={}
                for k, v in mg.attrs.iteritems():
                    d[k]=v
                xrddictlist+=[d]
    f.close()

    orderarray=numpy.abs(numpy.array([metadict['prevcoolrate_400C'] for metadict in metadictlist]))

    sortinds=numpy.argsort(orderarray)
    cols=['b', (160./256.,160./256.,0), 'r', 'g', 'c', 'm', 'k']


    ## plotting series of heat ramps
    mult=1.e6
    nplots=len(orderarray)
    pylab.figure(figsize=(8, 8))
    axl=[pylab.subplot(nplots, 1, nplots)]
    for i in range(1, nplots):
        #ax=pylab.subplot2grid((n, 3), (n-1-i, 0), colspan=2, sharex=axl[0], sharey=axl[0])
        #ax=pylab.subplot(nplots, 1, nplots-i, sharex=axl[0], sharey=axl[0])
        ax=pylab.subplot(nplots, 1, nplots-i, sharex=axl[0])
        pylab.setp(ax.get_xticklabels(), visible=False)
        axl+=[ax]

    namestack=[]
    for count, i in  enumerate(sortinds):
        hpsdl=allsegdict[i]
        metadict=metadictlist[i]
        namestack+=[metadict['name']]
    namestack='_'.join(namestack)

    ymin, ymax=None, None
    for count, i in  enumerate(sortinds):
        hpsdl=allsegdict[i]
        metadict=metadictlist[i]
        print metadict['name']
        T=hpsdl[metadict['heatseg']]['sampletemperature'][cycleindex]
        C=hpsdl[metadict['heatseg']]['sampleheatcapacity'][cycleindex]
        tp=hpsdl[metadict['heatseg']]['cyclepartition'][cycleindex]
        PdT=hpsdl[metadict['heatseg']]['samplepowerperrate'][cycleindex]
        if selectcell==10 and metadict['name']=='heat2':
            C=C[T<680]
            T=T[T<680]
            tp=tp[T<680]
            PdT=PdT[T<680]
        if selectcell==4 and metadict['name']=='heat1a':
            print T.max()
            C=C[T<615]
            T=T[T<615]
            tp=tp[T<615]
            PdT=PdT[T<615]
        if selectcell==4 and metadict['name']=='heat1b':
            print T.max()
            C=C[T<665]
            T=T[T<665]
            tp=tp[T<665]
            PdT=PdT[T<665]
        if selectcell==20 and metadict['name']=='heat1a':
            print T.max()
            C=C[T<665]
            T=T[T<665]
            tp=tp[T<665]
            PdT=PdT[T<665]
        if selectcell==20 and metadict['name']=='heat2':
            print T.max()
            C=C[T<635]
            T=T[T<635]
            tp=tp[T<635]
            PdT=PdT[T<635]
            
#        #Cp plots
#        axl[count].plot(T, mult*C, '-', color=cols[count], lw=1, label=metadict['name'])
#        
#        rxnindlist=[metadict['Cpregions_glassind'], metadict['Cpregions_xtalind'], metadict['Cpregions_meltind'], metadict['Cpregions_melt2ind']]
#        for regind, (enth, Tp, Cp, Tmean, i, j) in enumerate(zip(metadict['Cpregions_enthalpy'], metadict['Cpregions_T_Cmax'], metadict['Cpregions_Cmax'], metadict['Cpregions_Tweightedmean'], metadict['Cpregions_cycindstart'], metadict['Cpregions_cycindstop'])):
#            if regind in rxnindlist:
#                col=cols[count]
#                #col=['b', 'g', 'r'][rxnindlist.index(regind)]
#                hatch=['/', '\\', '+', '+'][rxnindlist.index(regind)]
#            else:
#                continue
#            axl[count].fill(T[i:j], mult*C[i:j], color=col, hatch=hatch, alpha=0.3)
#            #axl[count].plot(Tp, mult*Cp, 'kx')
#            #axl[count].plot(Tmean, 0, 'k*')
#    for ax in axl:
#        if selectcell==1:
#            ax.set_ylim(-1.2, 5.4)
#            ax.set_yticks([-1, 0, 2, 4])    
#        elif selectcell==2:
#            ax.set_ylim(-.9, 8.6)
#            ax.set_yticks([0, 2, 4, 6])    
#        elif selectcell==4:
#            ax.set_ylim(-2.5, 11.8)
#            ax.set_yticks([-2, 0, 2, 4, 6, 8])    
#        elif selectcell==5:
#            ax.set_ylim(-2.6, 8.7)
#            ax.set_yticks([-2, 0, 2, 4, 6]) 
#        elif selectcell==8:
#            ax.set_ylim(-1.1, 7.6)
#            ax.set_yticks([-1, 0, 2, 4, 6])    
#        elif selectcell==12:
#            ax.set_ylim(-1.8, 6.8)
#            ax.set_yticks([-1, 0, 2, 4, 6])    
#        elif selectcell==13:
#            ax.set_ylim(-1.6, 7.)
#            ax.set_yticks([-1, 0, 2, 4, 6])
#        elif selectcell==14:
#            ax.set_ylim(-1.8, 6.7)
#            ax.set_yticks([-1, 0, 2, 4, 6])   
#        elif selectcell==18:
#            ax.set_ylim(-1.1, 5.7)
#            ax.set_yticks([-1, 0, 2, 4])
#        elif selectcell==20:
#            ax.set_ylim(-2.1, 5.9)
#            ax.set_yticks([-2, 0, 2, 4])
#        else:
#            ax.set_ylim(-2.1, 4.9)
#            ax.set_yticks([-2, 0, 2, 4])
#        ax.set_xlim(plotTlim)
#    axl[2].set_ylabel(r'Heat Capacity ($\mu$J/K),  endothermic ->', fontsize=14)
#    axl[0].set_xlabel('Temperature (C)', fontsize=14)
#    headname='Cpstack'
##****    
        #Cp fit plots
        colors=['k']+['r', 'g', 'c', 'm', 'y', 'b']*5
        for i in range(-1, tp.max()+1):
            if i<0:
                al=0.5
            else:
                al=1.
            if numpy.any(tp==i):
                axl[count].plot(T[tp==i], mult*PdT[tp==i], '.', markersize=1, color=colors[i+1], alpha=al)
        axl[count].plot(T, mult*(PdT-C), '-', color='k', lw=1)
        temp=PdT[(T>plotTlim[0])&(T<plotTlim[1])]
        if ymin is None:
            ymin=temp.min()
            ymax=temp.max()
        else:
            ymin=min(ymin, temp.min())
            ymax=max(ymax, temp.max())
        axl[count].set_ylim(mult*temp.min(), mult*temp.max())

        
    ymin*=mult
    ymax*=mult
    for ax in axl:
        #ax.set_ylim(ymin, ymax)
#        temp=[yv for l in ax.get_lines() for xd, yd in l.get_data() for yv in yd]
#        ax.set_ylim(min(temp), max(temp))
        ax.set_xlim(plotTlim)

    headname='Cpfitstack'
    axl[2].set_ylabel(r'Power per heat rate ($\mu$J/K)', fontsize=14)
    axl[0].set_xlabel('Temperature (C)', fontsize=14)
#****        
    pylab.subplots_adjust(right=.95, top=0.95, hspace=0.01)
    pylab.savefig(os.path.join(savef, '%s_cell%02d_%s.png' %(headname, selectcell, namestack)))

##BELOW ARE THE EXTRA DATA PLOTS
    #only show xrd data if the prevname pnsc scan was the last performed before an xrd experiment
    phasecomps=[]
    for metadict in metadictlist:
        pcal=metadict['prevname'][:5]
        cal=metadict['name'][:5]
        c=numpy.zeros(3)
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
            for d in xrddictlist:
                if d['name']==pcal:#use xrd that happened after the prev scan
                    c=numpy.float32([d['amfrac'], d['othfrac'], d['fccfrac']])
                    c=c/c.sum()
                    print c
                    break
        phasecomps+=[c]

    ##barplot
    pylab.figure(figsize=(3, 8))
    listofarraysorderedwithmetadictlist=phasecomps
    axl=[pylab.subplot(nplots, 1, nplots)]
    for i in range(1, nplots):
        #ax=pylab.subplot2grid((n, 3), (n-1-i, 0), colspan=2, sharex=axl[0], sharey=axl[0])
        ax=pylab.subplot(nplots, 1, nplots-i, sharex=axl[0], sharey=axl[0])
        pylab.setp(ax.get_xticklabels(), visible=False)
        axl+=[ax]
    maxval=0.
    for count, i in  enumerate(sortinds):
        arr=listofarraysorderedwithmetadictlist[i]
        maxval=max(maxval,numpy.max(arr))
        axl[count].barh(range(len(arr)),arr,color=['y','g','r'],height=1)
    for ax in axl:
        ax.set_ylim(-.5, len(arr)+.5)
        ax.set_xlim((0., maxval*1.1))
        make_ticklabels_invisible(ax, x=False)
    axl[0].set_xlabel('phase fraction', fontsize=14)
    pylab.subplots_adjust(right=.95, top=0.95, hspace=0.01)

    pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'PhaseConcstack_cell%02d_%s.png' %(selectcell, metadict['name'])))
#
#    #hatchlegend
#    pylab.figure()
#    for hatch, lab in zip(['/', '\\', '+'], ['glass trans.', 'crystallization', 'melting']):
#        pylab.fill([0, 1, 0], [0, 1, 1], color='k', hatch=hatch, alpha=0.3, label=lab)
#    pylab.legend()
#    pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'CpRegionLegend_cell%02d_%s.png' %(selectcell, metadict['name'])))
#
#    #plot cooling rates
#    pylab.figure(figsize=(2, 8))
#    for count, i in enumerate(sortinds):
#        pylab.semilogx(numpy.abs(metadictlist[i]['prevcoolrate_400C']), count, 'o', color=cols[count], markersize=11)
#    make_ticklabels_invisible(pylab.gca(), x=False)
#    pylab.xlabel('cooling rate\nat 400C (K/s)', fontsize=14)
#    pylab.ylim(-.5, count+.5)
#    pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'Cool400Cstack_cell%02d_%s.png' %(selectcell, metadict['name'])))
#
#    ##plot cooling rates
#    pylab.figure(figsize=(2, 8))
#    for count, i in enumerate(sortinds):
#        pylab.semilogx(numpy.abs(metadictlist[i]['prevcoolrate_180C']), count, 'o', color=cols[count], markersize=11)
#    make_ticklabels_invisible(pylab.gca(), x=False)
#    pylab.xlabel('cooling rate\nat 180C (K/s)', fontsize=14)
#    pylab.ylim(-.5, count+.5)
#    pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'Cool180Cstack_cell%02d_%s.png' %(selectcell, metadict['name'])))
#
#    ##plot heat rates
#    pylab.figure(figsize=(2, 8))
#    for count, i in enumerate(sortinds):
#        pylab.semilogx(numpy.abs(metadictlist[i]['heatrate_170C500C']), count, 'o', color=cols[count], markersize=11)
#    make_ticklabels_invisible(pylab.gca(), x=False)
#    pylab.xlabel('heating rate (K/s)', fontsize=14)
#    pylab.ylim(-.5, count+.5)
#    pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell), 'HeatRatestack_cell%02d_%s.png' %(selectcell, metadict['name'])))


#pylab.show()
print 'done'
