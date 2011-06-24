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
    #p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'
    #p=mm.h5path
    #f=h5py.File(p, mode='r+')
    #f=h5py.File(p, mode='r')
    savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'

    plotTlim=(50., 700.)

    expdicts=[\
    dict([('name', 'heat1a'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat1b'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat1c'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat1d'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat2'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat3'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat4a'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat4b'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat4c'), ('heatseg', 2), ('coolseg', 3)]), \
    dict([('name', 'heat6a'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat6b'), ('heatseg', 2), ('coolseg', 4)]), \
    dict([('name', 'heat6c'), ('heatseg', 2), ('coolseg', 3)]), \
    dict([('name', 'heat7'), ('heatseg', 2), ('coolseg', 3)]), \
    dict([('name', 'heat8'), ('heatseg', 2), ('coolseg', 4)]), \
    ]

    metadictlist=[]
    allsegdict=[]
    for ed in expdicts:
        exp=ed['name']
        f, hppaths=experimenthppaths(p, exp)
        f.close()
        hpsdl=None
        saveh5hpname=None
        for hpp in hppaths:
            h5hpname=hpp.rpartition('/')[2]
            f, g=gethpgroup(p, exp, h5hpname=h5hpname)
            cell=g.attrs['CELLNUMBER']
            f.close()
            if cell!=selectcell:
                continue
            hpsdl=CreateHeatProgSegDictList(p, exp, h5hpname)
            saveh5hpname=h5hpname
        if not hpsdl is None:
            allsegdict+=[hpsdl]
            ed['h5hpname']=saveh5hpname
            metadictlist+=[copy.deepcopy(ed)]

    for i, (metadict, hpsdl) in enumerate(zip(metadictlist, allsegdict)):
        if hpsdl is None:
            continue
        prevhpsdl=None
        for tempmetadict, prevhpsdl in zip(metadictlist[:i][::-1], allsegdict[:i][::-1]):
            if not prevhpsdl is None:
                prevmetadict=tempmetadict
                #prevmetadict=metadictlist[allsegdict.index(prevhpsdl)]
                break
        if prevhpsdl is None:
            metadict['prevcoolrate_180C']=-1.e2
            metadict['prevcoolrate_250C']=-1.e2
            metadict['prevcoolrate_400C']=-1.e2
            metadict['prevname']='heat0'
        else:
            metadict['prevcoolrate_180C']=heatrate_T(prevhpsdl[prevmetadict['coolseg']], 180.)
            metadict['prevcoolrate_250C']=heatrate_T(prevhpsdl[prevmetadict['coolseg']], 250.)
            metadict['prevcoolrate_400C']=heatrate_T(prevhpsdl[prevmetadict['coolseg']], 400.)
            metadict['prevname']=prevmetadict['name']
        print metadict['name'], ', previous scan ', metadict['prevname']
        metadict['heatrate_170C500C']=heatrate_T(hpsdl[metadict['heatseg']], 335., Twin=165.)
        if not 'sampleheatcapacity' in hpsdl[metadict['heatseg']].keys():
            continue
        enthdlist=findenthalpyandpinacles(hpsdl[metadict['heatseg']], critenth=1.e-6, dTmin=10.)
        if len(enthdlist)>0:
            for k in enthdlist[0].keys():
                metadict['Cpregions_'+k]=[]
                for d in enthdlist:
                    metadict['Cpregions_'+k]+=[d[k]]
                metadict['Cpregions_'+k]=numpy.array(metadict['Cpregions_'+k])

        metadict['Cpregions_glassind']=-1
        metadict['Cpregions_xtalind']=-1
        metadict['Cpregions_meltind']=-1
        metadict['Cpregions_melt2ind']=-1
        ctemp=0.
        for enthind, d in enumerate(enthdlist):
            if d['Tweightedmean']<100. or d['Cmax']<0:
                continue
            if d['Tweightedmean']>250.:
                break
            if metadict['Cpregions_glassind']<0 or numpy.abs(d['enthalpy'])>ctemp:
                metadict['Cpregions_glassind']=enthind
                ctemp=numpy.abs(d['enthalpy'])
        ctemp=0.
        for enthind, d in enumerate(enthdlist):
            if d['Tweightedmean']<180. or d['Cmax']>0:
                continue
            #if metadict['Cpregions_glassind']<0 or d['Tweightedmean']>450.:
            if d['Tweightedmean']>450.:
                break
            if metadict['Cpregions_xtalind']<0 or numpy.abs(d['enthalpy'])>ctemp:
                metadict['Cpregions_xtalind']=enthind
                ctemp=numpy.abs(d['enthalpy'])
        ctemp=0.
        ctemp2=0.
        for enthind, d in enumerate(enthdlist):
            if d['Tweightedmean']<300. or d['Cmax']<0:
                continue
            if d['Tweightedmean']>600.:
                break
            if metadict['Cpregions_meltind']<0:
                metadict['Cpregions_meltind']=enthind
                ctemp=numpy.abs(d['enthalpy'])
            else:
                if numpy.abs(d['enthalpy'])>ctemp:
                    if metadict['Cpregions_meltind']<0 or ctemp>ctemp2:
                        metadict['Cpregions_melt2ind']=metadict['Cpregions_meltind']
                        ctemp2=ctemp
                    metadict['Cpregions_meltind']=enthind
                    ctemp=numpy.abs(d['enthalpy'])
                elif numpy.abs(d['enthalpy'])>ctemp2:
                    metadict['Cpregions_melt2ind']=enthind
                    ctemp2=numpy.abs(d['enthalpy'])
                
        if selectcell==1 and metadict['name']=='heat8':
            metadict['Cpregions_glassind']=-1
        if selectcell==7 and metadict['name']=='heat1a':
            metadict['Cpregions_xtalind']=-1
        if selectcell==12 and metadict['name']=='heat1a':
            metadict['Cpregions_glassind']=-1
            metadict['Cpregions_xtalind']=-1
        if selectcell==13 and metadict['name']=='heat1a':
            metadict['Cpregions_glassind']=-1
            metadict['Cpregions_xtalind']=-1
        if selectcell==14 and metadict['name']=='heat1a':
            metadict['Cpregions_glassind']=-1
            metadict['Cpregions_xtalind']=-1
        if selectcell==16 and metadict['name']=='heat1a':
            metadict['Cpregions_glassind']=-1
            metadict['Cpregions_xtalind']=-1
        if selectcell==17 and metadict['name']=='heat1a':
            metadict['Cpregions_glassind']=-1
            metadict['Cpregions_xtalind']=-1
        if selectcell==18 and metadict['name']=='heat8':
            metadict['Cpregions_glassind']=-1
        if selectcell==22 and metadict['name']=='heat1a':
            metadict['Cpregions_glassind']=-1
            metadict['Cpregions_xtalind']=-1   
        if selectcell==25 and metadict['name']=='heat8':
            metadict['Cpregions_glassind']=-1
    savebool=True
    if savebool:
        f=h5py.File(p, mode='r+')
        if 'calbycellmetadata' in f:
            g=f['calbycellmetadata']
        else:
            g=f.create_group('calbycellmetadata')
        if `selectcell` in g:
            cg=g[`selectcell`]
        else:
            cg=g.create_group(`selectcell`)
        for metadict in metadictlist:
            if selectcell==1 and metadict['name'].startswith('heat1'):
                continue
            if selectcell==6 and metadict['name']=='heat6a':
                continue
            if selectcell==16 and (metadict['name']=='heat2' or metadict['name']=='heat8'):
                continue
            if selectcell==19 and metadict['name'].startswith('heat1'):
                continue
            if selectcell==20 and metadict['name']=='heat4a':
                continue
            if selectcell==22 and metadict['name']=='heat4b':
                continue
            if selectcell==24 and metadict['name']=='heat6a':
                continue
            if selectcell==25 and (metadict['name']=='heat2' or metadict['name']=='heat3'):
                continue
            if metadict['name'] in cg:
                mg=cg[metadict['name']]
            else:
                mg=cg.create_group(metadict['name'])
            for k, v in metadict.iteritems():
                mg.attrs[k]=v
        f.close()
        
    plotcpregions=True
    if plotcpregions:
        #pylab.figure(figsize=(8, 6))
        for i, (metadict, hpsdl) in enumerate(zip(metadictlist, allsegdict)):
            if not 'Cpregions_enthalpy' in metadict.keys():
                continue
            T=hpsdl[metadict['heatseg']]['sampletemperature'][cycleindex]
            C=hpsdl[metadict['heatseg']]['sampleheatcapacity'][cycleindex]
            pylab.plot(T, C)
            pv=numpy.max(hpsdl[metadict['heatseg']]['sampleheatcapacity'][cycleindex])
            nv=numpy.min(hpsdl[metadict['heatseg']]['sampleheatcapacity'][cycleindex])
            lgen=lambda s:((s>0) and ([pv, pv], ) or ([nv, nv], ))[0]
            
            rxnindlist=[metadict['Cpregions_glassind'], metadict['Cpregions_xtalind'], metadict['Cpregions_meltind'], metadict['Cpregions_melt2ind']]
            for regind, (enth, Tp, Cp, Tmean, i, j) in enumerate(zip(metadict['Cpregions_enthalpy'], metadict['Cpregions_T_Cmax'], metadict['Cpregions_Cmax'], metadict['Cpregions_Tweightedmean'], metadict['Cpregions_cycindstart'], metadict['Cpregions_cycindstop'])):
                if regind in rxnindlist:
                    col=['b', 'g', 'r', 'r'][rxnindlist.index(regind)]
                else:
                    col=(.4, .4, .4)
                #pylab.plot([T[i], T[j]], lgen(numpy.sign(Cp)), 'k--')
                pylab.fill(T[i:j], C[i:j], color=col)
                pylab.plot(Tp, Cp, 'kx')
                pylab.plot(Tmean, 0, 'k*')
            pylab.xlim(plotTlim)
            Cinrange= C[(T>plotTlim[0])&(T<plotTlim[1])]
            a=Cinrange.min()
            b=Cinrange.max()
            pylab.ylim(a-0.02*(b-a), b+0.02*(b-a))
            pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell),'Cpregions_cell%02d_%s.png' %(selectcell, metadict['name'])))
            pylab.clf()




#orderarray=numpy.abs(numpy.array(cool400))

cols=['k', 'b', 'g', 'r', 'c', 'm']


### plotting series of heat ramps
#mult=1.e6
#nplots=len(orderarray)
#axl=[pylab.subplot(nplots, 1, nplots)]
#for i in range(1, nplots):
#    #ax=pylab.subplot2grid((n, 3), (n-1-i, 0), colspan=2, sharex=axl[0], sharey=axl[0])
#    ax=pylab.subplot(nplots, 1, nplots-i, sharex=axl[0], sharey=axl[0])
#    pylab.setp(ax.get_xticklabels(), visible=False)
#    axl+=[ax]

#for count, i in  enumerate(numpy.argsort(orderarray)):
#    hi, hseg, ci, cseg=heati_heatseg_prevcooli_prevcoolseg[i]
#    print hi, hseg, heatlist[hi], allsegdict[hi][hseg].keys()
#    axl[count].plot(allsegdict[hi][hseg]['sampletemperature'][cycleindex], mult*allsegdict[hi][hseg]['sampleheatcapacity'][cycleindex], cols[count]+'.', markersize=1, label=heatlist[i])

#for ax in axl:
#    ax.set_ylim(-2.1, 4.9)
#    ax.set_yticks([-2, 0, 2, 4])
#
#axl[2].set_ylabel(r'Heat Capacity ($\mu$J/K),  endothermic ->', fontsize=14)
#axl[0].set_xlabel('Temperature (C)', fontsize=14)
#pylab.subplots_adjust(right=.95, top=0.95, hspace=0.01)

###plot cooling rates
#pylab.figure(figsize=(1.5, 8))
#for count, x in enumerate(numpy.sort(orderarray)):
#    pylab.semilogx(numpy.abs(x), count, cols[count]+'o')
#make_ticklabels_invisible(pylab.gca(), x=False)
#pylab.xlabel('cooling rate at 400C (K/s)', fontsize=14)
#pylab.ylim(-.5, count+.5)



###extra stuff?
#pylab.ylim(t3, t4)
#        pylab.xlabel('T (C)')
#        pylab.ylabel('P / dT/dt')
#        pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
#    pylab.subplots_adjust(left=.1, right=.97, top=.93, bottom=.08, wspace=.25, hspace=.25)
##    pylab.show()
##    idialog=messageDialog(title='continue')
##    if not idialog.exec_():
##        break
##    break
#    pylab.savefig(os.path.join(savef,'SCcellplot_cell%02d' %(cellcount+1)+'.png'))
#    pylab.clf()


#pylab.show()



print 'done'
