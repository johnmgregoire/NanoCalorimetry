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

selectcell=1
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
    
    
    
nskip=100

cycleindex=0

savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'

cols=['k', 'b', 'g', 'r', 'c', 'm', 'y']*3

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

axl=get25pylabaxes(horizslowaxis=True)
taxl=[ax.twinx() for ax in axl]

metadictlist=[]
allsegdict=[]
for count, ed in enumerate(expdicts):
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
        hpsdl=CreateHeatProgSegDictList(p, exp, h5hpname)
        saveh5hpname=h5hpname
        if not hpsdl is None:
            ed['cell']=cell
            T=hpsdl[ed['heatseg']]['sampletemperature'][cycleindex][nskip:-1*nskip]
            maxT=T.max()
            T=numpy.concatenate([T, hpsdl[ed['coolseg']]['sampletemperature'][cycleindex][nskip:-1*nskip]])
            i=0
            dt=[]
            while len(dt)<2:
                try:
                    dt=hpsdl[i]['cycletime'][cycleindex][0:2]
                    i+=1
                except:
                    dt=[]
                    break
            if len(dt)<0:
                continue
            dt=dt[1]-dt[0]
            timeabove=len(numpy.where(T>800.)[0])*dt
            axl[cell-1].plot(count, maxT, '%sD' %cols[count])
            taxl[cell-1].plot(count, timeabove, '%sx' %cols[count])
            axl[cell-1].set_ylabel('cell %d, max Temp (C)' %cell)
            taxl[cell-1].set_ylabel('Time above 800C (s)(X)')
            ed['maxT']=maxT
            ed['secabove800C']=timeabove
            metadictlist+=[copy.deepcopy(ed)]
for ax in axl:
    ax.set_xlim(-0.5, len(expdicts)-0.5)
savebool=False
if savebool:
    f=h5py.File(p, mode='r+')
    g=f['calbycellmetadata']
    f.close()
#    if `selectcell` in g:
#        cg=g[`selectcell`]
#    else:
#        cg=g.create_group(`selectcell`)
#    for metadict in metadictlist:
#        if selectcell==1 and metadict['name'].startswith('heat1'):
#            continue
#        if selectcell==6 and metadict['name']=='heat6a':
#            continue
#        if selectcell==16 and (metadict['name']=='heat2' or metadict['name']=='heat8'):
#            continue
#        if selectcell==19 and metadict['name'].startswith('heat1'):
#            continue
#        if selectcell==20 and metadict['name']=='heat4a':
#            continue
#        if selectcell==22 and metadict['name']=='heat4b':
#            continue
#        if selectcell==24 and metadict['name']=='heat6a':
#            continue
#        if selectcell==25 and (metadict['name']=='heat2' or metadict['name']=='heat3'):
#            continue
#        if metadict['name'] in cg:
#            mg=cg[metadict['name']]
#        else:
#            mg=cg.create_group(metadict['name'])
#        for k, v in metadict.iteritems():
#            mg.attrs[k]=v
#    f.close()
#    
#plotcpregions=True
#if plotcpregions:
#    for i, (metadict, hpsdl) in enumerate(zip(metadictlist, allsegdict)):
#        if not 'Cpregions_enthalpy' in metadict.keys():
#            continue
#        T=hpsdl[metadict['heatseg']]['sampletemperature'][cycleindex]
#        C=hpsdl[metadict['heatseg']]['sampleheatcapacity'][cycleindex]
#        pylab.plot(T, C)
#        pv=numpy.max(hpsdl[metadict['heatseg']]['sampleheatcapacity'][cycleindex])
#        nv=numpy.min(hpsdl[metadict['heatseg']]['sampleheatcapacity'][cycleindex])
#        lgen=lambda s:((s>0) and ([pv, pv], ) or ([nv, nv], ))[0]
#        
#        rxnindlist=[metadict['Cpregions_glassind'], metadict['Cpregions_xtalind'], metadict['Cpregions_meltind'], metadict['Cpregions_melt2ind']]
#        for regind, (enth, Tp, Cp, Tmean, i, j) in enumerate(zip(metadict['Cpregions_enthalpy'], metadict['Cpregions_T_Cmax'], metadict['Cpregions_Cmax'], metadict['Cpregions_Tweightedmean'], metadict['Cpregions_cycindstart'], metadict['Cpregions_cycindstop'])):
#            if regind in rxnindlist:
#                col=['b', 'g', 'r', 'r'][rxnindlist.index(regind)]
#            else:
#                col=(.4, .4, .4)
#            #pylab.plot([T[i], T[j]], lgen(numpy.sign(Cp)), 'k--')
#            pylab.fill(T[i:j], C[i:j], color=col)
#            pylab.plot(Tp, Cp, 'kx')
#            pylab.plot(Tmean, 0, 'k*')
#        pylab.xlim(plotTlim)
#        Cinrange= C[(T>plotTlim[0])&(T<plotTlim[1])]
#        a=Cinrange.min()
#        b=Cinrange.max()
#        pylab.ylim(a-0.02*(b-a), b+0.02*(b-a))
#        pylab.savefig(os.path.join(os.path.join(savef, 'cell%02d' %selectcell),'Cpregions_cell%02d_%s.png' %(selectcell, metadict['name'])))
#        pylab.clf()

    


#orderarray=numpy.abs(numpy.array(cool400))



pylab.show()
print 'done'
