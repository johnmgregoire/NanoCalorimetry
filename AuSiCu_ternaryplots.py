import os
os.chdir('C:/Users/JohnnyG/Documents/PythonCode/ternaryplot')
from myternaryutility import TernaryPlot
import matplotlib.cm as cm
import numpy, operator
import pylab, h5py
from matplotlib.ticker import FuncFormatter

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

savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'


f=h5py.File(p, mode='r')
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
comp=f['CompThick/atfrac'][:, :]
nm=f['CompThick/nm'][:]


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

#following 7 lines are common to several belwo subroutines
z=[[] for cv in comp]
y=[[] for cv in comp]
x=[[] for cv in comp]
for i, (xv, yv, zv, cv, nmv, mdl, xdl) in enumerate(zip(x, y, z, comp, nm, metadictlists, xrddictlists)):
    for md in mdl:
        pcal=md['prevname'][:5]
        cal=md['name'][:5]
        if pcal!=cal:#if previous scan was in same heat# as scan then there was no xrd in between
        
        
##amfrac for a range of cool rate
#            if 'prevcoolrate_400C' in md.keys() and numpy.abs(md['prevcoolrate_400C'])>2.05e4 and numpy.abs(md['prevcoolrate_400C'])<2.25e4:
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[xd['amfrac']]
#                        yv+=[md['prevcoolrate_400C']]
#                        break
##if multiple qualifying values take the mean
#data=[(yv, numpy.mean(zv), cv) for yv, zv, cv in zip(y, z, comp) if len(zv)>0]
#yal=numpy.float32(map(operator.itemgetter(0), data))
#za=numpy.float32(map(operator.itemgetter(1), data))
#ca=map(operator.itemgetter(2), data)
#title='amorhpous fraction after ~2.1 10$^4$ K/s quench'

##Tg enthalpy per amfrac
#            if 'prevcoolrate_400C' in md.keys() and 'Cpregions_enthalpy' in md.keys() and numpy.abs(md['prevcoolrate_400C'])>1.3e4:
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[xd['amfrac']]
#                        yv+=[(md['Cpregions_glassind']<0 and (0., ) or (md['Cpregions_enthalpy'][md['Cpregions_glassind']], ))[0]]
#                        print 'cell ', i+1
#                        break
##take biggest Tg enthalpy
#data=[(cell, yv, zv, numpy.max(yv)/zv[numpy.argmax(yv)], cv) for cell, yv, zv, cv in zip(range(1,26),y, z, comp)) if len(zv)>0]
#cells=numpy.array(map(operator.itemgetter(0), data))
#amf=map(operator.itemgetter(1), data)
#Hg=map(operator.itemgetter(2), data)
#za=numpy.float32(map(operator.itemgetter(3), data))
#ca=map(operator.itemgetter(4), data)
#
#za*=1.e6
#title=r'Glass transition enthalpy per amorphous fraction ($\mu$J)'
#
##za/=(3.8*.8*nm[cells]*1.e-15) #puts it in inuts of J/m3
##za*=1.e-6 #J/cm3
##title=r'volumetric glass transition enthalpy (J/cm$^3$)'


##Tg weighted avergae temperature
#            if 'prevcoolrate_400C' in md.keys() and 'Cpregions_glassind' in md.keys() and md['Cpregions_glassind']>=0:
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[md['Cpregions_Tweightedmean'][md['Cpregions_glassind']]]
#                        yv+=[md['prevcoolrate_400C']]
#                        xv+=[md['Cpregions_enthalpy'][md['Cpregions_glassind']]]
#                        print 'cell ', i+1
#                        break
##Tg for fastest cool rate
##data=[(cell, xv,yv, zv, zv[numpy.argmax(yv)], cv) for cell,xv,yv, zv, cv in zip(range(1,26),x,y, z, comp) if len(zv)>0]
##Tg avergaed over multiple scans
#data=[(cell, xv, yv, zv, (numpy.array(zv)*numpy.array(xv)).sum()/numpy.array(xv).sum(), cv) for cell, xv, yv, zv, cv in zip(range(1,26),x, y, z, comp) if len(zv)>0]
#cells=numpy.array(map(operator.itemgetter(0), data))
#enth=map(operator.itemgetter(1), data)
#coolrate=map(operator.itemgetter(2), data)
#Tgave=map(operator.itemgetter(3), data)
#za=numpy.float32(map(operator.itemgetter(4), data))
#ca=map(operator.itemgetter(5), data)
#title=r'Glass transition characteristic Temp. (C)'

##Tx weighted avergae temperature
#            if 'prevcoolrate_400C' in md.keys() and 'Cpregions_xtalind' in md.keys() and md['Cpregions_xtalind']>=0:
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[md['Cpregions_Tweightedmean'][md['Cpregions_xtalind']]]
#                        yv+=[md['prevcoolrate_400C']]
#                        xv+=[md['Cpregions_enthalpy'][md['Cpregions_xtalind']]]
#                        print 'cell ', i+1
#                        break
##Tg for fastest cool rate
##data=[(cell, xv,yv, zv, zv[numpy.argmax(yv)], cv) for cell,xv,yv, zv, cv in zip(range(1,26),x,y, z, comp) if len(zv)>0]
##Tg avergaed over multiple scans
#data=[(cell, xv, yv, zv, (numpy.array(zv)*numpy.array(xv)).sum()/numpy.array(xv).sum(), cv) for cell, xv, yv, zv, cv in zip(range(1,26),x, y, z, comp) if len(zv)>0]
#cells=numpy.array(map(operator.itemgetter(0), data))
#enth=map(operator.itemgetter(1), data)
#coolrate=map(operator.itemgetter(2), data)
#Tgave=map(operator.itemgetter(3), data)
#za=numpy.float32(map(operator.itemgetter(4), data))
#ca=map(operator.itemgetter(5), data)
#title=r'Crystallization characteristic Temp. (C)'

##Tm weighted avergae temperature
#            if 'prevcoolrate_400C' in md.keys() and 'Cpregions_meltind' in md.keys() and md['Cpregions_meltind']>=0:
#                #no xrd required
#                zv+=[md['Cpregions_Tweightedmean'][md['Cpregions_meltind']]]
#                yv+=[md['prevcoolrate_400C']]
#                xv+=[md['Cpregions_enthalpy'][md['Cpregions_meltind']]]
#                print 'cell ', i+1
#                break
##Tg for fastest cool rate
##data=[(cell, xv,yv, zv, zv[numpy.argmax(yv)], cv) for cell,xv,yv, zv, cv in zip(range(1,26),x,y, z, comp) if len(zv)>0]
##Tg avergaed over multiple scans
#data=[(cell, xv, yv, zv, (numpy.array(zv)*numpy.array(xv)).sum()/numpy.array(xv).sum(), cv) for cell, xv, yv, zv, cv in zip(range(1,26),x, y, z, comp) if len(zv)>0]
#cells=numpy.array(map(operator.itemgetter(0), data))
#enth=map(operator.itemgetter(1), data)
#coolrate=map(operator.itemgetter(2), data)
#Tgave=map(operator.itemgetter(3), data)
#za=numpy.float32(map(operator.itemgetter(4), data))
#ca=map(operator.itemgetter(5), data)
#title=r'Melting characteristic Temp. (C)'


##melt enthalpy - xtal enthalpy / fcc, including when xtal enthalpy is zero
#            if 'prevcoolrate_400C' in md.keys() and 'Cpregions_xtalind' in md.keys() and 'Cpregions_meltind' in md.keys() and md['Cpregions_meltind']>=0:
#                for xd in xdl:
#                    if xd['name']==pcal:#use xrd that happened after the prev scan
#                        zv+=[md['Cpregions_enthalpy'][md['Cpregions_meltind']]]
#                        yv+=[xd['fccfrac']]
#                        xv+=[(md['Cpregions_xtalind']<0 and (0, ) or (md['Cpregions_enthalpy'][md['Cpregions_xtalind']], ))[0]]
#                        print 'cell ', i+1
#                        break
##Tg for fastest cool rate
##data=[(cell, xv,yv, zv, zv[numpy.argmax(yv)], cv) for cell,xv,yv, zv, cv in zip(range(1,26),x,y, z, comp) if len(zv)>0]
##Tg avergaed over multiple scans
#data=[(cell, xv, yv, zv, ((numpy.array(zv)-numpy.array(xv))/numpy.array(yv)*numpy.array(yv)).sum()/numpy.array(yv).sum(), cv) for cell, xv, yv, zv, cv in zip(range(1,26),x, y, z, comp) if len(zv)>0]
#cells=numpy.array(map(operator.itemgetter(0), data))
#enth=map(operator.itemgetter(1), data)
#coolrate=map(operator.itemgetter(2), data)
#Tgave=map(operator.itemgetter(3), data)
#za=numpy.float32(map(operator.itemgetter(4), data))
#ca=map(operator.itemgetter(5), data)
#title=r'melt enthalpy - xtal enthalpy // fcc frac'



#thickness for cells with cal and xrd data
            if 'prevcoolrate_400C' in md.keys() and 'Cpregions_xtalind' in md.keys() and 'Cpregions_meltind' in md.keys() and md['Cpregions_meltind']>=0:
                for xd in xdl:
                    if xd['name']==pcal:#use xrd that happened after the prev scan
                        zv+=[nmv]
                        yv+=[0.]
                        xv+=[0.]
                        print 'cell ', i+1
                        break
        if len(zv)>0.:
            break
#Tg for fastest cool rate
#data=[(cell, xv,yv, zv, zv[numpy.argmax(yv)], cv) for cell,xv,yv, zv, cv in zip(range(1,26),x,y, z, comp) if len(zv)>0]
#Tg avergaed over multiple scans
data=[(cell, xv, yv, zv,0., cv) for cell, xv, yv, zv, cv in zip(range(1,26),x, y, z, comp) if len(zv)>0]
cells=numpy.array(map(operator.itemgetter(0), data))
enth=map(operator.itemgetter(1), data)
coolrate=map(operator.itemgetter(2), data)
za=numpy.array(map(operator.itemgetter(3), data))
garb=numpy.float32(map(operator.itemgetter(4), data))
ca=numpy.array(map(operator.itemgetter(5), data))
title=r'Film thickness (nm)'

#********start ternary plotting
pylab.figure(figsize=(10, 6))
ax=pylab.subplot(111)
##stp = TernaryPlot(ax, ellabels=['Au', 'Si', 'Cu']) 
stp = TernaryPlot(ax, ellabels=['Au', 'Si', 'Cu'], minlist=[0.47, 0.12, .2])
stp.grid(nintervals=3)
stp.label(fontsize=14)




#plot za value from above by color
stp.scatter(ca, s=80, c=za, label='_', cmap=cm.jet, alpha=1., marker='o') #scatter doesn't work well with legend, so label="_" (hides it); cmap chooses the color scheme; alpha allows some transparency to see overlapping data
stp.colorbar(title, fontsize=14)
#stp.ax.set_title(title, fontsize=14)
#****

##plot cellnumbers
#celllist=[1, 2]+range(4, 21)+[22]+[24, 25]
#for selectcell in celllist:
#    stp.text(comp[selectcell-1], `selectcell`, ha='center', va='center', color='r', fontsize=12)

##plot points by color
#celllist=[1, 2]+range(4, 21)+[22]+[24, 25]
#minlist=[c.min() for c in comp[numpy.array(celllist)-1, :].T]
#rangelist=numpy.float32([[m, 1.-numpy.concatenate([minlist[:i], minlist[i+1:]]).sum()] for i, m in enumerate(minlist)])
#colors=stp.color_comp_calc(comp[numpy.array(celllist)-1, :], rangelist=rangelist)
#stp.colorcompplot(comp[numpy.array(celllist)-1, :], '.', colors=colors, markersize=20)
#stp.line([0.6123922371385736, 0.14728011617695427, 0.24032764668447215], [0.48297887710346687, 0.21946269934149762, 0.29755842355503548], fmt='k-')
#
#pylab.subplots_adjust(bottom=0., top=1., left=0.)#, right=1.)


pylab.show()
print 'done'
