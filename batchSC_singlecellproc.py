import time
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

#p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'
p=mm.h5path
#f=h5py.File(p, mode='r+')
#f=h5py.File(p, mode='r')
savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_May17'
allsegdict=[]
selectcell=11
#                            0                       1                      2           3                           4                 5
heatlist=['heat1a', 'heat2', 'heat3', 'heat4a', 'heat7', 'heat8']
for exp in heatlist:
    f, hppaths=experimenthppaths(p, exp)
    f.close()
    for hpp in hppaths:
        h5hpname=hpp.rpartition('/')[2]
        f, g=gethpgroup(p, exp, h5hpname=h5hpname)
        cell=g.attrs['CELLNUMBER']
        f.close()
        if cell!=selectcell:
            continue
        
        hpsdl=CreateHeatProgSegDictList(p, exp, h5hpname)
        allsegdict+=[hpsdl]


def heatrate_T(d, T, Twin=5.):
    #i=numpy.argmin((T-d['sampletemperature'])**2)
    Ta=d['sampletemperature'][0]
    x=numpy.where((Ta>=T-Twin)&(Ta<=T+Twin))[0]
    prev=numpy.array([not (t-1 in x) for t in x])
    previ=numpy.where(prev)[0]
    stopi=numpy.append(previ[1:],len(x))
    longestbunchind=numpy.argmax(stopi-previ)
    inds=x[previ[longestbunchind]:stopi[longestbunchind]]
    return d['sampleheatrate'][0][inds].mean()
    

cool180=[]
cool400=[]
heat200500=[]
heati_heatseg_prevcooli_prevcoolseg=[(0, 2, -1, -1), (1, 2, 0, 4), (2, 2, 1, 4), (3, 2, 2, 4), (5, 2, 4, 3)]
for hi, hseg, ci, cseg in heati_heatseg_prevcooli_prevcoolseg:
    if ci<0:
        cool180+=[-1.e2]
        cool400+=[-1.e2]
    else:
        cool180+=[heatrate_T(allsegdict[ci][cseg], 180.)]
        cool400+=[heatrate_T(allsegdict[ci][cseg], 400.)]
    heat200500+=[heatrate_T(allsegdict[hi][hseg], 350., Twin=150.)]

orderarray=numpy.abs(numpy.array(cool400))

cols=['k', 'b', 'g', 'r', 'c', 'm']



mult=1.e6
nplots=len(orderarray)
axl=[pylab.subplot(nplots, 1, nplots)]
for i in range(1, nplots):
    #ax=pylab.subplot2grid((n, 3), (n-1-i, 0), colspan=2, sharex=axl[0], sharey=axl[0])
    ax=pylab.subplot(nplots, 1, nplots-i, sharex=axl[0], sharey=axl[0])
    pylab.setp(ax.get_xticklabels(), visible=False)
    axl+=[ax]

for count, i in  enumerate(numpy.argsort(orderarray)):
    hi, hseg, ci, cseg=heati_heatseg_prevcooli_prevcoolseg[i]
    print hi, hseg, heatlist[hi], allsegdict[hi][hseg].keys()
    axl[count].plot(allsegdict[hi][hseg]['sampletemperature'][0], mult*allsegdict[hi][hseg]['sampleheatcapacity'][0], cols[count]+'.', markersize=1, label=heatlist[i])


for ax in axl:
    ax.set_ylim(-2.1, 4.9)
    ax.set_yticks([-2, 0, 2, 4])

axl[2].set_ylabel(r'Heat Capacity ($\mu$J/K),  endothermic ->', fontsize=14)
axl[0].set_xlabel('Temperature (C)', fontsize=14)
pylab.subplots_adjust(right=.95, top=0.95, hspace=0.01)

###plot cooling rates
#pylab.figure(figsize=(1.5, 8))
#for count, x in enumerate(numpy.sort(orderarray)):
#    pylab.semilogx(numpy.abs(x), count, cols[count]+'o')
#make_ticklabels_invisible(pylab.gca(), x=False)
#pylab.xlabel('cooling rate at 400C (K/s)', fontsize=14)
#pylab.ylim(-.5, count+.5)
#pylab.show()


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
