import numpy, h5py, pylab, copy
from PnSC_h5io import *
from PnSC_math import *

from matplotlib.ticker import FuncFormatter

def myexpformat(x, pos):
    for ndigs in range(5):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

ptspercyc=30.
n1wcyc=4
find_1w=4

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
f=h5py.File(p,mode='r')
dct=f['Calorimetry/Sn_10kHz_4e4Ks_meta/analysis/dc/sampletemperature'][0,:]

dcppr=f['Calorimetry/Sn_10kHz_4e4Ks_meta/analysis/dc/samplepowerperrate'][0,:]

act=f['Calorimetry/Sn_10kHz_4e4Ks_meta/analysis/ac1w/sampletemperature'][0,:]

acppr=f['Calorimetry/Sn_10kHz_4e4Ks_meta/analysis/ac1w/samplepowerperrate'][0,:]

dchpsd=CreateHeatProgSegDictList(p, 'Sn_10kHz_4e4Ks_meta', 'dc')[0]

achpsd=CreateHeatProgSegDictList(p, 'Sn_10kHz_4e4Ks_meta', 'ac1w')[0]
ac1wlia=achpsd['samplevoltage'][0]
time_s=achpsd['cycletime'][0]

hpsdl=CreateHeatProgSegDictList(p, 'Sn_10kHz_4e4Ks', 'cell7_57.75dc56.1ac_10kHz_12.6ms_1_of_1')
v=hpsdl[3]['samplevoltage'][0]

f.close()

if 0:
    pylab.figure()
    pylab.plot(dct,dcppr*1.e6,'k-',label='dc measurement')
    pylab.plot(act,acppr*1.e6,'r-',label='ac measurement')
    pylab.xlabel('Temperature (C)',fontsize=16)
    pylab.ylabel('power/heat rate ($\mu$J/K)',fontsize=16)

    t1=110
    t2=500
    pylab.xlim(t1, t2)
    i1=numpy.argmin((act-t1)**2)
    i2=numpy.argmin((act-t2)**2)

    pylab.ylim(1,2.6)

    leg=pylab.legend(loc=1)
    temp=[t.set_fontsize(14) for t in leg.texts]

    pylab.figure()
    pylab.plot(time_s[i1:i2]*1000., v[i1:i2], 'b-')
    pylab.xlim(time_s[i1]*1000., time_s[i2]*1000.)
    pylab.xlabel('time (ms)',fontsize=16)
    pylab.ylabel('device voltage (V)',fontsize=16)
    

if 1:
    skip=100
    skipe=150
    count=0
    for k in dchpsd.keys():
        if k in achpsd.keys() and isinstance(dchpsd[k], numpy.ndarray):
            x=dchpsd[k][0, skip:-1*skipe]
            y=achpsd[k][0, skip:-1*skipe]
            count+=1
            ax1=pylab.subplot(4, 2, count)
            pylab.plot(x, 'b', label='dc')
            pylab.plot(y, 'g', label='ac')
            ax2=ax1.twinx()
            ax2.plot((y-x)/x, 'r', label='(ac-dc)/dc')
            pylab.title(k)
    ax1.legend(loc=2)
    ax2.legend(loc=1)
pylab.show()
