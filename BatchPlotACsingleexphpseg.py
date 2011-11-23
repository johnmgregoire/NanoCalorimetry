import os
from PnSC_h5io import *

import time, copy
import os
import sys
import numpy
import h5py
import pylab
from matplotlib.ticker import FuncFormatter
import scipy.integrate
from scipy.interpolate import griddata
import matplotlib.cm as cm

from matplotlib.ticker import FuncFormatter

def myexpformat(x, pos):
    for ndigs in range(5):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        print lab, eval(lab), x
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

#enter the h5 path here and note that in Python you should only use "/" and not "\"
p='E:/pnscpythntraining/trainingcpy.h5'

#select the segment index that you wanted ploted - only one segment at a time
seg=4

#select the experiment group and heat program
exp='AC'
hp='cell24_20.4to71.8to12.2dc_20.4to71.8to12.2ac_120ppc_8s8s_again2_1_of_1'

#enter the time range for the plots in seconds
tlims=(0, 7.8)

#enter the plotting interval - the program will plot a subset of the data with this interval. If there are >100,000 data points you should probably use an interval to make plotting faster. (must be integer)
interv=6

#choose the temperatures in C for the tik labels on the temperature axis
Tmks=[100, 200, 300, 400, 500, 600]


f=h5py.File(p, mode='r')
cyc=0
d=CreateHeatProgSegDictList(p, exp, hp, expandmultdim=False)[seg]

for k in d.keys():
    if k!='cycletime' and isinstance(d[k], numpy.ndarray) and d[k].shape[:2]==d['cycletime'].shape:
        d[k]=d[k][cyc][::interv]
k='cycletime'
d[k]=d[k][cyc][::interv]
    
t=d['cycletime']
t-=t[0]
T=d['sampletemperature']


Tmklabs=['%d' %x for x in Tmks]
#the first plot is temperature vs time along with a polynomial fit for caclulatingf the temperature tick mark positions. if the fit doesn't look good try changing  these parameters.
fitinterv=200
Tfitpolyorder=4

t_Tfit=numpy.polyfit(T[::fitinterv], t[::fitinterv], Tfitpolyorder)
tmks=[numpy.polyval(t_Tfit, x) for x in Tmks]

ms=1
subadjdict={'left':.2, 'bottom':.17, 'right':.92, 'top':.85}
figsize=(5, 4)
if 1:
    pylab.figure(figsize=figsize)
    pylab.plot(t, T, 'b.', ms=ms)
    pylab.xlabel('elapsed time (s)')
    pylab.ylabel('Temperature (C)')
    pylab.plot(tmks, Tmks, 'ro')
    pylab.xlim(tlims)
    pylab.subplots_adjust(**subadjdict)
    if 0:
        pylab.show()    


def plot1(x=t, yk='sampletemperature',fmt='b.', ms=ms, mult=1., xlab='elapsed time (s)', ylab='Temperature (C)', \
          xlims=tlims, Ttopmks=(tmks, Tmklabs), subadjdict=subadjdict, figsize=figsize):

    pylab.figure(figsize=figsize)
    if isinstance(yk, str):
        pylab.plot(x, d[yk]*mult, fmt, ms=ms)
    else:
        pylab.plot(x, yk*mult, fmt, ms=ms)
    pylab.xlabel(xlab)
    pylab.ylabel(ylab)
    pylab.xlim(xlims)
    if not Ttopmks is None:
        ax2=pylab.twiny()
        ax2.set_xticks(Ttopmks[0])
        ax2.set_xticklabels(Ttopmks[1])
        ax2.set_xlabel('Temperature (C)')
        ax2.set_xlim(xlims)
    pylab.subplots_adjust(**subadjdict)

def plotamp(yk, freqind, **kwargs):
    plot1(yk=(d[yk][:, freqind, :]**2).sum(axis=1), **kwargs)
    

plot1(Ttopmks=None)
plot1(yk='sampleheatrate', ylab='Heat rate (K/s)')
pylab.show()
plot1(yk='samplepowerperrate', ylab='Power / dT/dt ($\mu$J/K)', mult=1.e6)
plotamp('WinFFT_voltage', 3, ylab='1$\omega$ amplitude (V)', mult=1.)
plotamp('WinFFT_voltage', 6, ylab='2$\omega$ amplitude (mV)', mult=1.e3)
plotamp('WinFFT_filteredvoltage', 6, ylab='filtered 2$\omega$ amplitude (mV)', mult=1.e3)
plot1(yk='acheatcapacity', ylab='mC$_p$, 2$\omega$ method ($\mu$J/K)', mult=1.e6)
plot1(yk='acheatcapacity_1', ylab='mC$_p$, filtered 2$\omega$ method ($\mu$J/K)', mult=1.e6)

pylab.show()    
