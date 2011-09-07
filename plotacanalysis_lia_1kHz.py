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


def lia_ampphaseTEST(x, ptspercyc, ncyclewin=1., returnphase=True, pad=True, phaseshift=0.):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=x*sinarr(ptspercyc, x, ph=phaseshift)
    c=x*sinarr(ptspercyc, x, ph=numpy.pi/2.+phaseshift)
    amp=(numpy.array([(numpy.abs(numpy.fft.fft(s[i:i+npts])[0]))**2+(numpy.abs(numpy.fft.fft(c[i:i+npts])[0]))**2 for i in numpy.arange(len(x)-npts)])**.5)*2./npts
    if returnphase:
        phase=numpy.array([numpy.arctan(numpy.abs(numpy.fft.fft(s[i:i+npts])[0])/numpy.abs(numpy.fft.fft(c[i:i+npts])[0])) for i in numpy.arange(len(x)-npts)])
    if pad:
        amp=numpy.concatenate([amp[:npts//2], amp, amp[-1*(len(x)-len(amp)-npts//2):]])
        if returnphase:
            phase=numpy.concatenate([phase[:npts//2], phase, phase[-1*(len(x)-len(phase)-npts//2):]])
    if returnphase:
        return amp, phase
    return amp
    
ptspercyc=300.
n1wcyc=4
if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')
    z=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplevoltage'][0,:]
    la=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_voltage'][0,:,:,0]
    lp=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_voltage'][0,:,:,1]
    f.close()

if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')
    z=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplefilteredvoltage'][0,:]
    la=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_filteredvoltage'][0,:,:,0]
    lp=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_filteredvoltage'][0,:,:,1]
    f.close()
    
if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')
    z=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplecurrent'][0,:]
    la=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_current'][0,:,:,0]
    lp=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_current'][0,:,:,1]
    f.close()
    
if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')
    z=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplevoltage'][0,:]/1000.
    f.close()
    hlist=[1, 2, 3]
    ans=numpy.empty((len(z), len(hlist), 2), dtype='float32')
    for j, h in enumerate(hlist):
        ans[:, j, 0], ans[:, j, 1]=lia_ampphaseTEST(z, ptspercyc/h, ncyclewin=n1wcyc*h, phaseshift=0.)
    la=ans[:, :, 0]
    lp=ans[:, :, 1]

if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')
    a=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplevoltage'][0,:]/1000.
    b=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplecurrent'][0,:]/1000.
    f.close()
    z=a/b
    hlist=[1, 2, 3]
    ans=numpy.empty((len(z), len(hlist), 2), dtype='float32')
    for j, h in enumerate(hlist):
        ans[:, j, 0], ans[:, j, 1]=lia_ampphase(z, ptspercyc/h, ncyclewin=n1wcyc*h, phaseshift=0.)
    la=ans[:, :, 0]
    lp=ans[:, :, 1]
    
if 1:
    z=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake_voltage.npy')
    hlist=[1, 2, 3]
    ans=numpy.empty((len(z), len(hlist), 2), dtype='float32')
    for j, h in enumerate(hlist):
        ans[:, j, 0], ans[:, j, 1]=lia_ampphase(z, ptspercyc/h, ncyclewin=n1wcyc*h, phaseshift=0.)
    la=ans[:, :, 0]
    lp=ans[:, :, 1]


if 0:
    z=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake_current.npy')
    hlist=[1, 2, 3]
    ans=numpy.empty((len(z), len(hlist), 2), dtype='float32')
    for j, h in enumerate(hlist):
        ans[:, j, 0], ans[:, j, 1]=lia_ampphase(z, ptspercyc/h, ncyclewin=n1wcyc*h, phaseshift=0.)
    la=ans[:, :, 0]
    lp=ans[:, :, 1]
    
if 0:
    z=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake2_smoothrcurrent.npy')
    hlist=[1, 2, 3]
    ans=numpy.empty((len(z), len(hlist), 2), dtype='float32')
    for j, h in enumerate(hlist):
        ans[:, j, 0], ans[:, j, 1]=lia_ampphase(z, ptspercyc/h, ncyclewin=n1wcyc*h, phaseshift=0.)
    la=ans[:, :, 0]
    lp=ans[:, :, 1]

def sinarr(nptspercycle, npts, ph=0.):
    if isinstance(npts, numpy.ndarray):
        npts=len(npts)
    return numpy.sin(numpy.arange(npts)*2.*numpy.pi/nptspercycle+ph)
    
if 0:
    z=3.*sinarr(ptspercyc, 500., .02)
    hlist=[1, 2, 3]
    ans=numpy.empty((len(z), len(hlist), 2), dtype='float32')
    for j, h in enumerate(hlist):
        ans[:, j, 0], ans[:, j, 1]=lia_ampphase(z, ptspercyc/h, ncyclewin=n1wcyc*h, phaseshift=0.)
    la=ans[:, :, 0]
    lp=ans[:, :, 1]
if 1:
    for j in range(3):
        for i in range(2):
            ax1=pylab.subplot(3, 2, j*2+1+i)
            pylab.plot(la[:, j], 'b')
            pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
            ax2 = ax1.twinx()
            ax2.plot(lp[:, j]*180./numpy.pi, 'r', alpha=.4)
            for tl in ax2.get_yticklabels():
                tl.set_color('r')

if 1:
    pylab.show()
