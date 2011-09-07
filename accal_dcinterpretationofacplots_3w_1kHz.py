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

fftfudgefactor=2

ptspercyc=300.
n1wcyc=4
find_1w=4
smoothpts=300
p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
f=h5py.File(p,mode='r')
v=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplevoltage'][0,:]/1000.
liava=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_voltage'][0,:,:,0]
liavp=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_voltage'][0,:,:,1]

liaca=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_current'][0,:,:,0]
liacp=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_current'][0,:,:,1]

fftva=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_voltage'][0,:,:,0]
fftca=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_current'][0,:,:,0]

T=f['Calorimetry/1kHzcrop/analysis/dc/sampletemperature'][0,:]
f.close()
vsmoothdc=savgolsmooth(fftva[:, 0], nptsoneside=smoothpts, order = 1)
csmoothdc=savgolsmooth(fftca[:, 0], nptsoneside=smoothpts, order = 1)

r=vsmoothdc/csmoothdc
    
c1wsmooth=savgolsmooth(fftca[:, find_1w], nptsoneside=smoothpts, order = 1)*fftfudgefactor
v3wsmooth=savgolsmooth(fftva[:, find_1w*3], nptsoneside=smoothpts, order = 1)*fftfudgefactor

w=1.e3*2.*numpy.pi
alpha=1.64e-3
k=alpha*r
mc=c1wsmooth**3*r*k/v3wsmooth/w/4.
mcmol=mc/7.5e-8


if 0:
    numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/mcfft_3w.npy', mc)

if 0:
    numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/v3wsmooth.npy', v3wsmooth)

    
if 0:
    fake2fft=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake2fft_smoothrcurrent.npy')[:, :, 0]
    fake2v2wsmooth=savgolsmooth(fake2fft[:, find_1w*2], nptsoneside=smoothpts, order = 1)*fftfudgefactor
    mcsub=c1wsmooth**2*csmoothdc*r*k/(v2wsmooth-fake2v2wsmooth)/w
    
    mcsub/=7.5e-8

for i, (data, nam) in enumerate(zip([vsmoothdc, csmoothdc, c1wsmooth, v3wsmooth,r, mc, mcmol], ['vsmoothdc', 'csmoothdc','c1wsmooth', 'v3wsmooth', 'R','$i^2 I R k / V_{2\omega}\omega$ (J/K)', ' 75nmol, mc (J/K/mol)'])):
    pylab.subplot(3, 3, i+1)
    pylab.plot(data)
    pylab.title(nam)

if 1:
    pylab.show()
