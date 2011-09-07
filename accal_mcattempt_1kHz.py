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

ptspercyc=300.
n1wcyc=4
find_1w=4

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
f=h5py.File(p,mode='r')
v=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplevoltage'][0,:]/1000.
liava=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_voltage'][0,:,:,0]
liavp=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_voltage'][0,:,:,1]

liaca=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_current'][0,:,:,0]
liacp=f['Calorimetry/1kHzcrop/analysis/again/LIAharmonics_current'][0,:,:,1]

fftva=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_voltage'][0,:,:,0]
fftca=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_current'][0,:,:,0]
f.close()
vsmoothdc=savgolsmooth(fftva[:, 0], nptsoneside=300, order = 1)
csmoothdc=savgolsmooth(fftca[:, 0], nptsoneside=300, order = 1)

r=vsmoothdc/csmoothdc
    
v1wsmooth=savgolsmooth(liava[:, 0], nptsoneside=300, order = 1)
c1wsmooth=savgolsmooth(liaca[:, 0], nptsoneside=300, order = 1)
v2wsmooth=savgolsmooth(liava[:, 1], nptsoneside=300, order = 1)

w=1.e3
alpha=1.64e-3
k=alpha*r
mc=c1wsmooth**2*csmoothdc*r*k/v2wsmooth/w
mcmol=mc/3.e-8


for i, (data, nam) in enumerate(zip([vsmoothdc, csmoothdc, v1wsmooth, c1wsmooth, v2wsmooth,r,  mc, mcmol], ['vsmoothdc', 'csmoothdc', 'v1wsmooth', 'c1wsmooth', 'v2wsmooth', 'R','$i^2 I R k / V_{2\omega}\omega$','mc (J/K/mol)'])):
    pylab.subplot(3, 3, i+1)
    pylab.plot(data)
    pylab.title(nam)

if 1:
    pylab.show()
