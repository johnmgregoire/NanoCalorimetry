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
dct=f['Calorimetry/1kHzcrop/analysis/dc/sampletemperature'][0,:]

dcppr=f['Calorimetry/1kHzcrop/analysis/dc/samplepowerperrate'][0,:]

act=f['Calorimetry/1kHzcrop/analysis/ac1w/sampletemperature'][0,:]

acppr=f['Calorimetry/1kHzcrop/analysis/ac1w/samplepowerperrate'][0,:]

fftv2w=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_filteredvoltage'][0,:,n1wcyc*2.,0]*2.

dchpsd=CreateHeatProgSegDictList(p, '1kHzcrop', 'dc')[0]

achpsd=CreateHeatProgSegDictList(p, '1kHzcrop', 'ac1w')[0]
ac1wlia=achpsd['samplevoltage'][0]
time_s=achpsd['cycletime'][0]

hpsdl=CreateHeatProgSegDictList(p, '1kHzcrop', 'again')
v=hpsdl[0]['samplevoltage'][0]
fv=hpsdl[0]['samplefilteredvoltage'][0]

f.close()

if 1:
    pylab.figure()
    if 0:
        pylab.plot(dct,dcppr*1.e6,'k-',label='analysis of dc-component')
        pylab.plot(act,acppr*1.e6,'r-',label='analysis of ac 1$\omega$-component')
        pylab.ylabel('power/heat rate ($\mu$J/K)',fontsize=16)
    else:
        pylab.plot(dct,dcppr/7.5e-8,'k-',label='analysis of dc-component', lw=2)
        pylab.plot(act,acppr/7.5e-8,'g-',label='analysis of ac 1$\omega$-component', lw=2)
        pylab.ylabel('$C_{Sn}$+$C_{device}$+d$H_{Sn}$/d$T$+heat loss (J/K/mol)',fontsize=16)
        pylab.ylim(10., 200.)
        pylab.title('dc calorimetry')
    if 1:
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/dcT.npy', dct)
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/acT.npy', act)
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/dcppr.npy', dcppr)
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/acppr.npy', acppr)
    pylab.xlabel('Temperature (C)',fontsize=16)

    t1=80.
    t2=280.
    pylab.xlim(t1, t2)
    i1=numpy.argmin((act-t1)**2)
    i2=numpy.argmin((act-t2)**2)

    leg=pylab.legend(loc=2)
    temp=[t.set_fontsize(14) for t in leg.texts]

    if 0:
        pylab.figure()
        pylab.plot(time_s[i1:i2]*1000., v[i1:i2], 'b-')
        pylab.xlim(time_s[i1]*1000., time_s[i2]*1000.)
        pylab.xlabel('time (ms)',fontsize=16)
        pylab.ylabel('device voltage (V)',fontsize=16)
        #pylab.xlim(50, 80)
        
    if 0:
        pylab.figure()
        pylab.plot(time_s[i1:i2]*1000., fv[i1:i2], 'b-')
        pylab.xlim(time_s[i1]*1000., time_s[i2]*1000.)
        pylab.xlabel('time (ms)',fontsize=16)
        pylab.ylabel('filtered voltage (V)',fontsize=16)
        pylab.xlim(50, 80)
    if 0:
        pylab.figure()
        if 1:
            fftv2w=savgolsmooth(fftv2w, nptsoneside=300., order = 1)
        pylab.plot(time_s[i1:i2]*1000., fftv2w[i1:i2]*1.e6, 'b-')
        pylab.xlim(time_s[i1]*1000., time_s[i2]*1000.)
        pylab.xlabel('time (ms)',fontsize=16)
        pylab.ylabel('device 2$\omega$ amplitude($\mu$V)',fontsize=16)
        
        
if 1:
    pylab.figure()
    skip=100
    skipe=300
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

if 1:
    pylab.figure()
    ax1=pylab.subplot(111)
    mc=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/mcfft.npy')
    mcmol=mc/7.5e-8
    if 0:
        pylab.plot(act, mcmol, 'b-', lw=2, label='ac-analysis of 2$\omega$-component')
    else:
        pylab.plot(dct, mcmol, 'b-', lw=2, label='ac-analysis of 2$\omega$-component')
    pylab.xlabel('Temperature (C)',fontsize=16)
    pylab.ylabel('$C_{Sn}$+$C_{device}$+background (J/K/mol)',fontsize=16)
    pylab.xlim(80, 280)
    pylab.ylim(8, 30)
    pylab.title('ac calorimetry')
    leg=pylab.legend(loc=2)
    temp=[t.set_fontsize(14) for t in leg.texts]
    if 0:
        ax2=ax1.twinx()
        Sn_bct=[298,505.08,-5855.135,65.443315,-15.961,-1.88702000E-02,3.12116700E-06,-6.19600000E+04,0]
        Sn_liq=[505.08,800,9496.31,-9.809114,-8.2590486,-1.68144290E-02,2.62313100E-06,-1.08124400E+06,0]
        
        def Cp_pars(p, dT=.1):
            T=numpy.linspace(p[0], p[1], (p[1]-p[0])//dT+1)
            return T, -1.*(p[4]+2*p[5]*T**1+6*p[6]*T**2+2*p[7]*T**-2+42*p[8]*T**6)
        Ts, Cps=Cp_pars(Sn_bct)
        Tl, Cpl=Cp_pars(Sn_liq)
        pylab.plot(Ts-298.15, Cps, 'r:', lw=2)
        pylab.plot(Tl-298.15, Cpl, 'r:', lw=2)
        pylab.ylabel('$C_{Sn}$, bulk (J/K/mol)',fontsize=16)
        pylab.xlabel('Temperature (C)')
        pylab.xlim(80, 280)
        pylab.ylim(27.5, 35)
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        
pylab.show()
