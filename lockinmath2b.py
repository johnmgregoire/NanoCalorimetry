import numpy, h5py, pylab
from PnSC_h5io import *
import scipy.optimize
from matplotlib.ticker import FuncFormatter

def myexpformat(x, pos):
    for ndigs in range(5):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110707_initAC.h5'
#f=h5py.File(p,mode='r')

exp='10Ohm'
skip=200
f, hpp=experimenthppaths(p, exp)
daqHz=f[hpp[0]].attrs['daqHz']
f.close()

#f['/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_1_of_1']
hp='cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_1_of_1'
hpsdl=CreateHeatProgSegDictList(p, exp, hp)
sampv=hpsdl[2]['samplevoltage'][0][skip:-1*skip]
diffv=hpsdl[2]['samplehighpassacvoltage'][0][skip:-1*skip]
t=hpsdl[2]['cycletime'][0][skip:-1*skip]
sampv-=sampv.mean()

def sinarr(nptspercycle, npts, ph=0.):
    if isinstance(npts, numpy.ndarray):
        npts=len(npts)
    return numpy.sin(numpy.arange(npts)*2.*numpy.pi/nptspercycle+ph)

ptspercyc=30.

sampv=numpy.concatenate([numpy.arange(600), numpy.arange(600)[::-1]/2.+300.])
sampv=sinarr(ptspercyc, sampv, ph=0.2)*sampv/600.


#sampv=sinarr(ptspercyc, 2000, ph=0.3)*10.+4.
#y=sinarr(ptspercyc/2., 2000, ph=2.)*2.+4.

#sampv=numpy.concatenate([numpy.arange(1000), numpy.arange(1000)[::-1]/2.+500.])
#sampv=sinarr(30, 2000, ph=.2)+sinarr(30/2., 2000, ph=.2)+sinarr(30/3., 2000, ph=.2)

minfcn=lambda n, x: 1./((x*sinarr(n, x, ph=0)).sum()**2+(x*sinarr(n, x, ph=numpy.pi/2.)).sum()**2)

def sfft2cfft2(x, ptspercyc, ncyclewin=1.):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=x*sinarr(ptspercyc, x, ph=0)
    c=x*sinarr(ptspercyc, x, ph=numpy.pi/2.)
    return numpy.array([(numpy.fft.fft(s[i:i+npts])[0])**2+(numpy.fft.fft(c[i:i+npts])[0])**2 for i in numpy.arange(len(x)-npts)])*4./npts**2

def lia_ampphase(x, ptspercyc, ncyclewin=1., returnphase=True, pad=True):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=x*sinarr(ptspercyc, x, ph=0)
    c=x*sinarr(ptspercyc, x, ph=numpy.pi/2.)
    amp=(numpy.array([(numpy.fft.fft(s[i:i+npts])[0])**2+(numpy.fft.fft(c[i:i+npts])[0])**2 for i in numpy.arange(len(x)-npts)])**.5)*2./npts
    if returnphase:
        phase=numpy.array([numpy.arctan(numpy.fft.fft(s[i:i+npts])[0]/numpy.fft.fft(c[i:i+npts])[0]) for i in numpy.arange(len(x)-npts)])
    if pad:
        amp=numpy.concatenate([amp[:npts//2], amp, amp[-1*(len(x)-len(amp)-npts//2):]])
        if returnphase:
            phase=numpy.concatenate([phase[:npts//2], phase, phase[-1*(len(x)-len(phase)-npts//2):]])
    if returnphase:
        return amp, phase
    return amp
    
def liaharmonic_relphase(x, phaserefdata, ptspercyc, ncyclewin_1w=2., ncyclewin_nw=6., harmonic=3., phaseshift=0., pad=True):
    #calculated phase is wrt a cosine reference
    # x is harmonic data, harmonic = 1 is ok, phaserefdata should be same length but if it is shorter it will be padded to make length x
    #ncyclewin_nw is number of harmonic cycles (as opposed to  1w cycles)
    npts=numpy.round(ptspercyc*ncyclewin_1w)
    s=phaserefdata*sinarr(ptspercyc, phaserefdata, ph=0)
    c=phaserefdata*sinarr(ptspercyc, phaserefdata, ph=numpy.pi/2.)
    ph1w=numpy.array([numpy.arctan(numpy.fft.fft(s[i:i+npts])[0]/numpy.fft.fft(c[i:i+npts])[0]) for i in numpy.arange(len(phaserefdata)-npts)])
    phnw=numpy.concatenate([ph1w[:npts//2.], ph1w, ph1w[-1*(len(x)-len(ph1w)-npts//2.):]])
    phnw-=phaseshift
    pylab.figure()
    pylab.plot(ph1w)
    pylab.plot(phnw)
    pylab.figure()
    hptspc=ptspercyc/harmonic
    nptsnw=numpy.round(hptspc*ncyclewin_nw)
    s=sinarr(hptspc, x, ph=0)
    c=sinarr(hptspc, x, ph=numpy.pi/2.)
    pylab.plot(numpy.array([(numpy.fft.fft(x[i:i+nptsnw]*(c[i:i+nptsnw]*numpy.cos(p)+s[i:i+nptsnw]*numpy.sin(p)))[0])**2 for i, p in zip(numpy.arange(len(x)-nptsnw), phnw[nptsnw//2:])])*4./nptsnw**2)
    pylab.plot(sfft2cfft2(x, hptspc, ncyclewin_nw), 'k--')
    pylab.show()
    amp=(numpy.array([(numpy.fft.fft(x[i:i+nptsnw]*(c[i:i+nptsnw]*numpy.cos(p)+s[i:i+nptsnw]*numpy.sin(p)))[0])**2 for i, p in zip(numpy.arange(len(x)-nptsnw), phnw[nptsnw//2:])])**0.5)*2./nptsnw
    if pad:
        amp=numpy.concatenate([amp[:nptsnw//2], amp, amp[-1*(len(x)-len(amp)-nptsnw//2):]])

#ans=liaharmonic_relphase(y, sampv, ptspercyc, ncyclewin_1w=2., ncyclewin_nw=4., harmonic=2., phaseshift=1.7)

#st=(sampv*sinarr(33.029, sampv, ph=0)).sum()
#ct=(sampv*sinarr(33.029, sampv, ph=numpy.pi/2.)).sum()
#print st**2+ct**2
#
#ptspercyc=scipy.optimize.fmin(minfcn,numpy.float64([30.]),(sampv, ))[0]
#print ptspercyc


cycs=numpy.arange(len(sampv))/ptspercyc
pylab.figure()
pylab.plot(cycs, sampv)
pylab.figure()
nplts=1
for i, h in enumerate((1, 2, 3)):
    for j, n in enumerate((1, 2, 10)):
        pylab.subplot(3, nplts, i*nplts+1)
        x=sfft2cfft2(sampv, ptspercyc/h, n*h)**.5
        offset=(len(cycs)-len(x))//2
        pylab.plot(cycs[offset:offset+len(x)], x, '.', markersize=1, label='%d$\omega$, %d cycle calc' %(h, n*h))

for i in range(nplts*3):
    pylab.subplot(3, nplts, i+1)
    pylab.legend(loc=1)
    pylab.rcParams.update({'legend.fontsize':10})
#pylab.figure(
pylab.show()
