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
sampv=sinarr(ptspercyc, sampv, ph=0.)*sampv/600.


#sampv=sinarr(ptspercyc, 2000, ph=.2)*10.+4.

#sampv=numpy.concatenate([numpy.arange(1000), numpy.arange(1000)[::-1]/2.+500.])
#sampv=sinarr(30, 2000, ph=.2)+sinarr(30/2., 2000, ph=.2)+sinarr(30/3., 2000, ph=.2)

minfcn=lambda n, x: 1./((x*sinarr(n, x, ph=0)).sum()**2+(x*sinarr(n, x, ph=numpy.pi/2.)).sum()**2)

def s2c2(x, ptspercyc, ncyclewin=1.):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=x*sinarr(ptspercyc, x, ph=0)
    c=x*sinarr(ptspercyc, x, ph=numpy.pi/2.)
    return numpy.array([(s[i:i+npts]).sum()**2+(c[i:i+npts]).sum()**2 for i in numpy.arange(len(x)-npts)])/npts**2

def sfft2cfft2(x, ptspercyc, ncyclewin=1.):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=sinarr(ptspercyc, npts, ph=0)
    c=sinarr(ptspercyc, npts, ph=numpy.pi/2.)
    return numpy.array([(numpy.fft.fft(x[i:i+npts]*s)[0])**2+(numpy.fft.fft(x[i:i+npts]*c)[0])**2 for i in numpy.arange(len(x)-npts)])*4./npts**2
    
#st=(sampv*sinarr(33.029, sampv, ph=0)).sum()
#ct=(sampv*sinarr(33.029, sampv, ph=numpy.pi/2.)).sum()
#print st**2+ct**2
#
#ptspercyc=scipy.optimize.fmin(minfcn,numpy.float64([30.]),(sampv, ))[0]
#print ptspercyc


t1=s2c2(sampv, ptspercyc, 1.)**.5
t2=s2c2(sampv, ptspercyc, 2.)**.5
t10=s2c2(sampv, ptspercyc, 10.)**.5

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
