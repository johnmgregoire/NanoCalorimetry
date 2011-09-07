import numpy, h5py, pylab
from PnSC_h5io import *
import scipy.optimize
from matplotlib.ticker import FuncFormatter

def myexpformat(x, pos):
    for ndigs in range(5):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)


def sinarr(nptspercycle, npts, ph=0.):
    if isinstance(npts, numpy.ndarray):
        npts=len(npts)
    return numpy.sin(numpy.arange(npts)*2.*numpy.pi/nptspercycle+ph)
    
def lia_ampphase(x, ptspercyc, ncyclewin=1., returnphase=True, pad=True, phaseshift=0.):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=x*sinarr(ptspercyc, x, ph=phaseshift)
    c=x*sinarr(ptspercyc, x, ph=numpy.pi/2.+phaseshift)
    amp=(numpy.array([(numpy.fft.fft(s[i:i+npts])[0].real)**2+(numpy.fft.fft(c[i:i+npts])[0].real)**2 for i in numpy.arange(len(x)-npts)])**.5)*2./npts
    if returnphase:
        phase=numpy.array([numpy.arctan(numpy.fft.fft(s[i:i+npts])[0].real/numpy.fft.fft(c[i:i+npts])[0].real) for i in numpy.arange(len(x)-npts)])
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
#    pylab.figure()
#    pylab.plot(ph1w)
#    pylab.plot(phnw)
#    pylab.figure()
    hptspc=ptspercyc/harmonic
    nptsnw=numpy.round(hptspc*ncyclewin_nw)
    s=sinarr(hptspc, x, ph=0)
    c=sinarr(hptspc, x, ph=numpy.pi/2.)
#    pylab.plot(numpy.array([(numpy.fft.fft(x[i:i+nptsnw]*(c[i:i+nptsnw]*numpy.cos(p)+s[i:i+nptsnw]*numpy.sin(p)))[0])**2 for i, p in zip(numpy.arange(len(x)-nptsnw), phnw[nptsnw//2:])])*4./nptsnw**2)
#    pylab.plot(sfft2cfft2(x, hptspc, ncyclewin_nw), 'k--')
#    pylab.show()
    amp=(numpy.array([(numpy.fft.fft(x[i:i+nptsnw]*(c[i:i+nptsnw]*numpy.cos(p)+s[i:i+nptsnw]*numpy.sin(p)))[0].real)**2 for i, p in zip(numpy.arange(len(x)-nptsnw), phnw[nptsnw//2:])])**0.5)*2./nptsnw
    if pad:
        amp=numpy.concatenate([amp[:nptsnw//2], amp, amp[-1*(len(x)-len(amp)-nptsnw//2):]])


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110707_initAC.h5'
#f=h5py.File(p,mode='r')

exp='10Ohm_1kHz'

skip=12000
maxpts=6000
f, hpp=experimenthppaths(p, exp)
daqHz=f[hpp[0]].attrs['daqHz']
f.close()

ptspercyc=300.

num1wcyclestoave=4.

#hpp=[\
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_nofilter_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_1_of_1'\
#]
#
#labs=['no filter','19-32kHz filter', '19-32kHz filter +20dB in 10dB out']

#hpp=[\
#'/Calorimetry/10Ohm_20kHz/measurement/HeatProgram/cell29_10ohm_10mAdc_9mAac_15ptscyc_gain8_nofilter_1_of_1', \
#'/Calorimetry/10Ohm_20kHz/measurement/HeatProgram/cell29_10ohm_10mAdc_9mAac_15ptscyc_gain1_SRS39kHz61kHz_20dBin10dBout_1_of_1'\
#]
#
#labs=['no filter','39-61kHz filter +20dB in 10dB out']


hpp=['/Calorimetry/10Ohm_1kHz/measurement/HeatProgram/cell29_1kHz_gain8_10mAdc9mAac_nofilter_1_of_1', '/Calorimetry/10Ohm_1kHz/measurement/HeatProgram/cell29_1kHz_10mAdc9mAac_1.9kHz3.2kHz_20dBin20dBout_1_of_1']
labs=['no filter','1.9-3.2kHz filter']

#hpp=['/Calorimetry/10Ohm_10kHz_100_90mA/measurement/HeatProgram/cell29_10kHz_gain2_100mAdc90mAac_nofilter_1_of_1', '/Calorimetry/10Ohm_10kHz_100_90mA/measurement/HeatProgram/cell29_10kHz_gain2_100mAdc90mAac_19kHz32kHz_10dBin20dBout_1_of_1']
#labs=['no filter','19-32kHz filter +20dB in 10dB out']

#ans=liaharmonic_relphase(y, sampv, ptspercyc, ncyclewin_1w=2., ncyclewin_nw=4., harmonic=2., phaseshift=1.7)
for hp, lab in zip(hpp, labs):
    print hp.rpartition('/')[2]
    hpsdl=CreateHeatProgSegDictList(p, exp, hp.rpartition('/')[2])
    fulllen=len(hpsdl[2]['samplecurrent'][0])
    sampc=hpsdl[2]['samplecurrent'][0][skip:skip+min(fulllen-2*skip, maxpts)]
    sampv=hpsdl[2]['samplevoltage'][0][skip:skip+min(fulllen-2*skip, maxpts)]
    diffv=hpsdl[2]['samplehighpassacvoltage'][0][skip:skip+min(fulllen-2*skip, maxpts)]
    t=hpsdl[2]['cycletime'][0][skip:skip+min(fulllen-2*skip, maxpts)]
    
#    sampc=sinarr(ptspercyc, 2000, ph=0.)+sinarr(ptspercyc/2., 2000, ph=0.)+sinarr(ptspercyc/3., 2000, ph=0.)
#    t=numpy.arange(2000)
    
    numcyc=numpy.arange(len(t))/ptspercyc
    pylab.figure(figsize=(19, 10))
    for i, (arr, nam) in enumerate(zip([sampc, sampv, diffv], ['current (A)', 'voltage (V)', 'filtered voltage (V)'])):
        pylab.subplot(3, 4, i*4+1)
        pylab.plot(t*1000., arr, 'b.', markersize=1)
        pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
        pylab.xlabel('time (ms)')
        pylab.ylabel(nam)
        if i==0:
            pylab.title(lab)
        for j, h in enumerate([1., 2., 3.]):
            if j==0:
                amp, phase=lia_ampphase(arr, ptspercyc/h, ncyclewin=num1wcyclestoave*h)
                phase1wave=phase.mean()
            amp, phase=lia_ampphase(arr, ptspercyc/h, ncyclewin=num1wcyclestoave*h, phaseshift=-1.*phase1wave)
            ax1=pylab.subplot(3, 4, i*4+j+2)
            pylab.title('%d$\omega$ lock-in with %.1f $\mu$s averaging' %(h, numpy.round(num1wcyclestoave*ptspercyc)/daqHz*1.e6))
            pylab.ylabel(nam)
            pylab.xlabel('num 1$\omega$ cycles')
            ax1.plot(numcyc, amp, 'k.', markersize=1)
            ax1.yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            ax1.yaxis.set_major_formatter(ExpTickLabels)
            ax2 = ax1.twinx()
            ax2.plot(numcyc, (phase)*180./numpy.pi, 'r.', markersize=1)
            ax2.set_ylabel('phase wrt 1$\omega$ ave. (deg)', color='r')
            ax2.yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            for tl in ax2.get_yticklabels():
                tl.set_color('r')
    pylab.subplots_adjust(left=0.08, right=0.92, wspace=0.85, hspace=.32, bottom=.05, top=.96)
    if True:
        pylab.savefig(os.path.join('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110708_initACtests/analysis', '_'.join(('LIA', exp, hp.rpartition('/')[2])))+'.png')
pylab.show()

