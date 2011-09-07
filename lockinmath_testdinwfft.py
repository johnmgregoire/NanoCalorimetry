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
    return amp

def windowfft_ampphase(x, npts_win, pad=True, ptspercyc=None):
    comp=numpy.array([numpy.fft.fft(x[i:i+npts_win])[:npts_win//2+1] for i in numpy.arange(len(x)-npts_win)])
    amp, phase=numpy.array([[numpy.abs(c), numpy.angle(c)] for c in comp]).swapaxes(0, 1)
#    if not ptspercyc is None:
#        ps=numpy.linspace(0, len(x)-npts_win, len(x)-npts_win, endpoint=False)*2.*numpy.pi/ptspercyc
#        phase=numpy.array([(p+numpy.pi+ps)%(2*numpy.pi)-numpy.pi for p in phase.T]).T
    if pad:
        amp=numpy.concatenate([amp[:npts_win//2], amp, amp[-1*(len(x)-len(amp)-npts_win//2):]])
        phase=numpy.concatenate([phase[:npts_win//2], phase, phase[-1*(len(x)-len(phase)-npts_win//2):]])
    return amp, phase
    
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
#f=h5py.File(p,mode='r')

seg=3
exp='Sn_10kHz_4e4Ks'
skip=50
skipe=20
f, hpp=experimenthppaths(p, exp)
daqHz=f[hpp[0]].attrs['daqHz']
f.close()

hpp=['/Calorimetry/Sn_10kHz_4e4Ks/measurement/HeatProgram/cell29_57.75dc56.1ac_10kHz_12.6ms_1_of_1', '/Calorimetry/Sn_10kHz_4e4Ks/measurement/HeatProgram/cell7_57.75dc56.1ac_10kHz_12.6ms_1_of_1']
labs=['10kHz, 10Ohm Res','fast ramp']


ptspercyc=30.

num1wcyclestoave=2.

#ans=liaharmonic_relphase(y, sampv, ptspercyc, ncyclewin_1w=2., ncyclewin_nw=4., harmonic=2., phaseshift=1.7)
for hp, lab in zip(hpp, labs):
    print hp.rpartition('/')[2]
    hpsdl=CreateHeatProgSegDictList(p, exp, hp.rpartition('/')[2])
    fulllen=len(hpsdl[seg]['samplecurrent'][0])
    sampc=hpsdl[seg]['samplecurrent'][0][skip:-1*skipe]
    sampv=hpsdl[seg]['samplevoltage'][0][skip:-1*skipe]
    diffv=hpsdl[seg]['samplehighpassacvoltage'][0][skip:-1*skipe]
    t=hpsdl[seg]['cycletime'][0][skip:-1*skipe]
    
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
        for j, (h, relphase) in enumerate(zip([1., 2., 3.], [0., 0., numpy.pi/2])):
            amp, phase=windowfft_ampphase(arr, ptspercyc*num1wcyclestoave, ptspercyc=ptspercyc)
            n=ptspercyc*num1wcyclestoave
            freqs=numpy.linspace(0,n//2)/n
            ffti=numpy.argmin((freqs-h/ptspercyc)**2)
            amp=amp[:, ffti]
            phase=phase[:, ffti]
            ax1=pylab.subplot(3, 4, i*4+j+2)
            pylab.title('%d$\omega$ lock-in with %.1f $\mu$s averaging' %(h, numpy.round(num1wcyclestoave*ptspercyc)/daqHz*1.e6))
            pylab.ylabel(nam)
            pylab.xlabel('num 1$\omega$ cycles')
            ax1.plot(numcyc, amp, 'k.', markersize=1)
            if i==1 or i==2:
                relphaseamp=liaharmonic_relphase(diffv, sampv, ptspercyc, ncyclewin_1w=num1wcyclestoave, ncyclewin_nw=num1wcyclestoave*h, harmonic=h, phaseshift=relphase)
                ax1.plot(numcyc, relphaseamp, 'g.', markersize=1, alpha=.5)
            ax1.yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            ax1.yaxis.set_major_formatter(ExpTickLabels)
            ax2 = ax1.twinx()
            ax2.plot(numcyc, (phase)*180./numpy.pi, 'r.', markersize=1, alpha=.4)
            ax2.set_ylabel('phase wrt 1$\omega$ ave. (deg)', color='r')
            ax2.yaxis.set_major_formatter(pylab.ScalarFormatter(useOffset=False))
            for tl in ax2.get_yticklabels():
                tl.set_color('r')
    pylab.subplots_adjust(left=0.08, right=0.92, wspace=0.85, hspace=.32, bottom=.05, top=.96)
    if False:
        pylab.savefig(os.path.join('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/', '_'.join(('LIA', exp, hp.rpartition('/')[2])))+'.png')
pylab.show()

