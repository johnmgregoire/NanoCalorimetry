import numpy, h5py, pylab
from PnSC_h5io import *

from matplotlib.ticker import FuncFormatter

def myexpformat(x, pos):
    for ndigs in range(5):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
#f=h5py.File(p,mode='r')

seg=3
exp='Sn_20kHz'
skip=200
skipe=100
f, hpp=experimenthppaths(p, exp)
daqHz=f[hpp[0]].attrs['daqHz']
f.close()


hpp=['/Calorimetry/Sn_20kHz/measurement/HeatProgram/cell29_17.5dc17ac_269ms_20kHz_1_of_1', '/Calorimetry/Sn_20kHz/measurement/HeatProgram/cell7_17.5dc17ac_269ms_20kHz_1_of_1', '/Calorimetry/Sn_20kHz/measurement/HeatProgram/cell7_17.5dc17ac_269ms_20kHzagain_1_of_1']
labs=['20kHz, 10Ohm Res','slow ramp, scan1', 'slow ramp, scan2']

targetf=2.e4


#labs=[hp.rpartition('/')[2] for hp in hpp]

nplots=4
pylab.figure(figsize=(20, 8))
for i, (hp, title) in enumerate(zip(hpp, labs)):
    hpsdl=CreateHeatProgSegDictList(p, exp, hp.rpartition('/')[2])
    sampv=hpsdl[seg]['samplevoltage'][0][skip:-1*skipe]
    diffv=hpsdl[seg]['samplehighpassacvoltage'][0][skip:-1*skipe]
    t=hpsdl[seg]['cycletime'][0][skip:-1*skipe]
    pylab.subplot(len(hpp), nplots, nplots*i+1)
    sy=sampv
    pylab.plot((t*1000.)[:4000], sy[:4000], 'g.', markersize=1)
    pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
    #pylab.ylim(-620, 620)
    pylab.title(title)
    pylab.ylabel('sample channel V')
    pylab.subplot(len(hpp), nplots, nplots*i+2)
    y=diffv
    pylab.plot((t*1000.)[:4000], y[:4000], 'r.', markersize=1)
    pylab.gca().yaxis.set_major_formatter(ExpTickLabels)
    #pylab.ylim(-620, 620)
    pylab.title(title)
    pylab.ylabel('filtered channel, V')
    pylab.subplot(len(hpp), nplots, nplots*i+3)
    fft=numpy.fft.fft(y)
    freq=numpy.fft.fftfreq(len(y))*daqHz
    pylab.loglog(freq[:len(freq)//2], numpy.abs(fft[:len(freq)//2]))
    pylab.ylabel('filtered channel fft mag.')
    pylab.subplot(len(hpp), nplots, nplots*i+4)
    pylab.loglog(freq[:len(freq)//2], numpy.abs(fft[:len(freq)//2]))
    pylab.xlim(.9*targetf, 4*targetf)
    pylab.xticks([targetf, 2.*targetf, 3.*targetf])
    pylab.ylabel('filtered channel fft mag.')
pylab.subplot(len(hpp), nplots, nplots*i+1)
pylab.xlabel('time (ms)')
pylab.subplot(len(hpp), nplots, nplots*i+2)
pylab.xlabel('time (ms)')
pylab.subplot(len(hpp), nplots, nplots*i+3)
pylab.xlabel('freq (Hz)')
pylab.subplot(len(hpp), nplots, nplots*i+4)
pylab.xlabel('freq (Hz)')

pylab.suptitle('response for 10mAdc+9mAac into 10$\Omega$')
pylab.subplots_adjust(left=.07, right=.97, wspace=.35, hspace=.25)
if True:
    pylab.savefig(os.path.join('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis', '_'.join(('FFT', exp)))+'.png')
pylab.show()
    
    



