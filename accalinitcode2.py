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

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110707_initAC.h5'
#f=h5py.File(p,mode='r')

exp='10Ohm_1kHz'
skip=12000
f, hpp=experimenthppaths(p, exp)
daqHz=f[hpp[0]].attrs['daqHz']
f.close()

#hpp=[\
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_wTtokeithley_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_9kHz11kHzfilter_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_9kHz11kHzfilter_wTtokeithley_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_nofilter_1_of_1'\
#]

#hpp=[\
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_nofilter_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_9kHz11kHzfilter_1_of_1'\
#]
#labs=['no filter', '9-11kHz filter']

#hpp=[\
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_wTtokeithley_1_of_1'\
#]
#labs=['19-32kHz filter without Keithley', '19-32kHz filter with Keithley']

#hpp=[\
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_nofilter_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_1_of_1', \
#'/Calorimetry/10Ohm/measurement/HeatProgram/cell29_10Ohm_10mAdc_9mA10kHz_19kHz32kHzfilter_gain1_20dBin10dBoutinout_1_of_1'\
#]
#
#labs=['no filter','19-32kHz filter', '19-32kHz filter +20dB in 10dB out']


#hpp=[\
#'/Calorimetry/10Ohm_20kHz/measurement/HeatProgram/cell29_10ohm_10mAdc_9mAac_15ptscyc_gain8_nofilter_1_of_1', \
#'/Calorimetry/10Ohm_20kHz/measurement/HeatProgram/cell29_10ohm_10mAdc_9mAac_15ptscyc_gain8_SRS19kHz21kHz_1_of_1', \
#'/Calorimetry/10Ohm_20kHz/measurement/HeatProgram/cell29_10ohm_10mAdc_9mAac_15ptscyc_gain1_SRS39kHz61kHz_20dBin10dBout_1_of_1'\
#]
#
#labs=['no filter','19-21kHz filter', '39-61kHz filter +20dB in 10dB out']

hpp=['/Calorimetry/10Ohm_1kHz/measurement/HeatProgram/cell29_1kHz_gain8_10mAdc9mAac_nofilter_1_of_1', '/Calorimetry/10Ohm_1kHz/measurement/HeatProgram/cell29_1kHz_10mAdc9mAac_1.9kHz3.2kHz_20dBin20dBout_1_of_1']
labs=['no filter','1.9-3.2kHz filter']

#hpp=['/Calorimetry/10Ohm_10kHz_100_90mA/measurement/HeatProgram/cell29_10kHz_gain2_100mAdc90mAac_nofilter_1_of_1', '/Calorimetry/10Ohm_10kHz_100_90mA/measurement/HeatProgram/cell29_10kHz_gain2_100mAdc90mAac_19kHz32kHz_10dBin20dBout_1_of_1']
#labs=['no filter','19-32kHz filter +20dB in 10dB out']

targetf=1.e3

#labs=[hp.rpartition('/')[2] for hp in hpp]
nplots=4
pylab.figure(figsize=(20, 8))
for i, (hp, title) in enumerate(zip(hpp, labs)):
    hpsdl=CreateHeatProgSegDictList(p, exp, hp.rpartition('/')[2])
    sampv=hpsdl[2]['samplevoltage'][0][skip:-1*skip]
    diffv=hpsdl[2]['samplehighpassacvoltage'][0][skip:-1*skip]
    t=hpsdl[2]['cycletime'][0][skip:-1*skip]
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
    pylab.savefig(os.path.join('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110708_initACtests/analysis', '_'.join(('FFT', exp)))+'.png')
pylab.show()
    
    



