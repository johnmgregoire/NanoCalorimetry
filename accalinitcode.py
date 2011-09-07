import numpy, h5py, pylab
from PnSC_h5io import *

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110707_initAC.h5'
#f=h5py.File(p,mode='r')

exp='10Ohmdc'
skip=200
f, hpp=experimenthppaths(p, exp)
daqHz=f[hpp[0]].attrs['daqHz']
f.close()

hpp=[\
'/Calorimetry/10Ohmdc/measurement/HeatProgram/cell29_10Ohm_10mA_gain1_noSRS_1_of_1', \
'/Calorimetry/10Ohmdc/measurement/HeatProgram/cell29_10Ohm_10mA_gain1_SRS0db_1_of_1',\
'/Calorimetry/10Ohmdc/measurement/HeatProgram/cell29_10Ohm_10mA_gain1_SRS20dbin_1_of_1',\
'/Calorimetry/10Ohmdc/measurement/HeatProgram/cell29_10Ohm_10mA_gain1_SRS20dbout_1_of_1'\
]

for i, (h, title) in enumerate(zip(hpp, ['no SRS', 'SRS 0dB', 'SRS 20dB in', 'SRS 20dB out'])):
    hpsdl=CreateHeatProgSegDictList(p, exp, h.rpartition('/')[2])
    sampv=hpsdl[2]['samplevoltage'][0][skip:-1*skip]
    diffv=hpsdl[2]['samplehighpassacvoltage'][0][skip:-1*skip]
    t=hpsdl[2]['cycletime'][0][skip:-1*skip]
    pylab.subplot(len(hpp), 2, 2*i+1)
    y=(diffv-diffv.mean())*1.e6
    pylab.plot(t*1000., y)
    pylab.ylim(-620, 620)
    pylab.title(title)
    pylab.ylabel('dev from mean, $\mu$V')
    pylab.subplot(len(hpp), 2, 2*i+2)
    fft=numpy.fft.fft(y)
    freq=numpy.fft.fftfreq(len(y))*daqHz
    pylab.loglog(freq[:len(freq)//2], numpy.abs(fft[:len(freq)//2]))
    pylab.xlabel('freq (Hz)')
    pylab.ylabel('Vsample fft mag.')
pylab.subplot(len(hpp), 2, 2*i+1)
pylab.xlabel('time (ms)')
    
pylab.suptitle('dc response for 10mA into 10$\Omega$')
pylab.show()
    
    



