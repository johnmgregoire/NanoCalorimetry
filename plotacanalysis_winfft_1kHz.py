import numpy, h5py, pylab, copy
from PnSC_h5io import *
from PnSC_math import *

if 0:
    p='/mnt/hgfs/HostDocuments/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    #p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')

    x=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_voltage'][0,:,:,0]
    f.close()

if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')

    x=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_filteredvoltage'][0,:,:,0]
    f.close()
    
if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')

    x=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_current'][0,:,:,0]
    f.close()
    
if 0:
    x=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fakefft_voltage.npy')[:, :, 0]
    
if 0:
    x=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fakefft_current.npy')[:, :, 0]

if 0:
    x=numpy.load('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake2fft_smoothrcurrent.npy')[:, :, 0]
if 0:
    pylab.imshow(numpy.log10(x),interpolation='nearest',origin='lower')
    pylab.ylim(20000,23000)
    pylab.gca().set_aspect(10.*x.shape[1]/x.shape[0])
    pylab.colorbar()


ptspercyc=300.
find_1w=4
if 0:
    for i in range(2):
        pylab.subplot(3, 2, 1+i)
        pylab.plot(x[:, 0], label='dc')
        pylab.plot(x[:, find_1w], label='1$\omega$')
        #pylab.gca().set_yscale('log')
        pylab.legend()

    for i in range(2):
        pylab.subplot(3, 2, 3+i)
        pylab.plot(x[:, find_1w*2], 'r', label='2$\omega$')
        pylab.plot(x[:, find_1w*2-1], 'b', label='2$\omega$-', alpha=.6)
        pylab.plot(x[:, find_1w*2+1], 'g', label='2$\omega$+', alpha=.6)
        #pylab.gca().set_yscale('log')
        pylab.legend()

    for i in range(2):
        pylab.subplot(3, 2, 5+i)
        pylab.plot(x[:, find_1w*3], 'r', label='3$\omega$')
        pylab.plot(x[:, find_1w*3-1], 'b', label='3$\omega$-', alpha=.6)
        pylab.plot(x[:, find_1w*3+1], 'g', label='3$\omega$+', alpha=.6)
        #pylab.gca().set_yscale('log')
        pylab.legend()

if 0:
    smooth1wfft=savgolsmooth(x[:, find_1w], nptsoneside=300, order = 1)
    smoothdc=savgolsmooth(x[:, 0], nptsoneside=300, order = 1)
    iarr=numpy.float32(range(x.shape[0]))
    fakedata=smoothdc+2.*smooth1wfft*numpy.sin(iarr/ptspercyc*numpy.pi*2.)

    if 1:
        pylab.plot(fakedata)

    fakefft=windowfft_ampphase(fakedata, 1200)
    fakefft=numpy.array(fakefft).swapaxes(0,2).swapaxes(0,1)
    if 1:
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake_voltage.npy', fakedata)
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fakefft_voltage.npy', fakefft[:, :22, :])
    if 0:
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake_current.npy', fakedata)
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fakefft_current.npy', fakefft[:, :22, :])

if 0:
    p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
    f=h5py.File(p,mode='r')
    x=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_voltage'][0,:,:,0]
    c=f['Calorimetry/1kHzcrop/analysis/again/WinFFT_current'][0,:,:,0]
    z=f['Calorimetry/1kHzcrop/measurement/HeatProgram/again/samplecurrent'][0,:]
    f.close()
    xsmoothdc=savgolsmooth(x[:, 0], nptsoneside=300, order = 1)
    csmoothdc=savgolsmooth(c[:, 0], nptsoneside=300, order = 1)
    r=xsmoothdc/csmoothdc/1000.
    rsmooth=savgolsmooth(r, nptsoneside=60, order = 1)

    iarr=numpy.float32(range(x.shape[0]))
    fakedata=rsmooth*z

    if 1:
        pylab.plot(fakedata)

    fakefft=windowfft_ampphase(fakedata, 1200)
    fakefft=numpy.array(fakefft).swapaxes(0,2).swapaxes(0,1)
    if 0:
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake2_smoothrcurrent.npy', fakedata)
        numpy.save('C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110714_Sn_analysis/Aug31kHzanalysis/fake2fft_smoothrcurrent.npy', fakefft[:, :22, :])
    
    if 0:
        x=fakefft[:, :, 0]
        pylab.imshow(numpy.log10(x[:, :-2]),interpolation='nearest',origin='lower')
        pylab.ylim(2000,2300)
        pylab.gca().set_aspect(10.*x.shape[1]/x.shape[0])
        pylab.colorbar()
if 1:

    iarr=numpy.float32(range(3600))
    z=1.+numpy.sin(iarr/ptspercyc*numpy.pi*2.)+numpy.sin(iarr/(ptspercyc/2.)*numpy.pi*2.)+numpy.sin(iarr/(ptspercyc/3.)*numpy.pi*2.)
    fft=windowfft_ampphase(z, 1200)
    fft=numpy.array(fft).swapaxes(0,2).swapaxes(0,1)[:, :22, :]
    if 1:
        x=fft[:, :, 0]
        pylab.imshow(numpy.log10(x[:, :]),interpolation='nearest',origin='lower')
        pylab.gca().set_aspect(50.*x.shape[1]/x.shape[0])
        pylab.colorbar()
if 1:
    pylab.show()
