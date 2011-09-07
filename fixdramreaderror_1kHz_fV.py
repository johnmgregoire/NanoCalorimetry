import numpy, h5py, pylab, copy
from PnSC_h5io import *
import scipy.optimize
from matplotlib.ticker import FuncFormatter

class fitfcns: #datatuples are x1,x2,...,y
    #.finalparams .sigmas .parnames useful, returns fitfcn(x)
    def genfit(self, fcn, initparams, datatuple, markstr='unspecified', parnames=[], interaction=0,  maxfev=2000, weights=None, optimizerfcn=None, full_output=1):
        self.maxfev=maxfev
        self.performfit=True
        self.initparams=initparams
        self.sigmas=scipy.zeros(len(initparams))
        self.parnames=parnames
        self.finalparams=initparams
        self.error=False
        if weights is None:
            def wts(x):
                return 1.
        elif weights=='parabolic':
            a=(datatuple[0][0]+datatuple[0][-1])/2.0
            b=(datatuple[0][-1]-datatuple[0][0])/2.0
            def wts(x):
                return 1.0+((x-a)/b)**2

        def res1(p, x1, y):
            return (y-fcn(p, x1))*wts(x1)

        def res2(p, x1,x2,y):
            return y-fcn(p, x1, x2)

        def res3(p, x1,x2,x3, y):
            return y-fcn(p, x1, x2, x3)

        def res4(p, x1,x2,x3,x4,  y):
            return y-fcn(p, x1, x2, x3, x4)

        resdic={1:res1,  2:res2,  3:res3,  4:res4}
        self.resfcn=resdic[len(datatuple)-1]
        i=0
        for arr in datatuple:  #if the numerical data is given as a list or tuple then convert to arrays. regardless convert to float64 because leastsq REQUIRES THIS
            datatuple=datatuple[0:i]+tuple([numpy.float64(arr)])+datatuple[i+1:]
            i=i+1
        while self.performfit:
            self.sigmas=scipy.zeros(len(self.finalparams))
            if not optimizerfcn is None:
                try:
                    self.finalparams=optimizerfcn(self.resfcn,self.initparams, args=datatuple, maxfun=self.maxfev, xtol=1.e-10, ftol=1.e-10)
                    self.error=0
                except:
                    self.error=1
            else:
                fitout = scipy.optimize.leastsq(self.resfcn,self.initparams, args=datatuple, maxfev=self.maxfev, full_output=1)#, warning=False)
                self.performfit=False
                self.finalparams=fitout[0]
                if full_output>0 and not fitout[4] in [1, 2]:
                    print 'Fitting Error ', fitout[4], ' at ', markstr,': ', fitout[3]
                    self.error=True
                else:
                    #self.finalparams=fitout[0]
                    self.covmat=fitout[1]
                    try:
                        self.sigmas=scipy.array([self.covmat[i, i] for i in range(len(self.sigmas))])
                    except:
                        pass

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
#f=h5py.File(p,mode='r')

seg=3
exp='Sn_1kHz'
skip=50
skipe=20
f, hpp=experimenthppaths(p, exp)
daqHz=f[hpp[0]].attrs['daqHz']
f.close()

hpp=['/Calorimetry/Sn_1kHz/measurement/HeatProgram/cell29_17.5dc17ac_269ms_1kHz_1_of_1', '/Calorimetry/Sn_1kHz/measurement/HeatProgram/cell7_17.5dc17ac_269ms_1kHz_1_of_1', '/Calorimetry/Sn_1kHz/measurement/HeatProgram/cell7_17.5dc17ac_269ms_1kHzagain_1_of_1']
labs=['1kHz, 10Ohm Res','slow ramp, scan1', 'slow ramp, scan2']

n=40
nptscyc=150
criterrv=0.00002
#criterrf=1.2
def fitfcn(p, x):
    return (p[0]+numpy.sin(x*2*numpy.pi/nptscyc-p[1]))*(p[2]+p[3]*x)
for hp, lab in zip(hpp, labs):
    f=h5py.File(p, mode='r')
    a=f[hp].attrs['fullscale_fields']
    f.close()
    #fullscale=[fullscale[temp] for temp in [0, 3, 4]]
    fullscale=[]
    a, b, c=a.partition('\t')
    a=a.lstrip('0').rstrip('.')
    a=(a=='' and (0,) or (eval(a),))[0]
    fullscale+=[a]
    a, b, c=c.partition('\t')
    a, b, c=c.partition('\t')
    a, b, c=c.partition('\t')
    a=a.lstrip('0').rstrip('.')
    a=(a=='' and (0,) or (eval(a),))[0]
    fullscale+=[a]
    a, b, c=c.partition('\t')
    a=a.lstrip('0').rstrip('.')
    a=(a=='' and (0,) or (eval(a),))[0]
    fullscale+=[a]
    
    
    print hp.rpartition('/')[2]
    hpsdl=CreateHeatProgSegDictList(p, exp, hp.rpartition('/')[2])
    fulllen=len(hpsdl[seg]['samplecurrent'][0])
    sampc=hpsdl[seg]['samplecurrent'][0]
    sampv=hpsdl[seg]['samplevoltage'][0]
    diffv=hpsdl[seg]['samplefilteredvoltage'][0]
    t=hpsdl[seg]['cycletime'][0]
    saveinds=hpsdl[seg]['segment_inds']
    pylab.figure(figsize=(19, 10))

    #for i, (arr, nam, fs) in enumerate(zip([sampc, sampv, diffv], ['current (A)', 'voltage (V)', 'filtered voltage (V)'], fullscale)):
    for i, (arr, nam, fs) in enumerate(zip([diffv], ['samplecurrent', 'samplevoltage','samplefilteredvoltage'], fullscale)):
        fixinds=numpy.uint32([])
        newarr=copy.copy(arr)
        starti=0
        while starti<len(arr):
            
            stopi=starti+n
            si=round(.2*n)-1
            if stopi+n>len(arr):
                stopi=len(arr)
            else:
                ei=stopi-starti
                ei=round(.8*n)+1
            Ya=arr[starti:stopi]
            
            Xa=numpy.float64(range(len(Ya)))

            ff=fitfcns()
            ff.genfit(fitfcn, [1., 0., Xa.mean(), 0.004], (Xa, Ya), full_output=0)
            fp=ff.finalparams
            z=fitfcn(fp, Xa)
#            pylab.plot(Ya)
#            pylab.plot(z)
#            pylab.show()

            diff=(Ya-z)
            #pow2diff=numpy.log((numpy.abs(diff)/fs*32768.))/numpy.log(2)
            #inds=numpy.where(pow2diff>4.5)[0]
            inds=numpy.where((numpy.abs(diff)>criterrv)&(Xa>=si)&(Xa<=ei))[0]
            if len(inds)>0:
                #ninds=numpy.where(pow2diff<=4.5)[0]
                ninds=numpy.where((numpy.abs(diff)<=criterrv)|(Xa<si)|(Xa>ei))[0]
                ff.genfit(fitfcn, fp, (Xa[ninds], Ya[ninds]), full_output=0)
                fp=ff.finalparams
                z=fitfcn(fp, Xa)
                diff=Ya-z
                inds=numpy.where((numpy.abs(diff)<=criterrv)&(Xa>=si)&(Xa<=ei))[0]
                diff=diff[inds]
                newarr[inds+starti]=numpy.float32(Ya[inds]-diff)
                fixinds=numpy.concatenate([fixinds, inds+starti])
                if len(inds)>10:
                    pylab.plot(Ya, 'k.')
                    pylab.plot(inds, z[inds], 'r.')
                    pylab.show()
            starti=stopi
        pylab.subplot(3, 1, i+1)
        print len(fixinds)
        pylab.plot(arr, 'k.', markersize=1)
        pylab.plot(fixinds, newarr[fixinds], 'r.', markersize=3)
        if 0:
            f=h5py.File(p, mode='r+')
            f[hp][nam][0, saveinds[0]:saveinds[1]]=newarr[:]*1000.#for V,A to mV,mA 
            f.close()

pylab.show()
            
