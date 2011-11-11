import numpy, scipy, scipy.optimize, scipy.interpolate
import h5py
import os, os.path, time, copy
from PnSC_ui import *

#def ndiv(x, nptsoneside=20):
#    nptsoneside=max(nptsoneside, 2.)
#    gapptsoneside=min(gapptsoneside, nptsoneside-2.)
#    for gap in range(gapptsoneside+1):
#        starti=numpy.uint32([max(i-(nptsoneside-gap), 0) for i in range(len(arr))])
#        stopi=numpy.uint32([min(i+(nptsoneside-gap)+1, len(arr)) for i in range(len(arr))])
#        #print [numpy.append(arr[i0:i], arr[i+1:i1]) for i, i0, i1 in zip(range(len(arr)), starti, stopi)][8]
#        #print [(((numpy.append(arr[i0:i], arr[i+1:i1]).mean()-arr[i]))**2, (numpy.append(arr[i0:i], arr[i+1:i1]).std()*nsig)**2) for i, i0, i1 in zip(range(len(arr)), starti, stopi)][8]
#        arr=numpy.array([(((numpy.append(arr[i0:i], arr[i+1:i1]).mean()-arr[i]))**2<(numpy.append(arr[i0:i], arr[i+1:i1]).std()*nsig)**2 and (arr[i],) or (numpy.append(arr[i0:i], arr[i+1:i1]).mean(),))[0] for i, i0, i1 in zip(range(len(arr)), starti, stopi)], dtype=arr.dtype)
#    return arr

def concat_extrap_ends(x, npts, polyorder=1, lowside=True, highside=True):
    i=numpy.arange(npts, dtype='float64')
    if lowside:
        ans=scipy.polyfit(-1*(i+1.), x[:npts], polyorder)
        x=numpy.concatenate([scipy.polyval(list(ans), i[::-1]), x])
    if highside:
        ans=scipy.polyfit(-1*(i[::-1]-1.), x[-1*npts:], polyorder)
        x=numpy.concatenate([x, scipy.polyval(list(ans), i)])
    return x    
    
def lininterpbetweenregularpoints(existy, interval):
    existy=numpy.array(existy)
    x=numpy.arange(interval,dtype='float32')/interval
    diff=existy[1:]-existy[:-1]
    o=numpy.outer(diff,x)
    return numpy.concatenate([arr+start for arr,start in zip(o,existy[:-1])]+[existy[-1:]])
    
def interpwithinarr(existind, existy, order=3, interpplotax=None, interpcols=['k', 'r']):
    if order==1:
        existind=numpy.array(existind)
        diff=existind[1:]-existind[:-1]
        if numpy.all(diff==diff[0]):
            return lininterpbetweenregularpoints(existy, diff[0])
    interind=sorted(list(set(numpy.arange(max(existind)+1))-set(existind)))
    yfull=numpy.zeros(max(existind)+1, existy.dtype)
    yfull[existind]=existy[:]
    yfull[interind]=scipy.interpolate.spline(existind, existy, interind, order=order)
    if not interpplotax is None:
        interpplotax.plot(existind, existy, interpcols[0])
        interpplotax.plot(interind,  yfull[interind], interpcols[1])
    return yfull
    
def savgolsmooth(x, nptsoneside=7, order = 4, dx=1.0, deriv=0, binprior=0): #based on scipy cookbook. x is 1-d array, window is the number of points used to smooth the data, order is the order of the smoothing polynomial, will return the smoothed "deriv"th derivative of x
    if binprior>1:
        origlen=len(x)
        x=numpy.array([x[i*binprior:(i+1)*binprior].mean() for i in range(origlen//binprior)])
        dx*=binprior
    side=numpy.uint16(max(nptsoneside, numpy.ceil(order/2.)))
    s=numpy.r_[2*x[0]-x[side:0:-1],x,2*x[-1]-x[-2:-1*side-2:-1]]
    # a second order polynomal has 3 coefficients
    b = numpy.mat([[k**i for i in range(order+1)] for k in range(-1*side, side+1)])
    m = numpy.linalg.pinv(b).A[deriv] #this gives the dth ? of the base array (.A) of the pseudoinverse of b

    # precompute the offset values for better performance
    offsets = range(-1*side, side+1)
    offset_data = zip(offsets, m)

    smooth_data=[numpy.array([(weight * s[i + offset]) for offset, weight in offset_data]).sum() for i in xrange(side, len(s) - side)]
    smooth_data=numpy.array(smooth_data)/(dx**deriv)
    
    if binprior>1:    
        ia=numpy.arange(binprior, dtype='float32')/binprior
        xr=numpy.concatenate([ia*(b-a)+b for a, b in zip(smooth_data[:-1], smooth_data[1:])])
        xr=numpy.concatenate([(smooth_data[1]-smooth_data[0])*ia[:binprior//2]+smooth_data[0], xr, (smooth_data[-1]-smooth_data[-2])*ia[:binprior//2]+smooth_data[-1]])
        smooth_data=numpy.concatenate([xr, (smooth_data[-1]-smooth_data[-2])*ia[:origlen-len(xr)]+smooth_data[-1]])


    return smooth_data

class fitfcns: #datatuples are x1,x2,...,y
    #.finalparams .sigmas .parnames useful, returns fitfcn(x)
    def genfit(self, fcn, initparams, datatuple, markstr='unspecified', parnames=[], interaction=0,  maxfev=2000, weights=None, optimizerfcn=None):
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
                if not fitout[4] in [1, 2]:
                    print 'Fitting Error ', fitout[4], ' at ', markstr,': ', fitout[3]
                    self.error=True
                else:
                    #self.finalparams=fitout[0]
                    self.covmat=fitout[1]
                    try:
                        self.sigmas=scipy.array([self.covmat[i, i] for i in range(len(self.sigmas))])
                    except:
                        pass

        def fitfcn(x):
            return fcn(self.finalparams, x)
        return fitfcn

    def poly(self, p, x):#both must be numpy arrays
        return numpy.array([p[i]*(x**i) for i in range(p.size)]).sum(0)

    def polyfit(self, datatuple, initparams, markstr='unspecified', interaction=0,  maxfev=2000, weights=None):
        #initparams can be an array of coefficients [constant,lin term, quad term,...] or an integer indicating the order of the polynomial
        if isinstance(initparams, int):
            initparams=numpy.ones(initparams+1)
        else:
            initparams=numpy.float64(initparams)
        parnames=[]
        i=0
        for par in initparams:
            parnames+=[''.join(('coef', `i`))]
            i+=1

        return self.genfit(self.poly, initparams, datatuple, markstr, parnames, interaction, maxfev, weights=weights)


    def gaussianfit(self, datatuple, initparams=scipy.array([1, 0, 1]), markstr='unspecified', interaction=0, showplot=True, maxfev=2000, weights=None):
        return self.genfit(self.gaussian, initparams, datatuple, markstr, parnames=['coef', 'center', 'sigma'], interaction=interaction, maxfev=maxfev, weights=weights)

    def gaussian(self, p, x):
        return p[0]*scipy.exp(-0.5*((x-p[1])/p[2])**2)

    def lorentzianfit(self, datatuple, initparams=scipy.array([1, 0, 1]), markstr='unspecified', interaction=0, showplot=True, maxfev=2000, weights=None):
        return self.genfit(self, self.lorentzian, initparams, datatuple, markstr, parnames=['coef', 'center', 'gamma'], interaction=interaction, maxfev=maxfev, weights=weights)

    def lorentzian(self, p, x):
        return (p[0]/scipy.pi)*p[2]/((x-p[1])**2+p[2]**2)

def Gaussian(pars, x):
    return pars[2]*numpy.exp(-0.5*((x-pars[0])/pars[1])**2) 

def Lorentzian(pars, x):#defined in nontraditional way so that pars[2] is the peak height
    return pars[2]/(1+((x-pars[0])/pars[1])**2)

def GaussLorentz(pars, x):
    gw=min(max(pars[3], 0.), 1.)
    return gw*Gaussian(pars, x)+(1.-gw)*Lorentzian(pars, x)
    
def GaussHalfLorentz(pars, x):
    return .5*Gaussian(pars, x)+.5*Lorentzian(pars, x)

PeakFcnLibrary={'Gaussian':Gaussian, 'Lorentzian':Lorentzian, 'GaussHalfLorentz':GaussHalfLorentz}

def fitpeakset(X, Y, initpars, peakfcn, negpeaks=True, optimizerfcn=None, nsigcut=3.):#peak function must be a function that accepts a list of 3 parameters (the reshape 3 needs to be changed if num params differs)
    numgauss=len(initpars)
    if numgauss==0:
        return (numpy.float32([]), numpy.float32([]), 0.)
    if nsigcut is None:
        imin=0
        imax=len(X)
    else:
        xmin=initpars[0][0]
        xmax=initpars[0][0]
        for p, w, h in initpars:
            xmin=min(xmin, p-w*nsigcut)
            xmax=max(xmax, p+w*nsigcut)
        imin=numpy.argmin((X-xmin)**2)
        imax=numpy.argmin((X-xmax)**2)
    
    zeroedpeakinds=[]
    repeatwithpkremoved=True #peaks are removed if their fitted height is <0. At the end, these peaks are added to the fit parameter list with 0 height and 0 error
    while repeatwithpkremoved:
        initparscpy=copy.copy(list(initpars))
        for pkind in reversed(zeroedpeakinds):#reverse so opo gets the right index
            initparscpy.pop(pkind)
        if len(initparscpy)==0:
            break
        initparsflat=numpy.float64(initparscpy).flatten()
        def fitfcn(p, x):
            allpars=numpy.reshape(p, (p.size//initpars.shape[1], initpars.shape[1]))
            if isinstance(x, numpy.ndarray):
                val=numpy.zeros(x.size, dtype='float32')
            else:
                val=0.0
            for pars in allpars:
                val+=peakfcn(pars, x)
            return val
#        def residfcn(p, x, y):
#            err=y-fitfcn(p, x)
#            return err
        Ya=numpy.float64(Y[imin:imax])
        Xa=numpy.float64(X[imin:imax])
        
        #if not optimizerfcn is None:
        ff=fitfcns()
        ff.genfit(fitfcn, initparsflat, (Xa, Ya), optimizerfcn=optimizerfcn)
        finalparams=ff.finalparams
#        else:
#            fitout=scipy.optimize.leastsq(residfcn, initparsflat, args=(X, Y), full_output=1)
#            if not (fitout[4] in [1, 2]):
#                print 'Fitting Error', fitout[4],': ', fitout[3]
#            finalparams=numpy.float32(fitout[0])
        finalparamsshaped=numpy.reshape(finalparams, (len(finalparams)//initpars.shape[1], initpars.shape[1]))
        if negpeaks:
            repeatwithpkremoved=False
        else:
            negpeakinds=numpy.where(finalparamsshaped[:, 2]<0)[0]
            zeroedpeakinds+=list(negpeakinds)
            zeroedpeakinds.sort()
            repeatwithpkremoved=len(negpeakinds)>0
#        print '^^^^^^^^^^^^^^^'
#        print initparsflat
#        print finalparamsshaped
#        pylab.plot(X, Y, 'b.')
#        pylab.show()
#    if not (fitout[1] is None):
#        covmat=fitout[1]
#        sigmas=numpy.float32([covmat[i, i] for i in range(len(finalparams))])
#    else:
#        print 'COVARIANCE NOT CALCULATED:', fitout[4],': ', fitout[3]
#        sigmas=numpy.zeros(len(finalparams), dtype='float32')
        sigmas=ff.sigmas
    finalresid=numpy.sqrt((ff.resfcn(finalparams, X, Y)**2).sum())
    #pylab.plot(X, Y, 'k.', X, fitfcn(finalparams, X), 'r-')

    sigmashaped=numpy.reshape(sigmas, (len(finalparams)//initpars.shape[1], initpars.shape[1]))
    for pkind in zeroedpeakinds:
        finalparamsshaped=list(finalparamsshaped)
        sigmashaped=list(sigmashaped)
        temp=copy.copy(initpars[pkind][:])
        temp[2]=0.#zero the height
        finalparamsshaped.insert(pkind, temp)
        sigmashaped.insert(pkind, numpy.zeros(initpars.shape[1], dtype='float64'))
        finalparamsshaped=numpy.float64(finalparamsshaped)
        sigmashaped=numpy.float64(sigmashaped)
    return (finalparamsshaped, sigmashaped, finalresid)
    
def arrayzeroind1d(arr, postoneg=False, negtopos=True):
    sarr=numpy.sign(arr)
    if postoneg:
        zeroind=numpy.where(sarr[:-1]>sarr[1:])[0]
        if negtopos:
            zeroind=numpy.append(zeroind, numpy.where(sarr[:-1]*sarr[1:]<=0)[0])
    else:#assume that if not postoneg then negtopos
        zeroind=numpy.where(sarr[:-1]*sarr[1:]<=0)[0]
    return (1.0*zeroind*arr[(zeroind+1,)]-(zeroind+1)*arr[(zeroind,)])/(arr[(zeroind+1,)]-arr[(zeroind,)]) #returns array of the floating point "index" linear interpolation between 2 indeces

def clustercoordsbymax1d(arr, pkind, critsepind):#results will be sorted. wherever there are peak indeces too close together. the peak index next to the peak index with highest arr value gets removed
    pkind.sort()
    indindslow=numpy.where((pkind[1:]-pkind[:-1])<critsepind)[0]
    indindshigh=indindslow+1
    while indindslow.size>0:
        maxindindindlow=numpy.nanargmax(arr[pkind[(indindslow,)]])
        maxindindindhigh=numpy.nanargmax(arr[pkind[(indindshigh,)]])
        if arr[pkind[indindslow[maxindindindlow]]]>arr[pkind[indindshigh[maxindindindhigh]]]:
            pkind=numpy.delete(pkind, indindshigh[maxindindindlow])
        else:
            pkind=numpy.delete(pkind, indindslow[maxindindindhigh])

        indindslow=numpy.where((pkind[1:]-pkind[:-1])<critsepind)[0]
        indindshigh=indindslow+1
    return pkind
    
def peaksearch1dSG(x, dx=1., critpeakheight=10, critsepind=5, critcurve=None, firstdernpts=7, firstderorder=1, secdernpts=14, secderorder=1, pospeaks=True, negpeaks=True):
    #dx is delta q for one index. zeros of the first derivative of inn are grouped together if within critsepind. only negative slope in the firstder is used so no secder is necessary unless specify a critical curvature in count nm^2
    if not (pospeaks or negpeaks):
        return numpy.float32([])
    ifirstder=savgolsmooth(x, nptsoneside=firstdernpts, order=firstderorder, dx=dx, deriv=1)
    fullpkind=numpy.float32([])
    if pospeaks:
        zeroind=arrayzeroind1d(ifirstder, postoneg=True, negtopos=False)
        temp=numpy.where(x[(numpy.uint32(numpy.round(zeroind)),)]>critpeakheight)
        fullpkind=numpy.append(fullpkind, zeroind[temp])
    if negpeaks:
        zeroind=arrayzeroind1d(ifirstder, postoneg=False, negtopos=True)
        temp=numpy.where(x[(numpy.uint32(numpy.round(zeroind)),)]<(-1*critpeakheight))
        fullpkind=numpy.append(fullpkind, zeroind[temp])
        
    if fullpkind.size==0:
        return fullpkind
    pkind=clustercoordsbymax1d(x, numpy.uint32(numpy.round(fullpkind)), critsepind)
    if critcurve is not None:
        isecder=savgolsmooth(x, nptsoneside=secdernpts, order=secderorder, dx=dx, deriv=2)
        temp=numpy.where(numpy.abs(isecder[(numpy.uint32(numpy.round(pkind)),)])>(critcurve))
        pkind=numpy.array(pkind)[temp]
#    pkind=list(pkind)
#    pkind.reverse()#highest to smallest for pairing below
    return numpy.array(pkind, dtype=numpy.float32)
    
def removeoutliers_meanstd(arr, nptsoneside, nsig, gapptsoneside=0): #avrages maximum of 2*nptoneside points and usees distance from mean scaled by std compared to nsig to determine if the value should be replaced by the mean. if gapptsoneside>0, will do this leaving a gap around the point in question and using nptsoneside-gaps points for the mean and std
    if nptsoneside==1 and gapptsoneside==0:
        return removesinglepixoutliers(arr, critratiotoneighbors=nsig)
    nsig=max(nsig, 1.)
    nptsoneside=max(nptsoneside, 2.)
    gapptsoneside=min(gapptsoneside, nptsoneside-2.)
    for gap in range(gapptsoneside+1):
        starti=numpy.uint32([max(i-(nptsoneside-gap), 0) for i in range(len(arr))])
        stopi=numpy.uint32([min(i+(nptsoneside-gap)+1, len(arr)) for i in range(len(arr))])
        #print [numpy.append(arr[i0:i], arr[i+1:i1]) for i, i0, i1 in zip(range(len(arr)), starti, stopi)][8]
        #print [(((numpy.append(arr[i0:i], arr[i+1:i1]).mean()-arr[i]))**2, (numpy.append(arr[i0:i], arr[i+1:i1]).std()*nsig)**2) for i, i0, i1 in zip(range(len(arr)), starti, stopi)][8]
        arr=numpy.array([(((numpy.append(arr[i0:i], arr[i+1:i1]).mean()-arr[i]))**2<(numpy.append(arr[i0:i], arr[i+1:i1]).std()*nsig)**2 and (arr[i],) or (numpy.append(arr[i0:i], arr[i+1:i1]).mean(),))[0] for i, i0, i1 in zip(range(len(arr)), starti, stopi)], dtype=arr.dtype)
    return arr

def removesinglepixoutliers(arr,critratiotoneighbors=1.5):
    c=numpy.where((arr[1:-1]>(critratiotoneighbors*arr[:-2]))*(arr[1:-1]>(critratiotoneighbors*arr[2:])))
    c0=c[0]+1
    #print len(c0), ' pixels being replaced'
    arr[c0]=(arr[c0-1]+arr[c0+1])/2
    return arr
    
def findlocmax(arr, critval=0.):
    inds=numpy.where((arr[1:-1]>arr[:-2]) & (arr[1:-1]>arr[2:]) & (arr[1:-1]>critval))
    return inds[0]+1
    

def CalcR0_segdict(d, critrelstd=0.02, AveBeforeDivision=True, dzero=None):
    amps=d['samplecurrent']#2d arrays of ncycles x npts
    volts=d['samplevoltage']
    if not dzero is None:
        i0=dzero['samplecurrent'].mean()
        v0=dzero['samplevoltage'].mean()
        print 'zero baseline of %.2e A, %.2e V subtracted' %(i0, v0)
        amps=amps-i0
        volts=volts-v0
    if AveBeforeDivision:
        ro_cycles=numpy.float32([v.mean()/a.mean() for a, v in zip(amps, volts)])
    else:
        inds=numpy.where(amps<=0.)
        amps=replacevalswithneighsin2nddim(amps, inds)
        volts=replacevalswithneighsin2nddim(volts, inds)
        ro_cycles=(volts/amps).mean(axis=1)
    if (ro_cycles.std()/ro_cycles.mean())>critrelstd:
        print 'The relative variation in the measured Ro from the %d cycles is unexpectedly high: %.3f' %(amps.shape[0], ro_cycles.std()/ro_cycles.mean()) #this is a nice way to do number->string formatting, see the internet for details. amps.shape is a tuple where each element givews the length of the array inthat dimension so amps.shape[0] is the number of cycles
        newro_cycles=[r for r in ro_cycles if abs((r-ro_cycles.mean()))<ro_cycles.std()]# this is a filtering for statement. read from left to right starting with "for". We see that we will iterate over newro_cycles and r will be the value in each iteration but the command in the for loop is only going to execute provided the if statement is True. The if statement checks if r is within 1 std dev of the mean. The "command in the for loop" is just to place the r value so the result is a list that does not include the outliers
        newro_cycles=ro_cycles[numpy.abs((ro_cycles-ro_cycles.mean()))<ro_cycles.std()]#this line does the same thing as the previous line - this is the numpy-style way. The expression inside the "[ ]" is a boolean array that is the same shape as ro_cycles and "[ ]" means index the array so numpy only returns the values at indeces where there is a True in the boolean array
        if len(newro_cycles)>0 and len(newro_cycles)<len(ro_cycles):#since this is an "and" statement, the 2nd boolean will not execute if the 1st is False. We have the 1st boolean becuase if all the Ro values were "outliers" then len(newro_cycles) will be zero and we have failed to do filtering. The 2nd boolean checks if we filtered anything.
            print 'Ro values from %d cycles were determined to be outliers and removed' %(len(ro_cycles)-len(newro_cycles))
            return numpy.float32(newro_cycles).mean()
    else:
        return volts.mean()/amps.mean()
def tcr(r1, r2, t1, t2):#either all a
    if isinstance(r1, numpy.ndarray) or isinstance(r2, numpy.ndarray):
        if not isinstance(r1, numpy.ndarray):
            r1=numpy.ones(len(r2), dtype='float32')*r1
        if not isinstance(r2, numpy.ndarray):
            r2=numpy.ones(len(r1), dtype='float32')*r2
        ratio=numpy.float32([max(rv1, rv2)/min(rv1, rv2) for rv1, rv2 in zip(r1, r2)])
    else:
        ratio=max(r1, r2)/min(r1, r2)
    return (ratio-1.)/numpy.abs(t2-t1)

def makearr_cyc(x, calcarr):
    if isinstance(x, numpy.ndarray):
        if len(x)==calcarr.shape[0]:
            return x
        else:
            return numpy.array([x[0]]*calcarr.shape[0])
    else:
        return numpy.array([x]*calcarr.shape[0])

def temp_res(R, R0, T0, alpha):
    print '$%^', R0
    return numpy.array([(Rv/R0v-1.)/alphav+T0v for Rv, R0v, T0v, alphav in zip(R, makearr_cyc(R0, R), makearr_cyc(T0, R), makearr_cyc(alpha, R))])

def dT_IVdIdV(I, V, dI, dV, R0, alpha):
    return numpy.array([(Iv*dVv-Vv*dIv)/alphav/R0v/Iv**2 for Iv, Vv, dIv, dVv, R0v, alphav in zip(I, V, dI, dV, makearr_cyc(R0, I), makearr_cyc(alpha, I))])

def D_IVdIdV(I, V, dI, dV, R0, alpha):
    return numpy.array([Vv*Iv**3*R0v*alphav/(Iv*dVv-Vv*dIv) for Iv, Vv, dIv, dVv, R0v, alphav in zip(I, V, dI, dV, makearr_cyc(R0, I), makearr_cyc(alpha, I))])
    
def replacevalswithneighsin2nddim(arr, inds):
    iall, jall=inds
    for n in range(arr.shape[0]):
        j=jall[iall==n]
        if len(j)==0:
            continue
        jgood=numpy.int64([jv for jv in range(arr.shape[1]) if not jv in j])
        juse=numpy.int64([jgood[numpy.argmin((jgood-jv)**2)] for jv in j])
        arr[n, j]=arr[n, juse]
    return arr

def replacevalswithneighs(arr, inds):
    jgood=numpy.int64([jv for jv in range(arr.shape[0]) if not jv in inds])
    juse=numpy.int64([jgood[numpy.argmin((jgood-jv)**2)] for jv in inds])
    arr[inds]=arr[juse]
    return arr

def timepartition(cycletime, timepartitionfcn='timepart_none', piecelist=[1], yvals=[]):
    if timepartitionfcn=='timepart_none':
        return numpy.zeros(cycletime.shape, dtype='float64')
    elif timepartitionfcn=='timepart_user':
        idialog=timepartDialog(None, cycletime, numpieces=len(piecelist), yvals=yvals)
        idialog.exec_()
        return idialog.timepart
    elif timepartitionfcn=='timepart_peakid':
        return numpy.zeros(cycletime.shape, dtype='float64')#not implemented yet
    else:
        print 'ERROR - ABOUT TO ABORT BECAUSE timepartitionfcn IS NOT VALID'
        return

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

def lia_xy(x, ptspercyc, ncyclewin=1., nptswinstartinterval=1, phaseshift=0., extrappolyorder=1, interporder=3, interpplotax=None):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=x*sinarr(ptspercyc, x, ph=numpy.pi+phaseshift)
    c=x*sinarr(ptspercyc, x, ph=numpy.pi/2.+phaseshift)
    if nptswinstartinterval>1:
        starti=range(0, len(x)-npts, nptswinstartinterval)
    else:
        starti=numpy.arange(len(x)-npts)
    liax=(numpy.array([numpy.fft.fft(c[i:i+npts])[0].real for i in starti]))*2./npts
    liay=(numpy.array([numpy.fft.fft(s[i:i+npts])[0].real for i in starti]))*2./npts
    
    if nptswinstartinterval>1:
        liax=interpwithinarr(starti, liax, order=interporder, interpplotax=interpplotax, interpcols=['k', 'r'])
        liay=interpwithinarr(starti, liay, order=interporder, interpplotax=interpplotax, interpcols=['g', 'y'])
        
    nptsextrap=npts//2
    liax=concat_extrap_ends(liax, nptsextrap, highside=False, polyorder=extrappolyorder)
    liay=concat_extrap_ends(liay, nptsextrap, highside=False, polyorder=extrappolyorder)
    liax=concat_extrap_ends(liax, len(x)-len(liax), lowside=False, polyorder=extrappolyorder)
    liay=concat_extrap_ends(liay, len(x)-len(liay), lowside=False, polyorder=extrappolyorder)
    return liax, liay
    
def lia_ampphase2(x, ptspercyc, ncyclewin=1., returnphase=True, pad=True, phaseshift=0.):
    npts=numpy.round(ptspercyc*ncyclewin)
    s=x*sinarr(ptspercyc, x, ph=phaseshift)
    c=x*sinarr(ptspercyc, x, ph=numpy.pi/2.+phaseshift)
    amp=(numpy.array([(s.sum())**2+(c.sum())**2 for i in numpy.arange(len(x)-npts)])**.5)*2./npts
    if returnphase:
        phase=numpy.array([numpy.arctan(s.sum()/c.sum()) for i in numpy.arange(len(x)-npts)])
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
    symmetryfactor=2.*numpy.ones(npts_win//2+1, dtype='float64')
    symmetryfactor[0]=1.
    comp=numpy.array([numpy.fft.fft(x[i:i+npts_win])[:npts_win//2+1]*symmetryfactor for i in numpy.arange(len(x)-npts_win)])
    amp, phase=numpy.array([[numpy.abs(c), numpy.angle(c)] for c in comp]).swapaxes(0, 1)
    if pad:
        amp=numpy.concatenate([amp[:npts_win//2], amp, amp[-1*(len(x)-len(amp)-npts_win//2):]])
        phase=numpy.concatenate([phase[:npts_win//2], phase, phase[-1*(len(x)-len(phase)-npts_win//2):]])
    return amp/npts_win, phase #frequencies are in units of daqHz and are  numpy.array(range(npts_win//2+1))/npts_win.

    
def windowfft_xy(x, npts_win, pad=True, ptspercyc=None, nptswinstartinterval=1, extrappolyorder=1, interporder=3, interpplotax=None, freqinds=None):
    if freqinds is None:
        freqinds=range(npts_win//2+1)
    symmetryfactor=2.*numpy.ones(len(freqinds), dtype='float64')
    symmetryfactor[numpy.where(numpy.array(freqinds)==0)]=1.
    
    if nptswinstartinterval>1:
        starti=range(0, len(x)-npts_win, nptswinstartinterval)
    else:
        starti=numpy.arange(len(x)-npts_win)
    print x.shape
    comp=numpy.array([numpy.fft.fft(x[i:i+npts_win])[freqinds]*symmetryfactor for i in starti])
    print comp.shape
    fftx, ffty=numpy.array([[numpy.real(c), numpy.imag(c)] for c in comp]).swapaxes(0, 1).swapaxes(1, 2) #swapaxes makes X,Y first index and then frequencies  and then data points like that of x
    print fftx.shape
    fftx/=npts_win
    ffty/=npts_win
    
    if nptswinstartinterval>1:
        fftx=numpy.array([interpwithinarr(starti, a, order=interporder, interpplotax=interpplotax, interpcols=['k', 'r']) for a in fftx])
        ffty=numpy.array([interpwithinarr(starti, a, order=interporder, interpplotax=interpplotax, interpcols=['g', 'y']) for a in ffty])
        
    nptsextrap=npts_win//2
    print fftx.shape
    fftx=numpy.array([concat_extrap_ends(a, nptsextrap, highside=False, polyorder=extrappolyorder) for a in fftx])
    print nptsextrap,  fftx.shape
    ffty=numpy.array([concat_extrap_ends(a, nptsextrap, highside=False, polyorder=extrappolyorder) for a in ffty])
    fftx=numpy.array([concat_extrap_ends(a, len(x)-len(a), lowside=False, polyorder=extrappolyorder) for a in fftx])
    print   len(x), len(x)-len(a),  fftx.shape
    ffty=numpy.array([concat_extrap_ends(a, len(x)-len(a), lowside=False, polyorder=extrappolyorder) for a in ffty])
    
    return fftx.T, ffty.T #frequencies are second index and are in units of daqHz and are  numpy.array(range(npts_win//2+1))/npts_win, first index is same as x
    
def performgenericfilter(arr, filterdict):#filterdict can contains unused key:val but it must contain all those necessary for a given filter step to be performed
    iterateoverlastdim=(arr.ndim>=3)
    if iterateoverlastdim:
        if arr.ndim>3:
            origshape=arr.shape
            arr=arr.reshape(origshape[:2]+(numpy.prod(origshape[2:]),))
        else:
            origshape=None
        arrupdim=copy.copy(arr).swapaxes(1,2).swapaxes(0,1)
    else:
        arrupdim=numpy.array([copy.copy(arr)])
    fcn_parname_fkey_eachcycle=[\
    (removeoutliers_meanstd, ['nptsoneside', 'nsig', 'gapptsoneside'], ['OLnpts', 'OLnsig', 'OLgappts'], [], [], True), \
    (savgolsmooth, ['nptsoneside', 'order', 'deriv'], ['SGnpts', 'SGorder', 'SGderiv'], ['binprior'], ['SGbin'], True), \
    #(timeintegrate, ['integwindow_s'], ['integwindow_s'], True)\
#    (timepartition, ['timepartitionfcn', 'piecelist', 'yvals'], ['timepartitionfcn', 'fitpars', 'yvals'], False)
    ]
    for f, nl, kl, nlopt, klopt, eachcycle in fcn_parname_fkey_eachcycle:
        parlist=[((not k in filterdict) or filterdict[k] is None) or (n, filterdict[k]) for n, k in zip(nl, kl)] 
        if True in parlist:
            continue
        parlist+=[(n, filterdict[k]) for n, k in zip(nlopt, klopt) if (k in filterdict and not filterdict[k] is None)] 
        print 'executing filter function ', f.func_name, dict(parlist)
        for i in range(len(arrupdim)):
            #print arrupdim[i].shape
            if eachcycle:
                arrupdim[i, :, :]=numpy.array([f(a, **dict(parlist)) for a in arrupdim[i]])
            else:
                arrupdim[i, :, :]=f(arrupdim[i], **dict(parlist))
        
    #arr2=removeoutliers_meanstd(arr2, nptsoneside=filterdict['OLnpts'], nsig=filterdict['OLnsig'], gapptsoneside=filterdict['OLgappts'])
    #savgolsmooth(arr2, nptsoneside=filterdict['SGnpts'], order=filterdict['SGorder'], deriv=filterdict['SGderiv'])
    if iterateoverlastdim:
        if origshape is None:
            return arrupdim.swapaxes(0,1).swapaxes(1,2)
        else:
            return arrupdim.swapaxes(0,1).swapaxes(1,2).reshape(origshape)
    else:
        return arrupdim[0]

#def polyorder4_T(polycoefs, T):
#    return numpy.array([polycoefs[i]*(x**i) for i in range(5)]).sum(axis=0)
#
#def T4_intdT_2pieceC(pars, T, intdT, fixedpardict={'startind_2pieceC':numpy.uint16([])}):#startind_2pieceC si the index after which pars[1] will be used as the heat capacity
#    c=numpy.ones(T.shape, dtype='float64')*pars[0]
#    c[fixedpardict['startind_2pieceC']:]=pars[1]
#    return c+pars[2]*T**4+pars[3]*intdT
    
#HeatLossFunctionLibrary={\ #key:function,list of segdict keys, list of (fitpar name, dflt), dictionary of fixed parameters
#                   'polyorder4_T':(polyorder4_T, ['sampletemperature'], 'a+bT+cT^2+dT^3+eT^4', [('a', 1.e-6), ('a', 1.e-8), ('a', 1.e-10), ('a', 1.e-12), ('a', 1.e-14)], None)\
#                    'T4_intdT_2pieceC':(T4_intdT_2pieceC, [])
#                   }

#x=numpy.linspace(10., 20., 40)
#x[1]=0.
#x[18]=666.
#x[19]=66.
##x[10]=44.
#x[-1]=0.
#y=removeoutliers_meanstd(x, 6, 1.5, 2)
#print '***', x-y

def evaluatefitfcn(fitd, segd, interp_kind=1, extrap_order=1, interppars=False):#for interp_kind, see scipy.interpolate.interp1d, this is for interpolating through the excluded parts of the time axis. if interppars then only the piecewise parameters are interpolated and the function is directly evaluated with the fit parameters. if there is only 1 piece to the piecewise function, this is a roundabout way of just evlauating the fit fcn regardless of the exlcuded region
    f=FitFcnLibrary[fitd['fcnname']]
    pns=fitd['parnames']
    ks=fitd['segdkeys']
    ans=[f(**dict([('p', p)]+[(pn, segd[k][i]) for pn, k in zip(pns, ks) if pn in f.func_code.co_varnames[:f.func_code.co_argcount]])) for i, p in enumerate(fitd['fitpars'])]
    pt=[segd[k] for pn, k in zip(pns, ks) if k=='cyclepartition' and pn in f.func_code.co_varnames[:f.func_code.co_argcount]]
    if len(pt)==0 or numpy.all(pt[0]>=0):#if the is no paritioning in this function that was just called or if the parititoning does not exclude regions then we're done
        return ans
    pt=pt[0]
    for i, (t, a, p) in enumerate(zip(pt, ans, fitd['fitpars'])):#iterates over cycles
        iarr=numpy.int32(range(len(a)))
        inds=numpy.where(t>=0.)
        inds2=numpy.where((t<0.)&(iarr>inds[0][0])&(iarr<inds[0][-1]))#use interpolations in the middle
        if len(inds2[0])>0:
            if interppars:
                fcn=scipy.interpolate.interp1d(iarr[inds], numpy.float32([p[int(round(j))] for j in t[inds]]), kind=interp_kind)
                intpar=numpy.float32(fcn(iarr[inds2]))
                #                                      take fitpars, replace 0th with intpar          use the args in fitfcn except make a cyclepartition that will access the 0th fitpar                              iteration to get the parameters needed by the fitfcn from the segd                                                             due to the intpar value change at every index, new set of fitpar and thus args for fitfcn at every array index
                argtuplist_ind2=[[('p', tuple([ip]+list(p[1:])))]+[(pn, numpy.array(k=='cyclepartition' and (0,) or (segd[k][i][j],))) for pn, k in zip(pns, ks) if pn in f.func_code.co_varnames[:f.func_code.co_argcount]] for j, ip in zip(inds2[0], intpar)]
                a[inds2]=numpy.float32([f(**dict(tuplist)) for tuplist in argtuplist_ind2])
            else:
                fcn=scipy.interpolate.interp1d(iarr[inds], a[inds], kind=interp_kind)
                a[inds2]=numpy.float32(fcn(iarr[inds2]))
        inds2=numpy.where((t<0.)&(iarr<inds[0][0]))
        if len(inds2[0])>0:#use spline on the ends and only use as many data points as you need to extrapolate over. interpolations happen with respect to index so that it is always monotonic
            if interppars:
                #fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][:len(inds2[0])]], numpy.float32([p[int(round(j))] for j in t[inds[0][:len(inds2[0])]]]), k=extrap_order)
                #intpar=numpy.float32(fcn(iarr[inds2]))
                #Univeritae spline has a problem on 27Aug2011 that was not seen previously, so replace this extrapolation with filling in the end value of the parameter, 4 places below ($$$)
                #$$$
                intpar=numpy.zeros(len(inds2[0]), dtype='float32')+p[int(round(t[inds[0][0]]))]
                
                #                                      take fitpars, replace 0th with intpar          use the args in fitfcn except make a cyclepartition that will access the 0th fitpar                              iteration to get the parameters needed by the fitfcn from the segd                                                             due to the intpar value change at every index, new set of fitpar and thus args for fitfcn at every array index
                argtuplist_ind2=[[('p', tuple([ip]+list(p[1:])))]+[(pn, numpy.array(k=='cyclepartition' and (0,) or (segd[k][i][j],))) for pn, k in zip(pns, ks) if pn in f.func_code.co_varnames[:f.func_code.co_argcount]] for j, ip in zip(inds2[0], intpar)]
                a[inds2]=numpy.float32([f(**dict(tuplist)) for tuplist in argtuplist_ind2])
            else:
                #fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][:len(inds2[0])]], a[inds[0][:len(inds2[0])]], k=extrap_order)
                #intpar=numpy.float32(fcn(iarr[inds2]))
                #$$$
                intpar=numpy.zeros(len(inds2[0]), dtype='float32')+p[int(round(t[inds[0][0]]))]
                a[inds2]=intpar
        inds2=numpy.where((t<0.)&(iarr>inds[0][-1]))
        if len(inds2[0])>0:
            starti=len(inds[0])-len(inds2[0])#use only that last len inds2 indeces of the inds 
            starti=max(0, starti)
            if interppars:
                #fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][starti:]], numpy.float32([p[int(round(j))] for j in t[inds[0][starti:]]]), k=extrap_order)
                #intpar=numpy.float32(fcn(iarr[inds2]))
                #$$$
                intpar=numpy.zeros(len(inds2[0]), dtype='float32')+p[int(round(t[inds[0][-1]]))]
                #                                      take fitpars, replace 0th with intpar          use the args in fitfcn except make a cyclepartition that will access the 0th fitpar                              iteration to get the parameters needed by the fitfcn from the segd                                                             due to the intpar value change at every index, new set of fitpar and thus args for fitfcn at every array index
                argtuplist_ind2=[[('p', tuple([ip]+list(p[1:])))]+[(pn, numpy.array(k=='cyclepartition' and (0,) or (segd[k][i][j],))) for pn, k in zip(pns, ks) if pn in f.func_code.co_varnames[:f.func_code.co_argcount]] for j, ip in zip(inds2[0], intpar)]
                a[inds2]=numpy.float32([f(**dict(tuplist)) for tuplist in argtuplist_ind2])
            else:
                #fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][starti:]], a[inds[0][starti:]], k=extrap_order)
                #intpar=numpy.float32(fcn(iarr[inds2]))
                #$$$
                intpar=numpy.zeros(len(inds2[0]), dtype='float32')+p[int(round(t[inds[0][-1]]))]
                a[inds2]=intpar
    return ans
#def piecewise(p, c):
#    ans=numpy.float64([p[int(round(i))] for i in c])
#    if numpy.all(c>=0):#get out as soon as possible for the instances where this is being called during a fit
#        return ans
#    a=numpy.where(c>=0)[0]
#    b=numpy.where(c<0)[0]
#    for i in b:# if in the middle of a no-fit region of a piecewiese, replaces by the average of the nearest points. if the no-fit region extends to the end, replace wit the nearest
#        x=[]
#        j=a[a>i]
#        if len(j)>0:
#            x+=[ans[j[0]]]
#        j=a[a<i]
#        if len(j)>0:
#            x+=[ans[j[-1]]]
#        ans[i]=numpy.mean(x)
#    return ans

#piecewise functions are required, if don't want piecewise make sure cyclepartion conatins nothinng above. the piecewise fitpars must come first in the list and there can only be one piecewise parameter
FitFcnLibrary=dict([\
    ('FIT_T4', lambda p, c, T: numpy.float64([p[int(round(i))] for i in c])+numpy.array([v*(T**(i+1)) for i, v in enumerate(p[-4:])]).sum(axis=0)),\
    ('FIT_t4', lambda p, c, t: numpy.float64([p[int(round(i))] for i in c])+numpy.array([v*(t**(i+1)) for i, v in enumerate(p[-4:])]).sum(axis=0)),\
    ('FIT_t5', lambda p, c, t: numpy.float64([p[int(round(i))] for i in c])+numpy.array([v*(t**(i+1)) for i, v in enumerate(p[-5:])]).sum(axis=0)),\
    ('FIT_T0124', lambda p, c, T: numpy.float64([p[int(round(i))] for i in c])+numpy.array([v*(T**(i)) for i, v in zip([1, 2, 4], p[-3:])]).sum(axis=0)),\
    #('pieceC_T4_intdT', lambda p, c, T, dT: numpy.float64([p[int(round(i))] for i in c])+p[-2]*T+p[-1]*dT),\
    ])


def calcRofromheatprogram(I1, V1, I2, V2, RoToAl, o_R2poly=1):
    d={}
    d['P1']=I1*V1
    d['P2']=I2*V2
    ff=fitfcns()
    inds=numpy.where(numpy.array([I1])==0.)#array() is to make it 2-d for the replace fcn
    if len(inds[0])>0:
        I1=replacevalswithneighsin2nddim(numpy.array([I1]), inds)[0]
        V1=replacevalswithneighsin2nddim(numpy.array([V1]), inds)[0]
    inds=numpy.where(I2==0.)
    if len(inds[0])>0:
        I2=replacevalswithneighsin2nddim(numpy.array([I2]), inds)[0]
        V2=replacevalswithneighsin2nddim(numpy.array([V2]), inds)[0]
    #R1=V1/I1
    R2=V2/I2
    coefs=[R2[0]]+[(R2[-1]-R2[0])/len(R2)**i for i in range(1, o_R2poly+1)]
    x=numpy.arange(len(R2), dtype='float64')
    fcn=ff.polyfit((x, numpy.float64(R2)), coefs)
    d['R2fit']=fcn(x)
    R12=d['R2fit'][0]
    #dR12=ff.finalparams[1]
    d['delT2']=(d['R2fit'][-1]-d['R2fit'][0])/(RoToAl[0]*RoToAl[2])
    
    
#    coefs=[T2[0]]+[(T2[-1]-T2[0])/len(T2)**i for i in range(1, o_T2poly+1)]
#    x=numpy.arange(len(T2), dtype='float64')
#    fcn=ff.polyfit((x, numpy.float64(T2)), coefs)
#    d['T2fit']=fcn(x)
#    delT2=d['T2fit'][-1]-d['T2fit'][0]
    
    d['c']=d['P2'].sum()/d['delT2']
    d['delT1']=(d['c']*(d['P1'].sum()))
    Ro=R12/(1.+RoToAl[2]*d['delT1'])
    d['calc_Ro']=Ro
    return Ro, d



#Fs = 1/SAMPLINGPERIOD;                    % Sampling frequency
# Fsmax=250e3; %maximum sampling rate: difference between neighbor I and V
#T = 1/Fs;                     % Sample time
#Ttotal=(length(U)-1)*T;            %total experiment time
#fq_base=Fs/POINTSPERCYCLE; 
#
#phi0(j)=-angle(i1(j));  % phi0
#
#Xadd2(j)=2*abs(i1(j))*R0(j)*lam*yita(j)/6/pi/fq_base*sin(phi0(j));
#Yadd2(j)=(3*abs(i0(j))*R0(j)+4*abs(i1(j))*R0(j)*cos(phi0(j)))*lam*yita(j)/3/2/pi/fq_base;
#Fv2(j)=v2(j)-Xadd2(j)-1i*Yadd2(j)-i2(j)*R0(j);
#
#mc2(j)=lam*abs(i1(j))^2*i0(j)*R0(j)^2*3/2/abs(Fv2(j))/2/pi/fq_base;
#
#Xadd3(j)=abs(i1(j))*R0(j)*lam*yita(j)/4/2/pi/fq_base*sin(phi0(j));
#Yadd3(j)=(8*abs(i0(j))*R0(j)+9*abs(i1(j))*R0(j)*cos(phi0(j)))*lam*yita(j)/12/2/pi/fq_base;
#Fv3(j)=v3f(j)-Xadd3(j)-1i*Yadd3(j);%-i3(j)*R0(j);
#mc3(j)=lam*abs(i1(j))^3*R0(j)^2/8/abs(Fv3(j))/2/pi/fq_base; 

def mCp_2w(VhX, VhY, I0X, I0Y, I1X, I1Y, IhX, IhY, dT, R, tcr, freq1w, Vsrcisfiltered=True,  returnall=False):
    V2=VhX+1j*VhY
    I1=I1X+1j*I1Y
    I0amp=numpy.sqrt(I0X**2+I0Y**2)
    angfreq1w=2.*numpy.pi*freq1w
    if Vsrcisfiltered:
        F2=V2
    else:
        I2=IhX+1j*IhY
        phi0=-1.*numpy.angle(I1)
        Xadd2=2.*numpy.abs(I1)*R*tcr*dT/3./angfreq1w*numpy.sin(phi0)
        Yadd2=(3.*I0amp*R+4.*numpy.abs(I1)*R*numpy.cos(phi0))*tcr*dT/3./angfreq1w
        F2=V2-Xadd2-1j*Yadd2-I2*R
    mc=tcr*numpy.abs(I1)**2*I0amp*R**2*1.5/numpy.abs(F2)/angfreq1w
    if returnall:
        if Vsrcisfiltered:
            return angfreq1w, I0amp, I1, V2, F2, mc
        else:
            return angfreq1w, phi0, I0amp, I1, I2, V2, Xadd2, Yadd2, F2, mc
    else:
        return mc
    
def mCp_3w(VhX, VhY, I0X, I0Y, I1X, I1Y, IhX, IhY, dT, R, tcr, freq1w, Vsrcisfiltered=True, returnall=False):
    V3=VhX+1j*VhY
    I1=I1X+1j*I1Y
    angfreq1w=2.*numpy.pi*freq1w
    
    if Vsrcisfiltered:
        F3=V3
    else:
        I3=IhX+1j*IhY
        I0amp=numpy.sqrt(I0X**2+I0Y**2)
        phi0=-1.*numpy.angle(I1)
        Xadd3=numpy.abs(I1)*R*tcr*dT/4./angfreq1w*numpy.sin(phi0)
        Yadd3=(8.*I0amp*R+9.*numpy.abs(I1)*R*numpy.cos(phi0))*tcr*dT/12./angfreq1w
        F3=V3-Xadd3-1j*Yadd3-I3*R

    mc=tcr*numpy.abs(I1)**3*R**2/8./numpy.abs(F3)/angfreq1w
    if returnall:
        if Vsrcisfiltered:
            return angfreq1w, I1, V3, F3, mc
        else:
            return angfreq1w, phi0, I0amp, I1, I3, V3, Xadd3, Yadd3, F3, mc
    else:
        return mc


