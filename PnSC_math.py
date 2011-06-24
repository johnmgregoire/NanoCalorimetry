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
    
def savgolsmooth(x, nptsoneside=7, order = 4, dx=1.0, deriv=0): #based on scipy cookbook. x is 1-d array, window is the number of points used to smooth the data, order is the order of the smoothing polynomial, will return the smoothed "deriv"th derivative of x
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

def temp_res(R, R0, T0, alpha):
    return (R/R0-1.)/alpha+T0

def dT_IVdIdV(I, V, dI, dV, R0, alpha):
    return (I*dV-V*dI)/alpha/R0/I**2

def D_IVdIdV(I, V, dI, dV, R0, alpha):
    return V*I**3*R0*alpha/(I*dV-V*dI)
    
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



def performgenericfilter(arr, filterdict):#filterdict can contains unused key:val but it must contain all those necessary for a given filter step to be performed
    arr=copy.copy(arr)
    fcn_parname_fkey_eachcycle=[\
    (removeoutliers_meanstd, ['nptsoneside', 'nsig', 'gapptsoneside'], ['OLnpts', 'OLnsig', 'OLgappts'], True), \
    (savgolsmooth, ['nptsoneside', 'order', 'deriv'], ['SGnpts', 'SGorder', 'SGderiv'], True), \
    #(timeintegrate, ['integwindow_s'], ['integwindow_s'], True)\
#    (timepartition, ['timepartitionfcn', 'piecelist', 'yvals'], ['timepartitionfcn', 'fitpars', 'yvals'], False)
    ]
    for f, nl, kl, eachcycle in fcn_parname_fkey_eachcycle:
        parlist=[((not k in filterdict) or filterdict[k] is None) or (n, filterdict[k]) for n, k in zip(nl, kl)] 
        if True in parlist:
            continue
        print 'executing filter function ', f.func_name, dict(parlist)
        print arr.shape
        if eachcycle:
            arr=numpy.array([f(a, **dict(parlist)) for a in arr])
        else:
            arr=f(arr, **dict(parlist))
    #arr=removeoutliers_meanstd(arr, nptsoneside=filterdict['OLnpts'], nsig=filterdict['OLnsig'], gapptsoneside=filterdict['OLgappts'])
    #savgolsmooth(arr, nptsoneside=filterdict['SGnpts'], order=filterdict['SGorder'], deriv=filterdict['SGderiv'])
    return arr

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
                fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][:len(inds2[0])]], numpy.float32([p[int(round(j))] for j in t[inds[0][:len(inds2[0])]]]), k=extrap_order)
                intpar=numpy.float32(fcn(iarr[inds2]))
                #                                      take fitpars, replace 0th with intpar          use the args in fitfcn except make a cyclepartition that will access the 0th fitpar                              iteration to get the parameters needed by the fitfcn from the segd                                                             due to the intpar value change at every index, new set of fitpar and thus args for fitfcn at every array index
                argtuplist_ind2=[[('p', tuple([ip]+list(p[1:])))]+[(pn, numpy.array(k=='cyclepartition' and (0,) or (segd[k][i][j],))) for pn, k in zip(pns, ks) if pn in f.func_code.co_varnames[:f.func_code.co_argcount]] for j, ip in zip(inds2[0], intpar)]
                a[inds2]=numpy.float32([f(**dict(tuplist)) for tuplist in argtuplist_ind2])
            else:
                fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][:len(inds2[0])]], a[inds[0][:len(inds2[0])]], k=extrap_order)
                a[inds2]=numpy.float32(fcn(iarr[inds2]))
        inds2=numpy.where((t<0.)&(iarr>inds[0][-1]))
        if len(inds2[0])>0:
            starti=len(inds[0])-len(inds2[0])#use only that last len inds2 indeces of the inds 
            starti=max(0, starti)
            if interppars:
                fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][starti:]], numpy.float32([p[int(round(j))] for j in t[inds[0][starti:]]]), k=extrap_order)
                intpar=numpy.float32(fcn(iarr[inds2]))
                #                                      take fitpars, replace 0th with intpar          use the args in fitfcn except make a cyclepartition that will access the 0th fitpar                              iteration to get the parameters needed by the fitfcn from the segd                                                             due to the intpar value change at every index, new set of fitpar and thus args for fitfcn at every array index
                argtuplist_ind2=[[('p', tuple([ip]+list(p[1:])))]+[(pn, numpy.array(k=='cyclepartition' and (0,) or (segd[k][i][j],))) for pn, k in zip(pns, ks) if pn in f.func_code.co_varnames[:f.func_code.co_argcount]] for j, ip in zip(inds2[0], intpar)]
                a[inds2]=numpy.float32([f(**dict(tuplist)) for tuplist in argtuplist_ind2])
            else:
                fcn=scipy.interpolate.UnivariateSpline(iarr[inds[0][starti:]], a[inds[0][starti:]], k=extrap_order)
                a[inds2]=numpy.float32(fcn(iarr[inds2]))
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
