import numpy
import h5py
import os, os.path, time, copy

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
    
def performgenericfilter(arr, filterdict):#filterdict can contains unused key:val but it must contain all those necessary for a given filter step to be performed
    arr=copy.copy(arr)
    fcn_parname_fkey=[\
    (removeoutliers_meanstd, ['nptsoneside', 'nsig', 'gapptsoneside'], ['OLnpts', 'OLnsig', 'OLgappts']), \
    (savgolsmooth, ['nptsoneside', 'order', 'deriv'], ['SGnpts', 'SGorder', 'SGderiv']), \
    ]
    for f, nl, kl in fcn_parname_fkey:
        parlist=[((not k in filterdict) or filterdict[k] is None) or (n, filterdict[k]) for n, k in zip(nl, kl)] 
        if True in parlist:
            continue
        print 'executing filter function ', f.func_name, dict(parlist)
        print arr.shape
        arr=numpy.array([f(a, **dict(parlist)) for a in arr])
    #arr=removeoutliers_meanstd(arr, nptsoneside=filterdict['OLnpts'], nsig=filterdict['OLnsig'], gapptsoneside=filterdict['OLgappts'])
    #savgolsmooth(arr, nptsoneside=filterdict['SGnpts'], order=filterdict['SGorder'], deriv=filterdict['SGderiv'])
    return arr

#x=numpy.linspace(10., 20., 40)
#x[1]=0.
#x[18]=666.
#x[19]=66.
##x[10]=44.
#x[-1]=0.
#y=removeoutliers_meanstd(x, 6, 1.5, 2)
#print '***', x-y
