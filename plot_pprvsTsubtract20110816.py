import numpy, h5py, os
from PnSC_main import *
from PnSC_h5io import *
from PnSC_math import *

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110816_Zr-Hf-B.h5'
h5f=h5py.File(p, mode='r')

ehl=[\
('quadlinheating2', 'pre_25mApluslinquad2_cell17_1_of_1', 'Zr-B, 1st'),\
('quadlinheating2', 'cell17_25malinquad2repeat_1_of_1', 'Zr-B, 2nd'),\
('quadlinheating2', 'cell17_25mAquadlincurrent_3rdtime_1_of_1', 'Zr-B, 3rd'), \
#('quadlinheating2', 'pre_25mApluslinquad2_cell16_1_of_1', 'Hf-B, nth'), \
#('quadlinheating2', 'cell11_25malinquad2_1_of_1', 'empty'), \
]
tarrs=[]
pprarrs=[]
for i, (e, h, l) in enumerate(ehl):
    hpsdl=CreateHeatProgSegDictList(p, e, h)
    T=hpsdl[2]['sampletemperature'][0, :]
    ppr=hpsdl[2]['samplepowerperrate'][0, :]
    if 0:
        pylab.plot(T, ppr*1.e6, label=l)
        pylab.xlabel('Temperature (C)')
        pylab.ylabel('power per rate ($\mu$J/K)')
        pylab.legend(loc=0)
    tarrs+=[T]
    pprarrs+=[ppr]

def extremesmooth(x, binn=70, SGpts=170, SGorder=3):
    xb=numpy.array([x[i*binn:(i+1)*binn].mean() for i in range(len(x)//binn)])
    xbf=savgolsmooth(xb, nptsoneside=SGpts, order =SGorder)
    ia=numpy.arange(binn, dtype='float32')/binn
    xr=numpy.concatenate([ia*(b-a)+b for a, b in zip(xbf[:-1], xbf[1:])])
    xr=numpy.concatenate([(xbf[1]-xbf[0])*ia[:binn//2]+xbf[0], xr, (xbf[-1]-xbf[-2])*ia[:binn//2]+xbf[-1]])
    xr=numpy.concatenate([xr, (xbf[-1]-xbf[-2])*ia[:len(x)-len(xr)]+xbf[-1]])
    return xr

if 1:
    x=extremesmooth(pprarrs[0])
    y=extremesmooth(pprarrs[1])
    z=extremesmooth(pprarrs[2])
    xt=tarrs[0]
    yt=tarrs[1]
    zt=tarrs[2]

        
    tmin=max([t.min() for t in [xt, yt, zt]])
    tmax=min([t.max() for t in [xt, yt, zt]])
    tinterp=numpy.linspace(tmin, tmax, 2000)

    xinterp=numpy.interp(tinterp, xt, x)
    yinterp=numpy.interp(tinterp, yt, y)
    zinterp=numpy.interp(tinterp, zt, z)

    pylab.figure()
    for i, (t, a, ai) in enumerate([(xt, x, xinterp), (yt, y, yinterp), (zt, z, zinterp)]):
        pylab.subplot(3, 1, i+1)
        pylab.plot(tinterp, ai)
        pylab.plot(t, a)
    pylab.figure()
    xsub=xinterp-(zinterp+yinterp)/2.
    
    for i, (a, l) in enumerate([(xinterp, '1st'), ((zinterp+yinterp)/2., 'subsequent')]):
        pylab.plot(tinterp, a*1.e6, label=l, lw=2)
    #pylab.legend(loc=2)
    pylab.xlabel('Temperature (C)', fontsize=14)
    pylab.ylabel('Calorimetric Signal ($\mu$J/K)', fontsize=14)
#    pylab.text(700, 14, '1st',color='b', ha='left', fontsize=14)
#    pylab.text(450, 14, 'subsequent',color='g', ha='right', fontsize=14)
    pylab.annotate('1st',(610, 14),xytext=(700, 14),fontsize=14,color='b',arrowprops={'arrowstyle':'->','color':'b'})
    pylab.annotate('subsequent',(555, 14),xytext=(450, 14),fontsize=14,color='g',arrowprops={'arrowstyle':'->','color':'g'}, ha='right')
    pylab.xlim(0, 1200)
    pylab.figure()    
    pylab.plot([0, 1200], [0, 0], 'k', lw=1)
    pylab.plot(tinterp, xsub*1.e6, 'r-', lw=2)
#    pylab.annotate(' ',(510, -2),xytext=(510, 0),color='k',arrowprops={'arrowstyle':'simple','color':'k'})
#    pylab.annotate(' ',(1010, -14),xytext=(1010, 0),color='k',arrowprops={'arrowstyle':'simple','color':'k'})
    #pylab.legend()
    pylab.xlabel('Temperature (C)', fontsize=14)
    pylab.ylabel('Differential signal ($\mu$J/K)', fontsize=14)
    pylab.xlim(0, 1200)
    pylab.subplots_adjust(right=.55, top=.5)
    print xsub[(tinterp>260)*(tinterp<670)].sum()*(tinterp[1]-tinterp[0])*1.e6
    print xsub[tinterp>670].sum()*(tinterp[1]-tinterp[0])*1.e6
pylab.show()

