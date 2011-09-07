import time
import os
import sys
import numpy
import h5py
#from PnSC_ui import *
#from PnSC_dataimport import *
from PnSC_SCui import *
#from PnSC_math import *
from PnSC_h5io import *
from PnSC_main import *


savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/checkRo'

#name, (cycles,daqHz)
exp_rec=[\
#('heat1a',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1b',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1d',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1c',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat2',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat2_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat3',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat3_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4a',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4b',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4c',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (3, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3')]),\
#('heat6a',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat6b',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat6c',[(2, 'IVOL8_RTPDSG10o1_SdtSG15o1'), (3, 'IVOL8_RTPDSG10o1_SdtSG15o1')]),\
#('heat6_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat7',[(2, 'IVOL8_RTPDSG10o1_SdtSG15o1'), (3, 'IVOL8_RTPDSG10o1_SdtSG15o1')]),\
#('heat8',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
]
p=mm.h5path
f=h5py.File(p, mode='r')
ff=fitfcns()

TTorange=[0., 10.]

#pd=simpleplotDialog(None, numpy.array([0, 0]))
#pd.show()
to=.3*100
twin=1.*100.#1ms 100 points per ms
for exp, recipelist in exp_rec:
    f=h5py.File(p, mode='r')
    h5pl=experimenthppaths(f, exp)
    temp=numpy.zeros(25, dtype='float32')
    RTC=numpy.float32([temp, temp+20.,  temp+.001]).T
    Rlinfit=numpy.float32([temp, temp]).T
    for h5p in h5pl:
        hp=h5p.strip('/').rpartition('/')[2]
        print exp, hp
        Ro, To, al=RoToAl_h5(p, exp, hp)
        segdl=CreateHeatProgSegDictList(p, exp, hp)
        xs=segdl[2]['cycletime'][0, :to]
        ys=segdl[2]['sampleresistance'][0, :to]
        x=segdl[2]['cycletime'][0, to:twin]
        y=segdl[2]['sampleresistance'][0, to:twin]
        b=(y[-1]-y[0])/(x[-1]-x[0])
        a=y[0]-b*x[0]
        x1=segdl[1]['cycletime'][0, 0]
        x2=segdl[2]['cycletime'][0,0]
        fcn=ff.polyfit((x, y), [a, b])
        a=ff.finalparams[0]
        b=ff.finalparams[1]
        R2=fcn(x2)
        R1=R2-.5*b*(x1**2/(x2-x1))-b*((.5*x2**2-x2*x1)/(x2-x1))
        x12=numpy.linspace(x1, x2, 100.)
        R12=R2+b*((.5*x12**2-x1*x12)/(x2-x1))-b*((.5*x2**2-x2*x1)/(x2-x1))
        Romin=R2/(TTorange[1]*al+1.)
        Romax=R2/(TTorange[0]*al+1.)
        pd=simpleplotDialog(None, numpy.array([0, 0]))
        pd.plotw.axes.cla()
        pd.plotw.axes.plot(xs, ys, 'c.', markersize=2)
        pd.plotw.axes.plot(x, y, 'b.', markersize=2)
        pd.plotw.axes.plot(x, fcn(x), 'g-', lw=2)
        pd.plotw.axes.plot(xs, fcn(xs), 'g-', lw=1)
        pd.plotw.axes.plot(x12, R12, 'g-', lw=1)
        pd.plotw.axes.plot([x1, x2], [R1, R2], 'ro', markersize=4)
        pd.plotw.axes.plot([x2, x2], [Romin, Romax], 'r-', lw=2)
        pd.plotw.axes.plot([x1, x2], [Ro, Ro], 'k-', markersize=4)
        pd.plotw.axes.set_ylim(min(Romin, Ro)*.99, y.max())
        f.close()
        f=h5py.File(p,mode='r')
        c=f[h5p].attrs['CELLNUMBER']
        RTC[c-1, 0]=R1
        Rlinfit[c-1, 0]=a
        Rlinfit[c-1, 1]=b
        pd.plotw.axes.set_xlabel('cycle time (s)')
        pd.plotw.axes.set_ylabel('R (Ohms)')
        pd.plotw.axes.set_title('Rocheck_%s_cell%02d' %(exp, c))
        pd.plotw.fig.canvas.draw()
        pd.plotw.fig.savefig(os.path.join(savef,'Rocheck_%s_cell%02d' %(exp, c)+'.png'))
        if not pd.exec_():
#        idialog=messageDialog(title='continue')
#        if not idialog.exec_():
            RTC=None
            break
            
        #pd.close()
        
    f.close()
    #RTC[0, 0]=10.511
    if RTC is None:
        break
    else:
        print 'saving ', exp
        f=h5py.File(p,mode='r+')
        g=f['Calorimetry/ResCal']
        if exp in g:
            del g[exp]
        g=g.create_group(exp)
        ds=g.create_dataset('Res_TempCal', data=RTC)
        f['Calorimetry/'+exp].attrs['Res_TempCalPath']=ds.name
        g.create_dataset('Rlinfit', data=Rlinfit)
        f.close()

f.close()

print 'batch done'
