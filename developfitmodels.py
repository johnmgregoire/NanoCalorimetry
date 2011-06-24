import time
import os
import sys
import numpy
import h5py
from scipy.special import erf as erf
from numpy import exp as exp
from numpy import pi as pi
from PnSC_ui import *
from PnSC_dataimport import *
from PnSC_SCui import *
from PnSC_math import *
from PnSC_h5io import *
from PnSC_main import *


def b5(t, b):
    coef1=[[1.], [.5, 1.], [-.75, 1., 1.], [15./8., -1.75, 1.5, 1.], [-105./16., 5., -3., 2., 1.]]
    coef2=[[.5, 1.], [-.75, 1., 1.], [.375, -.75, 1.5, 1.], [-15./16., 1.5, -1.5, 2., 1.], [105./32., -75./16., 15./4., -2.5, 2.5, 1.]]
    f1=numpy.float32([numpy.float32([c*t**(2.*j+1)*b**(j-i) for j, c in enumerate(l)]).sum() for i, l in enumerate(coef1)])
    f2=numpy.float32([numpy.float32([c*t**(j)*b**(j-i-.5) for j, c in enumerate(l)]).sum() for i, l in enumerate(coef2)])
    return f1*exp(-1.*t*b)/pi**.5+f2*erf((b*t)**.5)

def b52(t, b):#t needs to be an array but b should be a float
    coef1=[[1.], [.5, 1.], [-.75, 1., 1.], [15./8., -1.75, 1.5, 1.], [-105./16., 5., -3., 2., 1.]]
    coef2=[[.5, 1.], [-.75, 1., 1.], [.375, -.75, 1.5, 1.], [-15./16., 1.5, -1.5, 2., 1.], [105./32., -75./16., 15./4., -2.5, 2.5, 1.]]
    f1=numpy.float32([(exp(-1.*t*b)/pi**.5)*numpy.float32([c*t**(2.*j+1)*b**(j-i) for j, c in enumerate(l)]).sum(axis=0) for i, l in enumerate(coef1)])
    f2=numpy.float32([erf((b*t)**.5)*numpy.float32([c*t**(j)*b**(j-i-.5) for j, c in enumerate(l)]).sum(axis=0) for i, l in enumerate(coef2)])
    return f1+f2
    
def Dmodel(t, T_C, S, c, Thcoefs, condlist_coefbeta, radcoef, To_C):#dT should already have 0 removed
    if isinstance(c, float):
        c=numpy.ones(T_C.shape, dtype='float32')*c
    tapb=numpy.float32([tap*numpy.float32([b5(tv, beta) for tv in t]) for tap, beta in condlist_coefbeta]).sum(axis=0).T
    condterms=numpy.float32([a*tb for a, tb in zip(Thcoefs, tapb)]).sum(axis=0)
    return c+(condterms+radcoef*((T_C+273.15)**4-(To_C+273.15)**4))/S

def Dmodel2(t, T_C, S, c, Thcoefs, condlist_coefbeta, radcoef, To_C):#dT should already have 0 removed
    if isinstance(c, float):
        c=numpy.ones(T_C.shape, dtype='float32')*c
    tapb=numpy.float32([tap*b52(t, beta) for tap, beta in condlist_coefbeta]).sum(axis=0)
    condterms=numpy.float32([a*tb for a, tb in zip(Thcoefs, tapb)]).sum(axis=0)
    return c+(condterms+radcoef*((T_C+273.15)**4-(To_C+273.15)**4))/S
    
p='E:/CHESS2010PnSC/AuSiCu_pnsc_all.h5'
#p=mm.h5path
#f=h5py.File(p, mode='r')
segdlist=CreateHeatProgSegDictList(p, 'heat6_cell21_other', '2010Nov27_Cell21_0_90_5_0_mA_100_70_280_100_ms_1C_0_ mT')
d=segdlist[4]
S=d['sampleheatrate'][0]
T=d['sampletemperature'][0]
t=d['cycletime'][0]
D=d['samplepowerperrate'][0]
To=20.


ff=fitfcns()
fitfcn_T_t=ff.polyfit((t, T), 5)
fitpars_T_t=ff.finalparams
fitTo=fitpars_T_t[0]
fitpars_T_t=fitpars_T_t[1:]
#temp=time.time()
#Dm=Dmodel(t, T, S, .3, fitpars_T_t, [(1.E-2, .1), (3.E-2, .03)], 1.E-6, To)
#print 'Dmodeltime:', time.time()-temp
#
#temp=time.time()
#Dm2=Dmodel2(t, T, S, .3, fitpars_T_t, [(1.E-2, .1), (3.E-2, .03)], 1.E-6, To)
#print 'Dmodel2time:', time.time()-temp    


def Dfitfcn_2cond(p, t, T, S):#p is c, al,beta,al,beta,radcoef
    return Dmodel2(t, T, S, p[0], fitpars_T_t, [(p[1], p[2])], p[3], To)

ff=fitfcns()
fitfcn_T_t=ff.genfit(Dfitfcn_2cond, (1.e-5, .01, .0001, 1.e-10), (t, T, S, D))
fitp=ff.finalparams

#T_C=T
#c=.3
#Thcoefs=fitpars_T_t
#condlist_coefbeta=[(1.E-2, .1), (3.E-2, .03)]
#radcoef=1.E-6
#To_C=To
#
#
#if isinstance(c, float):
#    c=numpy.ones(T_C.shape, dtype='float32')*c
#tapb=numpy.float32([tap*numpy.float32([b5(tv, beta) for tv in t]) for tap, beta in condlist_coefbeta]).sum(axis=0).T
#condterms=numpy.float32([a*tb for a, tb in zip(Thcoefs, tapb)]).sum(axis=0)
#ans=c+condterms+radcoef*((T_C+273.15)**4-(To_C+273.15)**4)
#
#c=.3
#if isinstance(c, float):
#    c2=numpy.ones(T_C.shape, dtype='float32')*c
#tapb2=numpy.float32([tap*b52(t, beta) for tap, beta in condlist_coefbeta]).sum(axis=0)
#condterms2=numpy.float32([a*tb for a, tb in zip(Thcoefs, tapb)]).sum(axis=0)
#ans2=c2+condterms2+radcoef*((T_C+273.15)**4-(To_C+273.15)**4)
    
print 'calc done'
