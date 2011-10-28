import pylab 
import matplotlib.cm as cm
import numpy
import h5py
import os, os.path, time, copy, operator
import PnSC_ui
#from PnSC_ui import *
#from PnSC_SCui import *
#from PnSC_dataimport import *
from PnSC_math import *

v=numpy.linspace(-16,16,5)
xmm=numpy.float32([x for i in range(5) for x in v])
ymm=numpy.float32([x for x in v[::-1] for i in range(5)])

def myeval(c):
    if c=='None':
        c=None
    else:
        temp=c.lstrip('0')
        if (temp=='' or temp=='.') and '0' in c:
            c=0
        else:
            c=eval(temp)
    return c
def createh5file(h5path):
    if os.path.exists(h5path):
        mode='r+'
    else:
        mode='w'
    h5file=h5py.File(h5path, mode=mode)
    h5file.attrs['xmm']=xmm
    h5file.attrs['ymm']=ymm
    h5file.attrs['cells']=numpy.int32(range(1, 26))
    #node=h5file[h5groupstr]
    gstrlist=['Calorimetry']

    for gs in gstrlist:
        if not gs in h5file:
            h5file.create_group(gs)

    h5file.close()

def readh5pyarray(arrpoint):
    return eval('arrpoint'+('['+':,'*len(arrpoint.shape))[:-1]+']')
    
def create_exp_grp(h5pf, h5expname):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r+')
    else:
        h5file=h5pf
    if 'Calorimetry' in h5file:
        h5cal=h5file['Calorimetry']
    else:
        h5cal=h5file.create_group('Calorimetry')
    if h5expname in h5cal:
       del h5cal[h5expname] 
    h5exp=h5cal.create_group(h5expname)
    h5a=h5exp.create_group('analysis')
    h5m=h5exp.create_group('measurement')
    h5hp=h5m.create_group('HeatProgram')
    if isinstance(h5pf, str):
        h5file.close()
    
def writenewh5heatprogram(h5path, h5expname, grpname, AttrDict, DataSetDict, SegmentData):
    h5file=h5py.File(h5path, mode='r+')
    h5hp=h5file['Calorimetry'][h5expname]['measurement']['HeatProgram']
    if grpname in h5hp:
        del h5hp[grpname]
    h5hpg=h5hp.create_group(grpname)
    h5hpg.attrs['segment_ms']=SegmentData[0]
    h5hpg.attrs['segment_mA']=SegmentData[1]
    for k, v in AttrDict.iteritems():
        #print k, type(v), v
        h5hpg.attrs[k]=v
    for nam, (adict, data) in DataSetDict.iteritems():
        h5d=h5hpg.create_dataset(nam, data=data)
        for k, v in adict.iteritems():
            h5d.attrs[k]=v
    
    h5file.close()

def geth5attrs(h5pf, h5grppath):#closes the h5file, doesnt return it
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    d={}
    for k, v in h5file[h5grppath].attrs.iteritems():
        d[k]=(v=='None' and (None,) or (v,))[0]
    if isinstance(h5pf, str):
        h5file.close()
    return d
    
def getindex_cell(h5pf, cellnum):#closes the h5file, doesnt return it
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    #i=numpy.where(h5file.attrs['cells']==cellnum)[0][0]
    i=list(h5file.attrs['cells']).index(cellnum)
    if isinstance(h5pf, str):
        h5file.close()
    return i
    
def getcalanalysis(h5pf, h5expname):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5hp=h5file['Calorimetry'][h5expname]['analysis']
    if isinstance(h5pf, str):
        return h5file, h5hp
    return h5hp

def dt_h5(h5path, h5expname, h5hpname):
    h5file=h5py.File(h5path, mode='r')
    h5hp=h5file['Calorimetry'][h5expname]['measurement']['HeatProgram'][h5hpname]
    dt=1./h5hp.attrs['daqHz']
    h5file.close()
    return dt

def atm_h5(h5path, h5expname, h5hpname):
    h5file=h5py.File(h5path, mode='r')
    h5hp=h5file['Calorimetry'][h5expname]['measurement']['HeatProgram'][h5hpname]
    atm=h5hp.attrs['ambient_atmosphere']
    h5file.close()
    return atm
    
def pts_sincycle_h5(h5path, h5expname, h5hpname):
    h5file=h5py.File(h5path, mode='r')
    h5hp=h5file['Calorimetry'][h5expname]['measurement']['HeatProgram'][h5hpname]
    n=h5hp.attrs['pts_sincycle']
    h5file.close()
    return n
def gethpgroup(h5pf, h5expname, h5hpname=None):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5hp=h5file['Calorimetry'][h5expname]['measurement']['HeatProgram']
    if not h5hpname is None:
        h5hp=h5hp[h5hpname]
    if isinstance(h5pf, str):
        return h5file, h5hp
    return h5hp

def assign_segmsma(h5path, h5expname, segms, segmA):
    h5file=h5py.File(h5path, mode='r+')
    h5hp=gethpgroup(h5file, h5expname)
    for node in h5hp.values():
        if 'segment_ms' in node.attrs.keys():
            node.attrs['segment_ms']=segms
            node.attrs['segment_mA']=segmA
    h5file.close()

def experimenthppaths(h5pf, h5expname):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    p=[]
    h5hp=gethpgroup(h5file, h5expname)
    for pnt in h5hp.values():
        if isinstance(pnt, h5py.Group):
            p+=[pnt.name]
    if isinstance(h5pf, str):
        return h5file, p
    return p


def msarr_hpgrp(h5hpgrp, twod=False):
    for pnt in h5hpgrp.values():
        if isinstance(pnt, h5py.Dataset):
            ms=numpy.linspace(0., pnt.shape[1]-1., pnt.shape[1])
            arrshape=pnt.shape
            break
    ms/=(h5hpgrp.attrs['daqHz']/1000.)
    if twod:
        ms=numpy.float32([ms]*arrshape[0])
    return numpy.float32(ms)

def segtypes():
    return ['step', 'soak', 'ramp', 'zero']
def CreateHeatProgSegDictList(h5path, h5expname, h5hpname, critms_step=1., critmAperms_constmA=0.01, critdelmA_constmA=10., critmA_zero=0.1):
#the segment types are step, soak, ramp, zero
#the CreateHeatProgSegDictList function reads the data from the .h5 file and organizes in a way that will be useful for many types of analysis. 
#the function returns a list where there is one dict for each segment in the heat program. Each dict value that is an array is assumed to be data and all have the same shape
    h5file, h5hpgrp=gethpgroup(h5path, h5expname, h5hpname)
    ms=msarr_hpgrp(h5hpgrp, twod=True)
    segms=h5hpgrp.attrs['segment_ms'][:]
    segmA=h5hpgrp.attrs['segment_mA'][:]
    dlist=[]
    def indgen(t, ms1d=ms[0]):
        ind=numpy.where(ms1d<=t)[0][-1]
        if ind==len(ms1d)-1:
            ind=len(ms1d)
        return ind
    for count, (ms0, ms1, mA0, mA1) in enumerate(zip(segms[:-1], segms[1:], segmA[:-1], segmA[1:])):
        if (ms1-ms0)<=critms_step:
            d={'segmenttype':'step'}
        elif numpy.abs((mA1-mA0)/(ms1-ms0))<critmAperms_constmA and numpy.abs(mA1-mA0)<critdelmA_constmA:
            if (mA1+mA0)<2.*critmA_zero:
                d={'segmenttype':'zero'}
            else:
                d={'segmenttype':'soak'}
        else:
            d={'segmenttype':'ramp'}
            d['ramprate']=(mA1-mA0)/(ms1-ms0)
        iterpts=h5hpgrp.values()
        h5an=getcalanalysis(h5file, h5expname)
        #print h5hpname, h5hpname in h5an
        if h5hpname in h5an:
            iterpts+=h5an[h5hpname].values()
            #print h5an[h5hpname].items()
        for pnt in iterpts:
            if isinstance(pnt, h5py.Dataset):
                nam=pnt.name.rpartition('/')[2]
                arr=pnt[:, :][:, indgen(ms0):indgen(ms1)]
                #print nam, arr.shape, numpy.isnan(arr).sum()
                if numpy.any(numpy.isnan(arr)):
                    continue
                d[nam]=arr
                for key, val in pnt.attrs.iteritems():
                    if 'unit' in key:
                        d[nam]*=val
        d['cycletime']=ms[:, indgen(ms0):indgen(ms1)]/1000.
        d['segment_ms']=(ms0, ms1)
        d['segment_mA']=(mA0, mA1)
        d['segment_inds']=(indgen(ms0), indgen(ms1))
        d['segindex']=count
        dlist+=[copy.deepcopy(d)]
    h5hpgrp.file.close()
    return dlist

def extractcycle_SegDict(d, cycind):
    dc={}
    for k, v in d.iteritems():
        if isinstance(v, numpy.ndarray) and v.shape==d['cycletime'].shape:
            dc[k]=v[cycind:cycind+1, :]
        else:
            dc[k]=v
    return dc
    
def piecetogethersegments(arrlist):
    cycles=arrlist[0].shape[0]
    return numpy.array([numpy.concatenate([arr[c] for arr in arrlist]) for c in range(cycles)], dtype=arrlist[0].dtype)
    
def saveSCcalculations(h5path, h5expname, h5hpname, hpsegdlist, recname):
    h5file=h5py.File(h5path, mode='r+')
    h5an=getcalanalysis(h5file, h5expname)
    h5hp=gethpgroup(h5file, h5expname)
    if h5hpname in h5an:
        h5g=h5an[h5hpname]
    else:
        h5g=h5an.create_group(h5hpname)
    savekeys=set([(k, d[k].shape[2:], d[k].dtype) for d in hpsegdlist for k in d.keys() if not ('~' in k or k in h5hp[h5hpname] or k=='cycletime') and isinstance(d[k], numpy.ndarray) and d[k].shape[:2]==d['cycletime'].shape])
    #nansegssh=[numpy.ones(d['cycletime'].shape, dtype=d['cycletime'].dtype)*numpy.nan for d in hpsegdlist]
    nansegssh=[d['cycletime'].shape for d in hpsegdlist]
    mastershape=piecetogethersegments([d['cycletime'] for d in hpsegdlist]).shape
    for k, endshape, dtype in list(savekeys):
        savearr=piecetogethersegments([(k in d.keys() and (d[k],) or (numpy.ones(ns+endshape, dtype=dtype)*numpy.nan,))[0] for d, ns in zip(hpsegdlist, nansegssh)])
        if k in h5g:
            del h5g[k]
        ds=h5g.create_dataset(k, data=savearr)
        ds.attrs['recipe']=recname
        
    #now save fit results    
    savekeys=set([k for d in hpsegdlist for k in d.keys() if k.startswith('FITPARS_') and not ('~' in k or k in h5hp[h5hpname] or k=='cycletime') and isinstance(d[k], numpy.ndarray)])
    for k in list(savekeys):
        if k in h5g:
            h5fg=h5g[k]
        else:
            h5fg=h5g.create_group(k)
        for n, d in enumerate(hpsegdlist):
            if k in d.keys():
                if `n` in h5fg:
                    del h5fg[`n`]
                ds=h5fg.create_dataset(`n`, data=d[k])
                ds.attrs['recipe']=recname
    
    #now save profile analysis
    savekeys=set([k for d in hpsegdlist for k in d.keys() if k.startswith('PROFILEANALYSIS_') and not ('~' in k or k in h5hp[h5hpname] or k=='cycletime') and isinstance(d[k], dict)])
    for k in list(savekeys):
        if k in h5g:
            h5pa=h5g[k]
        else:
            h5pa=h5g.create_group(k)
        for n, d in enumerate(hpsegdlist):
            if k in d.keys():
                if `n` in h5pa:
                    h5pac=h5pa[`n`]
                else:
                    h5pac=h5pa.create_group(`n`)
                for arrk, arr in d[k].iteritems():
                    if isinstance(arr, numpy.ndarray):
                        if arrk in h5pac:
                            del h5pac[arrk]
                        ds=h5pac.create_dataset(arrk, data=arr)
                        ds.attrs['recipe']=recname
                
    h5file.close()
def getfitdictlist_hp(h5path, h5expname, h5hpname):
    hpsdl=CreateHeatProgSegDictList(h5path, h5expname, h5hpname)
    h5file=h5py.File(h5path, mode='r')
    fitdlist=[]
    for k in set([k for d in hpsdl for k in d.keys()]):
        for i in range(len(hpsdl)):
            temp=getfitdict_nameseg(h5file, h5expname, h5hpname, k, i)
            if not temp is None:
                fitdlist+=[temp]
    h5file.close()
    return fitdlist
def getfitdict_nameseg(h5pf, h5expname, h5hpname, dsname, seg):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    try:
        h5ang=getcalanalysis(h5file, h5expname)[h5hpname]
        nam='FITPARS_'+dsname
        ds=h5ang[nam][`seg`]
        h5rg=getSCrecipegrp(h5file, h5expname)[ds.attrs['recipe']]
        for fcnname in h5rg.attrs['fcns']:
            g=h5rg[fcnname]
            if g.attrs['savename']==nam:
                break
        d={}
        d['seg']=seg
        d['dsname']=dsname
        d['fitpars']=ds[:, :]
        d['fcnname']=fcnname
        for k in ['parnames', 'segdkeys', 'filters', 'postfilter']:
            temp=g.attrs[k]
            if isinstance(temp, numpy.ndarray):#these are all list of strings stored as arrays so return them to lists
                temp=list(temp)
            d[k]=temp
    except:
        d=None
    if isinstance(h5pf, str):
        return h5file, d
    return d
    
def writecellres(h5path, h5expname, h5hpname, R):
    h5file=h5py.File(h5path, mode='r+')
    h5hpgrp=gethpgroup(h5file, h5expname, h5hpname)
    h5calan=getcalanalysis(h5file, h5expname)
    if not 'CellResistance' in h5calan:
        temp=numpy.zeros(len(h5file.attrs['cells']), dtype='float32')
        h5res=h5calan.create_dataset('CellResistance', data=temp[:])
        h5res.attrs['ambient_tempC']=temp[:]
    h5res=h5calan['CellResistance']

    i=getindex_cell(h5file, h5hpgrp.attrs['CELLNUMBER'])
    h5res[i]=R
    #this doesn't work for unkonw reason: h5res.attrs['ambient_tempC'][i]=h5hpgrp.attrs['ambient_tempC']
    temp=h5res.attrs['ambient_tempC'][:]
    if 'ambient_tempC' in h5hpgrp.attrs.keys():
        temp[i]=h5hpgrp.attrs['ambient_tempC']
    else:
        temp[i]=21.
    h5res.attrs['ambient_tempC']=temp[:]
    h5file.close()

def writecellres_calc(h5path, h5expname, h5hpname, R):
    h5file=h5py.File(h5path, mode='r+')
    h5hpgrp=gethpgroup(h5file, h5expname, h5hpname)
    h5hpgrp.attrs['Ro']=R
    h5file.close()
    

def experimentgrppaths(h5pf):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    p=[]
    for pnt in h5file['Calorimetry'].values():
        if isinstance(pnt, h5py.Group) and 'measurement' in pnt and 'analysis' in pnt:
            p+=[pnt.name]
    if isinstance(h5pf, str):
        return h5file, p
    return p


#
#h5path=os.path.join(os.getcwd(), 'TestImport.h5')
#h5expname='experiment1'
#h5hpname='2010Nov27_Cell1_1mA_50ms_500ms_Ro_1C_a'
#critms_step=1.; critmAperms_constmA=0.0005; critmA_zero=0.05

def AddRes_CreateHeatProgSegDictList(hpsegdlist, SGnpts_curr=100, SGorder_curr=3, SGnpts_volt=100, SGorder_volt=3):
    for d in hpsegdlist:
        if d['segmenttype'] in ['ramp', 'soak']:
            if SGwindow_curr is None or SGorder_curr is None:
                c=copy.copy(d['samplecurrent'])
            else:
                c=numpy.array([savgolsmooth(copy.copy(x), window=SGnpts_curr, order=SGorder_curr) for x in d['samplecurrent']])
            if SGwindow_volt is None or SGorder_volt is None:
                v=copy.copy(d['samplevoltage'])
            else:    
                v=numpy.array([savgolsmooth(copy.copy(x), window=SGnpts_volt, order=SGorder_volt) for x in d['samplevoltage']])
            #print d['samplecurrent'].shape, d['samplevoltage'].shape
            inds=numpy.where(c<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value
            c=replacevalswithneighsin2nddim(c, inds)
            v=replacevalswithneighsin2nddim(v, inds)
            d['sampleresistance']=v/c
            
def RoToAl_h5(h5path, h5expname, h5hpname):#get from the array for the experiment group, an Ro attr in the hp will be used instead if it is available, in which case also use the local To if available
    #rtcpath=rescalpath_getorassign(h5path, h5expname)
    h5file=h5py.File(h5path, mode='r')
    rtcpath=rescalpath_getorassign(h5file, h5expname)
    if not h5file:#in case ti was closed in the fcn
        h5file=h5py.File(h5path, mode='r')
    h5hp=gethpgroup(h5file, h5expname, h5hpname)
    i=getindex_cell(h5file, h5hp.attrs['CELLNUMBER'])
    restempal=list(h5file[rtcpath][i])
    if 'Ro' in h5hp.attrs.keys():
        restempal[0]=h5hp.attrs['Ro']
        if 'ambient_tempC' in h5hp.attrs:
            restempal[1]=h5hp.attrs['ambient_tempC']
    h5file.close()
    return restempal[0], restempal[1], restempal[2]

def AddTemp_CreateHeatProgSegDictList(hpsegdlist, RoToAl, SGwindow_res=None, SGorder_res=3):
    for d in hpsegdlist:
        if 'sampleresistance' in d.keys():
            if SGwindow_res is None or SGorder_res is None:
                R=copy.copy(d['sampleresistance'])
            else:
                R=numpy.array([savgolsmooth(copy.copy(x), window=SGwindow_res, order=SGorder_res) for x in d['sampleresistance']])
            d['sampletemperature']=temp_res(R, RoToAl[0], RoToAl[1], RoToAl[2])
    
def tempvsms_heatprogram(h5path, h5expname, h5hpname, segind=None):
    hpsegdlist=CreateHeatProgSegDictList(h5path, h5expname, h5hpname)
    
    if segind==None:
        hpsegdlist=[d for d in hpsegdlist if d['segmenttype'] in ['ramp', 'soak']]
    else:
        hpsegdlist=[hpsegdlist[segind]]
    print len(hpsegdlist)
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    AddRes_CreateHeatProgSegDictList(hpsegdlist)
    AddTemp_CreateHeatProgSegDictList(hpsegdlist, RoToAl)
    print 'tempdone'
    return numpy.concatenate([d['cycletime'] for d in hpsegdlist], axis=1), numpy.concatenate([d['sampletemperature'] for d in hpsegdlist], axis=1)

def getfiltergrp(h5pf, h5expname):#will only create is passed 'r+'
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5an=getcalanalysis(h5file, h5expname)
    if 'filter' in h5an:
        h5filter=h5an['filter']
    elif h5file.mode=='r':
        return False
    else:
        h5filter=h5an.create_group('filter')
    if isinstance(h5pf, str):
        return h5file, h5filter
    return h5filter

def getfilter(h5pf, h5expname, filtername):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5filter=h5file['Calorimetry'][h5expname]['analysis']['filter']
    d={}
    for k, v in h5filter[filtername].attrs.iteritems():
        d[k]=(v=='None' and (None,) or (v,))[0]
        if isinstance(d[k], numpy.ndarray):
            d[k]=list(d[k])
    if isinstance(h5pf, str):
        return h5file, d
    return d

def getfilterdict(h5pf, h5expname):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5filter=getfiltergrp(h5file, h5expname)
    if not h5filter:
        h5file.close()
        return False
    filterd={}
    for pnt in h5filter.values():
        if isinstance(pnt, h5py.Group):
            d={}
            nam=pnt.name.rpartition('/')[2]
            for k, v in pnt.attrs.iteritems():
                d[k]=(v=='None' and (None,) or (v,))[0]
                if isinstance(d[k], numpy.ndarray):
                    d[k]=list(d[k])
            if len(d.keys())==0:
                continue
            filterd[nam]=d
    if isinstance(h5pf, str):
        return h5file, filterd
    return filterd

def getSCrecipegrp(h5pf, h5expname):#will create the grp only if pass an r+file
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5an=getcalanalysis(h5file, h5expname)
    if (not 'SCrecipe' in h5an) and h5file.mode=='r+':
        h5hp=h5an.create_group('SCrecipe')
    else:
        h5hp=h5an['SCrecipe']
    if isinstance(h5pf, str):
        return h5file, h5hp
    return h5hp
    
def savefilters(h5pf, h5expname, filterd):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r+')
    else:
        h5file=h5pf
    h5filter=getfiltergrp(h5file, h5expname)
    for nam, d in filterd.iteritems():
        if nam in h5filter:
            del h5filter[nam]
        h5g=h5filter.create_group(nam)
        for k, v in d.iteritems():
            h5g.attrs[k]=(v is None and ('None',) or (v,))[0]
    if isinstance(h5pf, str):
        h5file.close()

def saveSCrecipe(h5pf, h5expname, recname, fcns, recdlist):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r+')
    else:
        h5file=h5pf
    h5rec=getSCrecipegrp(h5file, h5expname)
    if recname in h5rec:
        del h5rec[recname]
    h5rg=h5rec.create_group(recname)
    h5rg.attrs['fcns']=fcns
    for f, d in zip(fcns, recdlist):
        h5g=h5rg.create_group(f)
        for k, v in d.iteritems():
            h5g.attrs[k]=(v is None and ('None',) or (v,))[0]
    if isinstance(h5pf, str):
        h5file.close()

def getSCrecipe(h5pf, h5expname, recname):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5rec=getSCrecipegrp(h5file, h5expname)
    h5filter=getfiltergrp(h5file, h5expname)
    h5rg=h5rec[recname]
    fcns=h5rg.attrs['fcns']
    f_saven_namsegkfilk_postfilk=[]
    for f in fcns:
        t=tuple([])
        attrs=h5rg[f].attrs
        t+=(f,)
        t+=(attrs['savename'],)
        t+=([(nam, (segk, filter)) for nam, segk, filter in zip(attrs['parnames'], attrs['segdkeys'], attrs['filters'])],)
        t+=(attrs['postfilter'],)
        f_saven_namsegkfilk_postfilk+=[t]

    if isinstance(h5pf, str):
        return h5file, f_saven_namsegkfilk_postfilk
    return f_saven_namsegkfilk_postfilk

def copySCrecipes(h5path, h5expname, h5expsource):
    h5file=h5py.File(h5path, mode='r+')
    h5srcrec=getSCrecipegrp(h5file, h5expsource)
    h5desrec=getSCrecipegrp(h5file, h5expname)
    filters=[]
    for pnt in h5srcrec.itervalues():
        if isinstance(pnt, h5py.Group):
            nam=pnt.name.rpartition('/')[2]
            if nam in h5desrec:
                del h5desrec[nam]
            h5file.copy(pnt, h5desrec)
            filters+=[fl for pnt2 in pnt.itervalues() if isinstance(pnt2, h5py.Group) and 'filters' in pnt2.attrs.keys() for fl in pnt2.attrs['filters']]
            filters+=[pnt2.attrs['postfilter'] for pnt2 in pnt.itervalues() if isinstance(pnt2, h5py.Group) and 'postfilter' in pnt2.attrs.keys()]
    
    filters=list(set(filters))
    h5srcfil=getfiltergrp(h5file, h5expsource)
    h5desfil=getfiltergrp(h5file, h5expname)
    for nam in filters:
        if nam in h5desfil:
            del h5desfil[nam]
        h5file.copy(h5srcfil[nam], h5desfil)
    h5file.close()
    
def rescalpath_getorassign(h5pf, h5expname, parent=None, forceassign=False, title='select the experiment whose R(T) calibration will be used'):
    openclose=isinstance(h5pf, str)
    if openclose:
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    if forceassign or not ('Res_TempCalPath' in h5file['Calorimetry'][h5expname].attrs.keys()):
        if parent is None:
            print 'Need to ask user for the Res Cal experiment group but UI parent not specified'
            return False
        p=experimentgrppaths(h5file)
        expname=[n.strip('/').rpartition('/')[2] for n in p if 'Res_TempCal' in getcalanalysis(h5file, n.strip('/').rpartition('/')[2])]
        if openclose:
            h5file.close()
        idialog=PnSC_ui.selectorDialog(parent, expname, title=title)#map(operator.itemgetter(1), pathname)
        if idialog.exec_():
            expname=expname[idialog.index]
            reopen=(not openclose) and h5file.mode=='r'
            if openclose or reopen:
                h5file.close()#this is an extra close for openclose
                h5file=h5py.File(h5pf, mode='r+')
            path=getcalanalysis(h5file, expname)['Res_TempCal'].name
            h5file['Calorimetry'][h5expname].attrs['Res_TempCalPath']=path
            if openclose or reopen:
                h5file.close()
            if reopen:
                print 'in rescalpath_getorassign,  file reference was passed but needed to close and there is not a way to reopen'
            return path
        else:
            return False
    else:
        path=h5file['Calorimetry'][h5expname].attrs['Res_TempCalPath']
        if openclose:
            h5file.close()
        return path

def performreferencesubtraction(segd, segkey, filter, h5path):
    ashape=segd[segkey].shape
    h5file=h5py.File(h5path, mode='r')
    refh5grp=h5file['REFh5path']
    if not 'REFalignment' in filter.keys() or not filter['REFalignment'] in segd.keys():
        ref_ind=[0 for temp in range(ashape[0])]
        print 'using start of scan as alignment'
    else:
        k=filter['REFalignment']
        arr=segd[k]
        ref=refh5grp[k][:, :]
        if ref.shape[0]==ashape[0]:
            pass
        elif ref.shape[0]<ashape[0]:
            ref=numpy.array([ref[0]]*ashape[0], dtype=ref.dtype) #use the first cycles over and over again
        else:
            ref=ref[:ashape[0], :]
        if ref.shape[1]==shape[1]:
            ref_ind=[0 for temp in range(ashape[0])]
        elif ref.shape[1]>shape[1]:#in this case, use the alignment dataset to find out which start index (=shift) provides minimum distance between datasets
            ref_ind=[numpy.argmin([((a-r[i:i+ashape[1]])**2).sum() for i in range(ref.shape[1]-ashape[1])]) for a, r in zip(arr, ref)]
        else:
            print 'ABORTING: THE REFERENCE DATA MUST BE AT LEAST AS LONG AS THE DATA'
            h5file.close()
            return
    arr=segd[segkey]
    ref=refh5grp[segkey][:, :]
    return numpy.array([a-r[i:i+ashape[1]] for i, a, r in zip(ref_ind, arr, ref)], dtype=arr.dtype)

def savetwopointres(h5path, resarr, savename):
    f=h5py.File(h5path, mode='r+')
    if 'TwoPointRes' in f:
        g=f['TwoPointRes']
    else:
        g=f.create_group('TwoPointRes')
    if savename in g:
        del g[savename]
    ds=g.create_dataset(savename, data=resarr)
    ds.attrs['doc']='numcellsx2 array, ordered by cells and then Res in Ohms for I-I and V-V'
    f.close()


def readtwopointres(h5path):
    f=h5py.File(h5path, mode='r')
    if 'TwoPointRes' in f:
        g=f['TwoPointRes']
    else:
        return []
    a=[]
    for ds in g.itervalues():
        if isinstance(ds, h5py.Dataset):
            a+=[(ds.name.rpartition('/')[2], ds[:, :])]
    f.close()
    return a


def calcRo_extraptoTo(h5path, h5expname, h5hpname, segind, inds_calcregion, nprevsegs=1, cycind=0, o_R2poly=1, iterations=1):
    hpsegdlist=CreateHeatProgSegDictList(h5path, h5expname, h5hpname)
    d={}
    d['I2']=hpsegdlist[segind]['samplecurrent'][cycind][inds_calcregion[0]:inds_calcregion[1]]
    d['V2']=hpsegdlist[segind]['samplevoltage'][cycind][inds_calcregion[0]:inds_calcregion[1]]
    #T2=hpsegdlist[segind]['sampletemperature'][cycind][inds_calcregion[0]:inds_calcregion[1]]
    I1=hpsegdlist[segind]['samplecurrent'][cycind][:inds_calcregion[0]]
    V1=hpsegdlist[segind]['samplevoltage'][cycind][:inds_calcregion[0]]
    for i in range(1, nprevsegs+1):
        I1=numpy.concatenate([I1, hpsegdlist[segind-i]['samplecurrent'][cycind][:inds_calcregion[0]]])
        V1=numpy.concatenate([V1, hpsegdlist[segind-i]['samplevoltage'][cycind][:inds_calcregion[0]]])
    d['I1']=I1
    d['V1']=V1
    d['RoToAl']=RoToAl_h5(h5path, h5expname, h5hpname)
    d['o_R2poly']=o_R2poly
    print d['RoToAl'][0]
    for i in range(iterations):
        f=calcRofromheatprogram
        Ro, d2=f(**dict([(k, v) for k, v in d.iteritems() if k in f.func_code.co_varnames[:f.func_code.co_argcount]]))
        for k, v in d2.iteritems():
            d[k]=v
        if i==0:
            d['RoToAl_original']=copy.copy(d['RoToAl'])
        d['RoToAl']=(Ro, d['RoToAl'][1], d['RoToAl'][2])
        print d['RoToAl'][0]
    return Ro, d



heatprogrammetadatafcns={\
'Cell Temperature':tempvsms_heatprogram, \
}#each must take h5path, h5expname, h5hpname as arguments


#p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110816_Zr-Hf-B.h5'
#e='quadlinheating2'
#h='cell17_25malinquad2repeat_1_of_1'
#
#ans, d=calcRo_extraptoTo(p, e, h, 2, (50, 1000), o_R2poly=2, iterations=3)
#print 'done'
#p='F:/CHESS2011_h5MAIN/2011Jun01b_AuSiCu.h5'
#e='AuSiCuheat1_Ro'
#h='cell02_ro_test_1_of_1'
#CreateHeatProgSegDictList(p, e, h, critms_step=1., critmAperms_constmA=0.01, critdelmA_constmA=10., critmA_zero=0.1)
    
