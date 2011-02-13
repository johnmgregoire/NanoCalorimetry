import pylab 
import matplotlib.cm as cm
import numpy
import h5py
import os, os.path, time, copy, operator
from PnSC_ui import *
#from PnSC_dataimport import *
#from PnSC_math import *

v=numpy.linspace(-16,16,5)
xmm=numpy.float32([x for i in range(5) for x in v])
ymm=numpy.float32([x for x in v[::-1] for i in range(5)])

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
    
def create_exp_grp(h5path, h5expname):
    h5file=h5py.File(h5path, mode='r+')
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
        print k, type(v), v
        h5hpg.attrs[k]=v
    for nam, (adict, data) in DataSetDict.iteritems():
        h5d=h5hpg.create_dataset(nam, data=data)
        for k, v in adict.iteritems():
            h5d.attrs[k]=v
    
    h5file.close()

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
    h5hp=h5file['Calorimetry'][h5expname]['measurement']['HeatProgram']
    dt=1./h5hp.attrs['daqHz']
    h5file.close()
    return dt

def gethpgroup(h5pf, h5expname, h5hpname):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5hp=h5file['Calorimetry'][h5expname]['measurement']['HeatProgram']
    if isinstance(h5pf, str):
        return h5file, h5hp[h5hpname]
    return h5hp[h5hpname]
    
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
    
def CreateHeatProgSegDictList(h5path, h5expname, h5hpname, critms_step=1., critmAperms_constmA=0.0005, critmA_zero=0.05):
#the segment types are step, soak, ramp, zero
#the CreateHeatProgSegDictList function reads the data from the .h5 file and organizes in a way that will be useful for many types of analysis. 
#the function returns a list where there is one element for each segment in the heat program
    h5file, h5hpgrp=gethpgroup(h5path, h5expname, h5hpname)
    ms=msarr_hpgrp(h5hpgrp, twod=True)
    segms=h5hpgrp.attrs['segment_ms'][:]
    segmA=h5hpgrp.attrs['segment_mA'][:]
    dlist=[]
    def indgen(t, ms1d=ms[0]):
        return numpy.where(ms1d<=t)[0][-1]
    for count, (ms0, ms1, mA0, mA1) in enumerate(zip(segms[:-1], segms[1:], segmA[:-1], segmA[1:])):
        if (ms1-ms0)<=critms_step:
            d={'segmenttype':'step'}
        elif (mA1-mA0)/(ms1-ms0)<critmAperms_constmA:
            if (mA1+mA0)<2.*critmA_zero:
                d={'segmenttype':'zero'}
            else:
                d={'segmenttype':'soak'}
        else:
            d={'segmenttype':'ramp'}
            d['ramprate']=(mA1-mA0)/(ms1-ms0)
        for pnt in h5hpgrp.values():
            if isinstance(pnt, h5py.Dataset):
                nam=pnt.name.rpartition('/')[2]
                d[nam]=pnt[:, :][:, indgen(ms0):indgen(ms1)]
                for key, val in pnt.attrs.iteritems():
                    if 'unit' in key:
                        d[nam]*=val
        d['ms']=ms[:, indgen(ms0):indgen(ms1)]
        d['segment_ms']=(ms0, ms1)
        d['segment_mA']=(mA0, mA1)
        d['segment_inds']=(indgen(ms0), indgen(ms1))
        d['segindex']=count
        dlist+=[copy.deepcopy(d)]
    h5hpgrp.file.close()
    return dlist

def writecellres(h5path, h5expname, h5hpname, R):
    h5file=h5py.File(h5path, mode='r+')
    h5hpgrp=gethpgroup(h5file, h5expname, h5hpname)
    h5calan=getcalanalysis(h5file, h5expname)
    if not 'CellResistance' in h5calan:
        temp=numpy.zeros(len(h5file.attrs['cells']), dtype='float32')
        h5res=h5calan.create_dataset('CellResistance', data=temp)
        h5res.attrs['ambient_tempC']=temp
    h5res=h5calan['CellResistance']

    i=getindex_cell(h5file, h5hpgrp.attrs['CELLNUMBER'])
    h5res[i]=R
    temp=h5hpgrp.attrs['ambient_tempC']
    temp[i]=h5hpgrp.attrs['ambient_tempC']
    h5res.attrs['ambient_tempC']=temp
    h5file.close()

def experimentgrppaths(h5pf):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    p=[]
    for pnt in h5file['Calorimetry'].values():
        if isinstance(pnt, h5py.Group):
            p+=[pnt.name]
    if isinstance(h5pf, str):
        return h5file, p
    return p
    
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
        idialog=selectorDialog(parent, expname, title=title)#map(operator.itemgetter(1), pathname)
        if idialog.exec_():
            expname=expname[idialog.index]
            reopen=(not openclose) and h5file.mode=='r'
            if openclose or reopen:
                h5file.close()#this is an extra close for openclose
                h5file=h5py.File(h5path, mode='r+')
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
            
def RoToAl_h5(h5path, h5expname, h5hpname):
    #rtcpath=rescalpath_getorassign(h5path, h5expname)
    h5file=h5py.File(h5path, mode='r')
    rtcpath=rescalpath_getorassign(h5file, h5expname)
    if not h5file:#in case ti was closed in the fcn
        h5file=h5py.File(h5path, mode='r')
    i=getindex_cell(h5file, gethpgroup(h5file, h5expname, h5hpname).attrs['CELLNUMBER'])
    restempal=h5file[rtcpath][i]
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
    return numpy.concatenate([d['ms'] for d in hpsegdlist], axis=1), numpy.concatenate([d['sampletemperature'] for d in hpsegdlist], axis=1)

def getfiltergrp(h5pf, h5expname):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    h5filter=h5file['Calorimetry'][h5expname]['analysis']['filter']
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
        d[k]=v
    if isinstance(h5pf, str):
        return h5file, d
    return d


def getSCrecipegrp(h5pf, h5expname):#will create the grp only if pass an r+file
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r')
    else:
        h5file=h5pf
    if (not 'filter' in h5file['Calorimetry'][h5expname]['analysis']) and h5file.mode=='r+':
        h5hp=h5file['Calorimetry'][h5expname]['analysis'].create_group('SCrecipe')
    else:
        h5hp=h5file['Calorimetry'][h5expname]['analysis']['SCrecipe']
    if isinstance(h5pf, str):
        return h5file, h5hp
    return h5hp
    
def savefilters(h5pf, h5expname, filterdlist):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r+')
    else:
        h5file=h5pf
    h5filter=getfiltergrp(h5file, h5expname)
    for d in filterdlist:
        if d['name'] in h5filter:
            del h5filter[d['name']]
        h5g=h5filter.create_group(d['name'])
        for k, v in d:
            if k=='name':
                continue
            h5g.attrs[k]=v
    if isinstance(h5pf, str):
        h5file.close()

def saveSCrecipe(h5pf, h5expname, recname, fcns, recdlist):
    if isinstance(h5pf, str):
        h5file=h5py.File(h5pf, mode='r+')
    else:
        h5file=h5pf
    h5rec=getSCrecipegrp(h5file, h5expname)
    if recname in h5recg:
        del h5rec[recname]
    h5rg=h5rec.create_group(recname)
    h5rg.attrs['fcns']=fcns
    for f, d in zip(fcns, recdlist):
        h5g=h5rg.create_group(f)
        for k, v in d:
            h5g.attrs[k]=v
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
    f_saven_segdns_partups_postfd=[]
    for f in fcns:
        fspp=tuple([])
        attrs=h5g[f].attrs
        fspp+=(eval(f),)
        fspp+=(attrs['savename'],)
        fspp+=(attrs['segkeys'],)
        fspp+=([('filterdict_'+nam, getfilter(filter)) for nam, filter in zip(attr['parnames'], attr['filters'])],)
        fspp+=(getfilter(attr['postfilter']),)
        f_saven_segdns_partups_postfd+=[fspp]

    if isinstance(h5pf, str):
        return h5file, f_saven_segdns_partups_postfd
    return f_saven_segdns_partups_postfd


heatprogrammetadatafcns={\
'Cell Temperature':tempvsms_heatprogram, \
}#each must take h5path, h5expname, h5hpname as arguments

#p='C:/Users/JohnnyG/Documents/HarvardWork/pharma/PnSC/Nanocopeia1_PnSC.h5'
#e='NoSampleRamps'
#h='20110203Nanocop_ramp_160ms_10mA_cell10'
#
#e='Res_nosample_150C'
#h='20110202nanocop_150C_1mA_cell2'
#ans=tempvsms_heatprogram(p, e, h)
#print 'done'
