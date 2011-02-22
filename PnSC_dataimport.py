import pylab 
import matplotlib.cm as cm
import numpy
import h5py
import os, os.path, time, copy

from PnSC_ui import *
from PnSC_math import *
from PnSC_h5io import *

def FileImport(parent, protocolname):
    if 'PatDAQ' in protocolname:
        fn='.dat'
        ms='Select .dat text file'
    elif 'JimDAQ' in protocolname:
        fn='.csv'
        ms='Select a .csv text file from the any cycle'
    p=mygetopenfile(parent=parent, markstr=ms, filename=fn)
    if p=='':
        return False
    print '***', FileFormatFunctionLibrary[protocolname]
    ans=FileFormatFunctionLibrary[protocolname](parent, p)
    if not ans:
        return ans
    AttrDict, DataSetDict, SegmentData=ans
    AttrDict['importpath']=str(p)
    AttrDict['protocolname']=protocolname
    if not importcheck(parent, AttrDict, title='Attributes of the Heat Program'):
        return False
    for nam, (d, arr) in DataSetDict.iteritems():
        if not importcheck(parent, d, arr=arr, title='Dataset %s' %nam):
            return False
    return AttrDict, DataSetDict, SegmentData

def importcheck(parent, AttrDict, arr=None, title=''):#AttrDict is a pointer to a dictionary that may be changed
    repeat=True
    count=0
    while repeat:
        if count==1:
            title='PLEASE CONFIRM CHANGES! '+title
        idialog=attreditorDialog(parent, AttrDict, arr=arr, title=title)
        if not idialog.exec_():
            print idialog
            return False
        if idialog.edited:
            for k, v in idialog.attrd.iteritems():
                AttrDict[k]=v
        else:
            repeat=False
        count+=1
    return True
    
    
    
def nanhandler(s):
    if ('NaN' in s):
       return numpy.nan
    try:
        return eval(s)
    except:
        return numpy.nan
        
def truncate_arrlist_shortest(arrlist):
    ln=numpy.uint16([arr.shape for arr in arrlist]).T
    ind=(numpy.uint16(range(l.min())) for l in ln)
    for i, l in enumerate(ln):
        l=l.min()
        arrlist=[arr.take(range(l), axis=i) for arr in arrlist]
    return numpy.array(arrlist)

def JimDAQ_getcell(filepath):
    folderpath, fn=os.path.split(filepath)
    while len(fn)>0 and not (fn[0].isdigit() and fn[0]!='0'): #skip any leading zeros to avoid eval problem
        fn=fn[1:]
    s=''
    while len(fn)>0 and fn[0].isdigit():
        s+=fn[0]
        fn=fn[1:]
    try:
        return eval(s)
    except:
        return 0
        
#The data reading functions should return 3 things, 0: a dictionary containing the attr for the heat program group, 1: a dictionary with key=datasetname, and val=tuple with 0th element an attr dict and 1st element the array, 2:(ms array, mA array) for segments
def JimDAQ_SC(parent, filepath):
    dlist, arrlist=JimDAQ_fileiterator(filepath)
    arr=truncate_arrlist_shortest(arrlist)
    print arr.shape
    d=dlist[0]
    d['daqHz']=100000.
    d['ncycles']=len(dlist)
    ds={}
    ds['samplecurrent']=({'Aunit':0.001}, arr[:, 0, :])
    ds['samplevoltage']=({'Vunit':0.001}, arr[:, 3, :])
    d['CELLNUMBER']=JimDAQ_getcell(filepath)
    d['ambient_atmosphere']='vacuum'
    d['ambient_tempC']=20.
    return d, ds, [[], []]
    
def JimDAQ_fileiterator(filepath):#filepath is a .dat file from any cycle
    folderpath, filename=os.path.split(filepath)
    a, c=os.path.splitext(filename)
    a, b, n=a.rpartition('_of_')
    a, atemp, i=a.rpartition('_')
    a+=atemp
    i=eval(i)
    n=eval(n)
    
    filelist=os.listdir(folderpath)
    dlist=[]
    arrlist=[]
    for cnum in range(1, n+1):
        fn='%s%d%s%d%s' %(a, cnum, b, n, c)
        if fn in filelist:
            p=os.path.join(folderpath, fn)
            print 'reading: ', p
            t1, t2=readdat_JimDAQ(p)
            dlist+=[t1]
            arrlist+=[t2]
    return dlist, arrlist

def readdat_JimDAQ(path):
    d={}
    f=open(path, mode='r')
    lines=f.readlines()
    f.close()
    
    a, b, c=lines[1].partition(':')
    d['epoch']=time.mktime(time.strptime(c.strip(),'%a, %m/%d/%Y %I:%M:%S %p'))
    
    a, b, c=lines[0].partition(':')
    d['operator']=c.strip()
    
    v=[]
    t=[]
    for i, l in enumerate(lines[11:]):
        t=[]
        l=l.strip()
        while len(l)>0:
            a, b, l=l.partition('\t')
            t+=[a]
    
        v+=[[eval(x) for x in t]]

    return d, numpy.float32(v).T

def PatDAQ_filenamedecode(filepath):
    folderpath, filename=os.path.split(filepath)
    fn=filename.lower()
    d={}
    d['CELLNUMBER']=0
    d['ambient_atmosphere']='vacuum'
    SegmentData=([], [])
    
    if 'cell' in fn:
        a, b, c=fn.partition('cell')
        n=''
        c=c.strip()
        while len(c)>0 and c[0].isdigit():
            n+=c[0]
            c=c[1:]
            c=c.strip()
        try:
            n=eval(n.lstrip('0'))
            d['CELLNUMBER']=n
        except:
            pass
        fn=a+c
    if 'mt' in fn:
        underscoreallowed=True
        a, b, c=fn.partition('mt')
        n=''
        a=a.strip()
        while len(a)>0 and (a[-1].isdigit() or (underscoreallowed and a[-1]=='_')):
            if a[-1].isdigit():
                n=a[-1]+n
            else:
                underscoreallowed=False #so only one underscore can be deleted
            a=a[:-1]
            a=a.strip()
        d['ambient_atmosphere']=str(n+'mT')
        for temp in ['He', 'N2', 'Ar', 'H2', 'CO2', 'O2']:
            if temp in filename:
                d['ambient_atmosphere']+=str(' '+temp)
        fn=a+c
    if 'c' in fn:#this is for the number of cycles, just remove this character and any neighboring numbers
        a, b, c=fn.partition('c')
        a=a.strip()
        while len(a)>0 and a[-1].isdigit():
            a=a[:-1]
            a=a.strip()
        c=c.strip()
        while len(c)>0 and c[0].isdigit():
            c=c[0]
            c=c.strip()
    fn=fn.partition('_')[2]
    ma=[]
    ms=[]
    if 'ma' in fn:
        a, b, c=fn.rpartition('ma')
        ma=PatDAQ_extractlist(a.replace('ma', ''))
    if 'ms' in fn:
        a, b, c=fn.rpartition('ms')
        ms=PatDAQ_extractlist(a.replace('ms', ''))
    if len(ma)>2 and len(ms)>2:
        sma=[]
        sms=[]
        totms=None
        for mav, msv in zip(ma, ms):
            sma+=[mav, mav]
            if totms is None:
                totms=0.
                sms+=[totms]
            else:
                sms+=[totms+0.01]
            totms+=msv
            sms+=[totms]
        SegmentData=(numpy.float32(sms), numpy.float32(sma))
        durationguess=SegmentData[0].max()
#    elif len(ms)==2 and len(ma)==1:
#        durationguess=max(ms)
    elif len(ms)>0:
        durationguess=max(ms)
    else:
        durationguess=None
    return d, SegmentData, durationguess

def PatDAQ_extractlist(s):
    nlist=[]
    c=''
    while len(s)>0 and (s[-1].isdigit() or s[-1]=='_'):
        c=s[-1]+c
        s=s[:-1]
        s=s.strip()
    while len(c)>0:
        a, b, c=c.partition('_')
        nlist+=[a]
    print nlist
    return [eval(n) for n in nlist if len(n)>0 and not (False in [nc.isdigit() for nc in n])]
#print PatDAQ_filenamedecode('c:/2010Nov27_Cell5_1mA_50ms_500ms_Ro_10C.dat')

def readdat_PatDAQ(path):
    f=open(path, mode='r')
    lines=f.readlines()
    f.close()
    
    
    v=[]
    t=[]
    for i, l in enumerate(lines):
        t=[]
        l=l.strip().strip('\t')
        while len(l)>0:
            a, b, l=l.partition('\t')
            t+=[a]
        try:
            v+=[[eval(x) for x in t]]
        except:
            print 'Error evaluating ', t

    return numpy.float32(v).T

def PatDAQ_SC(parent, filepath):
    d, SegmentData, durationguess=PatDAQ_filenamedecode(filepath)
    print '^', d, SegmentData, durationguess
    arr=readdat_PatDAQ(filepath)
    d['daqHz']=100000.
    d['operator']=''
    d['epoch']=os.path.getmtime(filepath)
    print durationguess, arr.shape
    if durationguess is None:
        durationguess=arr.shape[1]
    else:
        durationguess*=d['daqHz']/1000.
    idialog=PatDAQCycleEditor(parent, arr[0], durationguess, os.path.split(filepath)[1])
    if not idialog.exec_():
        return False
    temp=idialog.partition_array_triggers(arr[0])
    temp2=idialog.partition_array_triggers(arr[1])
    if temp is None or temp2 is None:
        QMessageBox.warning(self,"FAILED",  "ABORTED due to failure in partitioning data into cycles")
        return False
    mAcycles=temp[0]
    Vcycles=temp2[0]
    
    d['ncycles']=mAcycles.shape[0]
    ds={}
    ds['samplecurrent']=({'Aunit':0.01}, mAcycles[:, :])
    ds['samplevoltage']=({'Vunit':1.}, Vcycles[:, :])
    d['ambient_tempC']=20.
    return d, ds, SegmentData

FileFormatFunctionLibrary={\
                   'PatDAQ_SC':PatDAQ_SC, \
                   'JimDAQ_SC':JimDAQ_SC, \
#                   'JimDAQ_DSC':JimDAQ_DSC, \
#                   'JimDAQ_DTA':JimDAQ_DTA, \
                   }


#a, b, c=JimDAQ_SC('C:\Users\JohnnyG\Documents\HarvardWork\ExampleJimDAQ_SC\S11_85mA_60ms_1_of_11.csv')


