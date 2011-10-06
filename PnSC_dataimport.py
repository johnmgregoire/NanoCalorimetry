import pylab 
import matplotlib.cm as cm
import numpy
import h5py
import os, os.path, time, copy
import struct
from PnSC_ui import *
from PnSC_math import *
from PnSC_h5io import *

def getemptybatchattrdict():
    batchattrdict={}
    for k in ['grpname', 'protname','path','durSpinBox','nnoiseSpinBox', 'naboveSpinBox', 'nsigSpinBox', 'firstderptsSpinBox','secderptsSpinBox','secdervalSpinBox','savegrpname']:
        batchattrdict[k]=None
    return batchattrdict
def FileImport(parent, protocolname, batchattrdict=None):
    if 'PatDAQ' in protocolname:
        fn='.dat'
        ms='Select .dat text file'
    elif 'JimDAQ' in protocolname:
        fn='.csv'
        ms='Select a .csv text file from the any cycle'
    if batchattrdict is None:
        p=mygetopenfile(parent=parent, markstr=ms, filename=fn)
    else:
        p=batchattrdict['path']
        
    if p=='':
        return False
    ans=FileFormatFunctionLibrary[protocolname](parent, p, batchattrdict=batchattrdict)
    if not ans:
        return ans
    AttrDict, DataSetDict, SegmentData=ans
    AttrDict['importpath']=str(p)
    AttrDict['protocolname']=protocolname
    if batchattrdict is None:
        if not importcheck(parent, AttrDict, title='Attributes of the Heat Program'):
            return False
        for nam, (d, arr) in DataSetDict.iteritems():
            if not importcheck(parent, d, arr=arr, title='Dataset %s' %nam):
                return False
    else:
        for k in AttrDict.keys():
            if k in batchattrdict.keys():
                print 'replacing AttrDict key %s from %s to %s'  %(k, `AttrDict[k]`, `batchattrdict[k]`)
                AttrDict[k]=copy.copy(batchattrdict[k])
    return AttrDict, DataSetDict, SegmentData

def importcheck(parent, AttrDict, arr=None, title=''):#AttrDict is a pointer to a dictionary that may be changed
    repeat=True
    count=0
    while repeat:
        if count==1:
            title='PLEASE CONFIRM CHANGES! '+title
        idialog=attreditorDialog(parent, AttrDict, arr=arr, title=title)
        if not idialog.exec_():
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
    ln=numpy.uint32([arr.shape for arr in arrlist]).T
    ind=(numpy.uint32(range(l.min())) for l in ln)
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
def JimDAQ_SC(parent, filepath, batchattrdict=None):
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

def JimDAQ2011_SC(parent, filepath, batchattrdict=None):
    dlist, arrlist=JimDAQ2011_fileiterator(filepath)
    arr=truncate_arrlist_shortest(arrlist)
    d=JimDAQ2011_translateheader(dlist[0])
    d['ncycles']=len(dlist)
    ds={}
    ds['samplecurrent']=({'Aunit':0.001}, arr[:, :, 0])
    ds['samplevoltage']=({'Vunit':0.001}, arr[:, :, 3])
    d['CELLNUMBER']=JimDAQ_getcell(filepath)
    return d, ds, [[], []]

def JimDAQ2011_acSC(parent, filepath, batchattrdict=None):
    dlist, arrlist=JimDAQ2011_fileiterator(filepath)
    arr=truncate_arrlist_shortest(arrlist)
    print arr.shape
    d=JimDAQ2011_translateheader(dlist[0])
    d['ncycles']=len(dlist)
    if not 'pts_sincycle' in d.keys():
        d['pts_sincycle']=30.
    ds={}
    ds['samplecurrent']=({'Aunit':0.001}, arr[:, :, 0])
    ds['samplevoltage']=({'Vunit':0.001}, arr[:, :, 3])
    ds['samplefilteredvoltage']=({'Vunit':0.001}, arr[:, :, 4])
    d['CELLNUMBER']=JimDAQ_getcell(filepath)
    return d, ds, [[], []]

    
def JimDAQ2011_translateheader(d):
    dummyfcn=lambda x:x
    key_headerkey_fcn_dflt=[\
    ('operator', 'name', dummyfcn,''),\
    ('epoch', 'date', lambda c:time.mktime(time.strptime(c.strip(),'%a, %m/%d/%Y %I:%M:%S %p')), 0), \
    ('ambient_tempC', 'furnace temp (C)', dummyfcn, 20.),\
    ('ambient_atmosphere', 'atmosphere', dummyfcn, 'vacuum'),\
    ('daqHz', 'daqtime_us', lambda c:1.e6/c, 301204.8), \
    ]
    for k, hk, f, dflt in key_headerkey_fcn_dflt:
        if hk in d.keys():
            temp=d[hk]
            del d[hk]
            d[k]=f(temp)
        else:
            d[k]=dflt
    return d
    
def JimDAQ2011_fileiterator(filepath):#filepath is a .dat file from any cycle
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
            t1, t2=readdat_JimDAQ2011(p)
            dlist+=[t1]
            arrlist+=[t2]
    return dlist, arrlist

def CHESSDAQ2011(parent, filepath, batchattrdict=None):
    dlist, arrlist=CHESSDAQ2011_fileiterator(filepath)
    arr=truncate_arrlist_shortest(arrlist)
    print arr.shape
    d=CHESSDAQ2011_translateheader(dlist[0])
    d['ncycles']=len(dlist)
    ds={}
    ds['samplecurrent']=({'Aunit':0.001}, arr[:, 0, :])
    ds['samplevoltage']=({'Vunit':0.001}, arr[:, 1, :])
    if arr.shape[1]==3:
        ds['samplefilteredvoltage']=({'Vunit':0.001}, arr[:, 2, :])
        if not 'pts_sincycle' in d.keys():
            d['pts_sincycle']=30.
    d['CELLNUMBER']=JimDAQ_getcell(filepath)
    return d, ds, [[], []]

def readdat_CHESSDAQ2011(path, startofheader=':header_start:', endofheader=':header_end:\r\n', startofdata=':data_begin:', endofdata=':data_end:\r\n'): 
    #read all keyword attributes and any in the comments section
    f=open(path, 'rb')
    bdata=f.read()
    f.close()
    headstr, garb, datasection=bdata.partition(endofheader)

    
    
    def attemptnumericconversion(s):
        if (s.replace('.', '', 1).replace('e', '', 1).replace('+', '', 1).replace('-', '', 1)).isalnum():
            try:
                return eval(s)
            except:
                pass
        return s
    garb, garb, headstr=headstr.partition(startofheader)    
    d={}
    while len(headstr)>0:
        a, garb, headstr=headstr.partition('\n')
        a=a.strip()
        if a.startswith(':'):
            b, garb, c=a[1:].partition(':')
            d[b]=attemptnumericconversion(c.strip())
        if a.startswith(':coms:'):
            b, garb, c=a[1:].partition(':')
            i=0
            while i<len(c)-1:
                j=c.find(':', i)
                k=c.find(':', i+j+1)
                if j<0 or k<0:
                    break
                i=c.find(':', k+1)
                if i<0:
                    i=len(c)-1
                d[c[j+1:k]]=attemptnumericconversion((c[k+1:i]).strip())
    fullscalelist=listeval(d['fullscale_fields'])
    rangelist=listeval(d['NIai_mVrange'])
    nchan=len(fullscalelist)

    datastr, garb, datasection=datasection.partition(endofdata)
    if len(datastr)%nchan>0:
        datastr=datastr[:-1*(len(datastr)%nchan)]
    z=[struct.unpack('>f',datastr[i:i+4]) for i in range(0,len(datastr),4)] # '>d' is big-endian double float, Labview uses big-endian
    z=numpy.reshape(numpy.float32(z), nchan, (len(z)//nchan))
    
    z=numpy.float32([za/(rng*0.001)*fs for za, fs, rng in z, fullscalelist, rangelist])
    return d, z

def listeval(c):
    returnlist=[]
    if '\t' in c:
        delim='\t'
    else:
        delim=' '
    while len(c)>0:
        a, garb, c=c.partition(delim)
        a=a.strip()
        c=c.strip()
        a=a.lstrip('0').rstrip('.')
        a=(a=='' and (0,) or (eval(a),))[0]
        returnlist+=[a]
    return returnlist
    
def CHESSDAQ2011_translateheader(d):
    dummyfcn=lambda x:x
    key_headerkey_fcn_dflt=[\
    ('operator', 'name', dummyfcn,''),\
    ('epoch', 'date', lambda c:time.mktime(time.strptime(c.strip(),'%a, %m/%d/%Y %I:%M:%S %p')), 0), \
    ('ambient_atmosphere', 'atmosphere', dummyfcn, 'vacuum'),\
    ('pts_sincycle', 'writepts_sincycle', lambda n: n*d['daqHz']/d['writeHz']), \
    ]
    for k, hk, f, dflt in key_headerkey_fcn_dflt:
        if hk in d.keys():
            temp=d[hk]
            del d[hk]
            d[k]=f(temp)
        else:
            d[k]=dflt
    return d
    
def CHESSDAQ2011_fileiterator(filepath):#filepath is a .dat file from any cycle
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
            t1, t2=readdat_CHESSDAQ2011(p)
            dlist+=[t1]
            arrlist+=[t2]
    return dlist, arrlist
    

def uint16tofloat32(x, offsetbinary=True, posfullscale=1.):#a negative posfullscale value will invert
    if not offsetbinary:
        x+=32768
    return numpy.float32(posfullscale*(x/32768.-1.))

def readdat_JimDAQ2011(path, startofheader=':header_start:', endofheader=':header_end:\r\n'): 
    #startofheader does not have to be start of file but endofheader has to include the character that is just before binary uint16 data
    #read all keyword attributes and any in the comments section
    f=open(path, 'rb')
    bdata=f.read()
    f.close()
    headstr, garb, uintdata=bdata.partition(endofheader)
    if len(uintdata)%8>0:
        uintdata=uintdata[:-1*(len(uintdata)%8)]
    z=[struct.unpack('>H',uintdata[i:i+2]) for i in range(0,len(uintdata),2)] # '>8' is big-endian uint16, Labview uses big-endian
    z=numpy.reshape(numpy.uint16(z),(len(z)//8,8))
    
    
    def attemptnumericconversion(s):
        if (s.replace('.', '', 1).replace('e', '', 1).replace('+', '', 1).replace('-', '', 1)).isalnum():
            try:
                return eval(s)
            except:
                pass
        return s
    garb, garb, headstr=headstr.partition(startofheader)    
    d={}
    while len(headstr)>0:
        a, garb, headstr=headstr.partition('\n')
        a=a.strip()
        if a.startswith(':'):
            b, garb, c=a[1:].partition(':')
            d[b]=attemptnumericconversion(c.strip())
        if a.startswith(':coms:'):
            b, garb, c=a[1:].partition(':')
            i=0
            while i<len(c)-1:
                j=c.find(':', i)
                k=c.find(':', i+j+1)
                if j<0 or k<0:
                    break
                i=c.find(':', k+1)
                if i<0:
                    i=len(c)-1
                d[c[j+1:k]]=attemptnumericconversion((c[k+1:i]).strip())
    fullscalelist=[]
    c=d['fullscale_fields']
    if '\t' in c:
        delim='\t'
    else:
        delim=' '
    while len(c)>0:
        a, garb, c=c.partition(delim)
        a=a.strip()
        c=c.strip()
        a=a.lstrip('0').rstrip('.')
        a=(a=='' and (0,) or (eval(a),))[0]
        fullscalelist+=[a]
    
    z=numpy.float32([uint16tofloat32(z[:, i], offsetbinary=ziz, posfullscale=fullsc) for i, (ziz, fullsc) in enumerate(zip([0, 0, 0, 0, 0, 1, 1, 0], fullscalelist))]).T
    return d, z

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

def PatDAQ_SC(parent, filepath, batchattrdict=None):
    d, SegmentData, durationguess=PatDAQ_filenamedecode(filepath)
    arr=readdat_PatDAQ(filepath)
    d['daqHz']=100000.
    d['operator']=''
    d['epoch']=os.path.getmtime(filepath)

    if durationguess is None:
        durationguess=arr.shape[1]
    else:
        durationguess*=d['daqHz']/1000.#durationguess from filename is in ms
    if 'ncycles' in d.keys() and not (d['ncycles'] is None):
        durationguess=int(round(1.*arr.shape[1]/d['ncycles']))

    if batchattrdict is None:
        idialog=PatDAQCycleEditor(parent, arr[0], durationguess, os.path.split(filepath)[1])
        if not idialog.exec_():
            return False
    else:
        if 'ncycles' in batchattrdict.keys() and not (batchattrdict['ncycles'] is None):
            durationguess=int(round(1.*arr.shape[1]/batchattrdict['ncycles']))
        idialog=PatDAQCycleEditor(parent, arr[0], durationguess, os.path.split(filepath)[1])
        for sb, k in [(idialog.durSpinBox, 'durSpinBox'), (idialog.nnoiseSpinBox, 'nnoiseSpinBox'), (idialog.nsigSpinBox, 'nsigSpinBox'), (idialog.naboveSpinBox, 'naboveSpinBox')]:
            if k in batchattrdict and not (batchattrdict[k] is None):
                sb.setValue(batchattrdict[k])
        idialog.calccycles()
        idialog.ExitRoutine()

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

def readdat_PatDAQ2011(path):
    f=open(path, mode='r')
    lines=f.readlines()
    f.close()
    
    def evalrow(l):
        t=[]
        l=l.strip().strip('\t')
        while len(l)>0:
            a, b, l=l.partition('\t')
            t+=[a]
        try:
            v=numpy.array([eval(x) for x in t])
        except:
            v=t
        if len(v)==1:
            return v[0]
        else:
            return v
    v=[]
    d={}
    for i, l in enumerate(lines):
        if l.startswith('%'):
            k, garb, s=l.partition('%')[2].partition(':')
            if len(k)>0 and len(s)>0:
                d[k]=evalrow(s)
        else:
            l=l.strip().strip('\t')
            if len(l)>0:
                v+=[evalrow(l)]
    return d, numpy.float32(v).T

def PatDAQ2011_SC(parent, filepath, batchattrdict=None):
    durationguess=None
    #SegmentData=?
    d, arr=readdat_PatDAQ2011(filepath)
    if not 'daqHz' in d.keys():
        d['daqHz']=100000.
    if not 'operator' in d.keys():
        d['operator']=''
    d['epoch']=os.path.getmtime(filepath)

    if durationguess is None:
        durationguess=arr.shape[1]
    else:
        durationguess*=d['daqHz']/1000.#durationguess from filename is in ms
    if 'ncycles' in d.keys() and not (d['ncycles'] is None):
        durationguess=int(round(1.*arr.shape[1]/d['ncycles']))

    if batchattrdict is None:
        idialog=PatDAQCycleEditor(parent, arr[0], durationguess, os.path.split(filepath)[1])
        if not idialog.exec_():
            return False
    else:
        if 'ncycles' in batchattrdict.keys() and not (batchattrdict['ncycles'] is None):
            durationguess=int(round(1.*arr.shape[1]/batchattrdict['ncycles']))
        idialog=PatDAQCycleEditor(parent, arr[0], durationguess, os.path.split(filepath)[1])
        for sb, k in [(idialog.durSpinBox, 'durSpinBox'), (idialog.nnoiseSpinBox, 'nnoiseSpinBox'), (idialog.nsigSpinBox, 'nsigSpinBox'), (idialog.naboveSpinBox, 'naboveSpinBox')]:
            if k in batchattrdict and not (batchattrdict[k] is None):
                sb.setValue(batchattrdict[k])
        idialog.calccycles()
        idialog.ExitRoutine()

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
                    'CHESSDAQ2011':CHESSDAQ2011, \
                    'JimDAQ2011_acSC':JimDAQ2011_acSC, \
                   'PatDAQ_SC':PatDAQ_SC, \
                   'PatDAQ2011_SC':PatDAQ2011_SC, \
                   'JimDAQ_SC':JimDAQ_SC, \
                   'JimDAQ2011_SC':JimDAQ2011_SC, \
#                   'JimDAQ_DSC':JimDAQ_DSC, \
#                   'JimDAQ_DTA':JimDAQ_DTA, \
                   }


#p='C:/Users/JohnnyG/Documents/HarvardWork/ACcal/20110708_initACtests/cell29_10Ohm_10mAdc_9mA10kHz_9kHz11kHzfilter_wTtokeithley_1_of_1.dat'
#dlist, arrlist=JimDAQ2011_fileiterator(p)
