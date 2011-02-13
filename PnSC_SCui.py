import time
import os, os.path
import sys
import numpy
import h5py
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import operator
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import numpy.ma as ma
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import pylab
from PnSC_math import *
from PnSC_ui import *
from PnSC_h5io import *


class SCrecipeDialog(QDialog):
    def __init__(self, parent, h5path, h5expname, h5hpname):
        super(SCrecipeDialog, self).__init__(parent)
        self.parent=parent
        self.h5path=h5path
        self.h5expname=h5expname
        self.h5hpname=h5hpname
        self.hpsegdlist=CreateHeatProgSegDictList(h5path, h5expname, h5hpname)
        

        self.filterdlist=[self.Nonefilterdict()]
        self.filterdlist+=[self.dfltfilterdict(deriv=0)]
        temp=self.dfltfilterdict(deriv=1)
        temp['name']='dt_dflt'
        self.filterdlist+=[temp]
        
        importfiltersButton=QPushButton()
        importfiltersButton.setText("import filters/n(overwrites if same name)")
        QObject.connect(importfiltersButton, SIGNAL("pressed()"), self.importfilters)
        
        editfilterButton=QPushButton()
        editfilterButton.setText("open//edit filter/n(use deriv:1 for dt)")
        QObject.connect(editfilterButton, SIGNAL("pressed()"), self.editfilter)
        
        newfilterButton=QPushButton()
        newfilterButton.setText("create new filters/nstaring with:")
        QObject.connect(newfilterButton, SIGNAL("pressed()"), self.newfilter)
        
        self.filterComboBox=QComboBox()
        self.fillfilterComboBox()
        
        self.segComboBox=QComboBox()
        segLabel=QLabel()
        segLabel.setText('select segment(s)/nfor calculation')
        self.segcalcoptions=[]
        for count, d in enumerate(self.hpsegdlist):
            if d['segmenttype'] in ['ramp', 'soak']:
                self.segComboBox.insertItem(len(self.segcalcoptions), 'segment %d : %s' %(count, d['segmenttype']))
                self.segcalcoptions+=[numpy.array([count])]
        self.segComboBox.insertItem(len(self.segcalcoptions), 'all ramp+soak')
        self.segcalcoptions+=[numpy.concatenate(self.segcalcoptions)]
        
        calcallButton=QPushButton()
        calcallButton.setText("Calculate all quantities")
        QObject.connect(calcallButton, SIGNAL("pressed()"), self.calcall)
        
        savefiltersButton=QPushButton()
        savefiltersButton.setText("save all filters")
        QObject.connect(savefiltersButton, SIGNAL("pressed()"), self.savefilters)
        
        saveSCrecipeButton=QPushButton()
        saveSCrecipeButton.setText("save analysis recipe")
        QObject.connect(saveSCrecipeButton, SIGNAL("pressed()"), self.saveSCrecipe)
        
        parLayout = QVBoxLayout()
        self.pardlist=[]
        
        d=self.calclayoutgen('R', ['', ''])
        parLayout.addWidget(d['widget'])
        d['savename']='sampleresistance'
        d['fcns']=[R_IV]
        d['parnames']=[['I', 'V']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage']]
        d['postcombofcns']=[self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)
        self.slidermoved0()
        
        d=self.calclayoutgen('T', ['', ''])
        parLayout.addWidget(d['widget'])
        d['savename']='sampletemperature'
        d['fcns']=[T_IV]
        d['parnames']=[['I', 'V']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage']]
        d['postcombofcns']=[self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved1)
        self.slidermoved1()
        
        d=self.calclayoutgen('P', ['IV', 'I2R'])
        parLayout.addWidget(d['widget'])
        d['savename']='samplepower'
        d['fcns']=[P_IV, P_IR]
        d['parnames']=[['I', 'V'], ['I', 'R']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage'], ['samplecurrent', 'sampleresistance']]
        d['postcombofcns']=[self.filterfill, self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill], [self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved2)
        self.slidermoved2()
        
        d=self.calclayoutgen('S', ['(IdV-VdI)/dI2', 'dT/dt'])
        parLayout.addWidget(d['widget'])
        d['savename']='sampleheatrate'
        d['fcns']=[S_IV, S_T]
        d['parnames']=[['I', 'V', 'dIdt', 'dVdt'], ['dTdt']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage'], ['sampletemperature']]
        d['postcombofcns']=[self.filterfill, self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill, self.derfilterfill, self.derfilterfill], [self.derfilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved3)
        self.slidermoved3()
        
        d=self.calclayoutgen('Q', ['IVdI2/(IdV-VdI)', 'P/S'])
        parLayout.addWidget(d['widget'])
        d['savename']='sampleheatcapacity'
        d['fcns']=[Q_IV, Q_PS]
        d['parnames']=[['I', 'V', 'dIdt', 'dVdt'], ['P', 'S']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage'], ['samplepower', 'sampleheatrate']]
        d['postcombofcns']=[self.filterfill, self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill, self.derfilterfill, self.derfilterfill], [self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved4)
        self.slidermoved4()
        
        
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(155, 370, 341, 32))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        
        mainlayout=QGridLayout()

        warningLabel=QLabel()
        warningLabel.setText('If you run analysis on some segments, you must not change the\nfilters used in that analysis when analyzing other segments.\nThe filters are shared among all segments of\nall heat programs in the experiment Group.')
        mainlayout.addWidget(warningLabel, 0, 0, 1, 4)
        mainlayout.addWidget(importfiltersButton, 1, 0, 1, 2)
        mainlayout.addWidget(editfilterButton, 2, 0, 1, 2)
        mainlayout.addWidget(newfilterButton, 3, 0, 1, 2)
        mainlayout.addWidget(self.filterComboBox, 4, 0, 1, 2)
        
        mainlayout.addWidget(savefiltersButton, 1, 2, 1, 2)
        mainlayout.addWidget(saveSCrecipeButton, 2, 2, 1, 2)
        mainlayout.addWidget(calcallButton, 3, 2, 1, 2)
        mainlayout.addWidget(segLabel, 4, 2)
        mainlayout.addWidget(self.segComboBox, 4, 3)
        mainlayout.addLayout(parLayout, 5, 0, 20, 4)
        
        self.saveCheckBox=QCheckBox()
        self.saveCheckBox.setText('save calculations upon close')
        mainlayout.addWidget(self.saveCheckBox, 25, 0, 1, 2)
        mainlayout.addWidget(self.buttonBox, 25, 2, 1, 2)
        self.setLayout(mainlayout)
        
        self.plotdialog=None
        
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox, SIGNAL("rejected()"), self.reject)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.ExitRoutine)
        QMetaObject.connectSlotsByName(self)
    
    def ExitRoutine(self):
        if self.saveCheckBox.isChecked():
            self.savefilters()
            self.saveSCrecipe()
            saveSCcalculations(self.h5path, self.h5expname, self.hpsegdlist)
    def savefilters(self):
        savefilters(self.h5path, self.h5expname, self.filterdlist)
    
    def saveSCrecipe(self):
        savefilters(self.h5path, self.h5expname, self.filterdlist)
        fcns=[]
        recdlist=[]
        for d in self.pardlist:
            recd={}
            ind=d['slider'].sliderPosition()
            fcns+=[d['fcns'][ind].func_name]
            recd['segdkeys']=d['segdkeys'][ind]
            a, b=zip(*[[nam, str(cb.currentText())] for nam, cb in zip(d['parnames'][ind], d['parcombos'])])#names and cb put together to control which cb are used
            recd['parnames']=list(a)
            recd['filters']=list(b)
            recd['postfilter']=str(d['postcombo'].currentText())
            recd['savename']=d['savename']
            recdlist+=[recd]
        saveSCrecipe(fcns, recdlist)
        
    def calcall(self):
        for d in self.pardlist:
            ind=d['slider'].sliderPosition()
            fcnpars=[(a, b) for a, b in zip(['h5path', 'h5expname', 'h5hpname'], [self.h5path, self.h5expname, self.h5hpname])]
            fcnpars+=[('filterdict_'+nam, fcn(cb, fetchdictpoint=True)) for nam, cb, fcn in zip(d['parnames'][ind], d['parcombos'], d['parcombofcns'][ind])]
            postfilterdict=d['postcombofcns'][ind](d['postcombo'], fetchdictpoint=True)
            seginds=self.segcalcoptions[self.segComboBox.currentIndex()]
            for si in seginds:
                datafcnpars+=[(arrname, self.hpsegdlist[si][arrname]) for arrname in d['segdkeys'][ind]]
                print 'calculating ', d['savename']
                arr=d['fcns'][ind](*(fcnpars+datafcnpars))
                arr=performgenericfilter(arr, postfilterdict)
                self.hpsegdlist[si][d['savename']]=arr
        if self.plotdialog is None:
            self.plotdialog=analysisviewerDialog(self.parent, self.hpsegdlist)
        else:
            #self.plotdialog.hpsegdlist=self.hpsegdlist
            self.plotdialog.initwidgets()

    def slidermoved0(self):
        self.updateparwidgets(0)
    def slidermoved1(self):
        self.updateparwidgets(1)
    def slidermoved2(self):
        self.updateparwidgets(2)
    def slidermoved3(self):
        self.updateparwidgets(3)
    def slidermoved4(self):
        self.updateparwidgets(4)

    def updateparwidgets(self, widgetind):
        d=self.pardlist[widgetind]
        ind=d['slider'].sliderPosition()
        for count, (lb, s) in enumerate(zip(d['parlabels'], d['parnames'][ind])):
            lb.setText(s)
        for lb in d['parlabels'][count+1:]:
            lb.setText('')
        for count, (cb, fcn) in enumerate(zip(d['parcombos'], d['parcombofcns'][ind])):
            cb.setEnabled(True)
            fcn(cb)
        for cb in d['parcombos'][count+1:]:
            cb.clear()
            cb.setEnabled(False)
        d['postcombofcns'][ind](d['postcombo'])
    
    def importfilters(self):
        h5file, expp=experimentgrppaths(self.h5path)
        paths=[]
        names=[]
        for p in expp:
            if 'filter' in h5file[p]['analysis']:
                pnts+=[h5file[p]['analysis']['filter']]
                names+=[p.rpartition('/')[2]]
        idialog=selectorDialog(self, names, title='Select experiment whose filters will be imported')
        if not idialog.exec_():
            h5file.close()
            return
        fnlist=[]
        dlist=[]
        for filter in pnts[idialog.index].values():
            fn=filter.rpartition('/')[2]
            fnlist+=[fn]
            d={'name':fn}
            for k, v in filter.attrs.iteritems():
                d[k]=v
            dlist+=[d]
        self.filterdlist=[d for f in self.filterdlist if not d['name'] in fnlist]+dlist
        self.fillfilterComboBox()
    
    def editfilter(self, editfilter=True):
        i=self.filterComboBox.currentIndex()
        d=copy.deepcopy(self.filterdlist[i])
        edit=True
        while edit:
            if not filterattredit(self, d):
                return
            if not (isinstance(d['name'], str) and len(d['name'])>0):
                edit=True
                QMessageBox.warning(self,"ERROR",  "'name' is not valid")
            if d['name'] in [dd['name'] for count, dd in enumerate(self.filterdlist) if not (editfilter and count==i)]:
                edit=True
                QMessageBox.warning(self,"ERROR",  "'name' already exists")
            if d['deriv']<0 or d['deriv']>1:
                edit=True
                QMessageBox.warning(self,"ERROR",  "only 0th or 1st deriv supported")
        self.filterdicttypechanger(d)
        if editfilter:
            self.filterdlist[i]=d
        else:
            self.filterdlist+=d
        self.fillfilterComboBox()

    def newfilter(self):
        self.editfilter(editfilter=False)

    def fillfilterComboBox(self):
        self.filterComboBox.clear()
        for counter, d in enumerate(self.filterdlist):
            self.filterComboBox.insertItem(counter, d['name'])
        self.updateparwidgets()
        
    def filterfill(self, cb, fetchdictpoint=False):
        if fetchdictpoint:
            name=str(cb.currentText())
            for d in self.filterdlist:
                if d['name']==name:
                    return d
        i=cb.currentIndex()
        cb.clear()
        counter=0
        for d in self.filterdlist:
            if d['SGderiv']==0:
                cb.insertItem(counter, d['name'])
                counter+=1
        cb.setCurrentIndex(i)
    def derfilterfill(self, cb, fetchdictpoint=False):
        if fetchdictpoint:
            name=str(cb.currentText())
            for d in self.filterdlist:
                if d['name']==name:
                    return d
        i=cb.currentIndex()
        cb.clear()
        counter=0
        for d in self.filterdlist:
            if d['SGderiv']==1:
                cb.insertItem(counter, d['name'])
                counter+=1
        cb.setCurrentIndex(i)

    
    def dfltfilterdict(self, deriv=0):
        return {'name':'dflt',\
                'OLnpts':30, \
                'OLnsig':2., \
                'OLgappts':10, \
                'SGnpts':30, \
                'SGorder':3, \
                'SGderiv':deriv, \
                }
    def Nonefilterdict(self):
        return {'name':'None',\
                'OLnpts':None, \
                'OLnsig':None, \
                'OLgappts':None, \
                'SGnpts':None, \
                'SGorder':None, \
                'SGderiv':0, \
                }
    def filterdicttypechanger(self, d):
        dfltd=self.dfltfilterdict()
        for k, v in d.iteritems():
            if k in dfltd.keys():
                if v is None:
                    if k=='SGderiv':
                        d[k]=0
                else:
                    d[k]=type(dfltd[k])(v)
        
    def calclayoutgen(self, varname, eqnames):
        gridLayoutWidget = QWidget()
        gridLayoutWidget.setGeometry(QRect(55, 130, 430, 201))
        gridLayoutWidget.setObjectName("gridLayoutWidget")
        gridLayout = QGridLayout(gridLayoutWidget)
        gridLayout.setHorizontalSpacing(6)
        gridLayout.setObjectName("gridLayout")
        eqLabel0 = QLabel(gridLayoutWidget)
        font = QFont()
        font.setPointSize(14)
        eqLabel0.setFont(font)
        eqLabel0.setObjectName("eqLabel0")
        eqLabel0.setText(eqnames[0])
        gridLayout.addWidget(eqLabel0, 1, 4, 2, 1)
        filterComboBox3 = QComboBox(gridLayoutWidget)
        filterComboBox3.setObjectName("filterComboBox3")
        gridLayout.addWidget(filterComboBox3, 4, 7, 1, 1)
        filterComboBox2 = QComboBox(gridLayoutWidget)
        filterComboBox2.setObjectName("filterComboBox2")
        gridLayout.addWidget(filterComboBox2, 3, 7, 1, 1)
        filterComboBox1 = QComboBox(gridLayoutWidget)
        filterComboBox1.setObjectName("filterComboBox1")
        gridLayout.addWidget(filterComboBox1, 2, 7, 1, 1)
        filterComboBox0 = QComboBox(gridLayoutWidget)
        filterComboBox0.setObjectName("filterComboBox0")
        gridLayout.addWidget(filterComboBox0, 1, 7, 1, 1)
        eqLabel1 = QLabel(gridLayoutWidget)
        font = QFont()
        font.setPointSize(14)
        eqLabel1.setFont(font)
        eqLabel1.setObjectName("eqLabel1")
        eqLabel1.setText(eqnames[1])
        gridLayout.addWidget(eqLabel1, 3, 4, 2, 1)
        parLabel0 = QLabel(gridLayoutWidget)
        parLabel0.setObjectName("parLabel0")
        gridLayout.addWidget(parLabel0, 1, 6, 1, 1)
        parLabel3 = QLabel(gridLayoutWidget)
        parLabel3.setObjectName("parLabel3")
        gridLayout.addWidget(parLabel3, 4, 6, 1, 1)
        parLabel2 = QLabel(gridLayoutWidget)
        parLabel2.setObjectName("parLabel2")
        gridLayout.addWidget(parLabel2, 3, 6, 1, 1)
        postfilterComboBox = QComboBox(gridLayoutWidget)
        postfilterComboBox.setObjectName("postfilterComboBox")
        gridLayout.addWidget(postfilterComboBox, 4, 0, 1, 1)
        postfilterLabel = QLabel(gridLayoutWidget)
        postfilterLabel.setObjectName("postfilterLabel")
        gridLayout.addWidget(postfilterLabel, 3, 0, 1, 1)
        eqSlider = QSlider(gridLayoutWidget)
        eqSlider.setMaximum(1)
        eqSlider.setInvertedAppearance(True)
        eqSlider.setOrientation(Qt.Vertical)
        eqSlider.setTickPosition(QSlider.TicksRight)
        eqSlider.setObjectName("eqSlider")
        gridLayout.addWidget(eqSlider, 2, 3, 2, 1)
        parLabel1 = QLabel(gridLayoutWidget)
        parLabel1.setObjectName("parLabel1")
        gridLayout.addWidget(parLabel1, 2, 6, 1, 1)
        line_2 = QFrame(gridLayoutWidget)
        line_2.setFrameShape(QFrame.VLine)
        line_2.setFrameShadow(QFrame.Sunken)
        line_2.setObjectName("line_2")
        gridLayout.addWidget(line_2, 1, 5, 4, 1)
        varLabel = QLabel(gridLayoutWidget)
        font = QFont()
        font.setPointSize(30)
        varLabel.setFont(font)
        varLabel.setLayoutDirection(Qt.LeftToRight)
        varLabel.setAutoFillBackground(False)
        varLabel.setObjectName("varLabel")
        varLabel.setText(varname)
        gridLayout.addWidget(varLabel, 1, 0, 2, 1)
        line_3 = QFrame(gridLayoutWidget)
        line_3.setFrameShape(QFrame.HLine)
        line_3.setFrameShadow(QFrame.Sunken)
        line_3.setObjectName("line_3")
        gridLayout.addWidget(line_3, 5, 0, 1, 8)
        line = QFrame(gridLayoutWidget)
        line.setFrameShape(QFrame.VLine)
        line.setFrameShadow(QFrame.Sunken)
        line.setObjectName("line")
        gridLayout.addWidget(line, 1, 2, 4, 1)
        line_4 = QFrame(gridLayoutWidget)
        line_4.setFrameShape(QFrame.HLine)
        line_4.setFrameShadow(QFrame.Sunken)
        line_4.setObjectName("line_4")
        gridLayout.addWidget(line_4, 0, 0, 1, 8)
        d={'widget':gridLayoutWidget}
        d['parlabels']=[parLabel0, parLabel1, parLabel2, parLabel3]
        d['parcombos']=[filterComboBox0, filterComboBox1, filterComboBox2, filterComboBox3]
        d['slider']=eqSlider
        d['eqlabels']=[eqLabel0, eqLabel1]
        d['postcombo']=postfilterComboBox
        d['name']=varname
        d['eqnames']=eqnames
        return d


class
    def calcall(self):
        h5file=h5py.File(self.h5path, mode='r')
        f_saven_segdns_partups_postfd=getSCrecipe(h5file, self.h5expname, self.recname)
        h5file.close()    
        for f, saven, segdns, partups, postfd in f_segdns_partups_postfd:
            partups+=[(a, b) for a, b in zip(['h5path', 'h5expname', 'h5hpname'], [self.h5path, self.h5expname, self.h5hpname])]
            seginds=self.segcalcoptions[self.segComboBox.currentIndex()]
            for si in seginds:
                datafcnpars+=[(arrname, self.hpsegdlist[si][arrname]) for arrname in segdns]
                print 'calculating ', d['savename']
                arr=f(*(datafcnpars+partups))
                arr=performgenericfilter(arr, postfd)
                self.hpsegdlist[si][saven]]=arr

class analysisviewerDialog(QDialog):
    def __init__(self, parent, hpsegdlist):
        super(SegmentEditor, self).__init__(parent)
        self.hpsegdlist=hpsegdlist
        self.setWindowTitle('Provide the endpoint of segments in a single cycle')
        
        self.plotw=plotwidget(self)
        self.markersize=markersize
        self.cycledata=cycledata
        self.plotw.axes.plot(self.cycledata[0], self.cycledata[1], 'r.', ms=self.markersize)

        self.plotw.axes.set_xlabel('cycle time (ms)')
        self.plotw.axes.set_ylabel('applied current (mA)')
    
    def initwidgets(self):
        return
def R_IV(h5path, h5expname, h5hpname, samplecurrent, samplevoltage, filterdict_I, filterdict_V):
    I=performgenericfilter(samplecurrent, filterdict_I)
    V=performgenericfilter(samplevoltage, filterdict_V)
    inds=numpy.where(I<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value
    I=replacevalswithneighs(I, inds)
    V=replacevalswithneighs(V, inds)
    return V/I

def T_IV(h5path, h5expname, h5hpname, samplecurrent, samplevoltage, filterdict_I, filterdict_V):
    I=performgenericfilter(samplecurrent, filterdict_I)
    V=performgenericfilter(samplevoltage, filterdict_V)
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    inds=numpy.where(I<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value
    I=replacevalswithneighs(I, inds)
    V=replacevalswithneighs(V, inds)
    return temp_res(V/I, RoToAl[0], RoToAl[1], RoToAl[2])
    
def P_IV(h5path, h5expname, h5hpname, samplecurrent, samplevoltage, filterdict_I, filterdict_V):
    I=performgenericfilter(samplecurrent, filterdict_I)
    V=performgenericfilter(samplevoltage, filterdict_V)
    return V*I
def P_IR(h5path, h5expname, h5hpname, samplecurrent, sampleresistance, filterdict_I, filterdict_R):
    I=performgenericfilter(samplecurrent, filterdict_I)
    R=performgenericfilter(sampleresistance, filterdict_R)
    return R*I**2

def S_IV(h5path, h5expname, h5hpname, samplecurrent, samplevoltage, filterdict_I, filterdict_V, filterdict_dIdt, filterdict_dVdt):
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    dt=dt_h5(h5path, h5expname, h5hpname)
    I=performgenericfilter(samplecurrent, filterdict_I)
    V=performgenericfilter(samplevoltage, filterdict_V)
    dI=performgenericfilter(samplecurrent, filterdict_dIdt)/dt
    dV=performgenericfilter(samplevoltage, filterdict_dVdt)/dt
    inds=numpy.where(dI<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value
    dI=replacevalswithneighs(dI, inds)
    dV=replacevalswithneighs(dV, inds)
    I=replacevalswithneighs(I, inds)
    V=replacevalswithneighs(V, inds)
    return dT_IVdIdV(I, V, dI, dV, RoToAl[0], RoToAl[2])
    
def S_T(h5path, h5expname, h5hpname, sampletemperature, filterdict_dTdt):
    dt=dt_h5(h5path, h5expname, h5hpname)
    dT=performgenericfilter(sampletemperature, filterdict_dTdt)/dt
    return dT

def Q_IV(h5path, h5expname, h5hpname, samplecurrent, samplevoltage, filterdict_I, filterdict_V, filterdict_dIdt, filterdict_dVdt):
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    dt=dt_h5(h5path, h5expname, h5hpname)
    I=performgenericfilter(samplecurrent, filterdict_I)
    V=performgenericfilter(samplevoltage, filterdict_V)
    dI=performgenericfilter(samplecurrent, filterdict_dIdt)/dt
    dV=performgenericfilter(samplevoltage, filterdict_dVdt)/dt
    inds=numpy.where((I*dV-V*dI)<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value
    dI=replacevalswithneighs(dI, inds)
    dV=replacevalswithneighs(dV, inds)
    I=replacevalswithneighs(I, inds)
    V=replacevalswithneighs(V, inds)
    return Q_IVdIdV(I, V, dI, dV, RoToAl[0], RoToAl[2])
def Q_PS(h5path, h5expname, h5hpname, samplepower, sampleheatrate, filterdict_P, filterdict_S):
    P=performgenericfilter(samplepower, filterdict_P)
    S=performgenericfilter(sampleheatrate, filterdict_S)
    inds=numpy.where(I<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value
    P=replacevalswithneighs(P, inds)
    S=replacevalswithneighs(S, inds)
    return P/S



def filterattredit(parent, AttrDict, arr=None, title="Edit filter parameters. 'SG'= SavistkyGolay, 'OL'=OutLier"):#AttrDict is a pointer to a dictionary that may be changed
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
#    
#p='C:/Users/JohnnyG/Documents/HarvardWork/pharma/PnSC/Nanocopeia1_PnSC.h5'
#e='NoSampleRamps'
#h='20110203Nanocop_ramp_160ms_10mA_cell10'
#
#mainapp=QApplication(sys.argv)
#form=SCrecipeDialog(None, p, e, h)
#form.show()
#form.setFocus()
#global PARENT
#PARENT=form
#mainapp.exec_()
