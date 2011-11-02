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
    def __init__(self, parent, h5path, h5expname, h5hpname, calctype='RTPSD', par=None):
        super(SCrecipeDialog, self).__init__(parent)
        self.parent=parent
        self.h5path=h5path
        self.h5expname=h5expname
        self.h5hpname=h5hpname
        self.hpsegdlist=CreateHeatProgSegDictList(h5path, h5expname, h5hpname)
        
        self.calctype=calctype
        
        self.setWindowTitle('SC analysis recipe editor')
        
        self.parLayout = QVBoxLayout()
        self.pardlist=[]
        self.filterd={}
        self.filterd['None']=self.Nonefilterdict()
        
        if self.calctype=='RTPSD':
            ncalcs=self.RTPSDsetup()
            self.filterd['dflt']=self.dfltfilterdict(deriv=0)
            self.filterd['dt_dflt']=self.dfltfilterdict(deriv=1)
        elif self.calctype=='FitPS':
            ncalcs=self.FitPSsetup()
            for d in self.dfltfitfilterdicts():
                self.filterd[d['name']]=copy.deepcopy(d)
        elif self.calctype=='QUC':
            ncalcs=self.QUCsetup()
            self.filterd['fiteval']=self.fitevalfilterdflt()
            #self.filterd['Qdflt']=self.Qfilterdict(par)
            #self.filterd['refdflt']=self.reffilterdict(par)
        elif self.calctype=='CTpk':
            ncalcs=self.CTpksetup()
            self.filterd['gridata']=self.griddatafilterdflt()
            self.filterd['peaksearchfit']=self.peakfitfilterdflt()
        elif self.calctype=='AC':
            ncalcs=self.ACpksetup()
            self.filterd['ac_harmonics']=self.acfilterdflt()
        elif self.calctype=='AC_RTPSDsetup':
            ncalcs=self.AC_RTPSDsetup()
            self.filterd['ac_chooseharm']=self.acvfilterdflt()
        #*****************************************************
        importfiltersButton=QPushButton()
        importfiltersButton.setText("import filters\n(overwrites if same name)")
        QObject.connect(importfiltersButton, SIGNAL("pressed()"), self.importfilters)
        
        editfilterButton=QPushButton()
        editfilterButton.setText("open/edit filter\n(use deriv:1 for dt)")
        QObject.connect(editfilterButton, SIGNAL("pressed()"), self.editfilter)
        
        newfilterButton=QPushButton()
        newfilterButton.setText("create new filters\nstaring with:")
        QObject.connect(newfilterButton, SIGNAL("pressed()"), self.newfilter)
        
        self.filterComboBox=QComboBox()
        
        self.importfilters(h5exp=self.h5expname)
        self.fillfilterComboBox()
        
        self.segComboBox=QComboBox()
        segLabel=QLabel()
        segLabel.setText('select segment(s)\nfor calculation')
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
        
        recLabel=QLabel()
        recLabel.setText('recipe:')
        self.recLineEdit=QLineEdit()
        self.recLineEdit.setText('dflt')
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(155, 370, 341, 32))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        
        mainlayout=QGridLayout()

        warningLabel=QLabel()
        warningLabel.setText('If you run analysis on some segments, you must not edit filters used in that analysis\nwhen analyzing other segments. The filters are common to the experiment Group.')
        mainlayout.addWidget(warningLabel, 0, 0, 1, 4)
        leftlayout=QGridLayout()
        leftlayout.addWidget(importfiltersButton, 0, 0, 1, 2)
        leftlayout.addWidget(editfilterButton, 1, 0, 1, 2)
        leftlayout.addWidget(newfilterButton, 2, 0, 1, 2)
        leftlayout.addWidget(self.filterComboBox, 3, 0, 1, 2)
        mainlayout.addLayout(leftlayout, 1, 0, 4, 2)
        
        rightlayout=QGridLayout()
        rightlayout.addWidget(savefiltersButton, 2, 0, 1, 2)
        rightlayout.addWidget(saveSCrecipeButton, 3, 0, 1, 2)
        rightlayout.addWidget(calcallButton, 1, 0, 1, 2)
        rightlayout.addWidget(segLabel, 0, 0)
        rightlayout.addWidget(self.segComboBox, 0, 1)
        rightlayout.addWidget(recLabel, 4, 0)
        rightlayout.addWidget(self.recLineEdit, 4, 1)
        mainlayout.addLayout(rightlayout, 1, 2, 4, 2)
        
        mainlayout.addLayout(self.parLayout, 5, 0, ncalcs*4, 4)
        self.saveCheckBox=QCheckBox()
        self.saveCheckBox.setText('save calculations upon close')
        mainlayout.addWidget(self.saveCheckBox, 25, 0, 1, 2)
        mainlayout.addWidget(self.buttonBox, ncalcs*4+5, 2, 1, 2)
        self.setLayout(mainlayout)
        
        self.plotdialog=None
        
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox, SIGNAL("rejected()"), self.reject)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.ExitRoutine)
        QMetaObject.connectSlotsByName(self)
        
        self.updateparwidgets()
#        self.slidermoved0()
#        self.slidermoved1()
#        self.slidermoved2()
#        self.slidermoved3()
#        self.slidermoved4()
    
    def RTPSDsetup(self):
        d=self.calclayoutgen('R', [''])
        self.parLayout.addWidget(d['widget'])
        d['savename']='sampleresistance'
        d['fcns']=[R_IV]
        d['parnames']=[['I', 'V']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage']]
        d['postcombofcns']=[self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)
        
        d=self.calclayoutgen('T', [''])
        self.parLayout.addWidget(d['widget'])
        d['savename']='sampletemperature'
        d['fcns']=[T_IV]
        d['parnames']=[['I', 'V']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage']]
        d['postcombofcns']=[self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved1)
        
        d=self.calclayoutgen('P', ['IV', 'I2R'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='samplepower'
        d['fcns']=[P_IV, P_IR]
        d['parnames']=[['I', 'V'], ['I', 'R']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage'], ['samplecurrent', 'sampleresistance']]
        d['postcombofcns']=[self.filterfill, self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill], [self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved2)
        
        d=self.calclayoutgen('S', ['(IdV-VdI)/dI2', 'dT/dt'])#, 'avedT/dt'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='sampleheatrate'
        d['fcns']=[S_IV, S_T]#, S_Tavesl]
        d['parnames']=[['I', 'V', 'dIdt', 'dVdt'], ['dTdt']]#, ['dTdt']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage', 'samplecurrent', 'samplevoltage'], ['sampletemperature']]#, ['sampletemperature']]
        d['postcombofcns']=[self.filterfill, self.filterfill]#, self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill, self.derfilterfill, self.derfilterfill], [self.derfilterfill]]#, [self.derfilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved3)
        
        d=self.calclayoutgen('D', ['IVdI2/(IdV-VdI)', 'P/S'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='samplepowerperrate'
        d['fcns']=[D_IV, D_PS]
        d['parnames']=[['I', 'V', 'dIdt', 'dVdt'], ['P', 'S']]
        d['segdkeys']=[['samplecurrent', 'samplevoltage', 'samplecurrent', 'samplevoltage'], ['samplepower', 'sampleheatrate']]
        d['postcombofcns']=[self.filterfill, self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill, self.derfilterfill, self.derfilterfill], [self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved4)
        return 5

    def FitPSsetup(self):
        d=self.calclayoutgen('pt', ['none', 'use existing','user-def', 'autocalc'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='cyclepartition'
        d['fcns']=[pt_none, useexisting, pt_user, pt_calc]
        d['parnames']=[['t'], ['x'], ['t', 'D'], ['t', 'D']]
        d['segdkeys']=[['cycletime'], ['cyclepartition'], ['cycletime', 'samplepowerperrate'], ['cycletime', 'samplepowerperrate']]
        d['postcombofcns']=[self.nofilterfill, self.nofilterfill, self.nofilterfill, self.nofilterfill]
        d['parcombofcns']=[[self.nofilterfill], [self.nofilterfill], [self.timepartfilterfill, self.filterfill], [self.timepartfilterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)
        
        d=self.calclayoutgen('Dfit', ['c(t)+eT+fT^2+gT^3+hT^4', 'c(t)+eT+fT^2+hT^4','c(t)+et+ft^2+gt^3+ht^4', 'c(t)+et+ft^2+gt^3+ht^4+it^5'])#'c(t)+k<dT>+aT4'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='FITPARS_samplepowerperrate'
        d['fcns']=[FIT_T4, FIT_T0124, FIT_t4, FIT_t5]#pieceC_T4_intdT]
        d['parnames']=[['c','T', 'D'], ['c','T', 'D'], ['c','t', 'D'], ['c','t', 'D']]#['c','T', 'dT', 'D']]
        d['segdkeys']=[['cyclepartition','sampletemperature', 'samplepowerperrate'], ['cyclepartition','sampletemperature', 'samplepowerperrate'], ['cyclepartition','cycletime', 'samplepowerperrate'], ['cyclepartition','cycletime', 'samplepowerperrate']]#['cyclepartition', 'sampletemperature', 'sampleheatrate', 'samplepowerperrate']]
        d['postcombofcns']=[self.nofilterfill, self.nofilterfill, self.nofilterfill, self.nofilterfill]
        d['parcombofcns']=[[self.fitfilterfill, self.fitfilterfill, self.filterfill], [self.fitfilterfill, self.fitfilterfill, self.filterfill], [self.fitfilterfill, self.fitfilterfill, self.filterfill], [self.fitfilterfill, self.fitfilterfill, self.filterfill]]#[self.fitfilterfill, self.fitfilterfill, self.integfilterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved1)
        
        return 2
    
    def QUCsetup(self):
        d=self.calclayoutgen('C', ['D-Dfit'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='sampleheatcapacity'
        d['fcns']=[C_Dfit]
        d['parnames']=[['D', 'pt']]
        d['segdkeys']=[['samplepowerperrate', 'cyclepartition']]
        d['postcombofcns']=[self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.fitevalfilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)


#        d=self.calclayoutgen('Q', ['Qref','Qo+kA<T>+esA<T4>'])
#        self.parLayout.addWidget(d['widget'])
#        d['savename']='sampleheatloss'
#        d['fcns']=[Q_ref, Q_T]
#        d['parnames']=[['T'], ['T']]
#        d['segdkeys']=[['sampletemperature'], ['sampletemperature']]
#        d['postcombofcns']=[self.filterfill, self.filterfill]
#        d['parcombofcns']=[[self.refpathfilterfill], [self.Qmodelfilterfill]]#these filters need to be worked on
#        d['slider'].setMaximum(len(d['parnames'])-1)
#        self.pardlist+=[copy.copy(d)]
#        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)
#        
#        d=self.calclayoutgen('U', ['Uref','Q/S'])
#        self.parLayout.addWidget(d['widget'])
#        d['savename']='sampleheatlossperrate'
#        d['fcns']=[U_ref, U_QS]
#        d['parnames']=[['T'], ['Q', 'S']]
#        d['segdkeys']=[['sampletemperature'], ['sampleheatloss', 'sampleheatrate']]
#        d['postcombofcns']=[self.filterfill, self.filterfill]
#        d['parcombofcns']=[[self.refpathfilterfill], [self.filterfill, self.filterfill]]
#        d['slider'].setMaximum(len(d['parnames'])-1)
#        self.pardlist+=[copy.copy(d)]
#        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved1)
#        
#        d=self.calclayoutgen('C', ['D-U','(P-Q)/S'])
#        self.parLayout.addWidget(d['widget'])
#        d['savename']='sampleheatcapacity'
#        d['fcns']=[C_DU, C_PQS]
#        d['parnames']=[['D', 'U'], ['P', 'Q', 'S']]
#        d['segdkeys']=[['samplepowerperrate', 'sampleheatlossperrate'], ['samplepower', 'sampleheatloss', 'sampleheatrate']]
#        d['postcombofcns']=[self.filterfill, self.filterfill]
#        d['parcombofcns']=[[self.filterfill, self.filterfill], [self.filterfill, self.filterfill, self.filterfill]]
#        d['slider'].setMaximum(len(d['parnames'])-1)
#        self.pardlist+=[copy.copy(d)]
#        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved2)

        return 3
    
    def CTpksetup(self):
        d=self.calclayoutgen('Cpk', ['C(T) peaks'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='PROFILEANALYSIS_sampleheatcapacity_sampletemperature'
        d['fcns']=[CTpks_secder]
        d['parnames']=[['C', 'T']]
        d['segdkeys']=[['sampleheatcapacity', 'sampletemperature']]
        d['postcombofcns']=[self.nofilterfill]
        d['parcombofcns']=[[self.peakfitfilterfilterfill, self.griddatafilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)
        return 1

    def ACpksetup(self):
        d=self.calclayoutgen('FFTI', ['Windowed FFT'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='WinFFT_current'
        d['fcns']=[WinFFTI]
        d['parnames']=[['I']]
        d['segdkeys']=[['samplecurrent']]
        d['postcombofcns']=[self.nofilterfill]
        d['parcombofcns']=[[self.winfftfilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)
        
        d=self.calclayoutgen('FFTV', ['Windowed FFT'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='WinFFT_voltage'
        d['fcns']=[WinFFTV]
        d['parnames']=[['V']]
        d['segdkeys']=[['samplevoltage']]
        d['postcombofcns']=[self.nofilterfill]
        d['parcombofcns']=[[self.winfftfilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved1)
        
        d=self.calclayoutgen('FFTfV', ['Windowed FFT'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='WinFFT_filteredvoltage'
        d['fcns']=[WinFFTfV]
        d['parnames']=[['fV']]
        d['segdkeys']=[['samplefilteredvoltage']]
        d['postcombofcns']=[self.nofilterfill]
        d['parcombofcns']=[[self.winfftfilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved2)
        
        d=self.calclayoutgen('LIAI', ['LIA math'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='LIAharmonics_current'
        d['fcns']=[LIAI]
        d['parnames']=[['I']]
        d['segdkeys']=[['samplecurrent']]
        d['postcombofcns']=[self.nofilterfill]
        d['parcombofcns']=[[self.liafilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved3)
        
        d=self.calclayoutgen('LIAV', ['LIA math'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='LIAharmonics_voltage'
        d['fcns']=[LIAV]
        d['parnames']=[['V']]
        d['segdkeys']=[['samplevoltage']]
        d['postcombofcns']=[self.nofilterfill]
        d['parcombofcns']=[[self.liafilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved4)
        
        d=self.calclayoutgen('LIfV', ['LIA math'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='LIAharmonics_filteredvoltage'
        d['fcns']=[LIAfV]
        d['parnames']=[['fV']]
        d['segdkeys']=[['samplefilteredvoltage']]
        d['postcombofcns']=[self.nofilterfill]
        d['parcombofcns']=[[self.liafilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved5)
        return 6
        
    def AC_RTPSDsetup(self):
        d=self.calclayoutgen('R', [''])
        self.parLayout.addWidget(d['widget'])
        d['savename']='sampleresistance'
        d['fcns']=[R_fftIV, R_fftIfV, R_liaIV, R_liaIfV]
        d['parnames']=[['fftI', 'fftV'], ['fftI', 'fV'], ['liaI', 'liaV'], ['liaI', 'fV']]
        d['segdkeys']=[['WinFFT_current', 'WinFFT_voltage'], ['WinFFT_current', 'WinFFT_filteredvoltage'], ['LIAharmonics_current', 'LIAharmonics_voltage'], ['LIAharmonics_current', 'LIAharmonics_filteredvoltage']]
        d['postcombofcns']=[self.filterfill]*4
        d['parcombofcns']=[[self.acvfilterfill, self.filterfill]]*4
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved0)
        
        d=self.calclayoutgen('T', [''])
        self.parLayout.addWidget(d['widget'])
        d['savename']='sampletemperature'
        d['fcns']=[T_fftIV, T_fftIfV, T_liaIV, T_liaIfV]
        d['parnames']=[['fftI', 'fftV'], ['fftI', 'fV'], ['liaI', 'liaV'], ['liaI', 'fV']]
        d['segdkeys']=[['WinFFT_current', 'WinFFT_voltage'], ['WinFFT_current', 'WinFFT_filteredvoltage'], ['LIAharmonics_current', 'LIAharmonics_voltage'], ['LIAharmonics_current', 'LIAharmonics_filteredvoltage']]
        d['postcombofcns']=[self.filterfill]*4
        d['parcombofcns']=[[self.acvfilterfill, self.filterfill]]*4
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved1)
        
        d=self.calclayoutgen('P', ['IV', 'I2R'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='samplepower'
        d['fcns']=[P_fftIV, P_fftIfV, P_liaIV, P_liaIfV]
        d['parnames']=[['fftI', 'fftV'], ['fftI', 'fV'], ['liaI', 'liaV'], ['liaI', 'fV']]
        d['segdkeys']=[['WinFFT_current', 'WinFFT_voltage'], ['WinFFT_current', 'WinFFT_filteredvoltage'], ['LIAharmonics_current', 'LIAharmonics_voltage'], ['LIAharmonics_current', 'LIAharmonics_filteredvoltage']]
        d['postcombofcns']=[self.filterfill]*4
        d['parcombofcns']=[[self.acvfilterfill, self.filterfill]]*4
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved2)
        
        d=self.calclayoutgen('S', ['dT/dt'])#, 'avedT/dt'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='sampleheatrate'
        d['fcns']=[S_T]#, S_Tavesl]
        d['parnames']=[['T']]#, ['dTdt']]
        d['segdkeys']=[['sampletemperature']]#, ['sampletemperature']]
        d['postcombofcns']=[self.filterfill]#, self.filterfill]
        d['parcombofcns']=[[self.derfilterfill]]#, [self.derfilterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved3)
        
        d=self.calclayoutgen('D', ['P/S'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='samplepowerperrate'
        d['fcns']=[D_PS]
        d['parnames']=[['P', 'S']]
        d['segdkeys']=[['samplepower', 'sampleheatrate']]
        d['postcombofcns']=[self.filterfill]
        d['parcombofcns']=[[self.filterfill, self.filterfill]]
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved4)
        
        d=self.calclayoutgen('mCp', ['HeatCap'])
        self.parLayout.addWidget(d['widget'])
        d['savename']='acheatcapacity'
        d['fcns']=[mCp_fftVIPdT, mCp_fftfVIPdT, mCp_liaVIPdT, mCp_liafVIPdT]
        d['parnames']=[['fftV', 'fftI', 'P', 'dT'], ['ffV', 'fftI', 'P', 'dT'], ['liaV', 'liaI', 'P', 'dT'], ['lfV', 'liaI', 'P', 'dT']]
        d['segdkeys']=[['WinFFT_voltage', 'WinFFT_current', 'samplepower', 'sampleheatrate'], ['WinFFT_filteredvoltage', 'WinFFT_current', 'samplepower', 'sampleheatrate'], ['LIAharmonics_voltage', 'LIAharmonics_filteredcurrent', 'samplepower', 'sampleheatrate'], ['LIAharmonics_filteredvoltage', 'LIAharmonics_filteredcurrent', 'samplepower', 'sampleheatrate']]
        d['postcombofcns']=[self.filterfill]*4
        d['parcombofcns']=[[self.acvfilterfill, self.filterfill, self.filterfill, self.filterfill]]*4
        d['slider'].setMaximum(len(d['parnames'])-1)
        self.pardlist+=[copy.copy(d)]
        QObject.connect(d['slider'], SIGNAL("sliderReleased()"), self.slidermoved5)
        
        return 6
        
    def ExitRoutine(self):
        if self.saveCheckBox.isChecked():
            self.savefilters()
            if self.saveSCrecipe():
                saveSCcalculations(self.h5path, self.h5expname, self.h5hpname, self.hpsegdlist, self.recname)
            else:
                QMessageBox.warning(self,"ERROR",  "data not saved")

    def savefilters(self):
        savefilters(self.h5path, self.h5expname, self.filterd)
    
    def saveSCrecipe(self):
        self.recname=str(self.recLineEdit.text())
        if len(self.recname)==0:
            QMessageBox.warning(self,"ERROR",  "recipe name not valid")
            return False
        savefilters(self.h5path, self.h5expname, self.filterd)
        fcns=[]
        recdlist=[]
        for d in self.pardlist:
            recd={}
            ind=d['slider'].sliderPosition()
            if d['fcns'][ind] is None:
                continue
            fcns+=[d['fcns'][ind].func_name]
            recd['segdkeys']=d['segdkeys'][ind]
            a, b=zip(*[[nam, str(cb.currentText())] for nam, cb in zip(d['parnames'][ind], d['parcombos'])])#names and cb put together to control which cb are used
            recd['parnames']=list(a)
            recd['filters']=list(b)
            recd['postfilter']=str(d['postcombo'].currentText())
            recd['savename']=d['savename']
            recdlist+=[recd]
        saveSCrecipe(self.h5path, self.h5expname, self.recname, fcns, recdlist)
        return True
        
    def calcall(self):
        if self.calctype=='FitPS':
            xdata=[]
            fitans=[]
            fitdata=[]
        partups=[('fild', self.filterd)]
        partups+=[(a, b) for a, b in zip(['h5path', 'h5expname', 'h5hpname'], [self.h5path, self.h5expname, self.h5hpname])]
        seginds=self.segcalcoptions[self.segComboBox.currentIndex()]
        for si in seginds:
            #since filters may have changed, delete filtered data from segdict
            for k in self.hpsegdlist[si].keys():
                if '~' in k:
                    del self.hpsegdlist[si][k]
            partup_seg=[('segd', self.hpsegdlist[si])]
            for d in self.pardlist:
                ind=d['slider'].sliderPosition()
                f=d['fcns'][ind]
                if f is None:
                    continue
                saven=d['savename']
                namsegkfilk=[(nam, (segk, str(cb.currentText()))) for nam, segk, cb in zip(d['parnames'][ind], d['segdkeys'][ind], d['parcombos'])]
                postfilk=str(d['postcombo'].currentText())
                print 'calculating %s for segment %d' %(saven, si)
                arr=f(**dict(partups+partup_seg+namsegkfilk))
                arr=performgenericfilter(arr, self.filterd[postfilk])
                self.hpsegdlist[si][saven]=arr
                print self.hpsegdlist[si]['cycletime'].shape, arr.shape
            if self.calctype=='FitPS':
                segd=self.hpsegdlist[si]
                parttime=segd['cyclepartition']
                fitf=FitFcnLibrary[f.func_name]
                pns=d['parnames'][ind] # this assumes that the fitting was the last thing done in the above pardlist loop
                ks=d['segdkeys'][ind]
                xdata+=[ct[pt>=0] for ct, pt in zip(segd['cycletime'], parttime)]
                fitans+=[fitf(**dict([('p', p)]+[(pn, segd[k][i][pt>=0]) for pn, k in zip(pns, ks) if pn in fitf.func_code.co_varnames[:fitf.func_code.co_argcount]])) for i, (p, pt) in enumerate(zip(arr, parttime))]
                fitdata+=[y[pt>=0] for y, pt in zip(segd['samplepowerperrate'], parttime)]
        if self.calctype=='RTPSD' or self.calctype=='QUC':
            if self.plotdialog is None:
                self.plotdialog=analysisviewerDialog(None, self.hpsegdlist, hpname=self.h5hpname)
                self.plotdialog.drawall()
                self.plotdialog.show()
                self.plotdialog.activateWindow()
            else:
                self.plotdialog.hpsegdlist=self.hpsegdlist
                self.plotdialog.initwidgets()
                self.plotdialog.drawall()
                self.plotdialog.show()
                self.plotdialog.activateWindow()
        elif self.calctype=='FitPS':
            if not self.plotdialog is None:
                self.plotdialog.close()
            cols=['k', 'y', 'b', 'm', 'g', 'c', 'r']*10
            self.plotdialog=simpleplotDialog(None, fitans+fitdata, xdata=xdata+xdata, style=[c+'--' for c in cols[:len(xdata)]]+[c+'.' for c in cols[:len(xdata)]])
            self.plotdialog.show()
            self.plotdialog.activateWindow()
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
    def slidermoved5(self):
        self.updateparwidgets(5)
        
    def updateparwidgets(self, widgetind=None):
        if widgetind is None:
            widgetind=range(len(self.pardlist))
        else:
            if widgetind>=len(self.pardlist):
                return
            widgetind=[widgetind]
        for wind in widgetind:
            d=self.pardlist[wind]
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
    
    def importfilters(self, h5exp=None):
        if h5exp is None:
            h5file, expp=experimentgrppaths(self.h5path)
            paths=[]
            names=[]
            for p in expp:
                if 'filter' in h5file[p]['analysis']:
                    #pnts+=[h5file[p]['analysis']['filter']]
                    names+=[p.rpartition('/')[2]]
            idialog=selectorDialog(self, names, title='Select experiment whose filters will be imported')
            if not idialog.exec_():
                h5file.close()
                return
            newfilterd=getfilterdict(h5file, names[idialog.index])
            h5file.close()
        else:
            h5file=h5py.File(self.h5path, mode='r')
            newfilterd=getfilterdict(h5file, h5exp)
            h5file.close()
        if not newfilterd:
            print 'error importing filters from ', h5exp
            return
        for nam, d in newfilterd.iteritems():
            self.filterd[nam]=d
        self.fillfilterComboBox()

    def editfilter(self, editfilter=True):
        nam=str(self.filterComboBox.currentText())
        d=copy.deepcopy(self.filterd[nam])
        d['name']=nam
        edit=True
        while edit:
            edit=False
            if not filterattredit(self, d):
                return
            self.filterdicttypechanger(d)
            if not (isinstance(d['name'], str) and len(d['name'])>0):
                edit=True
                QMessageBox.warning(self,"ERROR",  "'name' is not valid")
            if d['name'] in self.filterd.keys() and not (editfilter and d['name']==nam):
                edit=True
                QMessageBox.warning(self,"ERROR",  "'name' already exists")
            if 'SGderiv' in d.keys() and (d['SGderiv']<0 or d['SGderiv']>1):
                edit=True
                QMessageBox.warning(self,"ERROR",  "only 0th or 1st deriv supported")
        if editfilter:
            del self.filterd[nam]
        nam=d['name']
        del d['name']
        self.filterd[nam]=d
        self.fillfilterComboBox()

    def newfilter(self):
        self.editfilter(editfilter=False)

    def fillfilterComboBox(self):
        self.filterComboBox.clear()
        for counter, nam in enumerate(self.filterd.keys()):
            self.filterComboBox.insertItem(counter, nam)
        self.updateparwidgets()
        
    def filterfill(self, cb, deriv=0, reqkeys=[]):
        dfltnam=str(cb.currentText())
        if dfltnam=='':
            dfltnam='None'
        i=0
        cb.clear()
        counter=0
        for nam, d in self.filterd.iteritems():
            if False in [k in d.keys() for k in reqkeys]:
                continue
            if ('SGderiv' in d.keys() and d['SGderiv']!=deriv) or (deriv>0 and not 'SGderiv' in d.keys()):
                continue
            cb.insertItem(counter, nam)
            if nam==dfltnam:
                i=counter
            counter+=1
        cb.setCurrentIndex(i)

    def derfilterfill(self, cb):
        self.filterfill(cb, deriv=1)
    
    def nofilterfill(self, cb):
        self.filterfill(cb, reqkeys=['None'])
    
    def fitfilterfill(self, cb):#fitpars is a list of stareting values for fitparameters applicable to the gicen variable. i.e. for a 4th order polynomial of T, fitpar should be a list of 4 coefficients. It is up to the user to get the number and order of the parameters correct.
        self.filterfill(cb, reqkeys=['fitpars'])
        
    def integfilterfill(self, cb):
        self.filterfill(cb, reqkeys=['integwindow_s', 'fitpars'])
    
    def fitevalfilterfill(self, cb):
        self.filterfill(cb, reqkeys=['interp_kind', 'extrap_order'])
    
    def griddatafilterfill(self, cb):
        self.filterfill(cb, reqkeys=['gridinterval'])
    
    def peakfitfilterfilterfill(self, cb):
        self.filterfill(cb, reqkeys=['useredit','pkfcn', 'critpeakheight', 'critsep', 'firstdernpts', 'firstderorder', 'secdernpts', 'secderorder', 'critcurve', 'pospeaks', 'negpeaks'])

    def timepartfilterfill(self, cb):
        self.filterfill(cb, reqkeys=['numpartitions'])#fitpars is going to determine the number of parititons so it is necessary. if there will not be a contant term in the fit, then will will still need to be icnluded here and excluded in the fitfcn defintion
        #time partition functions: timepart_user, timepart_none, timepart_peakid (peakid not developed as of April2011)
    def refpathfilterfill(self, cb):
        self.filterfill(cb, reqkeys=['REFh5path', 'REFalignment'])
    
    def liafilterfill(self, cb):
        self.filterfill(cb, reqkeys=['n1wcycles_window'])
    def winfftfilterfill(self, cb):
        self.filterfill(cb, reqkeys=['n1wcycles_window'])
    
    def acvfilterfill(self, cb):
        self.filterfill(cb, reqkeys=['harmonic'])
        
    def dfltfitfilterdicts(self, deriv=0):
        return [{'name':'timepart', 'numpartitions':1}, #contants and any other piecewise parameters have this many segments and these starting values\
                    {'name':'fit_T4','fitpars':[1.e-8, 1.e-10, 1.e-12, 1.e-14]}, \
#                    {'name':'fit_T4','fitpars':[1.e-14]}, \
#                    {'name':'fit_dT','fitpars':[1.e-8]}, \
                    {'name':'fit_c','fitpars':[1.e-6]}, \
                ]
    
    def dfltfilterdict(self, deriv=0):
        return {'OLnpts':1, \
                'OLnsig':1.5, \
                'OLgappts':0, \
                'SGnpts':40, \
                'SGorder':1, \
                'SGderiv':deriv, \
                'SGbin':0, \
                }
    
    def griddatafilterdflt(self):
        return {'gridinterval':0.2}
    
    def peakfitfilterdflt(self):
        return {'useredit':1,'pkfcn':'GaussHalfLorentz','critpeakheight':5.e-7, 'critsep':20., \
        'firstdernpts':10, 'firstderorder':1, 'secdernpts':20, 'secderorder':1, \
        'critcurve':None, 'pospeaks':1, 'negpeaks':1}
    
    def acfilterdflt(self):
        return {'n1wcycles_window':10}#, 'harmonics':[0, 1, 2, 3]}

    def acvfilterdflt(self):
        d=dfltfilterdict
        d['SGnpts']=0
        d['harmonic']=2
        return d
    
    def Nonefilterdict(self):
        return {'None':None, \
                }
    
    def fitevalfilterdflt(self):
        return {'interp_kind':'linear', 'extrap_order':1, 'interppars':0}
        
    def Qfilterdict(self, h5pathstr):
        attrs=geth5attrs(self.h5path, '/'.join((h5pathstr, 'sampleheatloss')))
        d=self.Nonefilterdict()
        d['kA']=attrs['kA']
        d['esA']=attrs['esA']
        return d
        
    def reffilterdict(self, h5pathstr):
        d=self.Nonefilterdict()
        d['REFh5path']=h5pathstr
        d['REFalignment']='sampletemperature'
        return d
        
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
        
        eqlablayout=QVBoxLayout()
        eqLabellist=[]
        if len(eqnames)>0:
            font = QFont()
            font.setPointSize(max(8, min(14, 29//len(eqnames))))
        for eqn in eqnames:
            eqLabel0 = QLabel(gridLayoutWidget)
            eqLabel0.setFont(font)
            #eqLabel0.setObjectName("eqLabel0")
            eqLabel0.setText(eqn)
            eqlablayout.addWidget(eqLabel0)
            eqLabellist+=[eqLabel0]
        gridLayout.addLayout(eqlablayout, 1, 4, 4, 1)
        
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
        font.setPointSize(24)
        varLabel.setFont(font)
        varLabel.setLayoutDirection(Qt.LeftToRight)
        varLabel.setAutoFillBackground(False)
        varLabel.setObjectName("varLabel")
        varLabel.setText(varname)
        gridLayout.addWidget(varLabel, 1, 0, 2, 1)
#        line_3 = QFrame(gridLayoutWidget)
#        line_3.setFrameShape(QFrame.HLine)
#        line_3.setFrameShadow(QFrame.Sunken)
#        line_3.setObjectName("line_3")
#        gridLayout.addWidget(line_3, 5, 0, 1, 8)
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
        d['eqlabels']=eqLabellist
        d['postcombo']=postfilterComboBox
        d['name']=varname
        d['eqnames']=eqnames
        return d


class SCanalysisDialog(QDialog):
    def __init__(self, parent, h5path, h5expname, h5hpdflt=None):
        super(SCanalysisDialog, self).__init__(parent)
        self.parent=parent
        self.h5path=h5path
        self.h5expname=h5expname

        reccopyButton=QPushButton()
        reccopyButton.setText("Copy recipe from other experiment\n(overwrite recipe and filters with same name)")
        QObject.connect(reccopyButton, SIGNAL("pressed()"), self.reccopy)
        
        self.recComboBox=QComboBox()
        
        self.hpComboBox=QComboBox()
        self.hpComboBox.insertItem(0, 'ALL heat programs')
        self.hpnamelist=[]
        h5file, hppl=experimenthppaths(self.h5path, self.h5expname)
        h5file.close()
        i=0
        for count, hp in enumerate([hpp.rpartition('/')[2] for hpp in hppl]):
            self.hpComboBox.insertItem(count+1, hp)
            self.hpnamelist+=[hp]
            if hp==h5hpdflt:
                i=count+1
        self.hpComboBox.setCurrentIndex(i)
        QObject.connect(self.hpComboBox,SIGNAL("activated(QString)"),self.hpchanged)
        
        self.segComboBox=QComboBox()
        segLabel=QLabel()
        segLabel.setText('segment(s)\nfor calculation')
        
        calcallButton=QPushButton()
        calcallButton.setText("Calculate all quantities")
        QObject.connect(calcallButton, SIGNAL("pressed()"), self.calcall)
        
        
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        
        
        mainlayout=QGridLayout()
        
        recLabel=QLabel()
        recLabel.setText('SC analysis\nrecipe:')
        hpLabel=QLabel()
        hpLabel.setText('heat\nprogram:')
        
        mainlayout.addWidget(hpLabel, 0, 0)
        mainlayout.addWidget(self.hpComboBox, 0, 1)
        mainlayout.addWidget(segLabel, 1, 0)
        mainlayout.addWidget(self.segComboBox, 1, 1)
        mainlayout.addWidget(recLabel, 2, 0)
        mainlayout.addWidget(self.recComboBox, 2, 1)
        mainlayout.addWidget(reccopyButton, 3, 0, 1, 2)
        mainlayout.addWidget(calcallButton, 4, 0, 1, 2)
        mainlayout.addWidget(self.buttonBox, 5, 1)
        self.setLayout(mainlayout)
        
        self.setWindowTitle('SC analysis recipe editor')
        self.plotdialog=None
        self.hpchanged()
        self.fillrecComboBox()
        
    
    def hpchanged(self):
        self.segComboBox.clear()
        hplist=self.readhpcb()
        self.segcalcnames=[]
        self.segcalcoptions=[]
        if len(hplist)==1:
            hpsegdlist=CreateHeatProgSegDictList(self.h5path, self.h5expname, hplist[0])
            ramp=[]
            soak=[]
            for count, d in enumerate(hpsegdlist):
                self.segcalcnames+=['segment %d : %s' %(count, d['segmenttype'])]
                self.segcalcoptions+=[numpy.array([count])]
                if d['segmenttype']=='ramp':
                    ramp+=[count]
                if d['segmenttype']=='soak':
                    soak+=[count]
            self.segcalcnames+=['all ramp']
            self.segcalcoptions+=[numpy.array(ramp)]
            self.segcalcnames+=['all soak']
            self.segcalcoptions+=[numpy.array(soak)]
            self.segcalcnames+=['all ramp+soak']
            self.segcalcoptions+=[numpy.array(ramp+soak)]
            self.segcalcnames+=['all segs']
            self.segcalcoptions+=[numpy.array(range(len(hpsegdlist)))]
        else:
            hpsegdlist=CreateHeatProgSegDictList(self.h5path, self.h5expname, hplist[0])
            for count, d in enumerate(hpsegdlist):
                self.segcalcnames+=['segment %d : %s' %(count, d['segmenttype'])]
                self.segcalcoptions+=[numpy.array([count])]
            self.segcalcnames+=['all ramp']
            self.segcalcoptions+=[['ramp']]
            self.segcalcnames+=['all soak']
            self.segcalcoptions+=[['soak']]
            self.segcalcnames+=['all ramp+soak']
            self.segcalcoptions+=[['ramp', 'soak']]
            self.segcalcnames+=['all segs']
            self.segcalcoptions+=[segtypes]
        for count, nam in enumerate(self.segcalcnames):
            self.segComboBox.insertItem(count, nam)
            
    def readhpcb(self):
        i=self.hpComboBox.currentIndex()
        if i==0:
            return self.hpnamelist
        else:
            return [self.hpnamelist[i-1]]
    def reccopy(self):
        h5file, expp=experimentgrppaths(self.h5path)
        paths=[]
        names=[]
        for p in expp:
            if 'SCrecipe' in h5file[p]['analysis']:
                #pnts+=[h5file[p]['analysis']['filter']]
                nam=p.rpartition('/')[2]
                if nam!=self.h5expname:
                    names+=[nam]
        h5file.close()
        idialog=selectorDialog(self, names, title='Select experiment whose recipes+filters will be imported')
        if not idialog.exec_():
            return
        copySCrecipes(self.h5path, self.h5expname, names[idialog.index])
        self.fillrecComboBox()
    
    def fillrecComboBox(self):
        self.recComboBox.clear()
        h5file=h5py.File(self.h5path, mode='r+')#r+ in case nSCrecipe needs to be created
        h5srcrec=getSCrecipegrp(h5file, self.h5expname)
        for pnt in h5srcrec.itervalues():
            if isinstance(pnt, h5py.Group):
                nam=pnt.name.rpartition('/')[2]
                self.recComboBox.insertItem(99, nam)
        h5file.close()
        
    def getseginds(self, hpsegdlist):
        segindsortypes=self.segcalcoptions[self.segComboBox.currentIndex()]
        if isinstance(segindsortypes[0], int):
            return segindsortypes
        return [count for count, d in enumerate(hpsegdlist) if d['segmenttype'] in segindsortypes]
        
    def calcall(self):
        print 'starting calcall'
        self.plotdialog=None
        self.recname=str(self.recComboBox.currentText())
        hplist=self.readhpcb()
        h5file=h5py.File(self.h5path, mode='r')
        filterd=getfilterdict(h5file, self.h5expname)
        f_saven_namsegkfilk_postfilk=getSCrecipe(h5file, self.h5expname, self.recname)
        h5file.close()
        partups=[('fild', filterd)]
        partups+=[(a, b) for a, b in zip(['h5path', 'h5expname'], [self.h5path, self.h5expname])]
        for h5hpname in hplist[0:]:#***
            hpsegdlist=CreateHeatProgSegDictList(self.h5path, self.h5expname, h5hpname)
            seginds=self.getseginds(hpsegdlist)
            hppartup=[('h5hpname', h5hpname)]
            xdata=[]#won't be used if not in fit recipe
            fitans=[]
            fitdata=[]
            for si in seginds:
                partup_seg=[('segd',hpsegdlist[si])]
                for f, saven, namsegkfilk, postfilk in f_saven_namsegkfilk_postfilk:
                    print 'calculating %s for segment %d in hp %s' %(saven, si, h5hpname)
                    f=eval(f)
                    arr=f(**dict(partups+hppartup+partup_seg+namsegkfilk))
                    arr=performgenericfilter(arr, filterd[postfilk])
                    hpsegdlist[si][saven]=arr
                if 'FIT' in f.func_name:
                    segd=hpsegdlist[si]
                    parttime=segd['cyclepartition']
                    fitf=FitFcnLibrary[f.func_name]
                    xdata+=[ct[pt>=0] for ct, pt in zip(segd['cycletime'], parttime)]
                    fitans+=[fitf(**dict([('p', p)]+[(nam, segd[segk][i][pt>=0]) for (nam, (segk, filter)) in namsegkfilk if nam in fitf.func_code.co_varnames[:fitf.func_code.co_argcount]])) for i, (p, pt) in enumerate(zip(arr, parttime))]
                    fitdata+=[y[pt>=0] for y, pt in zip(segd['samplepowerperrate'], parttime)]
            if 'FIT' in f.func_name:#this is not defitie, just up to function naming conventions
                if not self.plotdialog is None:
                    self.plotdialog.close()
                cols=['k', 'y', 'b', 'm', 'g', 'c', 'r']*10
                self.plotdialog=simpleplotDialog(None, fitans+fitdata, xdata=xdata+xdata, style=[c+'--' for c in cols[:len(xdata)]]+[c+'.' for c in cols[:len(xdata)]])
                self.plotdialog.show()
                self.plotdialog.activateWindow()
            elif 'CTpks_secder' in f.func_name:
                print 'post-calc plotting skipped'
            else:
                if self.plotdialog is None:
                    self.plotdialog=analysisviewerDialog(None, hpsegdlist, hpname=h5hpname)
                    self.plotdialog.drawall()
                    self.plotdialog.show()
                    self.plotdialog.activateWindow()
                else:
                    self.plotdialog.hpsegdlist=hpsegdlist
                    self.plotdialog.hpname=h5hpname
                    self.plotdialog.initwidgets()
                    self.plotdialog.drawall()
                    self.plotdialog.show()
                    self.plotdialog.activateWindow()
            saveSCcalculations(self.h5path, self.h5expname, h5hpname, hpsegdlist, self.recname)
#        if not newfilterd:
#            print 'error importing filters from ', h5exp
#            return
#        for nam, d in newfilterd.iteritems():
#            self.reciped[nam]=d
#        self.fillrecipeComboBox()

class analysisviewerDialog(QDialog):
    def __init__(self, parent, hpsegdlist, hpname='', markersize=2):
        super(analysisviewerDialog, self).__init__(parent)
        self.markersize=markersize
        self.hpsegdlist=hpsegdlist
        self.hpname=hpname
        
        self.plotsdlist=[]
        plotslayout=QGridLayout()
        rows=[0, 1]
        cols=[0, 1, 2]
        segcbfcns=[lambda garb=None: self.segcgchanged(ploti=0), lambda garb=None: self.segcgchanged(ploti=1), lambda garb=None: self.segcgchanged(ploti=2), lambda garb=None: self.segcgchanged(ploti=3), lambda garb=None: self.segcgchanged(ploti=4), lambda garb=None: self.segcgchanged(ploti=5)]
        self.drawfcns=[lambda : self.draw(ploti=0), lambda : self.draw(ploti=1), lambda : self.draw(ploti=2), lambda : self.draw(ploti=3), lambda : self.draw(ploti=4), lambda : self.draw(ploti=5)]
        for row in rows:
            for col in cols:
                playout=QGridLayout()
                d=self.plotmenuwidget()
                QObject.connect(d['draw'], SIGNAL("pressed()"), self.drawfcns[row*len(cols)+col])
                QObject.connect(d['segcb'],SIGNAL("activated(QString)"),segcbfcns[row*len(cols)+col])
                plotw=plotwidget(self)
                d['plotw']=plotw
                playout.addWidget(d['widget'], 0, 0)
                playout.addWidget(plotw, 1, 0, 3, 1)
                plotslayout.addLayout(playout, row, col)
                self.plotsdlist+=[d]
        self.setLayout(plotslayout)
        self.initwidgets()
    
    def drawall(self):
        for f in self.drawfcns:
            try:
                f()
            except:
                print 'error during auto plotting'

    def initwidgets(self):
        self.setWindowTitle('Multi-plot view of calorimetry calculations: %s' %self.hpname)
        self.segcalcoptions=[]
        self.segcalcnames=[]
        ramp=[]
        soak=[]
        for count, d in enumerate(self.hpsegdlist):
            self.segcalcnames+=['segment %d : %s' %(count, d['segmenttype'])]
            self.segcalcoptions+=[numpy.array([count])]
            if d['segmenttype']=='ramp':
                ramp+=[count]
            if d['segmenttype']=='soak':
                soak+=[count]
        self.segcalcnames+=['all ramp']
        self.segcalcoptions+=[numpy.array(ramp)]
        self.segcalcnames+=['all soak']
        self.segcalcoptions+=[numpy.array(soak)]
        self.segcalcnames+=['all ramp+soak']
        self.segcalcoptions+=[numpy.array(ramp+soak)]
        self.segcalcnames+=['all segs']
        self.segcalcoptions+=[numpy.array(range(len(self.hpsegdlist)))]
        ycbdflts=['samplecurrent', 'samplevoltage', 'sampletemperature', 'samplepower', 'sampleheatrate', 'samplepowerperrate']
        xcbdflts=['cycletime', 'cycletime', 'cycletime', 'cycletime', 'cycletime', 'sampletemperature']
        temp=len(self.segcalcoptions)-1
        segdflt=[temp, temp, temp-1, temp-1, temp-1, temp-1]
        for ploti, (d, yd, sd) in enumerate(zip(self.plotsdlist, ycbdflts, segdflt)):
            for count, nam in enumerate(self.segcalcnames):
                d['segcb'].insertItem(count, nam)
            d['segcb'].setCurrentIndex(sd)
            self.segcgchanged(ploti, ydflt=yd)
        return
    
    def segcgchanged(self, ploti, ydflt=None, xdflt='cycletime'):
        d=self.plotsdlist[ploti]
        d['xcb'].clear()
        d['ycb'].clear()
        seginds=self.segcalcoptions[d['segcb'].currentIndex()]
        if len(seginds)==0:
            d['shortkeyd']={}
            return
        segkl=list(set.intersection(*[set([k for k in self.hpsegdlist[si].keys() if isinstance(self.hpsegdlist[si][k], numpy.ndarray)]) for si in seginds]))
        shortkeyd=dict([((i<26 and (chr(65+i), ) or (chr(97+i%26)+chr(97+i%26+i//26-1), ))[0], k) for i, k in enumerate(segkl)])
        d['shortkeyd']=shortkeyd
        xind=0
        yind=0
        for count, (k, segk) in enumerate(d['shortkeyd'].iteritems()):
            if xdflt==segk:
                xind=count
            if ydflt==segk:
                yind=count
            d['xcb'].insertItem(count, ':'.join((k, segk)))
            d['ycb'].insertItem(count, ':'.join((k, segk)))
        d['xcb'].setCurrentIndex(xind)
        d['ycb'].setCurrentIndex(yind)

    def draw(self, ploti):
        d=self.plotsdlist[ploti]
        d['plotw'].axes.cla()
        seginds=self.segcalcoptions[d['segcb'].currentIndex()]
        xstr=str(d['xle'].text()).strip()
        if len(xstr)>0:
            for (k, segk) in d['shortkeyd'].iteritems():
                xstr=xstr.replace(k, "segd['%s']" %segk)
            xfcn=lambda segd:eval(xstr)
        else:
            xsegk=str(d['xcb'].currentText()).partition(':')[2]
            xfcn=lambda segd:segd[xsegk]
        ystr=str(d['yle'].text()).strip()
        if len(ystr)>0:
            for (k, segk) in d['shortkeyd'].iteritems():
                ystr=ystr.replace(k, "segd['%s']" %segk)
            yfcn=lambda segd:eval(ystr)
        else:
            ysegk=str(d['ycb'].currentText()).partition(':')[2]
            yfcn=lambda segd:segd[ysegk]
        temp=self.hpsegdlist[0]['cycletime']
        xcyc=[numpy.array([], dtype=temp.dtype) for i in range(temp.shape[0])]
        ycyc=[numpy.array([], dtype=temp.dtype) for i in range(temp.shape[0])]

        for si in seginds:
            segd=self.hpsegdlist[si]
            xarr=xfcn(segd)
            xcyc=[numpy.append(xc, x) for xc, x in zip(xcyc, xarr)]
            yarr=yfcn(segd)
            ycyc=[numpy.append(yc, y) for yc, y in zip(ycyc, yarr)]
#            print xarr==yarr
#        print [a==b for a, b in zip(xcyc, ycyc)]
        for xc, yc in zip(xcyc, ycyc):
            d['plotw'].axes.plot(xc, yc, '.', ms=self.markersize)
        d['plotw'].axes.set_xlabel(str(d['xlable'].text()))
        d['plotw'].axes.set_ylabel(str(d['ylable'].text()))
        d['plotw'].fig.canvas.draw()
        
    def plotmenuwidget(self):
        gridLayoutWidget = QWidget()
        gridLayoutWidget.setGeometry(QRect(25, 185, 476, 80))
        gridLayoutWidget.setObjectName("gridLayoutWidget")
        gridLayout = QGridLayout(gridLayoutWidget)
        gridLayout.setObjectName("gridLayout")
        segLabel = QLabel(gridLayoutWidget)
        segLabel.setObjectName("segLabel")
        gridLayout.addWidget(segLabel, 0, 0, 1, 1)
        segComboBox = QComboBox(gridLayoutWidget)
        segComboBox.setObjectName("segComboBox")
        gridLayout.addWidget(segComboBox, 1, 0, 1, 1)
        xLabel = QLabel(gridLayoutWidget)
        xLabel.setObjectName("xLabel")
        gridLayout.addWidget(xLabel, 0, 1, 1, 1)
        yLabel = QLabel(gridLayoutWidget)
        yLabel.setObjectName("yLabel")
        gridLayout.addWidget(yLabel, 1, 1, 1, 1)
        xComboBox = QComboBox(gridLayoutWidget)
        xComboBox.setObjectName("xComboBox")
        gridLayout.addWidget(xComboBox, 0, 2, 1, 1)
        yComboBox = QComboBox(gridLayoutWidget)
        yComboBox.setObjectName("yComboBox")
        gridLayout.addWidget(yComboBox, 1, 2, 1, 1)
        xLineEdit = QLineEdit(gridLayoutWidget)
        xLineEdit.setObjectName("xLineEdit")
        gridLayout.addWidget(xLineEdit, 0, 3, 1, 1)
        yLineEdit = QLineEdit(gridLayoutWidget)
        yLineEdit.setObjectName("yLineEdit")
        gridLayout.addWidget(yLineEdit, 1, 3, 1, 1)
        xlabLineEdit = QLineEdit(gridLayoutWidget)
        xlabLineEdit.setObjectName("xlabLineEdit")
        gridLayout.addWidget(xlabLineEdit, 0, 4, 1, 1)
        ylabLineEdit = QLineEdit(gridLayoutWidget)
        ylabLineEdit.setObjectName("ylabLineEdit")
        gridLayout.addWidget(ylabLineEdit, 1, 4, 1, 1)
        drawButton = QPushButton(gridLayoutWidget)
        drawButton.setObjectName("drawButton")
        gridLayout.addWidget(drawButton, 0, 5, 1, 1)
        xLabel.setText('x:')
        yLabel.setText('y:')
        segLabel.setText('segment(s)')
        drawButton.setText('draw')
        #QMetaObject.connectSlotsByName()
        
        d={'widget':gridLayoutWidget}
        d['segcb']=segComboBox
        d['xcb']=xComboBox
        d['ycb']=yComboBox
        d['xle']=xLineEdit
        d['yle']=yLineEdit
        d['xlable']=xlabLineEdit
        d['ylable']=ylabLineEdit
        d['draw']=drawButton
        return d
        
        
def R_IV(segd, fild, I, V, h5path=None, h5expname=None, h5hpname=None):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    for tup in [I, V]:
        (segkey, filkey)=tup[:2]
        if not '~'.join(tup[:2]) in segd.keys():
            if len(tup)==3:
                arr=tup[2]
            else:
                arr=segd[segkey]
            segd['~'.join(tup[:2])]=performgenericfilter(arr, fild[filkey])
    I_=segd['~'.join(I[:2])]
    V_=segd['~'.join(V[:2])]
    inds=numpy.where(I_<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value. These modification will no be reflected in the segd values for the source data
    if len(inds[0])>0:
        I_=replacevalswithneighsin2nddim(I_, inds)
        V_=replacevalswithneighsin2nddim(V_, inds)
    return V_/I_

def T_IV(segd, fild, I, V, h5path, h5expname, h5hpname):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    for tup in [I, V]:
        (segkey, filkey)=tup[:2]
        if not '~'.join(tup[:2]) in segd.keys():
            if len(tup)==3:
                arr=tup[2]
            else:
                arr=segd[segkey]
            segd['~'.join(tup[:2])]=performgenericfilter(arr, fild[filkey])
    I_=segd['~'.join(I[:2])]
    V_=segd['~'.join(V[:2])]
    inds=numpy.where(I_<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value. These modification will no be reflected in the segd values for the source data
    if len(inds[0])>0:
        I_=replacevalswithneighsin2nddim(I_, inds)
        V_=replacevalswithneighsin2nddim(V_, inds)
    return temp_res(V_/I_, RoToAl[0], RoToAl[1], RoToAl[2])
    

def P_IV(segd, fild, I, V, h5path=None, h5expname=None, h5hpname=None):
    for tup in [I, V]:
        (segkey, filkey)=tup[:2]
        if not '~'.join(tup[:2]) in segd.keys():
            if len(tup)==3:
                arr=tup[2]
            else:
                arr=segd[segkey]
            segd['~'.join(tup[:2])]=performgenericfilter(arr, fild[filkey])
    I_=segd['~'.join(I[:2])]
    V_=segd['~'.join(V[:2])]
    return I_*V_

def P_IR(segd, fild, I, R, h5path=None, h5expname=None, h5hpname=None):
    for tup in [I, R]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
        #    segd['~'.join(tup)]/=dt
    I_=segd['~'.join(I)]
    R_=segd['~'.join(R)]
    return I_**2*R_
    
def S_IV(segd, fild, I, V, dIdt, dVdt, h5path, h5expname, h5hpname):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    dt=dt_h5(h5path, h5expname, h5hpname)
    for tup in [I, V, dIdt, dVdt]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
            segd['~'.join(tup)]/=dt
    I_=segd['~'.join(I)]
    V_=segd['~'.join(V)]
    dIdt_=segd['~'.join(dIdt)]
    dVdt_=segd['~'.join(dVdt)]
    inds=numpy.where(I_==0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value. These modification will no be reflected in the segd values for the source data
    if len(inds[0])>0:
        I_=replacevalswithneighsin2nddim(I_, inds)
        V_=replacevalswithneighsin2nddim(V_, inds)
        dIdt_=replacevalswithneighsin2nddim(dIdt_, inds)
        dVdt_=replacevalswithneighsin2nddim(dVdt_, inds)
    return dT_IVdIdV(I_, V_, dIdt_, dVdt_, RoToAl[0], RoToAl[2])
    
def S_T(segd, fild, dTdt, h5path, h5expname, h5hpname):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    dt=dt_h5(h5path, h5expname, h5hpname)
    for tup in [dTdt]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
            segd['~'.join(tup)]/=dt
    dTdt_=segd['~'.join(dTdt)]
    return dTdt_

def S_Tavesl(segd, fild, dTdt, h5path, h5expname, h5hpname):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    dt=dt_h5(h5path, h5expname, h5hpname)
    for tup in [dTdt]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
            segd['~'.join(tup)]/=dt
    dTdt_=segd['~'.join(dTdt)]
    return dTdt_

def D_IV(segd, fild, I, V, dIdt, dVdt, h5path, h5expname, h5hpname):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    dt=dt_h5(h5path, h5expname, h5hpname)
    for tup in [I, V, dIdt, dVdt]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
            segd['~'.join(tup)]/=dt
    I_=segd['~'.join(I)]
    V_=segd['~'.join(V)]
    dIdt_=segd['~'.join(dIdt)]
    dVdt_=segd['~'.join(dVdt)]
    inds=numpy.where((I_*dVdt_-V_*dIdt_)==0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value. These modification will no be reflected in the segd values for the source data
    if len(inds[0])>0:
        I_=replacevalswithneighsin2nddim(I_, inds)
        V_=replacevalswithneighsin2nddim(V_, inds)
        dIdt_=replacevalswithneighsin2nddim(dIdt_, inds)
        dVdt_=replacevalswithneighsin2nddim(dVdt_, inds)
    return D_IVdIdV(I_, V_, dIdt_, dVdt_, RoToAl[0], RoToAl[2])

def D_PS(segd, fild, P, S, h5path=None, h5expname=None, h5hpname=None):
    for tup in [P, S]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
        #    segd['~'.join(tup)]/=dt
    P_=segd['~'.join(P)]
    S_=segd['~'.join(S)]
    inds=numpy.where(S_==0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value. These modification will no be reflected in the segd values for the source data
    if len(inds[0])>0:
        P_=replacevalswithneighsin2nddim(P_, inds)
        S_=replacevalswithneighsin2nddim(S_, inds)
    return P_/S_

def Q_T(segd, fild, T, h5path, h5expname, h5hpname):
    T0=T0_h5(h5path, h5expname, h5hpname)
    for tup in [T]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
        #    segd['~'.join(tup)]/=dt
    Q=segd['~'.join(T)]#here the filter does the calculation
    return Q
#these need to be written Q_ref, Q_T, U_ref, U_QS, C_DU, C_PQS
def U_QS(segd, fild, Q, S, h5path=None, h5expname=None, h5hpname=None):
    for tup in [Q, dTdt]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
        #    segd['~'.join(tup)]/=dt
    Q_=segd['~'.join(Q)]
    dTdt_=segd['~'.join(dTdt)]
    inds=numpy.where(dTdt_==0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value. These modification will no be reflected in the segd values for the source data
    if len(inds[0])>0:
        Q_=replacevalswithneighsin2nddim(Q_, inds)
        dTdt_=replacevalswithneighsin2nddim(dTdt_, inds)
    return Q_/dTdt_
    
def C_DU(segd, fild, D, U, h5path, h5expname=None, h5hpname=None):
    for tup in [D, U]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
            if True in ['REF' in k for k in fild[filkey].keys()]:
                segd['~'.join(tup)]=performreferencesubtraction(segd, segkey, fild[filkey], h5path)
        #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
        #    segd['~'.join(tup)]/=dt
    #next few lines takes care of reference - this could be automated/batch code but for now leave it as 1-offs because of unique math
    delD_=segd['~'.join(D)]
    delU_=segd['~'.join(U)]
    return delD-delU
    

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

def FIT_T0124(segd, fild, c, T, D, h5path=None, h5expname=None, h5hpname=None):
    fitpars=[]
    for tup in [c, T, D]:
        (segkey, filkey)=tup
        if 'fitpars' in fild[filkey].keys():
            fitpars+=list(fild[filkey]['fitpars'])
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
    
    paritionedtime=segd['~'.join(c)]
    T_=segd['~'.join(T)]
    D_=segd['~'.join(D)]
    
    fp=[]
    ff=fitfcns()
    for pt, cycT, cycD in zip(paritionedtime, T_, D_):
        ff.genfit(FitFcnLibrary[FIT_T0124.__name__], fitpars, (pt[pt>=0], cycT[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[:50].mean()]*len(fitpars)
            tempT=cycT[-50:].mean()
            tempD=cycD[-50:].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempT**(i+1)
            ff.genfit(FitFcnLibrary[FIT_T0124.__name__], fitpars, (pt[pt>=0], cycT[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[-50:].mean()]*len(fitpars)
            tempT=cycT[:50].mean()
            tempD=cycD[:50].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempT**(i+1)
            ff.genfit(FitFcnLibrary[FIT_T0124.__name__], fitpars, (pt[pt>=0], cycT[pt>=0], cycD[pt>=0]))
        fp+=[ff.finalparams]
    fp=numpy.float32(fp)
    return fp
    
def FIT_T4(segd, fild, c, T, D, h5path=None, h5expname=None, h5hpname=None):
    fitpars=[]
    for tup in [c, T, D]:
        (segkey, filkey)=tup
        if 'fitpars' in fild[filkey].keys():
            fitpars+=list(fild[filkey]['fitpars'])
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
    
    paritionedtime=segd['~'.join(c)]
    T_=segd['~'.join(T)]
    D_=segd['~'.join(D)]
    
    fp=[]
    ff=fitfcns()
    for pt, cycT, cycD in zip(paritionedtime, T_, D_):
        ff.genfit(FitFcnLibrary[FIT_T4.__name__], fitpars, (pt[pt>=0], cycT[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[:50].mean()]*len(fitpars)
            tempT=cycT[-50:].mean()
            tempD=cycD[-50:].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempT**(i+1)
            ff.genfit(FitFcnLibrary[FIT_T4.__name__], fitpars, (pt[pt>=0], cycT[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[-50:].mean()]*len(fitpars)
            tempT=cycT[:50].mean()
            tempD=cycD[:50].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempT**(i+1)
            ff.genfit(FitFcnLibrary[FIT_T4.__name__], fitpars, (pt[pt>=0], cycT[pt>=0], cycD[pt>=0]))
        fp+=[ff.finalparams]
    fp=numpy.float32(fp)
    return fp

def FIT_t4(segd, fild, c, t, D, h5path=None, h5expname=None, h5hpname=None):
    fitpars=[]
    for tup in [c, t, D]:
        (segkey, filkey)=tup
        if 'fitpars' in fild[filkey].keys():
            fitpars+=list(fild[filkey]['fitpars'])
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
    
    paritionedtime=segd['~'.join(c)]
    t_=segd['~'.join(t)]
    D_=segd['~'.join(D)]
    
    fp=[]
    ff=fitfcns()
    for pt, cyct, cycD in zip(paritionedtime, t_, D_):
        ff.genfit(FitFcnLibrary[FIT_t4.__name__], fitpars, (pt[pt>=0], cyct[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[:50].mean()]*len(fitpars)
            tempt=cyct[-50:].mean()
            tempD=cycD[-50:].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempt**(i+1)
            ff.genfit(FitFcnLibrary[FIT_t4.__name__], fitpars, (pt[pt>=0], cyct[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[-50:].mean()]*len(fitpars)
            tempt=cyct[:50].mean()
            tempD=cycD[:50].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempt**(i+1)
            ff.genfit(FitFcnLibrary[FIT_t4.__name__], fitpars, (pt[pt>=0], cyct[pt>=0], cycD[pt>=0]))
        fp+=[ff.finalparams]
    fp=numpy.float32(fp)
    return fp

def FIT_t5(segd, fild, c, t, D, h5path=None, h5expname=None, h5hpname=None):
    fitpars=[]
    for tup in [c, t, D]:
        (segkey, filkey)=tup
        if 'fitpars' in fild[filkey].keys():
            fitpars+=list(fild[filkey]['fitpars'])
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
    
    paritionedtime=segd['~'.join(c)]
    t_=segd['~'.join(t)]
    D_=segd['~'.join(D)]
    
    fp=[]
    ff=fitfcns()
    for pt, cyct, cycD in zip(paritionedtime, t_, D_):
        ff.genfit(FitFcnLibrary[FIT_t5.__name__], fitpars, (pt[pt>=0], cyct[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[:50].mean()]*len(fitpars)
            tempt=cyct[-50:].mean()
            tempD=cycD[-50:].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempt**(i+1)
            ff.genfit(FitFcnLibrary[FIT_t5.__name__], fitpars, (pt[pt>=0], cyct[pt>=0], cycD[pt>=0]))
        if ff.error:
            fitpars=[cycD[-50:].mean()]*len(fitpars)
            tempt=cyct[:50].mean()
            tempD=cycD[:50].mean()/4.
            for i in range(4):
                fitpars[-4+i]=tempD/tempt**(i+1)
            ff.genfit(FitFcnLibrary[FIT_t5.__name__], fitpars, (pt[pt>=0], cyct[pt>=0], cycD[pt>=0]))
        fp+=[ff.finalparams]
    fp=numpy.float32(fp)
    return fp


#def timeintegrate(arr, integwindow_s=1.)

#def FIT??pieceC_T4_intdT(segd, fild, c, T, dT, D, h5path=None, h5expname=None, h5hpname=None):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
#    fitpars=[]
#    for tup in [c, T, dT, D]:
#        (segkey, filkey)=tup
#        if 'fitpars' in fild[filkey].keys():
#            fitpars+=list(fild[filkey]['fitpars'])
#        if not '~'.join(tup) in segd.keys():
#            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
#    
#    paritionedtime=segd['~'.join(c)]
#    T_=segd['~'.join(T)]
#    dT_=segd['~'.join(dT)]
#    D_=segd['~'.join(D)]
#    
#    fp=[]
#    ff=fitfcns()
#    for pt, cycT, cycdT, cycD in zip(paritionedtime, T_, dT_, D_):
#        ff.genfit(FitFcnLibrary[pieceC_T4_intdT.__name__], fitpars, (pt[pt>=0], cycT[pt>=0], cycdT[pt>=0], cycD[pt>=0]))
#        fp+=[ff.finalparams]
#    fp=numpy.float32(fp)
#    return fp

def pt_none(segd, fild, t, h5path=None, h5expname=None, h5hpname=None):
    (segkey, filkey)=t
    return numpy.zeros(segd[segkey].shape, dtype='float32')
    
def useexisting(segd, fild, x, h5path=None, h5expname=None, h5hpname=None):
    (segkey, filkey)=x
    return segd[segkey]

def pt_user(segd, fild, t, D, h5path=None, h5expname=None, h5hpname=None):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    fitpars=[]
    for tup in [t, D]:
        (segkey, filkey)=tup
        if 'fitpars' in fild[filkey].keys():
            fitpars+=fild[filkey]['fitpars']
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
    
    D_=segd['~'.join(D)]
    (segkey, filkey)=t
    segd[segkey]
    idialog=timepartDialog(None, segd[segkey], numpieces=fild[filkey]['numpartitions'], yvals=D_)
    idialog.exec_()
    return idialog.timepart
    
def pt_calc(segd, fild, t, D, h5path=None, h5expname=None, h5hpname=None):#TODO: need tgo wirte auto time partiiooning. this is not implemented yet, so just using copy user-defined fcn
    fitpars=[]
    for tup in [t, D]:
        (segkey, filkey)=tup
        if 'fitpars' in fild[filkey].keys():
            fitpars+=fild[filkey]['fitpars']
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
    
    D_=segd['~'.join(D)]
    (segkey, filkey)=t
    segd[segkey]
    idialog=timepartDialog(None, segd[segkey], numpieces=fild[filkey]['numpartitions'], yvals=D_)
    idialog.exec_()
    return idialog.timepart
    
def C_Dfit(segd, fild, D, pt, h5path=None, h5expname=None, h5hpname=None):
    for tup in [D]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
        #    segd['~'.join(tup)]/=dt
    D_=segd['~'.join(D)]
    
    f, fitd=getfitdict_nameseg(h5path, h5expname, h5hpname, 'samplepowerperrate', segd['segindex'])
    f.close()
    (segkey, filkey)=pt
    Dfit=numpy.float32(evaluatefitfcn(fitd, segd, interp_kind=fild[filkey]['interp_kind'], extrap_order=fild[filkey]['extrap_order'], interppars=fild[filkey]['interppars']))
    return D_-Dfit
    
    
def CTpks_secder(segd, fild, C, T, h5path=None, h5expname=None, h5hpname=None):
    for tup in [C, T]:
        (segkey, filkey)=tup
        if not '~'.join(tup) in segd.keys():
            segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
        #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
        #    segd['~'.join(tup)]/=dt
    C_=segd['~'.join(C)]
    T_=segd['~'.join(T)]
    
    (segkey, filkey)=T
    Tfild=fild[filkey]
    dx=Tfild['gridinterval']
    (segkey, filkey)=C
    Cfild=fild[filkey]
    
    psSGdict=dict([(k, v) for k, v in Cfild.iteritems() if k in peaksearch1dSG.func_code.co_varnames[:peaksearch1dSG.func_code.co_argcount]])
    psSGdict['dx']=dx
    psSGdict['critsepind']=Cfild['critsep']/dx
    
    pkslist_cycles=[]
    siglist_cycles=[]
    residlist=[]
    for X, Y in zip(T_, C_):
        Xgrid=numpy.linspace(X.min(), X.max(), (X.max()-X.min())/dx+1)
        Ygrid=numpy.empty(Xgrid.shape, dtype='float32')
        gridind=[numpy.argmin((x-Xgrid)**2) for x in X]
        indsgot=numpy.sort(numpy.uint32(list(set(gridind))))
        indsinterp=numpy.sort(numpy.uint32(list(set(range(len(Xgrid)))-set(gridind))))
        gridind=numpy.uint32(gridind)
        for i in indsgot:
            Ygrid[i]=Y[gridind==i].mean()
        Ygrid[indsinterp]=numpy.float32(scipy.interpolate.interp1d(indsgot, Ygrid[indsgot])(indsinterp))
        psSGdict['x']=Ygrid
        pkind=peaksearch1dSG(**psSGdict)

        pkht=Ygrid[numpy.uint32(numpy.round(pkind))]
        pkposn=Xgrid[numpy.uint32(numpy.round(pkind))]
        iarr=numpy.uint32(range(len(Xgrid)))
        hwposns1=[(numpy.any((s*Ygrid<(s*h/2.))&(iarr>i)) and (numpy.where((s*Ygrid<(s*h/2.))&(iarr>i))[0][0]-i,) or (iarr[-1]/100.,))[0] for i, h, s in zip(pkind, pkht, numpy.sign(pkht))]
        hwposns0=[(numpy.any((s*Ygrid<(s*h/2.))&(iarr<i)) and (i-numpy.where((s*Ygrid<(s*h/2.))&(iarr<i))[0][-1],) or (iarr[-1]/100.,))[0] for i, h, s in zip(pkind, pkht, numpy.sign(pkht))]
        pkhw=dx*(numpy.float32(hwposns1)+numpy.float32(hwposns0))/2.


        pks=numpy.float32([pkposn, pkhw, pkht]).T#, numpy.ones(pkht.shape, dtype='float32')*.5]).T
        print '&&', pks
        #***
        fitpks=pks
        sigs=numpy.zeros(fitpks.shape, dtype='float32')
        resids=numpy.zeros(fitpks.shape[0], dtype='float32')
        #***
#        fitpks, sigs, resids=fitpeakset(X, Y, pks, PeakFcnLibrary[Cfild['pkfcn']])
#        print '**', fitpks
        fitpks=numpy.float32(fitpks)
        sigs=numpy.float32(sigs)
        resids=numpy.float32(resids)
        if Cfild['useredit']:
            idialog=peakfiteditDialog(None, X, Y, fitpks, PeakFcnLibrary[Cfild['pkfcn']])
            idialog.exec_()
            if idialog.sigs is None:#if fitting was not done then only update if user changed peaks. if fitting done, update everything
                if not (len(fitpks)==(idialog.fitpks) and numpy.all(fitpks==idialog.fitpks)):
                    fitpks=idialog.fitpks
                    sigs=numpy.zeros(fitpks.shape, dtype='float32')
                    resids=numpy.zeros(fitpks.shape[0], dtype='float32')
            else:
                fitpks=idialog.fitpks
                sigs=idialog.sigs
                resids=idialog.resids
        
        pkslist_cycles+=[fitpks]
        siglist_cycles+=[sigs]
        residlist+=[resids]
    
    PROFILEANALYSIS={}
    PROFILEANALYSIS['Xgrid']=Xgrid
    PROFILEANALYSIS['Ygrid']=Ygrid
    
    n=numpy.max([len(l) for l in pkslist_cycles])
    if n:
        parr=numpy.ones((len(pkslist_cycles), n, 3), dtype='float32')*numpy.nan
        sarr=numpy.ones((len(pkslist_cycles), n, 3), dtype='float32')*numpy.nan
        for pl, sl, p, s in zip(pkslist_cycles, siglist_cycles, parr, sarr):
            p[:len(pl), :]=pl[:, :]
            s[:len(sl), :]=sl[:, :]
        PROFILEANALYSIS['peaks']=parr
        PROFILEANALYSIS['sigmas_peaks']=sarr
        PROFILEANALYSIS['residuals_peaks']=numpy.float32(residlist)

    return PROFILEANALYSIS

def WinFFTI(segd, fild, I, h5path, h5expname, h5hpname):
    return WinFFT(segd, fild, I, h5path, h5expname, h5hpname)
def WinFFTV(segd, fild, V, h5path, h5expname, h5hpname):
    return WinFFT(segd, fild, V, h5path, h5expname, h5hpname)
def WinFFTfV(segd, fild, fV, h5path, h5expname, h5hpname):
    return WinFFT(segd, fild, fV, h5path, h5expname, h5hpname)
    
def WinFFT(segd, fild, X, h5path, h5expname, h5hpname):
    ptspercyc=pts_sincycle_h5(h5path, h5expname, h5hpname)

    (segkey, filkey)=X
    n=ptspercyc*fild[filkey]['n1wcycles_window']
    #for i in fild[filkey]['harmonics']:
    #ans=numpy.empty((segd[segkey].shape[0], segd[segkey].shape[1], n//2+1, 2), dtype='float32')
    ans=numpy.empty((segd[segkey].shape[0], segd[segkey].shape[1], 10, 2), dtype='float32')
    
    freq_1w=(numpy.array(range(n//2+1))/n)*ptspercyc
    i2w=numpy.argmin((1.-freq_1w)**2)
    i2w=numpy.argmin((2.-freq_1w)**2)
    i3w=numpy.argmin((3.-freq_1w)**2)
    saveinds=[0, i1w-1, i1w, i1w+1, i2w-1, i2w, i2w+1, i3w-1, i3w, i3w+1]
    for i, arr in enumerate(segd[segkey]):
        fftx, ffty=windowfft_xy(arr, n)
        ans[i, :, :, 0]=fftx[saveinds]
        ans[i, :, :, 1]=ffty[saveinds]
    return ans #makes array that is cycles x nptsinsegx 10 x2 , where the 10 are 0w,1w-,1w,1w+,2w-,2w,2w+,3w-,3w,3w+, the separation between neibhboring frequencies is daqHz/totnumptsinwindow or daqHz/(ptspercyc*n1wcycles_window).

def LIAI(segd, fild, I, h5path, h5expname, h5hpname):
    return LIA(segd, fild, I, h5path, h5expname, h5hpname)
def LIAV(segd, fild, V, h5path, h5expname, h5hpname):
    return LIA(segd, fild, V, h5path, h5expname, h5hpname)
def LIAfV(segd, fild, fV, h5path, h5expname, h5hpname):
    return LIA(segd, fild, fV, h5path, h5expname, h5hpname)

def LIA(segd, fild, X, h5path, h5expname, h5hpname):
    ptspercyc=pts_sincycle_h5(h5path, h5expname, h5hpname)
    (segkey, filkey)=X
    #hlist=fild[filkey]['harmonics']
    hlist=[1, 2, 3]
    ans=numpy.empty((segd[segkey].shape[0], segd[segkey].shape[1], len(hlist), 2), dtype='float32')
    for i, arr in enumerate(segd[segkey]):
        for j, h in enumerate(hlist):
            ans[i, :, j, 0], ans[i, :, j, 1]=lia_xy(arr, ptspercyc/h, ncyclewin=fild[filkey]['n1wcycles_window']*h)
            #ans[i, :, j, 0], ans[i, :, j, 1]=lia_ampphase(arr, ptspercyc/h, ncyclewin=fild[filkey]['n1wcycles_window']*h, phaseshift=0.)
    return ans #makes array that is cycles x nptsinsegxharmonic x 2for x,y

def R_fftIV(segd, fild, fftI, fftV, h5path=None, h5expname=None, h5hpname=None):
    h=fild[fftI[1]]['harmonic']
    freqind=h*3-1*(h>0)
    I=fftI+(segd[fftI[0]][:, :, freqind, 0],)
    V=fftV+(segd[fftV[0]][:, :, freqind, 0],)
    return R_IV(segd, fild, I, V, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def T_fftIV(segd, fild, fftI, fftV, h5path, h5expname, h5hpname):
    h=fild[fftI[1]]['harmonic']
    freqind=h*3-1*(h>0)
    I=fftI+(segd[fftI[0]][:, :, freqind, 0],)
    V=fftV+(segd[fftV[0]][:, :, freqind, 0],)
    return T_IV(segd, fild, I, V, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def P_fftIV(segd, fild, fftI, fftV, h5path, h5expname, h5hpname):
    h=fild[fftI[1]]['harmonic']
    freqind=h*3-1*(h>0)
    I=fftI+(segd[fftI[0]][:, :, freqind, 0],)
    V=fftV+(segd[fftV[0]][:, :, freqind, 0],)
    return P_IV(segd, fild, I, V, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def R_liaIV(segd, fild, liaI, liaV, h5path=None, h5expname=None, h5hpname=None):
    h=fild[liaI[1]]['harmonic']
    freqind=h-1
    I=liaI+(segd[liaI[0]][:, :, freqind, 0],)
    V=liaV+(segd[liaV[0]][:, :, freqind, 0],)
    return R_IV(segd, fild, I, V, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def T_liaIV(segd, fild, liaI, liaV, h5path=None, h5expname=None, h5hpname=None):
    h=fild[liaI[1]]['harmonic']
    freqind=h-1
    I=liaI+(segd[liaI[0]][:, :, freqind, 0],)
    V=liaV+(segd[liaV[0]][:, :, freqind, 0],)
    return T_IV(segd, fild, I, V, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def P_liaIV(segd, fild, liaI, liaV, h5path=None, h5expname=None, h5hpname=None):
    h=fild[liaI[1]]['harmonic']
    freqind=h-1
    I=liaI+(segd[liaI[0]][:, :, freqind, 0],)
    V=liaV+(segd[liaV[0]][:, :, freqind, 0],)
    return P_IV(segd, fild, I, V, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def R_fftIfV(segd, fild, fftI, fV, h5path=None, h5expname=None, h5hpname=None):
    return R_fftIV(segd, fild, fftI, fV, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)
def T_fftIfV(segd, fild, fftI, fV, h5path=None, h5expname=None, h5hpname=None):
    return T_fftIV(segd, fild, fftI, fV, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)
def P_fftIfV(segd, fild, fftI, fV, h5path=None, h5expname=None, h5hpname=None):
    return P_fftIV(segd, fild, fftI, fV, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def R_liaIfV(segd, fild, liaI, fV, h5path=None, h5expname=None, h5hpname=None):
    return R_liaIV(segd, fild, liaI, fV, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)
def T_liaIfV(segd, fild, liaI, fV, h5path=None, h5expname=None, h5hpname=None):
    return T_liaIV(segd, fild, liaI, fV, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)
def P_liaIfV(segd, fild, liaI, fV, h5path=None, h5expname=None, h5hpname=None):
    return P_liaIV(segd, fild, liaI, fV, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)

def mCp_fftVIPdT(segd, fild, fftV, fftI, P, dT, h5path, h5expname, h5hpname):
    h=fild[fftV[1]]['harmonic']
    harmind=h*3-1*(h>0)
    return mCp_gen(segd, fild, fftV, fftI, P, dT, h5path, h5expname, h5hpname, h, harmind)

def mCp_ffVIPdT(segd, fild, ffV, fftI, P, dT, h5path, h5expname, h5hpname):
    h=fild[ffV[1]]['harmonic']
    harmind=h*3-1*(h>0)
    return mCp_gen(segd, fild, ffV, fftI, P, dT, h5path, h5expname, h5hpname, h, harmind, Vsrcisfiltered=True)
    
def mCp_liaVIPdT(segd, fild, liaV, liaI, P, dT, h5path, h5expname, h5hpname):
    h=fild[liaV[1]]['harmonic']
    harmind=h-1
    return mCp_gen(segd, fild, liaV, liaI, P, dT, h5path, h5expname, h5hpname, h, harmind)

def mCp_lfVIPdT(segd, fild, lfV, liaI, P, dT, h5path, h5expname, h5hpname):
    h=fild[lfV[1]]['harmonic']
    harmind=h-1
    return mCp_gen(segd, fild, lfV, liaI, P, dT, h5path, h5expname, h5hpname, h, harmind, Vsrcisfiltered=True)
    
def mCp_gen(segd, fild, acV, acI, P, dT, h5path, h5expname, h5hpname, h, harmind, harm01inds, Vsrcisfiltered=False):#the I, dIdt, etc. should be tuples with a key for segd and then a key for fild
    RoToAl=RoToAl_h5(h5path, h5expname, h5hpname)
    ppc=pts_sincycle_h5(h5path, h5expname, h5hpname)
    dt=dt_h5(h5path, h5expname, h5hpname)
    freq1w=1./(ppc*dt)
    if h<2:
        return D_PS(segd, fild, P, dT, h5path=h5path, h5expname=h5expname, h5hpname=h5hpname)
    elif h==2:
        mCpfcn=mCp_2w
    elif h==3:
        mCpfcn=mCp_3w
    else:
        print 'invalid harmonic value - abort'
        return None
    acV+=(segd[acV[0]][:, :, harmind, :],)#the last ":" is for using X and Y, any smoothing filter will be applied to both
    acI+=(segd[acI[0]][:, :, harm01inds+[harmind], :],)
    for tup in [acV, acI, P, dT]:
        (segkey, filkey)=tup[:2]
        if not '~'.join(tup[:2]) in segd.keys():
            if len(tup)==3:
                arr=tup[2]
            else:
                arr=segd[segkey]
            segd['~'.join(tup[:2])]=performgenericfilter(arr, fild[filkey])
    acV_=segd['~'.join(acV[:2])]
    acI_=segd['~'.join(acI[:2])]
    #P_=segd['~'.join(P[:2])]
    dT_=segd['~'.join(dT[:2])]
    inds=numpy.where(acV_[:, :, 0]<=0.)# these 3 lines will effectively remove neg and inf Res and replace them with the res value of the nearest acceptable value. These modification will no be reflected in the segd values for the source data
    if len(inds[0])>0:
        acV_[:, :, 0]=replacevalswithneighsin2nddim(acV_[:, :, 0], inds)
        acV_[:, :, 1]=replacevalswithneighsin2nddim(acV_[:, :, 1], inds)
        acI_[:, :, 2, 0]=replacevalswithneighsin2nddim(acI_[:, :, 2, 0], inds)
        acI_[:, :, 2, 1]=replacevalswithneighsin2nddim(acI_[:, :, 2, 1], inds)
        #P_=replacevalswithneighsin2nddim(P_, inds)
        #dT_=replacevalswithneighsin2nddim(dT_, inds)
    mCpargdict={}
    mCpargdict['VX']=acV_[:, :, 0]
    mCpargdict['VY']=acV_[:, :, 1]
    mCpargdict['I0X']=acI_[:, :, 0, 0]
    mCpargdict['I0Y']=acI_[:, :, 0, 1]
    mCpargdict['I1X']=acI_[:, :, 1, 0]
    mCpargdict['I1Y']=acI_[:, :, 1, 1]
    mCpargdict['IX']=acI_[:, :, 2, 0]
    mCpargdict['IY']=acI_[:, :, 2, 1]
    #mCpargdict['P']=P_
    mCpargdict['dT']=dT_
    mCpargdict['Ro']=RoToAl[0]
    mCpargdict['tcr']=RoToAl[2]
    mCpargdict['freq1w']=freq1w
    mCpargdict['Vsrcisfiltered']=Vsrcisfiltered
    return mCpfcn([(k, v) for k, v in mCpargdict.itervalues() if k in f.func_code.co_varnames[:f.func_code.co_argcount]])


#p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/Nanocopeia1_PnSC.h5'
##p='C:/Users/JohnnyG/Documents/HarvardWork/pharma/PnSC/Nanocopeia1_PnSC.h5'
#e='NoSampleRamps'
#h='20110203Nanocop_ramp_160ms_10mA_cell10'
#
#mainapp=QApplication(sys.argv)
#form=SCrecipeDialog(None, p, e, h)
##form=analysisviewerDialog(None,CreateHeatProgSegDictList(p,e,h))
#form.show()
#form.setFocus()
#mainapp.exec_()


#p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110714_SnACcal.h5'
#e='Sn_10kHz_4e4Ks'
#h='cell29_57.75dc56.1ac_10kHz_12.6ms_1_of_1'
#
#mainapp=QApplication(sys.argv)
#form=SCrecipeDialog(None, p, e, h, calctype='AC')
##form=analysisviewerDialog(None,CreateHeatProgSegDictList(p,e,h))
#form.show()
#form.setFocus()
#mainapp.exec_()
