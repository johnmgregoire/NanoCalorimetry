import time
import os
import sys
import numpy
import h5py
from PnSC_ui import *
from PnSC_dataimport import *
from PnSC_SCui import *
from PnSC_math import *
from PnSC_h5io import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *

#class MainMenu(QMainWindow):
#    def __init__(self, TreeWidg):
#        super(MainMenu, self).__init__(None)
#
#        self.setObjectName("MainMenu")
#        self.bodywidget = QWidget(self)
#        self.bodywidget.setObjectName("bodywidget")
#        self.tasklistLabel = QLabel(self.bodywidget)
#        self.tasklistLabel.setGeometry(QRect(9, 10, 1006, 16))
#        self.tasklistLabel.setObjectName("tasklistLabel")
#        self.setCentralWidget(self.bodywidget)
#        self.main_menu_pulldown = QMenuBar(self)
#        self.main_menu_pulldown.setGeometry(QRect(0, 0, 1025, 27))
#        self.main_menu_pulldown.setObjectName("main_menu_pulldown")
#        self.menuExit = QMenu(self.main_menu_pulldown)
#        self.menuExit.setObjectName("menuExit")
#        self.menuExit.setTitle('EXIT')
#        self.setMenuBar(self.main_menu_pulldown)
#        self.statusbar = QStatusBar(self)
#        self.statusbar.setEnabled(False)
#        self.statusbar.setObjectName("statusbar")
#        self.setStatusBar(self.statusbar)
#        self.actionExit = QAction(self)
#        self.actionExit.setObjectName("actionExit")
#        self.actionExit.setText('exit')
#        self.menuExit.addAction(self.actionExit)
#        self.main_menu_pulldown.addAction(self.menuExit.menuAction())
#
#        QMetaObject.connectSlotsByName(self)
#
#    @pyqtSignature("")
#    def on_actionExit_triggered(self):
#        print 'init h5'
        
class MainMenu(QMainWindow):
    def __init__(self, previousmm):#, TreeWidg):
        super(MainMenu, self).__init__(None)
        #self.setupUi(self)
        self.setWindowTitle('Vlassak Group PnSC Analysis')
        
        #self.treeWidget=TreeWidg
        
        self.h5path="%s" % os.getcwd()
        
        self.bodywidget = QWidget(self)
        self.bodywidget.setObjectName("bodywidget")
        self.treeWidget=QTreeWidget(self.bodywidget)
        QObject.connect(self.treeWidget,SIGNAL("itemSelectionChanged()"),self.processtreeselection)
        self.setupmenu()
        self.setCentralWidget(self.bodywidget)
        
        self.statusdict={'h5open':False}
        self.actionenable()
        self.resize(820, 620)
        self.treeWidget.setGeometry(QRect(10, 10, 800, 520))
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        #sizePolicy.setHeightForWidth(self.treeWidget.sizePolicy().hasHeightForWidth())
        self.treeWidget.setSizePolicy(sizePolicy)

        self.redrawPushButton = QPushButton(self.bodywidget)
        self.redrawPushButton.setGeometry(QRect(10, 550, 200, 25))
        self.redrawPushButton.setText('Draw h5 Tree')
        QObject.connect(self.redrawPushButton, SIGNAL("pressed()"), self.redraw)
        self.expandPushButton = QPushButton(self.bodywidget)
        self.expandPushButton.setGeometry(QRect(210, 550, 200, 25))
        self.expandPushButton.setText('Expand h5 Tree')
        QObject.connect(self.expandPushButton, SIGNAL("pressed()"), self.expandtree)
        self.expandexceptPushButton = QPushButton(self.bodywidget)
        self.expandexceptPushButton.setGeometry(QRect(410, 550, 200, 25))
        self.expandexceptPushButton.setText('Expand Groups')
        QObject.connect(self.expandexceptPushButton, SIGNAL("pressed()"), self.expandgrouptree)
        
        if previousmm is None:
            self.on_action_openh5_triggered()
        else:
            oldselection=mm.geth5selectionpath(liststyle=True, removeformatting=False)
            self.h5path=previousmm.h5path
            h5file=h5py.File(self.h5path, mode='r')
            fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
            h5file.close()
            self.statusdict['h5open']=True
            self.actionenable()
        
    def setupmenu(self):
        self.setObjectName("MainMenu")
        self.main_menu_pulldown = QMenuBar(self)
        self.main_menu_pulldown.setObjectName("main_menu_pulldown")
        self.ActionDict={}
        
        #setup a menu section
        self.menufileio = QMenu(self.main_menu_pulldown)
        self.menufileio.setObjectName("menufileio")
        self.menufileio.setTitle('File IO')
        self.main_menu_pulldown.addAction(self.menufileio.menuAction())
        #end of menu head
        
        #setup a menu item in a menu section.    self.<NAME>=....,(self, <NAME>, <text>, <self.menufileio>, <list of tuples, tuple is name of requirement and list of acceptable values>, self.ActionDict), keep this last item the same
        self.action_openh5=MainMenuQAction(self,'action_openh5', 'open h5 file', self.menufileio, [], self.ActionDict)
        self.action_importscdata=MainMenuQAction(self,'action_importscdata', 'import calorimetry data', self.menufileio, [('h5open', [True])], self.ActionDict)
        self.action_batchimportscdata=MainMenuQAction(self,'action_batchimportscdata', 'batch import calorimetry data setup', self.menufileio, [('h5open', [True])], self.ActionDict)
        self.action_batchimportdatafixedmsma=MainMenuQAction(self,'action_batchimportdatafixedmsma', 'batch import calorimetry data using segment info from selected HeatProgram', self.menufileio, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_createh5=MainMenuQAction(self,'action_createh5', 'new h5 file', self.menufileio, [], self.ActionDict)
        self.action_createexpgrp=MainMenuQAction(self,'action_createexpgrp', 'new experiment group', self.menufileio, [('h5open', [True])], self.ActionDict)
        self.action_delh5grp=MainMenuQAction(self,'action_delh5grp', 'DELETE selected group', self.menufileio, [('h5open', [True]), ('selectiontype', ['Group'])], self.ActionDict)
        #self.action_delexpgrp=MainMenuQAction(self,'action_delexpgrp', 'DELETE experiment group', self.menufileio, [('h5open', [True])], self.ActionDict)
        self.action_editattrs=MainMenuQAction(self,'action_editattrs', 'Edit import attrs (select a heat program)', self.menufileio, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        
        #setup a menu section
        self.menuplot = QMenu(self.main_menu_pulldown)
        self.menuplot.setObjectName("menuplot")
        self.menuplot.setTitle('Visualization')
        self.main_menu_pulldown.addAction(self.menuplot.menuAction())
        #end of menu head
        
        #setup a menu item in a menu section.   
        self.action_plotraw=MainMenuQAction(self,'action_plotraw', 'plot Dataset values (select dataset)', self.menuplot, [('h5open', [True]), ('selectiontype', ['Dataset'])], self.ActionDict)
        self.action_printdata=MainMenuQAction(self,'action_printdata', 'print Dataset values (select dataset or attribute)', self.menuplot, [('h5open', [True]), ('selectiontype', ['Dataset', 'Attr'])], self.ActionDict)
        self.action_plotmetadata=MainMenuQAction(self,'action_plotmetadata', 'Plot Heat Program MetaData(select heat program)', self.menuplot, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_getsegd=MainMenuQAction(self,'action_getsegd', 'send SegDict to data (select a heat program)', self.menuplot, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_plotsegs=MainMenuQAction(self,'action_plotsegs', 'plot Segs by color (select a heat program)', self.menuplot, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_viewSCanalysis=MainMenuQAction(self,'action_viewSCanalysis', 'SC data viewer (select a heat program)', self.menuplot, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_viewFit=MainMenuQAction(self,'action_viewFit', 'Fit data viewer (select a heat program)', self.menuplot, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        
        
        #setup a menu section
        self.calprep = QMenu(self.main_menu_pulldown)
        self.calprep.setObjectName("calprep")
        self.calprep.setTitle('Calibration Prep')
        self.main_menu_pulldown.addAction(self.calprep.menuAction())
        #end of menu head
        
        #setup a menu item in a menu section.   
        self.action_calcresistance=MainMenuQAction(self,'action_calcresistance', 'Calc cell Res (select heat program or experiment)', self.calprep, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram', 'experiment'])], self.ActionDict)
        self.action_setuprescal=MainMenuQAction(self,'action_setuprescal', 'Setup R(T) cal', self.calprep, [('h5open', [True])], self.ActionDict)
        self.action_assignrescal=MainMenuQAction(self,'action_assignrescal', 'Assign R(T) cal (select experiment)', self.calprep, [('h5open', [True]), ('selectiongrouptype', ['experiment'])], self.ActionDict)
        self.action_calcresextraptoTo=MainMenuQAction(self,'action_calcresextraptoTo', 'Calc Res that gives To (select heat program or experiment)', self.calprep, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram', 'experiment'])], self.ActionDict)
        self.action_calcresbycycle=MainMenuQAction(self,'action_calcresbycycle', 'Calc Res for each cycle using first soak segment (select heat program or experiment)', self.calprep, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram', 'experiment'])], self.ActionDict)
        self.action_entertwopointres=MainMenuQAction(self,'action_entertwopointres', 'Enter list of 2 points R values', self.calprep, [('h5open', [True])], self.ActionDict)

        #end of actions
        
        #setup a menu section
        self.anmenu = QMenu(self.main_menu_pulldown)
        self.anmenu.setObjectName("anmenu")
        self.anmenu.setTitle('Calorimetry Analysis')
        self.main_menu_pulldown.addAction(self.anmenu.menuAction())
        #end of menu head
        
        #setup a menu item in a menu section. 
        self.action_delan=MainMenuQAction(self,'action_delan', 'Delete analysis Group (select analysis group)', self.anmenu, [('h5open', [True]),  ('selectiongrouptype', ['analysis'])], self.ActionDict)
        self.action_screcipe=MainMenuQAction(self,'action_screcipe', 'Build SC analysis recipe (select heat program)', self.anmenu, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_fitlossrecipe=MainMenuQAction(self,'action_fitlossrecipe', 'Build heat loss fit model recipe (select heat program)', self.anmenu, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram']), ('samplepowerperrateexists', [True])], self.ActionDict)
        self.action_heatcaprecipe=MainMenuQAction(self,'action_heatcaprecipe', 'Build heat capacity recipe (select heat program)', self.anmenu, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram']), ('samplepowerperrateexists', [True])], self.ActionDict)
        self.action_heatcappeaksrecipe=MainMenuQAction(self,'action_heatcappeaksrecipe', 'Build C(T) peak search+fit recipe (select heat program)', self.anmenu, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram']), ('samplepowerperrateexists', [True])], self.ActionDict)
        self.action_acrecipe=MainMenuQAction(self,'action_acrecipe', 'Build AC freq analysis recipe (select heat program)', self.anmenu, [('h5open', [True]),  ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_applyscrecipe=MainMenuQAction(self,'action_applyscrecipe', 'Apply SC analysis recipe (select experiment or heat program)', self.anmenu, [('h5open', [True]),  ('selectiongrouptype', ['experiment', 'heatprogram'])], self.ActionDict)
        
        self.setMenuBar(self.main_menu_pulldown)
        QMetaObject.connectSlotsByName(self)
    
    def redraw(self):
        if os.path.exists(self.h5path) and self.h5path.endswith('.h5'):
            h5file=h5py.File(self.h5path, mode='r')
            fillh5tree(self.treeWidget, h5file)
            h5file.close()
            self.statusdict['h5open']=True
            self.actionenable()
            

    def expandgrouptree(self):
        self.expandtree(groupsonly=True)
    def expandtree(self, groupsonly=False):
        def expandchildren(item):#recursive
            for i in range(item.childCount()):
                child=item.child(i)
                if not groupsonly or True in [not ('(' in child.child(j).text(0) or str(child.child(j).text(0)).startswith("'")) for j in range(child.childCount())]:
                    child.setExpanded(True)
                    expandchildren(child)
        for i in range(self.treeWidget.topLevelItemCount()):
            item=self.treeWidget.topLevelItem(i) 
            item.setExpanded(True)
            expandchildren(item)
    
    def settreeselection_list(self, selectionpathlist):
        item=self.treeWidget.topLevelItem(0)
        for itemname in selectionpathlist:
            chn=[item.child(i).text(0) for i in range(item.childCount())]
            if itemname in chn:
                item.setExpanded(True)
                item=item.child(chn.index(itemname))
            else:
                break
        self.treeWidget.setCurrentItem(item)
    def actionenable(self):
        for aname, ad in self.ActionDict.iteritems():
            ad['ref'].setDisabled(False in [(k in self.statusdict.keys()) and (self.statusdict[k] in vals) for k, vals in ad['enable_reqs']])
#            if aname=='action_calcresistance':
#                print [(k in self.statusdict.keys()) and (self.statusdict[k] in vals) for k, vals in ad['enable_reqs']]
#                print [(k, vals) for k, vals in ad['enable_reqs']]

    def h5nodename_treeitem(self, treeitem, removeformatting=True):
        if removeformatting:
            return ((str(treeitem.text(0)).partition(':')[0]).partition('(')[0]).strip("'")
        else:
            return str(treeitem.text(0))

    def geth5selectionpath(self, liststyle=False, removeformatting=True):
        try:
            treeitem=self.currenttreeitem
        except:
            return '/'
        attrname=None
        if self.statusdict['selectiontype']=='Attr':
            attrname=self.h5nodename_treeitem(treeitem, removeformatting=removeformatting)
            treeitem=treeitem.parent()
        s=[]
        while not treeitem.parent() is None:
            s=[self.h5nodename_treeitem(treeitem, removeformatting=removeformatting)]+s
            treeitem=treeitem.parent()
        if not liststyle:
            s='/'.join((s))
        if not attrname is None:
            return s, attrname
        return s
        

    def processtreeselection(self):
        treeitem=self.treeWidget.currentItem()
        self.currenttreeitem=treeitem
        print 'selection changed to ', treeitem.text(0)
        if treeitem.parent() is None:
            self.statusdict['selectiontype']='File'
        elif str(treeitem.text(0)).startswith("'"):
            self.statusdict['selectiontype']='Attr'
        elif '(' in treeitem.text(0):
            self.statusdict['selectiontype']='Dataset'
        else:
            self.statusdict['selectiontype']='Group'
        
        self.statusdict['selectionname']=self.h5nodename_treeitem(treeitem)
        if self.statusdict['selectiontype']=='File':
            self.statusdict['selectionparentname']=''
        else:
            self.statusdict['selectionparentname']=self.h5nodename_treeitem(treeitem.parent())
        
        self.statusdict['samplepowerperrateexists']=True#TODO: write code for checking on existence of analysis data arrays
        
        if self.statusdict['selectiontype']=='Group':
            if self.statusdict['selectionparentname']=='HeatProgram':
                self.statusdict['selectiongrouptype']='heatprogram'
            elif self.statusdict['selectionparentname']=='Calorimetry':
                self.statusdict['selectiongrouptype']='experiment'
            elif self.statusdict['selectionparentname']=='analysis':
                self.statusdict['selectiongrouptype']='analysis'
            else:
                self.statusdict['selectiongrouptype']='other'
        else:
            self.statusdict['selectiongrouptype']=''
        
        print self.statusdict
        self.actionenable()

    @pyqtSignature("")
    def on_action_createh5_triggered(self):
        temp=mygetsavefile(parent=self, xpath=self.h5path,markstr='Enter name of new h5 file', filename='.h5')
        if temp=='':
            return
        self.h5path=str(temp)
        createh5file(self.h5path)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
        self.statusdict['h5open']=True
        self.actionenable()
        
    @pyqtSignature("")
    def on_action_openh5_triggered(self):
        self.statusdict['h5open']=False
        self.actionenable()
        temp=mygetopenfile(parent=self, xpath=self.h5path, markstr='h5 file with calorimetry data', filename='.h5' )
        if temp=='':
            return
        self.h5path=str(temp)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
        self.statusdict['h5open']=True
        self.actionenable()

    @pyqtSignature("")
    def on_action_createexpgrp_triggered(self):
        idialog=lineeditDialog(self, title='Enter name for the new h5 experiment', deftext='')
        if not idialog.exec_():
            return
        h5expname=idialog.text
        
        h5file=h5py.File(self.h5path, mode='r')
        if h5expname in h5file['Calorimetry']:
            h5file.close()
            QMessageBox.warning(self,"FAILED",  "Experiment Group Exists - must first delete")
            return
        h5file.close()
        create_exp_grp(self.h5path, h5expname)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
       
       
    @pyqtSignature("")
    def on_action_editattrs_triggered(self):
        path=self.geth5selectionpath(liststyle=False)
        editattrs(self, self.h5path, path)

    @pyqtSignature("")
    def on_action_delh5grp_triggered(self):
        h5file=h5py.File(self.h5path, mode='r+')
        del h5file[self.geth5selectionpath(liststyle=False)]
        h5file.close()
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_delexpgrp_triggered(self):
        h5file=h5py.File(self.h5path, mode='r+')
        idialog=selectgroupDialog(self, h5file['Calorimetry'], title='Select h5 experiment group to DELETE')
        if not idialog:
            h5file.close()
            return
        del h5file[idialog.grp.name]
        h5file.close()
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()

    def batchrun_files(self, folder, startsendswith=('', ''), skiperrors=False):
        for fn in os.listdir(folder):
            if not (fn.startswith(startsendswith[0]) and fn.endswith(startsendswith[1])):
                continue
            p=os.path.join(folder, fn)
            self.batchattrdict['path']=p
            print 'running ', p
            if skiperrors:
                try:
                    self.batchfcn(batchattrdict=self.batchattrdict)
                except:
                    print 'ERROR IMPORTING ', p
            else:
                self.batchfcn(batchattrdict=self.batchattrdict)
    
    @pyqtSignature("")
    def on_action_batchimportdatafixedmsma_triggered(self):
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        pathlist=self.geth5selectionpath(liststyle=True)
        h5file, h5hpgrp=gethpgroup(self.h5path, pathlist[1], pathlist[4])
        segms=h5hpgrp.attrs['segment_ms'][:]
        segmA=h5hpgrp.attrs['segment_mA'][:]
        sgd=selectgroupDialog(self, h5file['Calorimetry'], title='Select h5 experiment group for import')
        if not sgd:
            h5file.close()
            return
        h5file.close()
        h5expname=sgd.grpname
        idialog=selectorDialog(self, FileFormatFunctionLibrary.keys(), title='Select data import protocol')
        if not idialog.exec_():
            return
        protname=idialog.name
        batchattrdict=getemptybatchattrdict()
        batchattrdict['grpname']=h5expname
        batchattrdict['protname']=protname
        plist=mygetopenfiles(parent=self, xpath=os.path.split(self.h5path)[0], markstr='Slsect ONE data file for each experiment to be imported', filename='.h5' )
        for p in plist:
            batchattrdict['path']=p
        
            ans=FileImport(self, protname, batchattrdict=batchattrdict)
            #print ans
            if not ans:
                continue
            AttrDict, DataSetDict, SegmentData=ans
            grpname=os.path.splitext(os.path.split(AttrDict['importpath'])[1])[0]
            writenewh5heatprogram(self.h5path, h5expname, grpname, AttrDict, DataSetDict, (segms, segmA))        
            
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_batchimportscdata_triggered(self, batchattrdict=None):
        self.batchfcn=self.on_action_importscdata_triggered
        self.batchattrdict=getemptybatchattrdict()
        #Not finished??
    @pyqtSignature("")
    def on_action_importscdata_triggered(self, batchattrdict=None):
        if batchattrdict is None:
            h5file=h5py.File(self.h5path, mode='r')
            sgd=selectgroupDialog(self, h5file['Calorimetry'], title='Select h5 experiment group for import')
            if not sgd:
                h5file.close()
                return
            h5file.close()
            h5expname=sgd.grpname
            idialog=selectorDialog(self, FileFormatFunctionLibrary.keys(), title='Select data import protocol')
            if not idialog.exec_():
                return
            protname=idialog.name
        else:
            h5expname=batchattrdict['grpname']
            protname=batchattrdict['protname']
        ans=FileImport(self, protname, batchattrdict=batchattrdict)
        #print ans
        if not ans:
            return
        AttrDict, DataSetDict, SegmentData=ans
        mA=DataSetDict['samplecurrent'][1][0]*DataSetDict['samplecurrent'][0]['Aunit']*1000.
        ms=1000.*numpy.float32(range(len(mA)))/AttrDict['daqHz']
        idialog=SegmentEditor(self, SegmentData, cycledata=(ms, mA))
        if batchattrdict is None:
            if not idialog.exec_():
                return
        else:
            for sb, k in [(idialog.firstderptsSpinBox, 'firstderptsSpinBox'), (idialog.secderptsSpinBox, 'secderptsSpinBox'), (idialog.secdervalSpinBox, 'secdervalSpinBox')]:
                if k in batchattrdict and not (batchattrdict[k] is None):
                    sb.setValue(batchattrdict[k])
            idialog.findsegs()
            idialog.ExitRoutine()
        SegmentData=idialog.SegmentData
        
        idialog=lineeditDialog(self, title='Enter name for h5 group', deftext=os.path.splitext(os.path.split(AttrDict['importpath'])[1])[0])
        if batchattrdict is None:
            if not idialog.exec_():
                return
            grpname=idialog.text
        else:
            grpname=os.path.splitext(os.path.split(AttrDict['importpath'])[1])[0]
            if 'savegrpname' in batchattrdict and not batchattrdict['savegrpname'] is None:
                grpname=batchattrdict['savegrpname']
        
        writenewh5heatprogram(self.h5path, h5expname, grpname, AttrDict, DataSetDict, SegmentData)
        
        oldselection=['Calorimetry', h5expname, 'measurement', 'HeatProgram', grpname]
        
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_calcresistance_triggered(self, critms_step=1., critmAperms_constmA=0.01, critdelmA_constmA=10., critmA_zero=0.1):
        pathlist=self.geth5selectionpath(liststyle=True)
        if self.statusdict['selectiongrouptype']=='experiment':
            h5file, hplist=experimenthppaths(self.h5path, pathlist[1])
            h5file.close()
            hplist=[hpp.rpartition('/')[2] for hpp in hplist]
        else:
            hplist=[pathlist[4]]
        for hp in hplist:
            #print hp
            dlist=CreateHeatProgSegDictList(self.h5path, pathlist[1], hp, critms_step=critms_step, critmAperms_constmA=critmAperms_constmA, critdelmA_constmA=critdelmA_constmA, critmA_zero=critmA_zero) 
            segtypelist=[d['segmenttype'] for d in dlist]
            if segtypelist.count('soak')==1:
                dsoak=dlist[segtypelist.index('soak')]
            elif segtypelist.count('soak')>1:
                print 'ERROR - MORE THAN ONE SOAK SEGMENT WAS FOUND - THIS IS UNEXPECTED FOR AN Ro HEAT PROGRAM'
                return
            else:
                print 'ERROR - NO SOAK SEGMENTS WERE FOUND - ONE IS REQUIRED FOR AN Ro HEAT PROGRAM'
                return
            if segtypelist.count('zero')>0:
                if segtypelist[segtypelist.index('soak')-1]=='zero':
                    dzero=dlist[segtypelist.index('soak')-1]#take the preceding zero if possible so that user can segment the initial transients separately and avoid issue
                else:
                    dzero=dlist[segtypelist.index('zero')]
            else:
                dzero=None
#            vals=[]
#            vals+=[CalcR0_segdict(dsoak, AveBeforeDivision=True, dzero=dzero)]
#            vals+=[CalcR0_segdict(dsoak, AveBeforeDivision=False, dzero=dzero)]
#            vals+=[(vals[0]+vals[1])/2.]
#            desc=['ratio of the means', 'mean of the ratios', 'Ave of these 2 values']
#            choices=['%.4f : %s' %(v, d) for v, d in zip(vals, desc)]
#            idialog=selectorDialog(self, choices, title='select value of R0 to use')
#            if not idialog.exec_():
#                return
#            R0=vals[idialog.index]
            R0=CalcR0_segdict(dsoak, AveBeforeDivision=True, dzero=dzero)
            writecellres(self.h5path, pathlist[1], hp, R0)
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
    

    @pyqtSignature("")
    def on_action_calcresextraptoTo_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        if self.statusdict['selectiongrouptype']=='experiment':
            h5file, hplist=experimenthppaths(self.h5path, pathlist[1])
            h5file.close()
            hplist=[hpp.rpartition('/')[2] for hpp in hplist]
        else:
            hplist=[pathlist[4]]
        pardict={}
        pardict['h5path']=self.h5path
        pardict['h5expname']=pathlist[1]
        
        for i, hp in enumerate(hplist):
            pardict['h5hpname']=hp
            if i==0:
                idialog=rescal_ExtraptoToDialog(self, pardict)
                idialog.exec_()
                Ro=idialog.Ro
                self.data=idialog.calcd
            else:
                Ro, self.data=calcRo_extraptoTo(**pardict)
            writecellres_calc(self.h5path, pathlist[1], hp, Ro)
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()

    @pyqtSignature("")
    def on_action_calcresbycycle_triggered(self, critms_step=1., critmAperms_constmA=0.01, critdelmA_constmA=10., critmA_zero=0.1):
        pathlist=self.geth5selectionpath(liststyle=True)
        if self.statusdict['selectiongrouptype']=='experiment':
            h5file, hplist=experimenthppaths(self.h5path, pathlist[1])
            h5file.close()
            hplist=[hpp.rpartition('/')[2] for hpp in hplist]
        else:
            hplist=[pathlist[4]]
        for hp in hplist:
            #print hp
            dlist=CreateHeatProgSegDictList(self.h5path, pathlist[1], hp, critms_step=critms_step, critmAperms_constmA=critmAperms_constmA, critdelmA_constmA=critdelmA_constmA, critmA_zero=critmA_zero) 
            segtypelist=[d['segmenttype'] for d in dlist]
            if 'soak' in segtypelist:
                dsoak=dlist[segtypelist.index('soak')]
            else:
                print 'ERROR - NO SOAK SEGMENTS WERE FOUND - ONE IS REQUIRED'
                return
            if segtypelist.count('zero')>0:
                if segtypelist[segtypelist.index('soak')-1]=='zero':
                    dzero=dlist[segtypelist.index('soak')-1]#take the preceding zero if possible so that user can segment the initial transients separately and avoid issue
                else:
                    dzero=dlist[segtypelist.index('zero')]
            else:
                dzero=None
            R0=numpy.array([CalcR0_segdict(extractcycle_SegDict(dsoak, i), AveBeforeDivision=True, dzero=extractcycle_SegDict(dzero, i)) for i in range(dsoak['cycletime'].shape[0])])
            writecellres_calc(self.h5path, pathlist[1], hp, R0)
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_setuprescal_triggered(self):
        idialog=rescalDialog(self, self.h5path)
        idialog.exec_()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()

    @pyqtSignature("")
    def on_action_assignrescal_triggered(self):
        rescalpath_getorassign(self.h5path, self.statusdict['selectionname'], parent=self, forceassign=True)
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        self.actionenable()

    @pyqtSignature("")
    def on_action_entertwopointres_triggered(self):
        idialog=TwoPointResTableDialog(self, self.h5path)
        idialog.exec_()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        

    @pyqtSignature("")
    def on_action_plotraw_triggered(self):
        h5file=h5py.File(self.h5path, mode='r')
        path=self.geth5selectionpath()
        self.data=readh5pyarray(h5file[path])
        h5file.close()
        if self.data.ndim==1:
            idialog=simpleplotDialog(self, self.data)
        else:
            plotdata=self.data.swapaxes(self.data.ndim-1,numpy.argmax(self.data.shape))#assume plot should be vs the longest dimmension
            plotdata=[plotdata[ind] for ind in numpy.ndindex(*plotdata.shape[:-1])]
            idialog=simpleplotDialog(self, plotdata)
        idialog.exec_()
    
    @pyqtSignature("")
    def on_action_plotmetadata_triggered(self):
        fcndict=heatprogrammetadatafcns
        idialog=selectorDialog(self, fcndict.keys(), title='select type of metadata to plot')
        if not idialog.exec_():
            return
        fcn=fcndict[idialog.name]
        
        pathlist=self.geth5selectionpath(liststyle=True)
        self.data=fcn(self.h5path, pathlist[1], pathlist[4])
        idialog=simpleplotDialog(self, self.data[1], xdata=self.data[0])
        idialog.exec_()

    @pyqtSignature("")
    def on_action_printdata_triggered(self):    
        h5file=h5py.File(self.h5path, mode='r')
        if self.statusdict['selectiontype']=='Attr':
            path, attrname=self.geth5selectionpath()
            self.data=h5file[path].attrs[attrname]
            print attrname, ': ', self.data
        elif self.statusdict['selectiontype']=='Dataset':
            path=self.geth5selectionpath()
            self.data=readh5pyarray(h5file[path])
            print path.rpartition('/')[2], ': ', self.data
        h5file.close()
    
    @pyqtSignature("")
    def on_action_getsegd_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        self.data=CreateHeatProgSegDictList(self.h5path, pathlist[1], pathlist[4])

    @pyqtSignature("")
    def on_action_plotsegs_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        self.data=CreateHeatProgSegDictList(self.h5path, pathlist[1], pathlist[4])
        idialog=SegmentCyclePlot(self, self.data)
        idialog.show()
        return idialog
        
    @pyqtSignature("")
    def on_action_viewSCanalysis_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        self.data=CreateHeatProgSegDictList(self.h5path, pathlist[1], pathlist[4])
        idialog=analysisviewerDialog(self, self.data, pathlist[1])
        idialog.show()
    
    @pyqtSignature("")
    def on_action_viewFit_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        self.data=getfitdictlist_hp(self.h5path, pathlist[1], pathlist[4])
        hpsdl=CreateHeatProgSegDictList(self.h5path, pathlist[1], pathlist[4])
        h5file, filterdict=getfilterdict(self.h5path, pathlist[1])
        h5file.close()
        fitviewer(self, hpsdl, self.data, filterdict)

    @pyqtSignature("")
    def on_action_delan_triggered(self):
        path=self.geth5selectionpath(liststyle=False)
        g, garb, p=path.strip('/').rpartition('/')
        h5file=h5py.File(self.h5path, mode='r+')
        h5g=h5file[g]
        del h5g[p]
        h5file.close()
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
        self.actionenable()
        
    @pyqtSignature("")
    def on_action_screcipe_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        idialog=SCrecipeDialog(self, self.h5path, pathlist[1], pathlist[4])
        idialog.show()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()

    @pyqtSignature("")
    def on_action_fitlossrecipe_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        idialog=SCrecipeDialog(self, self.h5path, pathlist[1], pathlist[4], calctype='FitPS')
        idialog.show()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_heatcaprecipe_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        idialog=SCrecipeDialog(self, self.h5path, pathlist[1], pathlist[4], calctype='QUC')
        idialog.show()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
    
    @pyqtSignature("")
    def on_action_heatcappeaksrecipe_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        idialog=SCrecipeDialog(self, self.h5path, pathlist[1], pathlist[4], calctype='CTpk')
        idialog.show()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_acrecipe_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        idialog=SCrecipeDialog(self, self.h5path, pathlist[1], pathlist[4], calctype='AC')
        idialog.show()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_applyscrecipe_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        if self.statusdict['selectiongrouptype']=='heatprogram':
            h5hpdflt=pathlist[4]
        else:
            h5hpdflt=None
        idialog=SCanalysisDialog(self, self.h5path, pathlist[1], h5hpdflt=h5hpdflt)
        idialog.show()
        oldselection=self.geth5selectionpath(liststyle=True, removeformatting=False)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection)
        h5file.close()

class MainMenuQAction(QAction):
    def __init__(self, parent, actionname, actiontext, hostmenu, reqs, adict):
        super(MainMenuQAction, self).__init__(parent)
        self.setObjectName(actionname)
        self.setText(actiontext)
        hostmenu.addAction(self)
        adict[actionname]={'ref':self, 'enable_reqs':reqs}

        
#class TreeWindow(QDialog):#Creates a side window for displaying the contents of an .h5 file in tree format
#    def __init__(self):
#        super(TreeWindow, self).__init__(None)
#        self.setWindowTitle('h5 File contents')
#        self.treeWidget=QTreeWidget()
#        mainlayout=QGridLayout()
#        mainlayout.addWidget(self.treeWidget, 0, 0)
#        self.setLayout(mainlayout)

        
def start(previousmm=None):
    mainapp=QApplication(sys.argv)

    #TW=TreeWindow()
    #TW.show()

    #form=MainMenu(TW.treeWidget)
    form=MainMenu(previousmm)
    form.show()
    form.setFocus()
    global PARENT
    PARENT=form
    mainapp.exec_()
    return form
mm=None
mm=start()
print 'done'

#fitd=getfitdictlist_hp(mm.h5path, 'heat1a','2010Nov27_Cell2_61mA_50ms_500ms_cool_1C')[0]
#segd=CreateHeatProgSegDictList(mm.h5path, 'heat1a', '2010Nov27_Cell2_61mA_50ms_500ms_cool_1C')[2]
#
#
#fild={}
#
#fild['reggrid']={'gridinterval':0.2}
#fild['peaksearch']={'pkfcn':'GaussHalfLorentz', 'critpeakheight':5.e-7, 'critsep':20., 'firstdernpts':10, 'firstderorder':1, 'secdernpts':20, 'secderorder':1, 'critcurve':None, 'pospeaks':1, 'negpeaks':1}
#def Cpk_secder(segd, fild, C, T, h5path=None, h5expname=None, h5hpname=None):
#    C=('sampleheatcapacity', 'peaksearch')
#    T=('sampletemperature', 'reggrid')
#for tup in [C, T]:
#    (segkey, filkey)=tup
#    if not '~'.join(tup) in segd.keys():
#        segd['~'.join(tup)]=performgenericfilter(segd[segkey], fild[filkey])
#    #if True in [v>0 for k, v in fild[filkey].iteritems() if 'deriv' in k]: #this handle deriv filters other than SG but if the deriv is not wrt dt something needs to change
#    #    segd['~'.join(tup)]/=dt
#C_=segd['~'.join(C)]
#T_=segd['~'.join(T)]
#(segkey, filkey)=T
#dx=fild[filkey]['gridinterval']
#(segkey, filkey)=C
#ch=fild[filkey]['critpeakheight']
#cs=fild[filkey]['critsep']/dx
#fp=fild[filkey]['firstdernpts']
#fo=fild[filkey]['firstderorder']
#sp=fild[filkey]['secdernpts']
#so=fild[filkey]['secderorder']
#
#X=T_[0]
#Y=C_[0]
#Xgrid=numpy.linspace(X.min(), X.max(), (X.max()-X.min())/dx+1)
#Ygrid=numpy.empty(Xgrid.shape, dtype='float32')
#gridind=[numpy.argmin((x-Xgrid)**2) for x in X]
#indsgot=numpy.sort(numpy.uint32(list(set(gridind))))
#indsinterp=numpy.sort(numpy.uint32(list(set(range(len(Xgrid)))-set(gridind))))
#gridind=numpy.uint32(gridind)
#for i in indsgot:
#    Ygrid[i]=Y[gridind==i].mean()
#Ygrid[indsinterp]=numpy.float32(scipy.interpolate.interp1d(indsgot, Ygrid[indsgot])(indsinterp))
#
#import pylab
#pylab.plot(X, Y, 'b.', markersize=1)
#pylab.plot(Xgrid, Ygrid, 'k-', lw=1)
#
##pkind=peaksearch1dSG(Ygrid, dx=dx, critcounts=ch, critsepind=cs, critcurve=None, max_withincritsep=False, firstdernpts=fp, firstderorder=fo, secdernpts=sp, secderorder=so)
#x=Ygrid
#dx=dx
#critcounts=ch
#critsepind=cs
#critcurve=None
#max_withincritsep=False
#firstdernpts=fp
#firstderorder=fo
#secdernpts=sp
#secderorder=so
#ifirstder=savgolsmooth(x, nptsoneside=firstdernpts, order=firstderorder, dx=dx, deriv=1)
#zeroind=arrayzeroind1d(ifirstder, postoneg=True)
#temp=numpy.where(x[(numpy.uint32(numpy.round(zeroind)),)]>critcounts)
#fullpkind=zeroind[temp]
#if fullpkind.size==0:
#    print '#$%^#$%^#%&$%^&'
#pkind=clustercoordsbymax1d(x, numpy.uint32(numpy.round(fullpkind)), critsepind)
#if critcurve is not None:
#    isecder=savgolsmooth(x, nptsoneside=secdernpts, order=secderorder, dx=dx, deriv=2)
#    temp=numpy.where(isecder[(numpy.uint32(numpy.round(pkind)),)]<(-1*critcurve))
#    pkind=numpy.array(pkind)[temp]
##    pkind=list(pkind)
##    pkind.reverse()#highest to smallest for pairing below
#pkind=numpy.array(pkind, dtype=numpy.float32)
#pkht=Ygrid[numpy.uint32(numpy.round(pkind))]
#pkposn=Xgrid[numpy.uint32(numpy.round(pkind))]
#iarr=numpy.uint32(range(len(Xgrid)))
#hwposns1=[(numpy.any((Ygrid<(h/2.))&(iarr>i)) and (numpy.where((Ygrid<(h/2.))&(iarr>i))[0][0],) or (i,))[0] for i, h in zip(pkind, pkht)]
#hwposns0=[(numpy.any((Ygrid<(h/2.))&(iarr<i)) and (numpy.where((Ygrid<(h/2.))&(iarr<i))[0][0],) or (i,))[0] for i, h in zip(pkind, pkht)]
#pkhw=dx*(numpy.float32(hwposns1)+numpy.float32(hwposns0))/2.
#pylab.plot(Xgrid[numpy.uint32(numpy.round(pkind))], Ygrid[numpy.uint32(numpy.round(pkind))], 'ro')
#
#pks=numpy.float32([pkposn, pkhw, pkht]).T#, numpy.ones(pkht.shape, dtype='float32')*.5]).T
#
#pkfcn=PeakFcnLibrary[fild[filkey]['pkfcn']]
#
#pars, sigs, resid=fitpeakset(X, Y, pks, GaussHalfLorentz)
#
#
#
#fitY=numpy.float32([GaussHalfLorentz(p, X) for p in pars]).sum(axis=0)
#gridfitY=numpy.float32([GaussHalfLorentz(p, Xgrid) for p in pars]).sum(axis=0)
#igridfitY=numpy.float32([gridfitY[:i+1].sum() for i in range(len(Xgrid))])*dx
#print igridfitY[-1]
#critval0_iY=igridfitY[-1]/25.
#critval1_iY=igridfitY[-1]/5.
#critval2_iY=igridfitY[-1]/2.
#print Xgrid[igridfitY>critval0_iY][0]
#print Xgrid[igridfitY>critval1_iY][0]
#print Xgrid[igridfitY>critval2_iY][0]
#
#
#pylab.plot(X, fitY, 'r-')
#pylab.show()



