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
    def __init__(self):#, TreeWidg):
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
        
        self.action_createh5=MainMenuQAction(self,'action_createh5', 'new h5 file', self.menufileio, [], self.ActionDict)
        self.action_createexpgrp=MainMenuQAction(self,'action_createexpgrp', 'new experiment group', self.menufileio, [('h5open', [True])], self.ActionDict)
        self.action_delexpgrp=MainMenuQAction(self,'action_delexpgrp', 'DELETE experiment group', self.menufileio, [('h5open', [True])], self.ActionDict)
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
        def expandchildren(item):
            for i in range(item.childCount()):
                child=item.child(i)
                if not groupsonly or True in [not ('(' in child.child(j).text(0) or str(child.child(j).text(0)).startswith("'")) for j in range(child.childCount())]:
                    child.setExpanded(True)
                    expandchildren(child)
        for i in range(self.treeWidget.topLevelItemCount()):
            item=self.treeWidget.topLevelItem(i) 
            item.setExpanded(True)
            expandchildren(item)
       
    def actionenable(self):
        for aname, ad in self.ActionDict.iteritems():
            ad['ref'].setDisabled(False in [(k in self.statusdict.keys()) and (self.statusdict[k] in vals) for k, vals in ad['enable_reqs']])
#            if aname=='action_calcresistance':
#                print [(k in self.statusdict.keys()) and (self.statusdict[k] in vals) for k, vals in ad['enable_reqs']]
#                print [(k, vals) for k, vals in ad['enable_reqs']]

    def h5nodename_treeitem(self, treeitem):
        return ((str(treeitem.text(0)).partition(':')[0]).partition('(')[0]).strip("'")

    def geth5selectionpath(self, liststyle=False):
        treeitem=self.currenttreeitem
        attrname=None
        if self.statusdict['selectiontype']=='Attr':
            attrname=self.h5nodename_treeitem(treeitem)
            treeitem=treeitem.parent()
        s=[]
        while not treeitem.parent() is None:
            s=[self.h5nodename_treeitem(treeitem)]+s
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

    @pyqtSignature("")
    def on_action_importscdata_triggered(self):
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
        ans=FileImport(self,protname)
        print ans
        if not ans:
            return
        AttrDict, DataSetDict, SegmentData=ans
        mA=DataSetDict['samplecurrent'][1][0]*DataSetDict['samplecurrent'][0]['Aunit']*1000.
        ms=1000.*numpy.float32(range(len(mA)))/AttrDict['daqHz']
        idialog=SegmentEditor(self, SegmentData, cycledata=(ms, mA))
        if not idialog.exec_():
            return
        SegmentData=idialog.SegmentData
        
        idialog=lineeditDialog(self, title='Enter name for h5 group', deftext=os.path.splitext(os.path.split(AttrDict['importpath'])[1])[0])
        if not idialog.exec_():
            return
        grpname=idialog.text
        
        writenewh5heatprogram(self.h5path, h5expname, grpname, AttrDict, DataSetDict, SegmentData)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
        
    @pyqtSignature("")
    def on_action_calcresistance_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        if self.statusdict['selectiongrouptype']=='experiment':
            h5file, hplist=experimenthppaths(self.h5path, pathlist[1])
            h5file.close()
            hplist=[hpp.rpartition('/')[2] for hpp in hplist]
        else:
            hplist=[pathlist[4]]
        for hp in hplist:
            print hp
            dlist=CreateHeatProgSegDictList(self.h5path, pathlist[1], hp) 
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
                dzero=dlist[segtypelist.index('zero')]
            else:
                dzero=None
            vals=[]
            vals+=[CalcR0_segdict(dsoak, AveBeforeDivision=True, dzero=dzero)]
            vals+=[CalcR0_segdict(dsoak, AveBeforeDivision=False, dzero=dzero)]
            vals+=[(vals[0]+vals[1])/2.]
            desc=['ratio of the means', 'mean of the ratios', 'Ave of these 2 values']
            choices=['%.4f : %s' %(v, d) for v, d in zip(vals, desc)]
            idialog=selectorDialog(self, choices, title='select value of R0 to use')
            if not idialog.exec_():
                return
            R0=vals[idialog.index]
            writecellres(self.h5path, pathlist[1], hp, R0)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
    
    @pyqtSignature("")
    def on_action_setuprescal_triggered(self):
        idialog=rescalDialog(self, self.h5path)
        idialog.exec_()

    @pyqtSignature("")
    def on_action_assignrescal_triggered(self):
        rescalpath_getorassign(self.h5path, self.statusdict['selectionname'], parent=self, forceassign=True)
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file)
        h5file.close()
        self.actionenable()


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
        
    @pyqtSignature("")
    def on_action_viewSCanalysis_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        self.data=CreateHeatProgSegDictList(self.h5path, pathlist[1], pathlist[4])
        idialog=analysisviewerDialog(self, self.data, pathlist[1])
        idialog.show()
    
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

    @pyqtSignature("")
    def on_action_applyscrecipe_triggered(self):
        pathlist=self.geth5selectionpath(liststyle=True)
        if self.statusdict['selectiongrouptype']=='heatprogram':
            h5hpdflt=pathlist[4]
        else:
            h5hpdflt=None
        idialog=SCanalysisDialog(self, self.h5path, pathlist[1], h5hpdflt=h5hpdflt)
        idialog.show()

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



def start():
    mainapp=QApplication(sys.argv)

    #TW=TreeWindow()
    #TW.show()

    #form=MainMenu(TW.treeWidget)
    form=MainMenu()
    form.show()
    form.setFocus()
    global PARENT
    PARENT=form
    mainapp.exec_()
    return form

mm=start()
print 'done'
