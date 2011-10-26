import time
import os
import sys
import numpy
import h5py
from PnSC_ui import *
#from PnSC_dataimport import *
#from PnSC_SCui import *
#from PnSC_math import *
from PnSC_h5io import *
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class MainMenu(QMainWindow):
    def __init__(self, previousmm):#, TreeWidg):
        super(MainMenu, self).__init__(None)
        #self.setupUi(self)
        self.setWindowTitle('Vlassak Group PnSC Analysis')
        
        #self.treeWidget=TreeWidg
        
        self.h5path="%s" % os.getcwd()
        self.readh5path="%s" % os.getcwd()
        
        self.bodywidget = QWidget(self)
        self.bodywidget.setObjectName("bodywidget")
        self.treeWidget=QTreeWidget(self.bodywidget)
        QObject.connect(self.treeWidget,SIGNAL("itemSelectionChanged()"),self.processtreeselection)
        
        self.readh5treeWidget=QTreeWidget(self.bodywidget)
        QObject.connect(self.readh5treeWidget,SIGNAL("itemSelectionChanged()"),self.processtreeselectionreadh5)
#        self.setupmenu()
#        self.setCentralWidget(self.bodywidget)
        
        self.setupmenu()
        self.setCentralWidget(self.bodywidget)
        
        self.statusdict={'readh5open':False}
        self.actionenable()
        self.resize(1630, 620)
        self.readh5treeWidget.setGeometry(QRect(820, 10, 800, 520))
        
        self.statusdict={'h5open':False}
        self.actionenable()
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
        self.sortattrLineEdit = QLineEdit(self.bodywidget)
        self.sortattrLineEdit.setGeometry(QRect(610, 550, 200, 25))
        self.sortattrLineEdit.setText('epoch')
        self.openreadh5PushButton = QPushButton(self.bodywidget)
        self.openreadh5PushButton.setGeometry(QRect(910, 550, 300, 25))
        self.openreadh5PushButton.setText('Open source h5 file')   
        QObject.connect(self.openreadh5PushButton, SIGNAL("pressed()"), self.openreadh5)
        self.openreadh5PushButton = QPushButton(self.bodywidget)
        self.openreadh5PushButton.setGeometry(QRect(1210, 550, 300, 25))
        self.openreadh5PushButton.setText('Auto copy spec scans to h5 root')   
        QObject.connect(self.openreadh5PushButton, SIGNAL("pressed()"), self.autocopyspec)
        
        if previousmm is None:
            self.on_action_openh5_triggered()
        else:
            oldselection=mm.geth5selectionpath(liststyle=True, removeformatting=False)
            self.h5path=previousmm.h5path
            h5file=h5py.File(self.h5path, mode='r')
            fillh5tree(self.treeWidget, h5file, selectionpathlist=oldselection, hpsortattr=str(self.sortattrLineEdit.text()))
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
        self.action_createh5=MainMenuQAction(self,'action_createh5', 'new h5 file', self.menufileio, [], self.ActionDict)
        self.action_createexpgrp=MainMenuQAction(self,'action_createexpgrp', 'new experiment group', self.menufileio, [('h5open', [True])], self.ActionDict)
        self.action_delh5grp=MainMenuQAction(self,'action_delh5grp', 'DELETE selected group', self.menufileio, [('h5open', [True]), ('selectiontype', ['Group'])], self.ActionDict)
        #self.action_delexpgrp=MainMenuQAction(self,'action_delexpgrp', 'DELETE experiment group', self.menufileio, [('h5open', [True])], self.ActionDict)
        self.action_editattrs=MainMenuQAction(self,'action_editattrs', 'Edit import attrs (select a heat program)', self.menufileio, [('h5open', [True]), ('selectiongrouptype', ['heatprogram'])], self.ActionDict)
        self.action_copynode=MainMenuQAction(self,'action_copynode', 'Copy read-only selected group or dataset to selected group', self.menufileio, [('h5open', [True]), ('selectiontype', ['Group', 'File']), ('readh5selectiontype', ['Group', 'Dataset'])], self.ActionDict)
        self.action_copygroup=MainMenuQAction(self,'action_copygroup', 'Copy contents of read-only selected group to selected group', self.menufileio, [('h5open', [True]), ('selectiontype', ['Group', 'File']), ('readh5selectiontype', ['Group'])], self.ActionDict)
  
        self.setMenuBar(self.main_menu_pulldown)
        QMetaObject.connectSlotsByName(self)
    
    def redraw(self):
        if os.path.exists(self.h5path) and self.h5path.endswith('.h5'):
            h5file=h5py.File(self.h5path, mode='r')
            fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
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

    def geth5selectionpath(self, liststyle=False, removeformatting=True, treeitem=None):
        try:
            if treeitem is None:
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
            if not '/' in s:
                s='/'+s
        if not attrname is None:
            return s, attrname
        return s
        
    def processtreeselectionreadh5(self):
        self.processtreeselection(readh5=True)
        
    def processtreeselection(self, readh5=False):
        if readh5:
            treeitem=self.readh5treeWidget.currentItem()
            readstr='readh5'
        else:
            readstr=''
            treeitem=self.treeWidget.currentItem()
            self.currenttreeitem=treeitem
        print 'selection changed to ', treeitem.text(0)
        if treeitem.parent() is None:
            self.statusdict[readstr+'selectiontype']='File'
        elif str(treeitem.text(0)).startswith("'"):
            self.statusdict[readstr+'selectiontype']='Attr'
        elif '(' in treeitem.text(0):
            self.statusdict[readstr+'selectiontype']='Dataset'
        else:
            self.statusdict[readstr+'selectiontype']='Group'
        
        self.statusdict[readstr+'selectionname']=self.h5nodename_treeitem(treeitem)
        if self.statusdict[readstr+'selectiontype']=='File':
            self.statusdict[readstr+'selectionparentname']=''
        else:
            self.statusdict[readstr+'selectionparentname']=self.h5nodename_treeitem(treeitem.parent())
        
        self.statusdict[readstr+'samplepowerperrateexists']=True#TODO: write code for checking on existence of analysis data arrays
        
        if self.statusdict[readstr+'selectiontype']=='Group':
            if self.statusdict[readstr+'selectionparentname']=='HeatProgram':
                self.statusdict[readstr+'selectiongrouptype']='heatprogram'
            elif self.statusdict[readstr+'selectionparentname']=='Calorimetry':
                self.statusdict[readstr+'selectiongrouptype']='experiment'
            elif self.statusdict[readstr+'selectionparentname']=='analysis':
                self.statusdict[readstr+'selectiongrouptype']='analysis'
            else:
                self.statusdict[readstr+'selectiongrouptype']='other'
        else:
            self.statusdict[readstr+'selectiongrouptype']=''
        
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
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
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
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
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
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
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
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
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
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
        h5file.close()

    def openreadh5(self):
        self.statusdict['readh5open']=False
        self.actionenable()
        temp=mygetopenfile(parent=self, xpath=self.readh5path, markstr='h5 file with calorimetry data', filename='.h5' )
        if temp=='':
            return
        self.readh5path=str(temp)
        h5file=h5py.File(self.readh5path, mode='r')
        fillh5tree(self.readh5treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
        h5file.close()
        self.statusdict['readh5open']=True
        self.actionenable()

    @pyqtSignature("")
    def on_action_copynode_triggered(self):
        writepath=self.geth5selectionpath(liststyle=False)
        readpath=self.geth5selectionpath(liststyle=False, treeitem=self.readh5treeWidget.currentItem())
        
        idialog=lineeditDialog(self, title='Enter new name', deftext=readpath.rpartition('/')[2])
        if not idialog.exec_():
            return
        newname=idialog.text
            
        h5file=h5py.File(self.h5path, mode='r+')
        samepaths=are_paths_equivalent(self.readh5path, self.h5path)
        if samepaths:
            readh5file=h5file
        else:
            readh5file=h5py.File(self.readh5path, mode='r')
        
        g=h5file[writepath]
        myh5copy(h5file, g, readh5file[readpath], newname=newname)

        if not samepaths:
            readh5file.close()
        h5file.close()
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
        h5file.close()

    @pyqtSignature("")
    def on_action_copygroup_triggered(self):
        writepath=self.geth5selectionpath(liststyle=False)
        readpath=self.geth5selectionpath(liststyle=False, treeitem=self.readh5treeWidget.currentItem())

        h5file=h5py.File(self.h5path, mode='r+')
        samepaths=are_paths_equivalent(self.readh5path, self.h5path)
        if samepaths:
            readh5file=h5file
        else:
            readh5file=h5py.File(self.readh5path, mode='r')
        
        g=h5file[writepath]
        for node in readh5file[readpath].itervalues():
            myh5copy(h5file, g, node)
        for k, v in readh5file[readpath].attrs.iteritems():
            g.attrs[k]=v

        if not samepaths:
            readh5file.close()
        h5file.close()
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
        h5file.close()

    def autocopyspec(self):
        h5file=h5py.File(self.h5path, mode='r')
        scannames=[]
        for node in h5file['Calorimetry'].itervalues():
            if not 'measurement/HeatProgram' in node:
                continue
            for node2 in node['measurement/HeatProgram'].itervalues():
                for k in ['specscan', 'postspecscan', 'prespecscan']:
                    if k in node2.attrs.keys():
                        v=node2.attrs[k]
                        if not isinstance(v, str):
                            v='%d' %v
                        scannames+=[v]
        h5file.close()
        scannames=set(scannames)
        print 'preparing to copy ', scannames
        
        h5file=h5py.File(self.h5path, mode='r+')
        samepaths=are_paths_equivalent(self.readh5path, self.h5path)
        if samepaths:
            readh5file=h5file
        else:
            readh5file=h5py.File(self.readh5path, mode='r')
        
        for v in scannames:
            print 'copying ', v
            myh5copy(h5file, h5file['/'], readh5file[v], overwrite=False)
        if not samepaths:
            readh5file.close()
        h5file.close()
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
        h5file.close()
    
def myh5copy(h5file, destgrp, srcnode, newname=None, overwrite='ask'):
    if newname is None:
        newname=srcnode.name.rpartition('/')[2]
    if newname in destgrp:
        if overwrite=='ask':
            idialog=messageDialog(self, 'press OK to overwrite ' +newname)
            overwrite=idialog.exec_()
        if not overwrite:
            return
        del destgrp[newname]
    if isinstance(srcnode, h5py.Dataset):
        ds=destgrp.create_dataset(newname, data=readh5pyarray(srcnode))
        for k, v in srcnode.attrs.iteritems():
            ds.attrs[k]=v
    else:
        h5file.copy(srcnode, destgrp, name=newname)
        
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
