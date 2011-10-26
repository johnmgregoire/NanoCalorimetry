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


are_paths_equivalent=lambda path1, path2:os.path.normcase(os.path.abspath(path1))==os.path.normcase(os.path.abspath(path2))
    
class MainMenu(QMainWindow):
    def __init__(self, previousmm):#, TreeWidg):
        super(MainMenu, self).__init__(None)
        #self.setupUi(self)
        self.setWindowTitle('h5 editor')
        
        #self.treeWidget=TreeWidg
        
        self.h5path="%s" % os.getcwd()
        self.readh5path="%s" % os.getcwd()
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
        
        
        
        self.readh5treeWidget=QTreeWidget(self.bodywidget)
        QObject.connect(self.readh5treeWidget,SIGNAL("itemSelectionChanged()"),self.processtreeselectionreadh5)
#        self.setupmenu()
#        self.setCentralWidget(self.bodywidget)
        
        self.statusdict={'readh5open':False}
        self.actionenable()
        self.resize(1630, 620)
        self.readh5treeWidget.setGeometry(QRect(820, 10, 800, 520))
        
        
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
        self.openreadh5PushButton.setGeometry(QRect(1010, 550, 400, 25))
        self.openreadh5PushButton.setText('Open source h5 file')   
        QObject.connect(self.openreadh5PushButton, SIGNAL("pressed()"), self.openreadh5)
        
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
        self.setupmenu()
        self.setCentralWidget(self.bodywidget)
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
        self.action_batchimportdatadfltmsma=MainMenuQAction(self,'action_batchimportdatadfltmsma', 'batch import calorimetry data with no segment info', self.menufileio, [('h5open', [True])], self.ActionDict)
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
        myh5copy(h5file, g, readh5file[readpath])

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
        for node in readh5file[readpath].iteritems():
            myh5copy(h5file, g, node)
        for k, v in readh5file[readpath].attrs.iteritems():
            g.attrs[k]=v

        if not samepaths:
            readh5file.close()
        h5file.close()
        h5file=h5py.File(self.h5path, mode='r')
        fillh5tree(self.treeWidget, h5file, hpsortattr=str(self.sortattrLineEdit.text()))
        h5file.close()
        
def myh5copy(h5file, destgrp, srcnode):
    if isinstance(srcnode, h5py.Dataset):
        ds=destgrp.create_dataset(newname, data=readh5pyarray(srcnode))
        for k, v in srcnode.attrs.iteritems():
            ds.attrs[k]=v
    else:
        h5file.copy(srcnode, destgrp, name=newname)
        
class messageDialog(QDialog):
    def __init__(self, parent=None, title=''):
        super(messageDialog, self).__init__(parent)
        self.setWindowTitle(title)
        mainlayout=QGridLayout()
  
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox, SIGNAL("rejected()"), self.reject)
        mainlayout.addWidget(self.buttonBox, 0, 0)
         
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
    def ExitRoutine(self):
        return
        
def mygetopenfile(parent=None, xpath="%s" % os.getcwd(),markstr='', filename='' ):
    if parent is None:
        xapp = QApplication(sys.argv)
        xparent = QWidget()
        returnfn = unicode(QFileDialog.getOpenFileName(xparent,''.join(['Select file to open:', markstr]),os.path.join(xpath, filename).replace('\\','/')))
        xparent.destroy()
        xapp.quit()
        return returnfn
    return unicode(QFileDialog.getOpenFileName(parent,''.join(['Select file to open: ', markstr]),os.path.join(xpath, filename).replace('\\','/')))

def mygetopenfiles(parent=None, xpath="%s" % os.getcwd(),markstr='', filename='' ):
    if parent is None:
        xapp = QApplication(sys.argv)
        xparent = QWidget()
        returnfns=QFileDialog.getOpenFileNames(xparent,''.join(['Select file to open:', markstr]),os.path.join(xpath, filename).replace('\\','/'))
        xparent.destroy()
        xapp.quit()
    else:
        returnfns=QFileDialog.getOpenFileNames(parent,''.join(['Select file to open: ', markstr]),os.path.join(xpath, filename).replace('\\','/'))
    return [str(s) for s in returnfns]

def mygetsavefile(parent=None, xpath="%s" % os.getcwd(),markstr='', filename='' ):
    if parent is None:
        xapp = QApplication(sys.argv)
        xparent = QWidget()
        returnfn = unicode(QFileDialog.getSaveFileName(xparent,''.join(['Select file for save: ', markstr]),os.path.join(xpath, filename).replace('\\','/')))
        xparent.destroy()
        xapp.quit()
        return returnfn
    return unicode(QFileDialog.getSaveFileName(parent,''.join(['Select file for save: ', markstr]),os.path.join(xpath, filename).replace('\\','/')))

def mygetdir(parent=None, xpath="%s" % os.getcwd(),markstr='' ):
    if parent is None:
        xapp = QApplication(sys.argv)
        xparent = QWidget()
        returnfn = unicode(QFileDialog.getExistingDirectory(xparent,''.join(['Select directory:', markstr]), xpath))
        xparent.destroy()
        xapp.quit()
        return returnfn
    return unicode(QFileDialog.getExistingDirectory(parent,''.join(['Select directory:', markstr]), xpath))

class fillh5tree():
    def __init__(self, tree, h5file, showattrs=True, selectionpathlist=None, hpsortattr='epoch'):
        self.treeWidget=tree
        self.treeWidget.clear()
        self.showattrs=showattrs
        if hpsortattr=='':
            hpsortattr=None
        self.hpsortattr=hpsortattr
        mainitem=QTreeWidgetItem([os.path.split(h5file.filename)[1]],  0)
        self.treeWidget.addTopLevelItem(mainitem)
        self.createTree(h5file, mainitem)
        
        if not selectionpathlist is None:
            item=mainitem
            for itemname in selectionpathlist:
                chn=[item.child(i).text(0) for i in range(item.childCount())]
                if itemname in chn:
                    item.setExpanded(True)
                    item=item.child(chn.index(itemname))
                else:
                    break
            tree.setCurrentItem(item)


    def createTree(self, startnode, parentitem):
        #print startnode.name
        #print startnode.listobjects()
        sortbool=(not self.hpsortattr is None) and 'HeatProgram' in parentitem.text(0)
        if sortbool:
            sortvals=[]
            items=[]
        for node in startnode.iterobjects():
            if isinstance(node, h5py.Dataset):
                item=QTreeWidgetItem([node.name.rpartition('/')[2]+`node.shape`],  0)
                if self.showattrs:
                    for attrname, attrval in node.attrs.iteritems():
                        attritem=QTreeWidgetItem([self.attrstring(attrname, attrval)],  0)
                        item.addChild(attritem)
            elif isinstance(node, h5py.Group):
                item=QTreeWidgetItem([node.name.rpartition('/')[2]],  0)
                self.createTree(node, item)
                if self.showattrs:
                    for attrname, attrval in node.attrs.iteritems():
                        attritem=QTreeWidgetItem([self.attrstring(attrname, attrval)],  0)
                        item.addChild(attritem)
            
            if sortbool:
                items+=[item]
                sortvals+=[(self.hpsortattr in node.attrs.keys() and (node.attrs[self.hpsortattr],) or (0.,))[0]]
            else:
                parentitem.addChild(item)
        if sortbool:
            for i in numpy.argsort(sortvals):
                parentitem.addChild(items[i])

    def attrstring(self, attrname, attrval):
        s="'"+attrname+"':"
        try:
            if isinstance(attrval, str):
                if len(attrval)>100:
                    s+=attrval[:20]+' ... '+attrval[-20:]
                else:
                    s+=attrval
            elif isinstance(attrval, int) or isinstance(attrval, float):
                s+=self.numfmt(attrval)
            elif isinstance(attrval, list) or isinstance(attrval, numpy.ndarray):
                temp=attrval
                temp2=attrval
                ndim=0
                while isinstance(temp, list) or isinstance(temp, numpy.ndarray):
                    if len(temp)==0 or len(temp2)==0:
                        s+='contains empty list'
                        return s
                    temp=temp[0]
                    temp2=temp2[-1]
                    ndim+=1

                    if isinstance(temp, str):
                        attrvalstr=`attrval`
                        attrvalstr=attrvalstr.partition('(')[2].rpartition(',')[0]
                        if len(attrvalstr)>100:
                            s+=attrvalstr[:20]+' ... '+attrvalstr[-20:]
                        else:
                            s+=attrvalstr
                        return s
                if ndim==1:
                    if len(attrval)<10:
                        s+='['+','.join([self.numfmt(attrel) for attrel in attrval])+']'
                    else:
                       s+= '['+',...,'.join([self.numfmt(attrel) for attrel in [temp, temp2]])+']'
                elif ndim==2:
                    s+= '[]'+',..][..,'.join([self.numfmt(attrel) for attrel in [temp, temp2]])+']]'
                else:
                    s+='%d' %ndim +' dimmension structure with first value of '+self.numfmt(temp)
            else:
                raise
        except:
            s+='type is '+`type(attrval)`
        return s

    def numfmt(self, num):
        if isinstance(num, int):
            s='%d' %num
        elif num==0.:
            s='0.0'
        elif numpy.abs(num)<100 and numpy.abs(num)>=.01:
            s='%.4f' %num
        else:
            s=myexpformat(num)
        return s
def myexpformat(x, pos=0):
    for ndigs in range(6):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','')
        if eval(lab)==x:
            return lab
    return lab
    
class lineeditDialog(QDialog):
    def __init__(self, parent, title='', deftext=''):
        super(lineeditDialog, self).__init__(parent)

        self.setWindowTitle(title)
        self.le=QLineEdit()
        self.le.setText(deftext)

        mainlayout=QGridLayout()
        mainlayout.addWidget(self.le, 0, 0)
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        mainlayout.addWidget(self.buttonBox, 1, 0)
         
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        self.setLayout(mainlayout)
    def ExitRoutine(self):
        self.text=str(self.le.text())

def readh5pyarray(arrpoint):
    return eval('arrpoint'+('['+':,'*len(arrpoint.shape))[:-1]+']')
    
class attreditorDialog(QDialog):
    def __init__(self, parent, attrd, arr=None, title=''):
        super(attreditorDialog, self).__init__(parent)

        arrbool=not arr is None
        self.setWindowTitle(title)
        self.lw=QListWidget()
        self.lw.setObjectName("lw")
        #self.lw.setEditTriggers(QAbstractItemView.AllEditTriggers)
        #QObject.connect(self.lw, SIGNAL("itemDoubleClicked(QListWidgetItem)"), self.edititem(QListWidgetItem))
        #QObject.connect(self.lw, SIGNAL("itemDoubleClicked()"), self.edititem)
 


        mainlayout=QGridLayout()
        
        if arrbool:
            plotw=plotwidget(self)
            mainlayout.addWidget(plotw, 0, 0, 1, 2)
            for a in arr:
                plotw.axes.plot(a, '.')
            plotw.axes.set_xlabel('array index')
            plotw.axes.set_ylabel('array value')
            mainlayout.addWidget(self.lw, 0, 2, 1, 1)
        else:
            mainlayout.addWidget(self.lw, 0, 0, 1, 1)    
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox, SIGNAL("rejected()"), self.reject)
        mainlayout.addWidget(self.buttonBox, 1, 0)
         
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        self.setLayout(mainlayout)
        
        self.attrd=attrd
        self.strlist=[]
        for k, v in attrd.iteritems():
            s="'%s':%s" %(k, `v`)
            self.strlist+=[s]
            it=self.lw.addItem(s)
        for ind in range(self.lw.count()):
            it=self.lw.item(ind)
            it.setFlags(it.flags() | Qt.ItemIsEditable)
        #self.lw.setEditTriggers(QAbstractItemView.AllEditTriggers)
        #QObject.connect(self.lw,SIGNAL("itemSelectionChanged()"),self.updateattrd)
        
        QMetaObject.connectSlotsByName(self)
#    def edititem(self):
#        print '*****', #self.lw.currentItem().text()
#        self.lw.editItem(self.lw.currentItem())    
#    def updateattrd(self):
#        s=str(self.lw.currentItem().text())
#        a, b, c=s.partition(':')
#        try:
#            c=eval(c)
#        except:
#            pass
#        self.attrd[a]=c
#        print a,  ' updated to ',  `c`

    def ExitRoutine(self):
        self.edited=False
        for ind, origs in enumerate(self.strlist):
            it=self.lw.item(ind)
            s=str(it.text())
            if origs!=s:
                a, b, c=s.partition(':')
                try:
                    c=myeval(c)
                except:
                    pass
                a=a.strip("'").strip('"')
                print a, a in self.attrd.keys()
                self.attrd[a]=c
                print a,  ' updated to ',  `c`
                self.edited=True
def editattrs(parent, h5path, path):
    h5file=h5py.File(h5path, mode='r')
    ad=dict([(k, v) for k, v in h5file[path].attrs.iteritems()])
    h5file.close()
    repeat=True
    edited=True
    while repeat and edited:
        idialog=attreditorDialog(parent, ad)
        repeat=idialog.exec_()
        edited=idialog.edited
    if repeat:
        h5file=h5py.File(h5path, mode='r+')
        for k, v in ad.iteritems():
            h5file[path].attrs[k]=v
        h5file.close()
    
class MainMenuQAction(QAction):
    def __init__(self, parent, actionname, actiontext, hostmenu, reqs, adict):
        super(MainMenuQAction, self).__init__(parent)
        self.setObjectName(actionname)
        self.setText(actiontext)
        hostmenu.addAction(self)
        adict[actionname]={'ref':self, 'enable_reqs':reqs}

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
