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
from PnSC_h5io import *
import PnSC_h5io as io

from matplotlib.ticker import FuncFormatter

def myexpformat(x, pos):
    for ndigs in range(5):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

are_paths_equivalent=lambda path1, path2:os.path.normcase(os.path.abspath(path1))==os.path.normcase(os.path.abspath(path2))
    
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
                nam=node.name.rpartition('/')[2]
#                if nam=='epoch' or nam=='Epoch':
#                    item=QTreeWidgetItem([nam+(`node.shape`[:-1])+!!!],  0)
#                else:
                item=QTreeWidgetItem([nam+`node.shape`],  0)
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
                if attrname=='epoch' or attrname=='Epoch':
                    s+='%.3f' %attrval
                else:
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

class plotwidget(FigureCanvas):
    def __init__(self, parent, width=12, height=6, dpi=72):

        #plotdata can be 2d array for image plot or list of 2 1d arrays for x-y plot or 2d array for image plot or list of lists of 2 1D arrays
        
        self.fig=Figure(figsize=(width, height), dpi=dpi)
        self.axes=self.fig.add_subplot(111, navigate=True)
        
        
        self.axes.hold(True)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        #self.parent=parent
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        #NavigationToolbar(self, parent)
        NavigationToolbar(self, self)
        
        self.mpl_connect('button_press_event', self.myclick)
        self.clicklist=[]
    
    def myclick(self, event):
        if not (event.xdata is None or event.ydata is None):
            arrayxy=[event.xdata, event.ydata]
            print 'clicked on image: array indeces ', arrayxy
            self.clicklist+=[arrayxy]
            self.emit(SIGNAL("genericclickonplot"), [event.xdata, event.ydata])

class attreditorDialog(QDialog):
    def __init__(self, parent, attrd, arr=None, title=''):
        super(attreditorDialog, self).__init__(parent)
        self.edited=False
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

class SegmentCyclePlot(QDialog):
    def __init__(self, parent, SegmentData, markersize=2):
        super(SegmentCyclePlot, self).__init__(parent)

        self.setWindowTitle('segments by color and cycles by color')
        
        cols=['k', 'y', 'b', 'm', 'g', 'c', 'r']*10#just to be sure there's enough

        self.plotw=plotwidget(self)
        self.fig=self.plotw.fig
        self.fig.clf()
        self.plotwlist=[[self.fig.add_subplot(2, 3, i*3+j+1) for j in range(3)] for i in range(2)]
        self.fig.subplots_adjust(left=.09, bottom=.1, right=.95, top=.92, wspace=.4, hspace=.3)
        self.markersize=markersize
        
        mainlayout=QGridLayout()
        mainlayout.addWidget(self.plotw, 0, 0)
        for i, lab1 in enumerate(['samplecurrent', 'samplevoltage']):
            for j, lab2 in enumerate(['1st cycle, %d segs' %len(SegmentData), 'all cycles', '1st 6 and last (red) cycles']):
                self.plotwlist[i][j].set_xlabel('cycle time (s)')
                self.plotwlist[i][j].set_ylabel(lab1)
                self.plotwlist[i][j].set_title(lab2)
        for d, col in zip(SegmentData, cols):
            first=True
            for count, (i, v) in enumerate(zip(d['samplecurrent'], d['samplevoltage'])):
                if first:
                    first=False
                    self.plotwlist[0][0].plot(d['cycletime'][0], i, col+'.', markersize=markersize)
                    self.plotwlist[1][0].plot(d['cycletime'][0], v, col+'.', markersize=markersize)
                self.plotwlist[0][1].plot(d['cycletime'][0], i, col+'.', markersize=markersize)
                self.plotwlist[1][1].plot(d['cycletime'][0], v, col+'.', markersize=markersize)
                if count<6 or count==(d['samplecurrent'].shape[0]-1):
                    if count<6:
                        col2=cols[count]
                    else:
                        col2=cols[6]
                    self.plotwlist[0][2].plot(d['cycletime'][0], i, col2+'.', markersize=markersize)
                    self.plotwlist[1][2].plot(d['cycletime'][0], v, col2+'.', markersize=markersize)
        self.setLayout(mainlayout)
        
#class SegmentCyclePlot(QDialog):
#    def __init__(self, parent, SegmentData, markersize=2):
#        super(SegmentCyclePlot, self).__init__(parent)
#
#        self.setWindowTitle('segments by color and cycles by color')
#        
#        cols=['k', 'y', 'b', 'm', 'g', 'c', 'r']*10#just to be sure there's enough
#
#        self.plotwlist=[[plotwidget(self) for j in range(3)] for i in range(2)]
#
#        self.markersize=markersize
#        
#        mainlayout=QGridLayout()
#        for i, lab1 in enumerate(['samplecurrent', 'samplevoltage']):
#            for j, lab2 in enumerate(['first cycle', 'all cycles', '1st 6 and last (red) cycles']):
#                mainlayout.addWidget(self.plotwlist[i][j], i, j)
#                self.plotwlist[i][j].axes.set_xlabel('cycle time (ms)')
#                self.plotwlist[i][j].axes.set_ylabel(lab1)
#                self.plotwlist[i][j].axes.set_title(lab2)
#        for d, col in zip(SegmentData, cols):
#            first=True
#            for count, (i, v) in enumerate(zip(d['samplecurrent'], d['samplevoltage'])):
#                if first:
#                    first=False
#                    self.plotwlist[0][0].axes.plot(d['cycletime'][0], i, col+'.', markersize=markersize)
#                    self.plotwlist[1][0].axes.plot(d['cycletime'][0], v, col+'.', markersize=markersize)
#                self.plotwlist[0][1].axes.plot(d['cycletime'][0], i, col+'.', markersize=markersize)
#                self.plotwlist[1][1].axes.plot(d['cycletime'][0], v, col+'.', markersize=markersize)
#                if count<6 or count==(d['samplecurrent'].shape[0]-1):
#                    if count<6:
#                        col2=cols[count]
#                    else:
#                        col2=cols[6]
#                    self.plotwlist[0][2].axes.plot(d['cycletime'][0], i, col2+'.', markersize=markersize)
#                    self.plotwlist[1][2].axes.plot(d['cycletime'][0], v, col2+'.', markersize=markersize)
#        self.setLayout(mainlayout)

class SegmentEditor(QDialog):
    def __init__(self, parent, SegmentData, cycledata, maxpts=15, markersize=2):
        super(SegmentEditor, self).__init__(parent)

        self.setWindowTitle('Provide the endpoint of segments in a single cycle')
        
        self.plotw=plotwidget(self)
        self.markersize=markersize
        self.cycledata=cycledata
        self.plotw.axes.plot(self.cycledata[0], self.cycledata[1], 'r.', ms=self.markersize)

        self.plotw.axes.set_xlabel('cycle time (ms)')
        self.plotw.axes.set_ylabel('applied current (mA)')
            
        mainlayout=QGridLayout()

        pwid=4
        nrows_ctrls=6
        mainlayout.addWidget(self.plotw, 0, 0, nrows_ctrls+maxpts, pwid)
        
        firstderptsLabel=QLabel()
        firstderptsLabel.setText('num pts/n1st der')
        self.firstderptsSpinBox=QSpinBox()
        self.firstderptsSpinBox.setRange(5, 101)
        self.firstderptsSpinBox.setValue(15)

        secderptsLabel=QLabel()
        secderptsLabel.setText('num pts/n2nd der')
        self.secderptsSpinBox=QSpinBox()
        self.secderptsSpinBox.setRange(5, 101)
        self.secderptsSpinBox.setValue(30)

        secdervalLabel=QLabel()
        secdervalLabel.setText('crit val\n2nd der')
        self.secdervalSpinBox=QDoubleSpinBox()
        self.secdervalSpinBox.setDecimals(8)
        self.secdervalSpinBox.setValue(0.003)
        
#        segcalcrowLabel=QLabel()
#        secdervalLabel.setText('')
#        self.segcalcrowComboBox=QComboBox()
#        self.segcalcrowComboBox.clear()
#        self.segcalcrowComboBox.insertItem(-1, 'ave')
        
        findsegsButton=QPushButton()
        findsegsButton.setText("find segs using derivatives")
        QObject.connect(findsegsButton, SIGNAL("pressed()"), self.findsegs)
        
        mainlayout.addWidget(firstderptsLabel, 0, pwid, 1, 1)
        mainlayout.addWidget(self.firstderptsSpinBox, 0, pwid+1, 1, 1)
        mainlayout.addWidget(secderptsLabel, 1, pwid, 1, 1)
        mainlayout.addWidget(self.secderptsSpinBox, 1, pwid+1, 1, 1)
        
        mainlayout.addWidget(secdervalLabel, 2, pwid, 1, 1)
        mainlayout.addWidget(self.secdervalSpinBox, 2, pwid+1, 1, 1)
        mainlayout.addWidget(findsegsButton, 3, pwid, 1, 2)
        
        delButton=QPushButton()
        delButton.setText("del row")
        QObject.connect(delButton, SIGNAL("pressed()"), self.delrow)

        addButton=QPushButton()
        addButton.setText("add row")
        QObject.connect(addButton, SIGNAL("pressed()"), self.addrow)

        self.rowComboBox=QComboBox()
        self.rowComboBox.clear()
        
        mainlayout.addWidget(delButton, 4, pwid, 1, 1)
        mainlayout.addWidget(addButton, 5, pwid, 1, 1)
        mainlayout.addWidget(self.rowComboBox, 4, pwid+1, 2, 1)
        
        tLabel=QLabel()
        tLabel.setText('cycletime')
        vLabel=QLabel()
        vLabel.setText('mA')
        mainlayout.addWidget(tLabel, nrows_ctrls, pwid, 1, 1)
        mainlayout.addWidget(vLabel, nrows_ctrls, pwid+1, 1, 1)
        
        if not cycledata is None:
            self.tmin, self.tmax=(cycledata[0].min(), cycledata[0].max())
            self.cmin, self.cmax=(cycledata[1].min(), cycledata[1].max())
        else:
            tmin, tmax=(0, 999999.)
            cmin, cmax=(0, 999999.)
        self.dflts=[self.tmax, 0]
        self.tSpinBoxList=[]
        self.cSpinBoxList=[]
        for i in range(maxpts):
            self.rowComboBox.insertItem(i, `i`)
            tsb=QDoubleSpinBox()
            csb=QDoubleSpinBox()
            tsb.setRange(cycledata[0][0], cycledata[0][-1])
            csb.setRange(cycledata[0][0], cycledata[0][-1])
            tsb.setDecimals(3)
            csb.setDecimals(3)
            self.fillspinboxlist(SegmentData)
            mainlayout.addWidget(tsb, i+1+nrows_ctrls, pwid, 1, 1)
            mainlayout.addWidget(csb, i+1+nrows_ctrls, pwid+1, 1, 1)
#            QObject.connect(tsb, SIGNAL("valueChanged(double)"), self.plotsegments)
#            QObject.connect(csb, SIGNAL("valueChanged(double)"), self.plotsegments)
            self.tSpinBoxList+=[tsb]
            self.cSpinBoxList+=[csb]
        
        self.rowComboBox.setCurrentIndex(0)
        
        plotButton=QPushButton()
        plotButton.setText('plot segments')
        QObject.connect(plotButton, SIGNAL("pressed()"), self.plotsegments)
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        mainlayout.addWidget(plotButton, nrows_ctrls+maxpts+1, pwid, 1, 2)
        mainlayout.addWidget(self.buttonBox, nrows_ctrls+maxpts+1, 0, 1, pwid)
        
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        self.setLayout(mainlayout)
    
    def fillspinboxlist(self, SegmentData):
        for i in range(len(self.tSpinBoxList)):
            if i<len(SegmentData[0]):
                self.tSpinBoxList[i].setValue(SegmentData[0][i])
            else:
                self.tSpinBoxList[i].setValue(self.dflts[0])
            if i<len(SegmentData[1]):
                self.cSpinBoxList[i].setValue(SegmentData[1][i])
            else:
                self.cSpinBoxList[i].setValue(self.dflts[1])
                
    def addrow(self):
        ind=self.rowComboBox.currentIndex()
        for i in range(ind+1, len(self.tSpinBoxList))[::-1]:
            self.tSpinBoxList[i].setValue(self.tSpinBoxList[i-1].value())
            self.cSpinBoxList[i].setValue(self.cSpinBoxList[i-1].value())

    def delrow(self):
        ind=self.rowComboBox.currentIndex()
        if ind==(len(self.tSpinBoxList)-1):
            return
        for i in range(ind, len(self.tSpinBoxList)-1):
            self.tSpinBoxList[i].setValue(self.tSpinBoxList[i+1].value())
            self.cSpinBoxList[i].setValue(self.cSpinBoxList[i+1].value())

        self.tSpinBoxList[-1].setValue(self.dflts[0])
        self.cSpinBoxList[-1].setValue(self.dflts[1])
    
    def findsegs(self):#the analysis is done with normalization of the array so it is insensitive to 'Aunit'. Also, the data time interval is not used in the derivative calculation
        fd_nptsoneside=self.firstderptsSpinBox.value()
        sd_nptsoneside=self.secderptsSpinBox.value()
        critval_sd=self.secdervalSpinBox.value()
        arr=self.cycledata[1]/self.cycledata[1].max()

        fd=savgolsmooth(arr, nptsoneside=sd_nptsoneside, order=3, deriv=1)
        sd=savgolsmooth(arr, nptsoneside=fd_nptsoneside, order=3, deriv=2)

        inds=findlocmax(numpy.abs(sd), critval=critval_sd)
        avelen=10
        l=avelen
        segvals=numpy.array([(abs(fd[max(0, i-l):i].mean())<abs(fd[i:min(len(fd), i+l)].mean())) and arr[max(0, i-l):i].mean() or arr[i:min(len(fd), i+l)].mean() for i in inds])
        self.fillspinboxlist((numpy.append(0., self.cycledata[0][inds]), numpy.append(0., segvals*self.cycledata[1].max())))
        
    def plotsegments(self):
        self.readsegdata()
        self.plotw.axes.cla()
        self.plotw.axes.plot(self.cycledata[0], self.cycledata[1], 'r.', ms=self.markersize)
        carr=numpy.zeros(self.cycledata[1].shape, dtype=self.cycledata[1].dtype)
        
        for counter, (t, c) in enumerate(zip(self.SegmentData[0], self.SegmentData[1])):
            if counter==0:
                iprev=0
                cprev=c
                continue
            temp=numpy.where(self.cycledata[0]<t)[0]
            if len(temp)==0 or temp[-1]<=iprev-1:#if len(temp)==0 or temp[-1]==iprev:
                continue
            i=temp[-1]
            if i+1<len(self.cycledata[0]):
                i+=1
            carr[iprev:i]=(c-cprev)/(i-iprev)*numpy.float32(range(i-iprev))+cprev
            iprev=i
            cprev=c
        self.plotw.axes.plot(self.cycledata[0], carr, 'k.', ms=self.markersize)
        self.plotw.fig.canvas.draw()

    def readsegdata(self):
        tprev=-1
        tvl=[]
        cvl=[]
        for tsb, csb in zip(self.tSpinBoxList, self.cSpinBoxList):
            tv=tsb.value()
            cv=csb.value()
            if tv==tprev or tv==self.tmax:
                break
            tvl+=[tv]
            cvl+=[cv]
            tprev=tv
        self.SegmentData=(numpy.float32(tvl), numpy.float32(cvl))
    def ExitRoutine(self):
        self.readsegdata()

class PatDAQCycleEditor(QDialog):
    def __init__(self, parent, mA, durationguess, filename='', markersize=2):
        super(PatDAQCycleEditor, self).__init__(parent)

        self.setWindowTitle('Ensure that the period and phase of each cycle is correct')
        
        self.mA=mA
        self.markersize=markersize
        self.plotw=plotwidget(self)

        self.plotw.axes.set_xlabel('array index')
        self.plotw.axes.set_ylabel('applied current (arb)')
        
        fnLabel=QLabel()
        fnLabel.setText('filename: '+filename)
        durLabel=QLabel()
        durLabel.setText('num indeces in cycle (e.g. 50000 is 500ms at 100kHz)')
        nnoiseLabel=QLabel()
        nnoiseLabel.setText('num indeces to establish noise level')
        nsigLabel=QLabel()
        nsigLabel.setText('num sigma to define trigger value')
        naboveLabel=QLabel()
        naboveLabel.setText('critical num indeces above trigger value')
        
        self.durSpinBox=QSpinBox()
        self.durSpinBox.setRange(1, len(self.mA))
        self.durSpinBox.setValue(int(durationguess))
        self.nnoiseSpinBox=QSpinBox()
        self.nnoiseSpinBox.setRange(1, len(self.mA))
        self.nsigSpinBox=QDoubleSpinBox()
        self.nsigSpinBox.setDecimals(3)
        self.nsigSpinBox.setRange(.001, 100.)
        self.nsigSpinBox.setValue(2.)
        self.naboveSpinBox=QSpinBox()
        self.naboveSpinBox.setRange(1, len(self.mA))
        
        self.ncycLabel=QLabel()
        self.ncycLabel.setText('')
        self.noverlapLabel=QLabel()
        self.noverlapLabel.setText('')
        
        plotindLabel=QLabel()
        plotindLabel.setText('cycles to plot:')
        plotindLabel.setAlignment(Qt.AlignRight)
        cyclestartLabel=QLabel()
        cyclestartLabel.setText('cycle start indeces:')
        cyclestartLabel.setAlignment(Qt.AlignRight)
        self.plotindLineEdit=QLineEdit()
        self.cyclestartLineEdit=QLineEdit()

        plotButton=QPushButton()
        plotButton.setText("plot cycles and \n'save' answer")
        QObject.connect(plotButton, SIGNAL("pressed()"), self.plotcycles)
        
        calcButton=QPushButton()
        calcButton.setText('determine cycles\nfrom trigger\nsettings')
        QObject.connect(calcButton, SIGNAL("pressed()"), self.calccycles)
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        mainlayout=QGridLayout()
        prow=10
        spinrow=4;spincol=2;slabcol=5;bcol=2;
        ansrow=6;
        
        mainlayout.addWidget(self.plotw, 0, 0, 10, prow)
        
        mainlayout.addWidget(fnLabel, prow, spincol, 1, slabcol)
        mainlayout.addWidget(durLabel, prow+1, spincol, 1, slabcol)
        mainlayout.addWidget(nnoiseLabel, prow+2, spincol, 1, slabcol)
        mainlayout.addWidget(nsigLabel, prow+3, spincol, 1, slabcol)
        mainlayout.addWidget(naboveLabel, prow+4, spincol, 1, slabcol)
        
        mainlayout.addWidget(self.durSpinBox, prow+1, 0, 1, spincol)
        mainlayout.addWidget(self.nnoiseSpinBox, prow+2, 0, 1, spincol)
        mainlayout.addWidget(self.nsigSpinBox, prow+3, 0, 1, spincol)
        mainlayout.addWidget(self.naboveSpinBox, prow+4, 0, 1, spincol)
        
        mainlayout.addWidget(calcButton, prow, spincol+slabcol, spinrow//2, bcol)
        mainlayout.addWidget(plotButton, prow+spinrow//2, spincol+slabcol, spinrow-spinrow//2, bcol)
        
        mainlayout.addWidget(self.ncycLabel, prow+spinrow+1, 0, 1, 5)
        mainlayout.addWidget(self.noverlapLabel, prow+spinrow+1, 5, 1, 5)
        
        mainlayout.addWidget(plotindLabel, prow+spinrow+2, 0, 1, 10)
        mainlayout.addWidget(cyclestartLabel, prow+spinrow+4, 0, 1, 10)
        mainlayout.addWidget(self.plotindLineEdit, prow+spinrow+3, 0, 1, 10)
        mainlayout.addWidget(self.cyclestartLineEdit, prow+spinrow+5, 0, 1, 10)
        
        mainlayout.addWidget(self.buttonBox, prow+spinrow+ansrow, 3, 1, 3)
        
        
        self.setLayout(mainlayout)
    
    def calccycles(self):
        self.dur=self.durSpinBox.value()
        self.nnoise=self.nnoiseSpinBox.value()
        self.nsig=self.nsigSpinBox.value()
        self.nabove=self.naboveSpinBox.value()
        
        critfracforlastcycle=.97
        ncyc=len(self.mA)//self.dur
        if ncyc%self.dur>self.dur*critfracforlastcycle:
            ncyc+=1
        self.cycind=numpy.uint32(range(ncyc))*self.dur
        self.triggind=[]
        self.overlaplength=self.dur
        self.preoverlaplength=self.dur
        for counter, (i0, i1) in enumerate(zip(self.cycind, numpy.append(self.cycind[1:], len(self.mA)))):
            noisevals=self.mA[i0:i0+self.nnoise]
            v0=noisevals.mean()
            v1=noisevals.std()*self.nsig
            abovebool=(self.mA[i0:i1]-v0)>v1
            if len(self.cycind)==1: #if only one cycle dont' search for trigger
                i=0
            else:
                trig=False
                for i in numpy.where(abovebool)[0]:
                    trig=(i+self.nabove)<len(abovebool) and numpy.all(abovebool[i:i+self.nabove])
                    if trig:
                        break
                if not trig:
                    QMessageBox.warning(self,"FAILED",  "ABORTED because no trigger found in cycle %d" %counter)
                    return
            self.triggind+=[i0+i]
            self.preoverlaplength=min(self.preoverlaplength, i)
            self.overlaplength=min(self.overlaplength, i1-(i0+i))
        self.ncycLabel.setText('%d cycles identified' %len(self.cycind))
        self.noverlapLabel.setText('cycle overlap is %d indeces' %(self.overlaplength+self.preoverlaplength))
        self.plotindLineEdit.setText(','.join(tuple([`temp` for temp in range(len(self.triggind))])))
        self.cyclestartLineEdit.setText(','.join(tuple([`temp` for temp in self.triggind])))

    def readuintlineedit(self, le):#designed only for nonnegativeintegers
        t=str(le.text())
        l=[]
        while len(t)>0:
            a, b, t=t.partition(',')
            l+=[a]
        l=[eval(v.strip()) for v in l if v.strip().isdigit()]
        return numpy.uint32(l)

    def partition_array_triggers(self, arr): #cyci is the indeces of arr at which each cycle begins and triggi is the indeces of the trigger in each cycle and these are used to produce a 2D array indexed by cycle num and containing the parts of arr that correspond to overlapped regions of the cycles. a list of the "excluded" segments is also included
        temp=min([i1-it for it, i1 in zip(self.triggind, numpy.append(self.cycind[1:], len(arr)))])
        temp2=min([it-i0 for it, i0 in zip(self.triggind, self.cycind)])
        if temp<1 or temp2<0:
            QMessageBox.warning(self,"FAILED",  "ABORTED because at least one trigger was not within the bounds of its cycle")
            return
        self.overlaplength=temp
        self.preoverlaplength=temp2
        self.noverlapLabel.setText('cycle overlap is %d indeces' %(self.overlaplength+self.preoverlaplength))
        temp=[[arr[i0:it-self.preoverlaplength], arr[it-self.preoverlaplength:it+self.overlaplength], arr[it+self.overlaplength:i1]] for i0, it, i1 in zip(self.cycind, self.triggind, numpy.append(self.cycind[1:], len(arr)))]
        cyclearr=numpy.array(map(operator.itemgetter(1),temp))
        cutsegs=map(operator.itemgetter(0, 2),temp) #a list the same length as cyclearr but each elemnt is tuple of the pre-cycle and post-cycle segments, which may be empty
        return cyclearr, cutsegs
        
    def plotcycles(self):
        self.plotind=self.readuintlineedit(self.plotindLineEdit)
        temp=self.readuintlineedit(self.cyclestartLineEdit)
        if len(temp)!=len(self.triggind):
            QMessageBox.warning(self,"FAILED",  "ABORTED because number of values in 'cycle start indeces' was changed")
            return
        self.triggind=temp
        temp=self.partition_array_triggers(self.mA)
        if temp is None:
            return
        cyclemA, cutsegsmA=temp
        
        self.plotw.axes.cla()
        maxval=numpy.max(cyclemA[self.plotind])
        minval=numpy.min(cyclemA[self.plotind])
        cols=['b', 'g', 'r', 'c', 'm', 'y', 'k']
        for counter, (arr, segs) in enumerate(zip(cyclemA, cutsegsmA)):
            if not counter in self.plotind:
                continue
            c=cols[counter%len(cols)]
            for counter2, s in enumerate(segs):
                if len(s)>0:
                    segi=numpy.int32(range(len(s)))
                    if counter2==0:
                        segi=(1+segi)*-1
                    else:
                        segi+=len(arr)
                    self.plotw.axes.plot(segi, s, c+'x', ms=self.markersize)
                    maxval=max([maxval, max(s)])
                    minval=min([minval, min(s)])
            self.plotw.axes.plot(range(len(arr)), arr, c+'.', ms=self.markersize)
        self.plotw.axes.plot([.5, .5], [minval, maxval], 'k-')
        self.plotw.axes.plot([len(arr)-.5, len(arr)-.5], [minval, maxval], 'k-')
        self.plotw.axes.set_ylim(minval, maxval)
        self.plotw.fig.canvas.draw()


    def ExitRoutine(self):
        self.plotcycles()

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

class selectorDialog(QDialog):
    def __init__(self, parent, strlist, title=''):
        super(selectorDialog, self).__init__(parent)

        self.setWindowTitle(title)
        
        self.ComboBox=QComboBox()
        self.ComboBox.clear()
        count=0
        for count, s in enumerate(strlist):
            self.ComboBox.insertItem(count, s)
        self.ComboBox.setCurrentIndex(0)

        mainlayout=QGridLayout()
        mainlayout.addWidget(self.ComboBox, 0, 0)
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        mainlayout.addWidget(self.buttonBox, 1, 0)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        self.setLayout(mainlayout)
    def ExitRoutine(self):
        self.name=str(self.ComboBox.currentText())
        self.index=self.ComboBox.currentIndex()
        
class selectgroupDialog():
    def __init__(self, parent, h5grp, title=''):
        strlist=[]
        for grp in h5grp.iterobjects():
            if isinstance(grp, h5py.Group):
                strlist+=[(grp.name).rpartition('/')[2]]
        
        idialog=selectorDialog(parent, strlist, title=title)
        if not idialog.exec_():
            return
        self.grpname=idialog.name
        self.grp=h5grp[self.grpname]

class rescalDialog(QDialog):
    def __init__(self, parent, h5path):
        super(rescalDialog, self).__init__(parent)
        
        self.arrComboBoxlist=[QComboBox() for i in range(4)]
        self.dfltSpinBoxlist=[QDoubleSpinBox() for i in range(3)]
        self.useaveCheckBoxlist=[QCheckBox() for i in range(2)]
        for sb, n in zip(self.dfltSpinBoxlist, [2, 2, 5]):
            sb.setDecimals(n)
        for cb, s in zip(self.useaveCheckBoxlist, ['Ro,To', 'dR/dT']):
            cb.setText('Use measured ave of %s \nas value for other cells' %s)
        
        l0=QLabel()
        l0.setText('Use Ro,To from dataset or enter const value')
        l1=QLabel()
        l1.setText('Calc dR/dT from 2 datasets or enter const value')
        l2=QLabel()
        l2.setText('Use const Ro,To\nfor all cells:')
        l3=QLabel()
        l3.setText('Use const dR/dT\nfor all cells:')
        l4=QLabel()
        l4.setText('Ave Ro with:')
        
        mainlayout=QGridLayout()
        mainlayout.addWidget(l0, 0, 0, 4, 1)
        mainlayout.addWidget(l1, 4, 0, 4, 1)
        
        for w, n, n2 in zip(self.arrComboBoxlist, [0, 1, 4, 6], [1, 3, 1, 1]):
            mainlayout.addWidget(w, n, n2, 1, 1)
        mainlayout.addWidget(l2, 0, 2, 2, 1)
        mainlayout.addWidget(l3, 4, 2, 2, 1)
        mainlayout.addWidget(l4, 0, 3, 1, 1)
        for w, n in zip(self.dfltSpinBoxlist, [2, 3, 6]):
            mainlayout.addWidget(w, n, 2, 1, 1)
        for w, n in zip(self.useaveCheckBoxlist, [2, 6]):
            mainlayout.addWidget(w, n, 3, 2, 1)

        self.buttonBox = QDialogButtonBox(self)
        #self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.calcansave)
        QMetaObject.connectSlotsByName(self)
        mainlayout.addWidget(self.buttonBox, 4, 4, 4, 1)
        
        self.setLayout(mainlayout)
        
        self.h5path=h5path
        self.h5file=h5py.File(self.h5path, mode='r')
        h5cal=self.h5file['Calorimetry']

        self.arrpathlist=[]
        #self.imagenamelist=[]
        counter=0
        for grp in h5cal.iterobjects():
            if 'analysis' in grp and 'CellResistance' in grp['analysis']:
                #self.imagenamelist+=[(grp.name).rpartition('/')[2]]
                for cb in self.arrComboBoxlist:
                    cb.insertItem(counter, (grp.name).rpartition('/')[2])
                self.arrpathlist+=[grp['analysis/CellResistance'].name]
                counter+=1
        self.arrComboBoxlist[1].insertItem(counter, 'None')
        self.arrComboBoxlist[1].setCurrentIndex(counter)
        self.h5file.close()
    def calcansave(self):
        self.h5file=h5py.File(self.h5path, mode='r+')
        sbvals=[sb.value() for sb in self.dfltSpinBoxlist]
        Rdfltbool=sbvals[0]>0. or len(self.arrpathlist)==0
        Adfltbool=sbvals[2]>0. or len(self.arrpathlist)==0
        
        pnts=[(cb.currentIndex()<len(self.arrpathlist) and (self.h5file[self.arrpathlist[cb.currentIndex()]],) or (None,))[0] for cb in self.arrComboBoxlist if len(self.arrpathlist)>0]
        
        savearr=numpy.zeros((len(self.h5file.attrs['cells']), 3), dtype='float32')
        if Rdfltbool:
            savearr[:, 0]=sbvals[0]
            savearr[:, 1]=sbvals[1]
            names=[p.rpartition('/')[2] for p in experimentgrppaths(self.h5file)]
            idialog=selectorDialog(self, names, title='Select experiment to save calibration')
            if not idialog.exec_():
                self.h5file.close()
                return
            grp=getcalanalysis(self.h5file, names[idialog.index])
        else:
            grp=pnts[0].parent
            if pnts[1] is None:
                r=pnts[0][:]
            else:
                r0=pnts[0][:]
                r1=pnts[1][:]
                b0=numpy.float32(r0>0.)
                b1=numpy.float32(r1>0.)
                r=numpy.zeros(len(r0), dtype='float32')
                r[(b0+b1)>0]=((r0+r1)/(b0+b1))[(b0+b1)>0]
            savearr[:, 0]=r[:]
            savearr[:, 1]=pnts[0].attrs['ambient_tempC']
            if self.useaveCheckBoxlist[0].isChecked():
                inds=numpy.where(savearr[:, 0]<=0.)
                notinds=numpy.where(savearr[:, 0]>0.)
                savearr[inds, 0]=numpy.mean(savearr[notinds, 0])
                savearr[inds, 1]=numpy.mean(savearr[notinds, 1])
        if Adfltbool:
            savearr[:, 2]=sbvals[2]
        else:
            r1=pnts[2][:]
            t1=pnts[2].attrs['ambient_tempC']
            r2=pnts[3][:]
            t2=pnts[3].attrs['ambient_tempC']
            
            inds=numpy.where((r1>0.) & (r2>0.) & (t2!=t1))
            savearr[inds, 2]=tcr(r2[inds], r1[inds], t2[inds], t1[inds])
            
            if self.useaveCheckBoxlist[1].isChecked():
                inds=numpy.where(savearr[:, 2]<=0.)
                notinds=numpy.where(savearr[:, 2]>0.)
                savearr[inds, 2]=numpy.mean(savearr[notinds, 2])
        
        if 'Res_TempCal' in grp:
            del grp['Res_TempCal']
        ds=grp.create_dataset('Res_TempCal', data=savearr)
        ds.attrs['doc']='numpts x 3 array where the 3 are R0, T0 and alpha'
        if Rdfltbool:
            ds.attrs['Ro']='user-entered'
        else:
            ds.attrs['Ro']=pnts[0].name
            if not pnts[1] is None:
                ds.attrs['Rovaluesaveragedwith']=pnts[1].name
        self.h5file.close()


def fitviewer(parent, hpsdl, fitdlist, filterdict):
    strlist=['%s, seg %d' %(d['dsname'], d['seg']) for d in fitdlist]
    idialog=selectorDialog(parent, strlist, title='select from the available fit data')
    if not idialog.exec_():
        return
    fitd=fitdlist[idialog.index]
    
    strlist=['cycletime']+fitd['segdkeys']
    idialog=selectorDialog(parent, strlist, title='select an x-coordinate')
    if not idialog.exec_():
        return
    xdata=hpsdl[fitd['seg']][strlist[idialog.index]]
    ydata=hpsdl[fitd['seg']][fitd['dsname']]
    
    strlist=['None']+[filkey for filkey, fild in filterdict.iteritems() if not (False in [k in fild.keys() for k in ['interp_kind', 'extrap_order', 'interppars']])]
    idialog=selectorDialog(parent, strlist, title='select a filter for fit fcn evaluate')
    if not idialog.exec_():
        return
    if idialog.index==0:
        yfitdata=evaluatefitfcn(fitd, hpsdl[fitd['seg']])
    else:
        fild=filterdict[strlist[idialog.index]]
        yfitdata=evaluatefitfcn(fitd, hpsdl[fitd['seg']], interp_kind=fild['interp_kind'], extrap_order=fild['extrap_order'], interppars=fild['interppars'])
    #return yfitdata
    idialog=fitplotDialog(parent, xdata, ydata, yfitdata, hpsdl[fitd['seg']]['cyclepartition'])
    idialog.show()

class fitplotDialog(QDialog):
    def __init__(self, parent, xdata, ydata, yfitdata, cyclepartition):#if xdata not provided just plots data, if provided must be same format as data. if data is 2-d or list or arrays, iterates over 1st d and plots vs second
        super(fitplotDialog, self).__init__(parent)

        self.xdata=xdata
        self.ydata=ydata
        self.yfitdata=yfitdata
        self.cyclepartition=cyclepartition
        
        mainlayout=QGridLayout()
        self.cycleComboBox=QComboBox()
        cycLabel=QLabel()
        cycLabel.setText('select cycle(s)\nfor calculation')

        self.cycleComboBox.clear()
        for i in range(xdata.shape[0]):
            self.cycleComboBox.insertItem(i, `i`)
        self.cycleComboBox.setCurrentIndex(0)
        
        QObject.connect(self.cycleComboBox,SIGNAL("activated(QString)"),self.plot)
        
        self.plotw=plotwidget(self)
        mainlayout.addWidget(cycLabel, 0, 0)
        mainlayout.addWidget(self.cycleComboBox, 0, 1)
        mainlayout.addWidget(self.plotw, 1, 0, 5, 3)
        self.setLayout(mainlayout)
        self.plot()
    
    def plot(self):
        self.plotw.axes.cla()
        i=self.cycleComboBox.currentIndex()
        (x, y, yf, tp) =(self.xdata[i], self.ydata[i], self.yfitdata[i], self.cyclepartition[i])
        if numpy.any(tp>=0):
            self.plotw.axes.plot(x[tp>=0], y[tp>=0], '.', markersize=1, color='b', label='data')
        if numpy.any(tp<0):
            self.plotw.axes.plot(x[tp<0], y[tp<0], '.', markersize=1, color='k', alpha=.4, label='excluded')
        self.plotw.axes.plot(x, yf, 'r-', lw=1, label='fit')
        self.plotw.axes.legend(loc=2)
        self.plotw.fig.canvas.draw()
        
class simpleplotDialog(QDialog):
    def __init__(self, parent, data, xdata=None, style='.'):#if xdata not provided just plots data, if provided must be same format as data. if data is 2-d or list or arrays, iterates over 1st d and plots vs second
        super(simpleplotDialog, self).__init__(parent)
        if isinstance(data, numpy.ndarray) and data.ndim==1:
            data=[data]
            xdata=((xdata is None) and (None,) or ([xdata],))[0]
#        elif data.ndim==2:
#            if data.shape[0]>data.shape[1]:
#                data=data.T
#            arrlist=[dt for dt in data]
        mainlayout=QGridLayout()
        
        self.plotw=plotwidget(self)
        mainlayout.addWidget(self.plotw, 0, 0, 1, 2)
        if isinstance(style, str):
            style=[style]*len(data)
        for count, (arr, sty) in enumerate(zip(data, style)):
            if xdata is None:
                self.plotw.axes.plot(arr, sty, markersize=1)
            else:
                self.plotw.axes.plot(xdata[count], arr, sty, markersize=1)
        self.plotw.axes.set_xlabel('array index')
        self.plotw.axes.set_ylabel('array value')

        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        mainlayout.addWidget(self.buttonBox, 1, 0)
        self.setLayout(mainlayout)

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
            if isinstance(v,str) and 'array' in v:
                try:
                    v=eval(v.replace('array', 'numpy.array'))
                except:
                    print 'attribute %s with value %s was not converted to an array' %(k, v)
            try:
                h5file[path].attrs[k]=v
            except:
                print 'error in saveing key "', k, '" with value ', v
        h5file.close()

def get25pylabaxes(horizslowaxis=True, oneontop=True, oneonleft=True, figsize=(9., 9.)):

    
    pylab.figure(figsize=figsize)
    fig=pylab.gcf()

#    xstart=list(numpy.linspace(0, 1., 5))*5
#    xwidth=[.2]*25
#    borderpadx1=[.02]*25
#    borderpadx2=[.02]*25
#
#    ystart=[.9-i*.2 for i in range(1, 6) for j in range(5)]
#    ywidth=[.2]*25
#    borderpady1=[.02]*25
#    borderpady2=[.02]*25
#
#    axl=[]
#    for x, xw, bx1, bx2, y, yw, by1, by2 in zip(xstart, xwidth, borderpadx1, borderpadx2, ystart, ywidth, borderpady1, borderpady2):
#        axl+=[fig.add_axes([x+bx1, y+by1, xw-bx1-bx2, yw-by1-by2])]
    axl=[pylab.subplot(5, 5, i*5+j+1) for i in range(5) for j in range(5)]
    
    
    if horizslowaxis:
        if oneontop:
            if oneonleft:
                inds=range(25)
            else:
                inds=[i*5+j for i in range(5) for j in range(4,-1,-1)]
        else:
            if oneonleft:
                inds=[i*5+j for i in range(4,-1,-1) for j in range(5)]
            else:
                inds=range(25)[::-1]
    else:
        if oneontop:
            if oneonleft:
                inds=[i+j*5 for i in range(5) for j in range(5)]
            else:
                inds=[i+j*5 for i in range(4,-1,-1) for j in range(5)]
        else:
            if oneonleft:
                inds=[i+j*5 for i in range(5) for j in range(4,-1,-1)]
            else:
                inds=[i+j*5 for i in range(4, -1, -1) for j in range(4,-1,-1)]
    return [axl[i] for i in inds]

def plotsegmentsongrid(h5path, h5expname, axl, xkey, ykey, seginds, cycinds=[0], cellnums=range(1, 26), plotstyle='b-', samexlim=True):
    h5file, nams=experimenthppaths(h5path, h5expname)
    hpcells=[h5file[nam].attrs['CELLNUMBER'] for nam in nams]
    print hpcells
    cellsh5hpnames=[[cell, nams[hpcells.index(cell)].rpartition('/')[2]] for cell in cellnums if cell in hpcells]
    h5file.close()
    xmin=None
    for cell, h5hpname in cellsh5hpnames:
        print 'plotting cell %d' %cell
        hpsdl=CreateHeatProgSegDictList(h5path, h5expname, h5hpname)
        if isinstance(plotstyle, list):
            ps=plotstyle[cellnums.index(cell)]
        else:
            ps=plotstyle
        if not samexlim:
            xmin=None
        for ci in cycinds:
            xvals=numpy.concatenate([hpsdl[i][xkey][ci] for i in seginds])
            yvals=numpy.concatenate([hpsdl[i][ykey][ci] for i in seginds])
            axl[cell-1].plot(xvals, yvals, ps)
            if xmin is None:
                xmin=xvals.min()
                xmax=xvals.max()
                ymin=yvals.min()
                ymax=yvals.max()
            else:
                xmin=min(xmin, xvals.min())
                xmax=max(xmax, xvals.max())
                ymin=min(ymin, yvals.min())
                ymax=max(ymax, yvals.max())
        if not samexlim:
            axl[cell-1].set_xlim(xmin, xmax)
            axl[cell-1].set_ylim(ymin, ymax)
        print xvals.min(), xvals.max(), yvals.min(), yvals.max(), numpy.isnan(yvals).sum()
    if samexlim:
        for cell, h5hpname in cellsh5hpnames:
            axl[cell-1].set_xlim(xmin, xmax)
            axl[cell-1].set_ylim(ymin, ymax)
    
def TEMP():
    h5path='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/2010Nov27_AuSiCu_pnsc.h5'
    h5expname='heat1d'
    axl=get25pylabaxes()
    #axl=[pylab.subplot(111)]
    print '****'
    cellnums=sorted(list(set(range(1, 26))-set([3, 15, 21])))
    #cellnums=[1]
    plotsegmentsongrid(h5path, h5expname, axl, 'sampletemperature', 'samplepowerperrate', [2], cellnums=cellnums)
    for ax in axl:
        ax.set_yticks([])
    pylab.show()
    
    print 'done'
    

#class fitpowerperheatrateDialog(QDialog):
#    def __init__(self, parent, segdlist):
#        super(fitpowerperheatrateDialog, self).__init__(parent)
#        self.segdlist=segdlist
#        mainlayout=QGridLayout()
#        
#        self.plotw=plotwidget(self)
#        QObject.connect(self.plotw, SIGNAL("genericclickonplot"), self.clickprocess)
#
#        self.segComboBox=QComboBox()
#        segLabel=QLabel()
#        segLabel.setText('select segment\nfor calculation')
#        self.cycleComboBox=QComboBox()
#        cycLabel=QLabel()
#        cycLabel.setText('select cycle(s)\nfor calculation')
#        self.keyComboBox=QComboBox()
#        keyLabel=QLabel()
#        keyLabel.setText('select dataset\nfor calculation')
#        self.fitfcnComboBox=QComboBox()
#        fitcbLabel=QLabel()
#        fitcbLabel.setText('select function\nfor fitting')
#        self.fitfcnLabel=QLabel()
#        
#        self.segcalcoptions=[]
#        for count, d in enumerate(self.hpsegdlist):
#            if d['segmenttype'] in ['ramp', 'soak']:
#                self.segComboBox.insertItem(len(self.segcalcoptions), 'segment %d : %s' %(count, d['segmenttype']))
#                self.segcalcoptions+=[count]
#        
#        QObject.connect(self.segComboBox,SIGNAL("activated(QString)"),self.segcgchanged)
#        QObject.connect(self.cycleComboBox,SIGNAL("activated(QString)"),self.plot)
#        QObject.connect(self.keyComboBox,SIGNAL("activated(QString)"),self.plot)
#        QObject.connect(self.fitfcnComboBox,SIGNAL("activated(QString)"),self.fitfcnchanged)
#        self.segcgchanged()
#        
#        self.clickforexcludeCheckBox=QCheckBox()
#        self.clickforexcludeCheckBox.setText('click to define excluded\nregions (x-coord only)\nclick low then high')
#        self.lastclickx=None
#        
#        clearremoveindsButton=QPushButton()
#        QObject.connect(clearremoveindsButton, SIGNAL("pressed()"), clearremoveindslist)
#        self.clearremoveindslist()
#        
#        mainlayout.addWidget(delButton, 4, pwid, 1, 1)
#        mainlayout.addWidget(addButton, 5, pwid, 1, 1)
#        mainlayout.addWidget(self.rowComboBox, 4, pwid+1, 2, 1)
#
#
#        mainlayout.addWidget(plotw, 0, 0, 1, 2)
#        
#        self.buttonBox = QDialogButtonBox(self)
#        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
#        self.buttonBox.setOrientation(Qt.Horizontal)
#        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
#        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
#        mainlayout.addWidget(self.buttonBox, 1, 0)
#        self.setLayout(mainlayout)
#        
#    
#    def segcgchanged(self, kdflt='samplepowerperheatrate'):
#        self.segind=self.segcalcoptions[self.segComboBox.currentIndex()]
#        d=self.segdlist[self.segind]
#        self.cycletime0=d['cycletime'][0]
#        self.cycinds=numpy.array(range(d['cycletime'].shape[1]))
#        self.allcycinds=numpy.array(range(d['cycletime'].shape[1]))
#        self.cycleComboBox.clear()
#        self.cycleComboBox.insertItem(0, 'all')
#        for i in range(d['cycletime'].shape[0]):
#            self.cycleComboBox.insertItem(i+1, `i`)
#        self.cycleComboBox.setCurrentIndex(0)
#        
#        self.keyComboBox.clear()
#        temp=0
#        i=0
#        for i in range(d['cycletime'].shape[0]):
#            if v.shape==d['cycletime'].shape:
#                self.keyComboBox.insertItem(i, k)
#                if k==kdflt:
#                    temp=i
#        self.keyComboBox.setCurrentIndex(temp)
#        
#        self.heatlossfcnsnames=[(k, fcnlabel, pard) for k, (fcn, reqdkeys, fcnlabel) in HeatLossFunctionLibrary.iteritems() if numpy.all([segk in d.keys() for segk in reqdkeys])]
#        self.fitfcnComboBox.clear()
#        for i, (nam, l) in enumerate(self.heatlossfcnsnames):
#            self.fitfcnComboBox.insertItem(i, nam)
#    
#    def fitfcnchanged(self):
#        fcnlabel=self.heatlossfcnsnames[self.fitfcnComboBox.currentIndex()][1]
#        self.fitfcnLabel.setText(fcnlabel)
#        
#        
#    def clearremoveindslist(self):
#        self.removeindslist=[]
#    def clickprocess(self, coords):
#        if self.clickforexcludeCheckBox.isChecked():
#            if self.lastclickx is None:
#                self.lastclickx=coords[0]
#            else:
#                self.removeindslist+=[(min(coords[0], self.lastclickx), max(coords[0], self.lastclickx))]
#                self.lastclickx=None
#                self.plot()
#        else:
#            self.lastclickx=None
#    
#    def removeinds(self):
#        self.cycinds=set(self.cycinds)
#        for t0, t1 in self.removeindslist:
#            self.cycinds-=set(numpy.where((self.cycletime0>t0)&(self.cycletime1<t1))[0])
#        self.cycinds=numpy.array(sorted(list(self.cycinds)))
#        
#    def plot(self):
#        self.key=str(self.keyComboBox.currentText())
#        d=self.segdlist[self.segind]
#        self.cycles=[self.cycleComboBox.currentIndex()]
#        if self.cycles[0]<0:
#            self.cycles=range(d['cycletime'].shape[0])
#        self.removeinds()
#        self.plotw.axes.cla()
#        nofitinds=numpy.array(sorted(set(self.allcycinds)-set(self.cycinds)))
#        for cyc in self.cycles:
#            if len(self.cycinds)>0:
#                self.plotw.axes.plot(self.cycletime0[self.cycinds], d[self.key][cyc][self.cycinds], 'k.', markersize=1)
#            if len(nofitinds)>0:
#                self.plotw.axes.plot(self.cycletime0[nofitinds], d[self.key][cyc][nofitinds], '.', markersize=1, color=(.6, .6, .6))

class timepartDialog(QDialog):
    def __init__(self, parent, cycletime, yvals=None, numpieces=1):
        if yvals is None:
            self.yvals=cycletime
        else:
            self.yvals=yvals
        super(timepartDialog, self).__init__(parent)
        self.cycletime=cycletime
        self.timepart=numpy.zeros(cycletime.shape, dtype='int32')
        self.numpieces=numpieces
        mainlayout=QGridLayout()

        self.piecebndrylist=[[] for i in range(self.numpieces-1)]
        self.removeindslist=[[] for i in range(self.cycletime.shape[0])]
        
        self.plotw=plotwidget(self, width=16, height=16)
        QObject.connect(self.plotw, SIGNAL("genericclickonplot"), self.clickprocess)

        self.cycleComboBox=QComboBox()
        cycLabel=QLabel()
        cycLabel.setText('select cycle(s)\nfor calculation')

        self.cycleComboBox.clear()
        self.cycleComboBox.insertItem(0, 'all')
        for i in range(self.cycletime.shape[0]):
            self.cycleComboBox.insertItem(i+1, `i`)
        self.cycleComboBox.setCurrentIndex(0)
        
        QObject.connect(self.cycleComboBox,SIGNAL("activated(QString)"),self.plot)

        self.clickforexcludeCheckBox=QCheckBox()
        self.clickforexcludeCheckBox.setText('click to define excluded\nregions (x-coord only)\nclick low then high')
        self.clickforexcludeCheckBox.setChecked(True)
        self.lastclickx=None
        
        clearremoveindsButton=QPushButton()
        clearremoveindsButton.setText('clear the list\nof excluded regions')
        QObject.connect(clearremoveindsButton, SIGNAL("pressed()"), self.clearremoveindslist)
        self.clearremoveindslist()
  
        autopiecebndryButton=QPushButton()
        autopiecebndryButton.setText('Use exclude regions to\ndefine piece-wise fcn boundaries')
        QObject.connect(autopiecebndryButton, SIGNAL("pressed()"), self.autopiecebndry)
        self.clearremoveindslist()
        
        self.clickforpieceCheckBox=QCheckBox()
        self.clickforpieceCheckBox.setText('click to define excluded\nregions (x-coord only)\nclick low then high')
        self.clickforpieceCheckBox.setChecked(False)

        clearpiecebndryButton=QPushButton()
        clearpiecebndryButton.setText('clear the list of\npiecewise boundaries')
        QObject.connect(clearpiecebndryButton, SIGNAL("pressed()"), self.clearpiecebndrylist)
        self.clearpiecebndrylist()
        
        mainlayout.addWidget(cycLabel,0,0)
        mainlayout.addWidget(self.cycleComboBox,0,1)
        mainlayout.addWidget(self.clickforexcludeCheckBox, 1,0, 1, 2)
        mainlayout.addWidget(clearremoveindsButton, 2,0, 1, 2)
        mainlayout.addWidget(autopiecebndryButton, 3,0, 1, 2)
        mainlayout.addWidget(self.clickforpieceCheckBox, 4,0, 1, 2)
        mainlayout.addWidget(clearpiecebndryButton, 5, 0, 1, 2)

        mainlayout.addWidget(self.plotw, 0, 2, 6, 6)
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        mainlayout.addWidget(self.buttonBox, 6, 2, 1, 6)
        self.setLayout(mainlayout)

    
    def clearremoveindslist(self):
        i=self.cycleComboBox.currentIndex()
        if i==0:
            self.removeindslist=[[] for i in range(self.cycletime.shape[0])]
        else:
            self.removeindslist[i-1]=[]
        self.plot()
    
    def clearpiecebndrylist(self):
        i=self.cycleComboBox.currentIndex()
        if i==0:
            self.piecebndrylist=[[] for i in range(self.numpieces-1)]
        else:
            self.piecebndrylist[i-1]=[]
        self.plot()
            
    def clickprocess(self, coords):
        if self.clickforexcludeCheckBox.isChecked():
            if self.lastclickx is None:
                self.lastclickx=coords[0]
            else:
                i=self.cycleComboBox.currentIndex()
                if i==0:
                    for l in self.removeindslist:
                        l+=[(min(coords[0], self.lastclickx), max(coords[0], self.lastclickx))]
                else:
                    self.removeindslist[i-1]+=[(min(coords[0], self.lastclickx), max(coords[0], self.lastclickx))]
                self.lastclickx=None
                self.plot()
        else:
            self.lastclickx=None
            
        if self.clickforpieceCheckBox.isChecked():
            i=self.cycleComboBox.currentIndex()
            if i==0:
                for l in self.piecebndrylist:
                    if len(l)<self.numpieces:
                        l+=[coords[0]]
                    else:
                        print 'replacing last piece boundary because enough boundaries have already been defined'
                        l[-1]=coords[0]
            else:
                if len(self.piecebndrylist[i-1])<self.numpieces:
                    self.piecebndrylist[i-1]+=[coords[0]]
                else:
                    print 'replacing last piece boundary because enough boundaries have already been defined'
                    self.piecebndrylist[i-1][-1]=coords[0]
            self.plot()
            
    def autopiecebndry(self):
        cycles=[self.cycleComboBox.currentIndex()-1]
        if cycles[0]<0:
            cycles=range(self.cycletime.shape[0])
        for cyc in cycles:
            self.piecebndrylist[cyc]=[]
            for i, rl in zip(range(self.numpieces-1), self.removeindslist[cyc]):#if there are more parititons than piece bndrys then use as many as necessary
                self.piecebndrylist[cyc]+=[(rl[0]+rl[1])/2.]#put the bndry in the middle of the exlcude region
        self.plot()
    def calctimepart(self):
        self.timepart=numpy.zeros(self.cycletime.shape, dtype='int32')
        for tp, time, rem in zip(self.timepart, self.cycletime, self.removeindslist):
            for (t0, t1) in rem:
                tp[(time>t0)&(time<t1)]=-1
        for tp, time, bdry in zip(self.timepart, self.cycletime, self.piecebndrylist):
            bdry=numpy.sort(bdry)
            for i, t0 in enumerate(bdry):
                tp[(time>=t0)&(tp>=0)]=i+1
    def plot(self):
        colors=['k']+['r', 'g', 'c', 'm', 'y', 'b']*5
        cycles=[self.cycleComboBox.currentIndex()-1]
        if cycles[0]<0:
            cycles=range(self.cycletime.shape[0])
        self.calctimepart()
        self.plotw.axes.cla()
        for cyc in cycles:
            tp=self.timepart[cyc]
            for i in range(-1, self.numpieces):
                if numpy.any(tp==i):
                    self.plotw.axes.plot(self.cycletime[cyc][tp==i], self.yvals[cyc][tp==i], '.', markersize=1, color=colors[i+1])
        self.plotw.fig.canvas.draw()


class peakfiteditDialog(QDialog):
    def __init__(self, parent, X, Y, fitpks, fitfcn, maxpts=9, markersize=1):#only works on one cycle, ie X and Y are 2-d
        super(peakfiteditDialog, self).__init__(parent)

        self.X=X
        self.Y=Y
        self.fitfcn=fitfcn
        self.fitpks=copy.copy(fitpks)
        self.sigs=None
        self.resids=None
        self.setWindowTitle('Provide the endpoint of segments in a single cycle')
        
        self.plotw=plotwidget(self)
        self.markersize=markersize
        self.plotw.axes.plot(X, Y, 'k.', ms=self.markersize)
        QObject.connect(self.plotw, SIGNAL("genericclickonplot"), self.clickprocess)            

#        self.plotw.axes.set_xlabel('cycle time (ms)')
#        self.plotw.axes.set_ylabel('applied current (mA)')
            
        mainlayout=QGridLayout()

        pwid=4
        nrows_ctrls=5
        mainlayout.addWidget(self.plotw, 0, 0, nrows_ctrls+maxpts, pwid)
        
        fitButton=QPushButton()
        fitButton.setText("fit for all parameters")
        QObject.connect(fitButton, SIGNAL("pressed()"), self.fit)
        

        delButton=QPushButton()
        delButton.setText("del row")
        QObject.connect(delButton, SIGNAL("pressed()"), self.delrow)

        addButton=QPushButton()
        addButton.setText("add row")
        QObject.connect(addButton, SIGNAL("pressed()"), self.addrow)

        self.rowComboBox=QComboBox()
        self.rowComboBox.clear()
        
        oflabel=QLabel()
        oflabel.setText('scipy.optimize.fcn:')
        self.optimizefcnLineEdit=QLineEdit()
        
        self.clickforpeakCheckBox=QCheckBox()
        self.clickforpeakCheckBox.setText('click to add peak')
        self.clickforpeakCheckBox.setChecked(False)
        
        mainlayout.addWidget(oflabel, 0, pwid, 1, 2)
        mainlayout.addWidget(self.optimizefcnLineEdit, 0, pwid+2, 1, 1)
        mainlayout.addWidget(self.clickforpeakCheckBox, 1, pwid, 1, 3)
        
        mainlayout.addWidget(fitButton, 2, pwid, 1, 3)
        mainlayout.addWidget(delButton, 3, pwid, 1, 2)
        mainlayout.addWidget(addButton, 4, pwid, 1, 2)
        mainlayout.addWidget(self.rowComboBox, 3, pwid+2, 2, 1)
        
        Label0=QLabel()
        Label0.setText('posn')
        Label1=QLabel()
        Label1.setText('halfw')
        Label2=QLabel()
        Label2.setText('height')
        mainlayout.addWidget(Label0, nrows_ctrls, pwid, 1, 1)
        mainlayout.addWidget(Label1, nrows_ctrls, pwid+1, 1, 1)
        mainlayout.addWidget(Label2, nrows_ctrls, pwid+2, 1, 1)
        
        self.range0=(X.min(), X.max())
        self.range1=(0., X.max()/2.)
        self.range2=(min(0., Y.min()), max(0., Y.max()))
        self.LineEditList0=[]
        self.LineEditList1=[]
        self.LineEditList2=[]
        for i in range(maxpts):
            self.rowComboBox.insertItem(i, `i`)
            sb0=QLineEdit()
            sb1=QLineEdit()
            sb2=QLineEdit()
#            sb0.setRange(X.min(), X.max())
#            sb0.setRange(0., X.max()/2.)
#            sb2.setRange(min(0., Y.min()), max(0., Y.max()))
#            tsb.setDecimals(3)
#            csb.setDecimals(3)
            self.filllineedits()
            mainlayout.addWidget(sb0, i+1+nrows_ctrls, pwid, 1, 1)
            mainlayout.addWidget(sb1, i+1+nrows_ctrls, pwid+1, 1, 1)
            mainlayout.addWidget(sb2, i+1+nrows_ctrls, pwid+2, 1, 1)
#            QObject.connect(tsb, SIGNAL("valueChanged(double)"), self.plotsegments)
#            QObject.connect(csb, SIGNAL("valueChanged(double)"), self.plotsegments)
            self.LineEditList0+=[sb0]
            self.LineEditList1+=[sb1]
            self.LineEditList2+=[sb2]
        
        self.rowComboBox.setCurrentIndex(0)
        
        plotButton=QPushButton()
        plotButton.setText('plot segments')
        QObject.connect(plotButton, SIGNAL("pressed()"), self.plotall)
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        mainlayout.addWidget(plotButton, nrows_ctrls+maxpts+1, pwid, 1, 3)
        mainlayout.addWidget(self.buttonBox, nrows_ctrls+maxpts+1, 0, 1, pwid)
        
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        self.setLayout(mainlayout)
    
    def clickprocess(self, coords):
        if self.clickforpeakCheckBox.isChecked():
            self.readdata()
            j=len(self.fitpks)
            if j<len(self.LineEditList0):
                self.LineEditList0[j].setText('%.3e' %coords[0])
                self.LineEditList1[j].setText('10')
                self.LineEditList2[j].setText('%.3e' %coords[1])
    def filllineedits(self):
        for i in range(len(self.LineEditList0)):
            if i<len(self.fitpks): 
                self.LineEditList0[i].setText('%.3e' %self.fitpks[i][0])
                self.LineEditList1[i].setText('%.3e' %self.fitpks[i][1])
                self.LineEditList2[i].setText('%.3e' %self.fitpks[i][2])
            else:
                self.LineEditList0[i].setText('')
                self.LineEditList1[i].setText('')
                self.LineEditList2[i].setText('')

    def addrow(self):
        ind=self.rowComboBox.currentIndex()
        for i in range(ind+1, len(self.LineEditList0))[::-1]:
            self.LineEditList0[i].setText(self.LineEditList0[i-1].text())
            self.LineEditList1[i].setText(self.LineEditList1[i-1].text())
            self.LineEditList2[i].setText(self.LineEditList2[i-1].text())

    def delrow(self):
        ind=self.rowComboBox.currentIndex()
        if ind==(len(self.LineEditList0)-1):
            return
        for i in range(ind, len(self.LineEditList0)-1):
            self.LineEditList0[i].setText(self.LineEditList0[i+1].text())
            self.LineEditList1[i].setText(self.LineEditList1[i+1].text())
            self.LineEditList2[i].setText(self.LineEditList2[i+1].text())
        self.LineEditList0[-1].setText('')
        self.LineEditList1[-1].setText('')
        self.LineEditList2[-1].setText('')

    
    def fit(self):#the analysis is done with normalization of the array so it is insensitive to 'Aunit'. Also, the data time interval is not used in the derivative calculation
        self.readdata()
        try:
            s=str(self.optimizefcnLineEdit.text())
            optimizerfcn=eval(s)
        except:
            if s!='':
                print 'optimize function ', s, ' not understood'
            optimizerfcn=None
        self.fitpks, self.sigs, self.resid=fitpeakset(self.X, self.Y, self.fitpks, self.fitfcn, optimizerfcn=optimizerfcn)
        self.filllineedits()
        
    def plotall(self):
        self.readdata()
        self.plotw.axes.cla()
        self.plotw.axes.plot(self.X, self.Y, 'k.', ms=self.markersize)
        fitYarr=numpy.float32([self.fitfcn(p, self.X) for p in self.fitpks])
        cols=['k', 'y', 'b', 'm', 'g', 'c', 'r']*10#just to be sure there's enough
        for count, arr in enumerate(fitYarr):
            self.plotw.axes.plot(self.X, arr, cols[i]+'--', lw=self.markersize)
        fitYall=fitYarr.sum(axis=0)
        self.plotw.axes.plot(self.X, fitYall, 'r', lw=self.markersize)
        self.resid=(((fitYall-self.Y)**2)**.5).sum()
        self.plotw.fig.canvas.draw()

    def readdata(self):
        j=0
        allvals=[]
        for i in range(len(self.LineEditList0)):
            strs=[str(le.text()) for le in [self.LineEditList0[j], self.LineEditList1[j], self.LineEditList2[j]]]
            if numpy.any([len(s)==0 for s in strs]):
                continue
            try:
                vals=[myeval(s) for s in strs]
                j+=1
            except:
                self.rowComboBox.setCurrentIndex(j)
                self.delrow()
            vals=[min(h, max(v, l)) for v, (l, h) in zip(vals, [self.range0, self.range1, self.range2])]
            allvals+=[vals]
        temp=numpy.float32(allvals)
        if not (len(temp)==len(self.fitpks) and numpy.all(temp==self.fitpks)):
            sigs=numpy.zeros(temp.shape, dtype='float32')
        self.fitpks=temp
        
    def ExitRoutine(self):
        #self.readdata()
        self.plotall()
        
    def savefig(self, p, dpi=300):
        self.plotw.fig.savefig(p, dpi=dpi)
            
class rescal_ExtraptoToDialog(QDialog):
    def __init__(self, parent, pardict, Rinit=0.):
        super(rescal_ExtraptoToDialog, self).__init__(parent)
        self.Rinit=Rinit
        self.pardict=pardict
        self.hpsdl=CreateHeatProgSegDictList(self.pardict['h5path'], self.pardict['h5expname'], self.pardict['h5hpname'])
        segtypelist=[d['segmenttype'] for d in self.hpsdl]
        ncyc=self.hpsdl[0]['cycletime'].shape[0]
        segi=2
        if segtypelist.count('soak')>=1:
            segi=segtypelist.index('soak')
            nprev=0
            for i in range(1, segi):
                if segtypelist[segi-i]=='zero':
                    break
                nprev+=1
            dsoak=self.hpsdl[segi]
        self.key_dflt_lab_min_max=[\
        ('segind', segi, 'cycle to use for extrapolation', 0, len(self.hpsdl)), \
        ('o_R2poly', 1, 'R(t) fit polynomial order', 0, 4), \
        ('iterations', 1, 'number iterations', 1, 10), \
        ('cycind', 0, 'cycle to use for extrapolation', 0, ncyc), \
        ('nprevsegs', nprev, 'number of previous cycles to extrap over', 0, 100), \
        (None, 0, 'starting index of extrap cycle', 0, 10000), \
        (None, 1000, 'stoping index of extrap cycle', 0, 10000), \
        ]
        
        mainlayout=QGridLayout()
        
        nrows=len(self.key_dflt_lab_min_max)+2
        self.spinboxlist=[]
        for i, (k, v, l, mi, ma) in enumerate(self.key_dflt_lab_min_max):
            lb=QLabel()
            lb.setText(l)
            sb=QSpinBox()
            sb.setRange(mi, ma)
            sb.setValue(v)
            self.spinboxlist+=[sb]
            mainlayout.addWidget(lb, i, 0)
            mainlayout.addWidget(sb, i, 1)
        
        calcButton=QPushButton()
        calcButton.setText("Calc Ro")
        QObject.connect(calcButton, SIGNAL("pressed()"), self.calc)
        
        limsButton=QPushButton()
        limsButton.setText("update limits")
        QObject.connect(limsButton, SIGNAL("pressed()"), self.updatelims)
        
        mainlayout.addWidget(limsButton, nrows-2, 0)
        mainlayout.addWidget(calcButton, nrows-2, 1)
        
        self.resultLabel=QLabel()
        mainlayout.addWidget(self.resultLabel, 0, 2, 2, 4)
        self.plotw=plotwidget(self)
        mainlayout.addWidget(self.plotw, 2, 2, nrows-2, 4)

        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox, SIGNAL("rejected()"), self.reject)
        mainlayout.addWidget(self.buttonBox, nrows-1, 0)
         
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        self.setLayout(mainlayout)

        #self.lw.setEditTriggers(QAbstractItemView.AllEditTriggers)
        #QObject.connect(self.lw,SIGNAL("itemSelectionChanged()"),self.updateattrd)
        
        QMetaObject.connectSlotsByName(self)

    def updatelims(self):
        print 'update limits for indeces and previous num of segs to avoid error in calc - not implemented'
    
    def calc(self):
        self.ExitRoutine()
        self.Ro, self.calcd=calcRo_extraptoTo(**self.pardict)
        s='Calculated Ro is %.3f.' %self.Ro
        if self.Rinit>0:
            s+='The prev Ro values was %.3f.' %self.Rinit
        s+='\nThe %d points were used to fit R in the range %.3f to %.3f.' %(len(self.calcd['P2']), self.calcd['R2fit'][0], self.calcd['R2fit'][-1])
        s+='\ndelTemp was %.1f for the fit, %.1f in the extrap range of %d points' %(self.calcd['delT2'], self.calcd['delT1'], len(self.calcd['P1']))
        self.resultLabel.setText(s)
    def ExitRoutine(self):
        nokeyvals=[]
        for i, (sb, (k, v, l, mi, ma)) in enumerate(zip(self.spinboxlist, self.key_dflt_lab_min_max)):
            if k is None:
                nokeyvals+=[sb.value()]
            else:
                self.pardict[k]=sb.value()
        self.pardict['inds_calcregion']=tuple(nokeyvals)

class TwoPointResTableDialog(QDialog):
    def __init__(self, parent, h5path):
        super(TwoPointResTableDialog, self).__init__(parent)
        self.h5path=h5path
        self.setWindowTitle('Select pinout and enter resistances')
        mainlayout=QGridLayout()
    
        namL=QLabel()
        namL.setText('save name')
        pinL=QLabel()
        pinL.setText('select pinout')
        
        self.namLE=QLineEdit()
        self.namLE.setText('PreCondition2PointRes')
        
        self.pinCB=QComboBox()
        self.pinCB.clear()
        for i, (k, v) in enumerate(PinoutLibrary.iteritems()):
            self.pinCB.insertItem(i, k)
        self.pinCB.setCurrentIndex(0)
        
            
        mainlayout.addWidget(namL, 0, 0, 1, 2)
        mainlayout.addWidget(self.namLE, 0, 2, 1, 3)
        
        mainlayout.addWidget(pinL, 1, 0, 1, 2)
        mainlayout.addWidget(self.pinCB, 1, 2, 1, 3)
        
        cellL=QLabel()
        cellL.setText('cell')
        l0=QLabel()
        l0.setText('I,I pads')
        l1=QLabel()
        l1.setText('I,I R(Ohms)')
        l2=QLabel()
        l2.setText('V,V pads')
        l3=QLabel()
        l3.setText('V,V R(Ohms)')
        
        mainlayout.addWidget(cellL, 2, 0)
        mainlayout.addWidget(l0, 2, 1)
        mainlayout.addWidget(l1, 2, 2)
        mainlayout.addWidget(l2, 2, 3)
        mainlayout.addWidget(l3, 2, 4)
        
        self.secdervalSpinBox=QDoubleSpinBox()
        self.secdervalSpinBox.setDecimals(1)
        self.secdervalSpinBox.setValue(-1)

    
        rowoff=3
        self.iLabelList=[]
        self.vLabelList=[]
        self.iSpinBoxList=[]
        self.vSpinBoxList=[]
        for i, j in zip(range(rowoff, rowoff+25), range(1, 26)):
            isb=QDoubleSpinBox()
            isb.setDecimals(1)
            isb.setRange(-1, 100000)
            isb.setValue(-1)
            vsb=QDoubleSpinBox()
            vsb.setDecimals(1)
            vsb.setRange(-1, 100000)
            vsb.setValue(-1)
            self.iSpinBoxList+=[isb]
            self.vSpinBoxList+=[vsb]
            cl=QLabel()
            cl.setText(`j`)
            il=QLabel()
            vl=QLabel()
            self.iLabelList+=[il]
            self.vLabelList+=[vl]
            mainlayout.addWidget(cl, i, 0)
            mainlayout.addWidget(il, i, 1)
            mainlayout.addWidget(isb, i, 2)
            mainlayout.addWidget(vl, i, 3)
            mainlayout.addWidget(vsb, i, 4)
        self.filllabels()
        
        QObject.connect(self.pinCB,SIGNAL("activated(QString)"),self.filllabels)
        
        readsavedButton=QPushButton()
        readsavedButton.setText("read saved\nresistances")
        QObject.connect(readsavedButton, SIGNAL("pressed()"), self.readsavedvalues)
        mainlayout.addWidget(readsavedButton, i+1, 0, 1, 2)
        
        self.buttonBox = QDialogButtonBox(self)
        self.buttonBox.setGeometry(QRect(520, 195, 160, 26))
        self.buttonBox.setOrientation(Qt.Horizontal)
        self.buttonBox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        QObject.connect(self.buttonBox, SIGNAL("accepted()"), self.accept)
        QObject.connect(self.buttonBox, SIGNAL("rejected()"), self.reject)
        mainlayout.addWidget(self.buttonBox, i+1, 2, 1, 3)
         
        QObject.connect(self.buttonBox,SIGNAL("accepted()"),self.ExitRoutine)
        
        self.setLayout(mainlayout)
        QMetaObject.connectSlotsByName(self)
        
    def ExitRoutine(self):
        arr=[]
        for isb, vsb in zip(self.iSpinBoxList, self.vSpinBoxList):
            arr+=[[isb.value(), vsb.value()]]
        self.resarr=numpy.float32(arr)
        savetwopointres(self.h5path, self.resarr, str(self.namLE.text()))
    def filllabels(self):
        k=str(self.pinCB.currentText())
        arr=PinoutLibrary[k]
        for il, vl, (ia, va) in zip(self.iLabelList,  self.vLabelList, arr):
            il.setText('%d,%d' %tuple(ia))
            vl.setText('%d,%d' %tuple(va))
    def readsavedvalues(self):
        nam_resarrlist=readtwopointres(self.h5path)
        if len(nam_resarrlist)==0:
            print 'no saved 2-point resistance arrays found'
            return
        elif len(nam_resarrlist)==1:
            nam, resarr=nam_resarrlist[0]
        else:
            naml=[nam for nam, resarr in nam_resarrlist]
            idialog=selectorDialog(self, naml, 'select saved 2point  res array')
            if not idialog.exec_():
                return
            nam=idialog.name
            resarr=nam_resarrlist[naml.index(nam)][1]
        self.namLE.setText(nam)
        for ra, isb, vsb in zip(resarr, self.iSpinBoxList, self.vSpinBoxList):
            isb.setValue(ra[0])
            vsb.setValue(ra[1])

class acharmonicsDialog(QDialog):
    def __init__(self, parent, h5path, h5expname, h5hpname, title='AC harmonics viewer', markersize=1):
        super(acharmonicsDialog, self).__init__(parent)
        
        mainlayout=QGridLayout()
        plotwlist=[]
        for i in range(3):
            for j in range(4):
                plotw=plotwidget(self)
                plotw.ax=plotw.axes
                plotw.ax2=None
                plotwlist+=[plotw]
                mainlayout.addWidget(plotw, i, j)
        
        
        hpsegdlist=CreateHeatProgSegDictList(h5path, h5expname, h5hpname)
        
        segindlist=[count for count, d in enumerate(hpsegdlist) if 'WinFFT_current' in d.keys()]
        strlist=['segment %d : %s' %(count, hpsegdlist[count]['segmenttype']) for count in segindlist]
        idialog=selectorDialog(self, strlist, title='Select segment to plot')
        if not idialog.exec_():
            return
        d=hpsegdlist[segindlist[idialog.index]]
        
        strlist=[k.partition('WinFFT_')[2] for k in d.keys() if k.startswith('WinFFT_')]
        idialog=selectorDialog(self, strlist, title='Select data type to plot')
        if not idialog.exec_():
            return
        segk=idialog.name
        fftsegk='WinFFT_'+segk
        liasegk='LIAharmonics_'+segk
        segk='sample'+segk
        ppc=pts_sincycle_h5(h5path, h5expname, h5hpname)
        x=numpy.array([0, .25, .5, .75, 1., 1.5, 2., 3., 4., 6., 8., 10., 15., 20., 50., 100.])*ppc
        x=numpy.round(x)
        x[x==0]=1
        strlist=['%d' %v for v in x]
        idialog=selectorDialog(self, strlist, title='Select plot point interval. %d pts in segment, %dppc' %(d['cycletime'].shape[1], ppc))
        if not idialog.exec_():
            return
        interval=x[idialog.index]
        indlist=range(0, d['cycletime'].shape[1], interval)
        
        if d['cycletime'].shape[0]>1:
            strlist=['%d' %i for i in  range(d['cycletime'].shape[0])]
            idialog=selectorDialog(self, strlist, title='Select cycle')
            if not idialog.exec_():
                return
            cycind=idialog.index
        else:
            cycind=0
            
        fftlab=['0', '0+', '1w-', '1w', '1w+', '2w-', '2w', '2w+', '3w-', '3w', '3w+']
        fftxyplot=[0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
        fftphplot=[-1, -1, -1, 5, -1, -1, 6, -1, -1, 7, -1]
        fftampplot=[-1, -1, 5, 5, 5, 6, 6, 6, 7, 7, 7]
        fftcol=['r', 'g', 'b', 'r', 'g', 'b', 'r', 'g', 'b', 'r', 'g']
        lialab=['lia1w', 'lia2w', 'lia3w']
        liaxyplot=[9, 10, 11]
        liaphplot=[-1, -1, -1]
        liaampplot=[5, 6, 7]
        liacol=['m', 'm', 'm']
        
        plotw=plotwlist[4]
        plotw.ax.plot(indlist, d[segk][cycind][indlist], 'k-', label=segk, markersize=markersize)
        
        fftarr=d[fftsegk][cycind][indlist].swapaxes(0, 1)
        for arr, l, xyp, php, ampp, c in zip(fftarr, fftlab, fftxyplot, fftphplot, fftampplot, fftcol):
            if xyp>=0:
                plotw=plotwlist[xyp]
                plotw.ax.plot(indlist, arr[:, 0], c+'-', label=l+'X', markersize=markersize)
                if plotw.ax2 is None:
                    plotw.ax2=plotw.ax.twinx()
                plotw.ax2.plot(indlist, arr[:, 1], c+':', label=l+'Y', markersize=markersize)

            if php>=0:
                plotw=plotwlist[php]
                if plotw.ax2 is None:
                    plotw.ax2=plotw.ax.twinx()
                plotw.ax2.plot(indlist, numpy.arctan(arr[:, 1]/arr[:, 0]), c+':', label=l+r'$\phi$', markersize=markersize)
            
            if ampp>=0:
                plotw=plotwlist[ampp]
                plotw.ax.plot(indlist, numpy.sqrt(arr[:, 1]**2+arr[:, 0]**2), c+'-', label=l+'A', markersize=markersize)

        liaarr=d[liasegk][cycind][indlist].swapaxes(0, 1)
        for arr, l, xyp, php, ampp, c in zip(liaarr, lialab, liaxyplot, liaphplot, liaampplot, liacol):
            if xyp>=0:
                plotw=plotwlist[xyp]
                plotw.ax.plot(indlist, arr[:, 0], c+'-', label=l+'X', markersize=markersize)
                if plotw.ax2 is None:
                    plotw.ax2=plotw.ax.twinx()
                plotw.ax2.plot(indlist, arr[:, 1], c+':', label=l+'Y', markersize=markersize)

            if php>=0:
                plotw=plotwlist[php]
                if plotw.ax2 is None:
                    plotw.ax2=plotw.ax.twinx()
                plotw.ax2.plot(indlist, numpy.arctan(arr[:, 1]/arr[:, 0]), c+':', label=l+r'$\phi$', markersize=markersize)
            
            if ampp>=0:
                plotw=plotwlist[ampp]
                plotw.ax.plot(indlist, numpy.sqrt(arr[:, 1]**2+arr[:, 0]**2), c+'--', label=l+'A', markersize=markersize)
        
        for count, plotw in enumerate(plotwlist):
            try:
                plotw.ax.yaxis.set_major_formatter(ExpTickLabels)
                leg=plotw.ax.legend(loc=2, frameon=False)
                #temp=[t.set_fontsize(12) for t in leg.texts]
            except:
                pass
            if not plotw.ax2 is None:
                try:
                    leg=plotw.ax2.legend(loc=1, frameon=False)
                    temp=[t.set_fontsize(12) for t in leg.texts]
                    plotw.ax2.yaxis.set_major_formatter(ExpTickLabels)
                except:
                    pass

        self.setLayout(mainlayout)

        QMetaObject.connectSlotsByName(self)

class textexportDialog(QDialog):
    def __init__(self, parent, h5path, h5expname, h5hpname, title='Export arrays as tab-delim text'):
        super(textexportDialog, self).__init__(parent)
        
        self.hpsdl=CreateHeatProgSegDictList(h5path, h5expname, h5hpname, expandmultdim=True)
        self.ncycs=self.hpsdl[0]['cycletime'].shape[0]
        for d in self.hpsdl:
            d['index']=numpy.float32([numpy.arange(d['segment_inds'][1]-d['segment_inds'][0])+d['segment_inds'][0] for count in range(self.ncycs)])
        
        self.keys=set([k for d in self.hpsdl for k, v in d.iteritems() if isinstance(v, numpy.ndarray) and v.shape==d['cycletime'].shape])
        self.keys=sorted(list(self.keys))
        
        ncols=int(numpy.ceil(len(self.keys)/10.))
        ninitrows=6
        
        seglayout=QGridLayout()
        lab=QLabel()
        s='Enter comma-delim list of segments to export. The length of each segment is\n'+','.join(['%d' %d['cycletime'].shape[1] for d in self.hpsdl])
        lab.setText(s)
        self.segLineEdit=QLineEdit()
        self.segLineEdit.setText(','.join(['%d' %count for count in range(len(self.hpsdl))]))
        seglayout.addWidget(lab, 0, 0)
        seglayout.addWidget(self.segLineEdit, 1, 0)
        
        cyclayout=QGridLayout()
        lab=QLabel()
        lab.setText('select cycle(s) to export')
        self.cycleComboBox=QComboBox()
        for i in range(self.ncycs):
            self.cycleComboBox.insertItem(i, `i`)
        self.cycleComboBox.insertItem(i+1, 'all')
        cyclayout.addWidget(lab, 0, 0)
        cyclayout.addWidget(self.cycleComboBox, 0, 1)
    
        fmtlayout=QGridLayout()
        lab=QLabel()
        lab.setText('fmt string:')
        self.fmtLineEdit=QLineEdit()
        self.fmtLineEdit.setText('%.4e')
        
        lab2=QLabel()
        lab2.setText('index interval:')
        self.intervSpinBox=QSpinBox()
        self.intervSpinBox.setRange(1, 1000000)
        self.intervSpinBox.setValue(1)
        
        fmtlayout.addWidget(lab, 0, 0)
        fmtlayout.addWidget(self.fmtLineEdit, 0, 1)
        fmtlayout.addWidget(lab2, 0, 2)
        fmtlayout.addWidget(self.intervSpinBox, 0, 3)
        
        
        saveButton=QPushButton()
        saveButton.setText("Save .txt file (make all selections first)")
        QObject.connect(saveButton, SIGNAL("pressed()"), self.save)
        
        doclab=QLabel()
        doclab.setText('below is a list of all array quantities and the segments over which they are available.\n'\
                    +'Check the boxes of those you want exported')
        mainlayout=QGridLayout()
        mainlayout.addLayout(seglayout, 0, 0, 2, ncols)
        mainlayout.addLayout(cyclayout, 2, 0, 1, ncols)
        mainlayout.addLayout(fmtlayout, 3, 0, 1, ncols)
        mainlayout.addWidget(saveButton, 4, 0, 1, ncols)
        mainlayout.addWidget(doclab, 5, 0, 1, ncols)
        
        
        self.cblist=[]
        for i, k in enumerate(self.keys):
            cb=QCheckBox()
            s=k+':'
            s+=','.join(['%d' %count for count, d in enumerate(self.hpsdl) if k in d.keys() and not numpy.any(numpy.isnan(d[k]))])
            cb.setText(s)
            mainlayout.addWidget(cb, ninitrows+i%10, i//10)
            self.cblist+=[cb]


        self.setLayout(mainlayout)

        QMetaObject.connectSlotsByName(self)

    def save(self):
        try:
            segs=eval('numpy.int16(['+str(self.segLineEdit.text())+'])')
        except:
            QMessageBox.warning(self,"FAILED",  "ABORTED because could not interpret segment list")
            return
        try:
            fmt=str(self.fmtLineEdit.text())
            fmt %1.1
        except:
            QMessageBox.warning(self,"FAILED",  "ABORTED because could not interpret fmt string")
            return
        
        fmtfcn=lambda x:(numpy.isnan(x) and ('NaN',) or (fmt %x,))[0]
        
        indinterv=self.intervSpinBox.value()
        
        p=mygetsavefile(parent=self, markstr='filename for text data export - WILL OVERWRITE IF EXISTS')
        if len(p)==0:
            return
        cycoptind=self.cycleComboBox.currentIndex()
        if cycoptind>=self.ncycs:
            cycl=range(self.ncycs)
            pstart, pext=os.path.splitext(p)
            getp=lambda cyci:''.join((pstart, '_%d_of_%d' %(cyci+1, self.ncycs), pext))
        else:
            cycl=[cycoptind]
            getp=lambda cyci:p
        keyl=[k for k, cb in zip(self.keys, self.cblist) if cb.isChecked()]
        
        getarr=lambda k, segi, cyci, indinterv:(k in self.hpsdl[segi] and (self.hpsdl[segi][k][cyci][::indinterv],) or (self.hpsdl[segi]['cycletime'][cyci][::indinterv]+numpy.nan,))[0]
        for cyci in cycl:
            path=getp(cyci)
            arr=numpy.array([numpy.concatenate([getarr(k, segi, cyci, indinterv) for segi in segs]) for k in keyl])
        s='\n'.join(['\t'.join([fmtfcn(x) for x in a]) for a in arr.T])
        f=open(path, mode='w')
        f.write('\t'.join(keyl)+'\n'+s)
        f.close()

PinoutLibrary={\
'xiaodong2010':numpy.array(\
[[[1,4],[2,3]],\
[[4,7],[5,6]],\
[[10,13],[11,12]],\
[[18,21],[19,20]],\
[[21,25],[22,24]],\
[[93,97],[94,95]],\
[[97,100],[98,99]],\
[[14,17],[15,16]],\
[[28,31],[29,30]],\
[[31,34],[32,33]],\
[[87,89],[86,88]],\
[[89,91],[90,92]],\
[[8,58],[9,59]],\
[[35, 37],[36, 41]],\
[[37, 39],[38, 40]],\
[[80,82],[79,81]],\
[[82,85],[83,84]],\
[[64,66],[65,67]],\
[[42,47],[43,48]],\
[[45,47],[44,46]],\
[[71,74],[72,75]],\
[[68,71],[69,70]],\
[[61,64],[62,63]],\
[[52,58],[55,57]],\
[[49,52],[50,51]]])\
}
