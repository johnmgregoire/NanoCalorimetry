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


def mygetopenfile(parent=None, xpath="%s" % os.getcwd(),markstr='', filename='' ):
    if parent is None:
        xapp = QApplication(sys.argv)
        xparent = QWidget()
        returnfn = unicode(QFileDialog.getOpenFileName(xparent,''.join(['Select file to open:', markstr]),os.path.join(xpath, filename).replace('\\','/')))
        xparent.destroy()
        xapp.quit()
        return returnfn
    return unicode(QFileDialog.getOpenFileName(parent,''.join(['Select file to open: ', markstr]),os.path.join(xpath, filename).replace('\\','/')))

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
    def __init__(self, tree, h5file, showattrs=True):
        self.treeWidget=tree
        self.treeWidget.clear()
        
        self.showattrs=showattrs
        
        mainitem=QTreeWidgetItem([os.path.split(h5file.filename)[1]],  0)
        self.treeWidget.addTopLevelItem(mainitem)
        self.createTree(h5file, mainitem)

    def createTree(self, startnode, parentitem):
        #print startnode.name
        #print startnode.listobjects()
        for node in startnode.iterobjects():
            if isinstance(node, h5py.Dataset):
                item=QTreeWidgetItem([node.name.rpartition('/')[2]+`node.shape`],  0)
                parentitem.addChild(item)
                if self.showattrs:
                    for attrname, attrval in node.attrs.iteritems():
                        attritem=QTreeWidgetItem([self.attrstring(attrname, attrval)],  0)
                        item.addChild(attritem)
            elif isinstance(node, h5py.Group):
                item=QTreeWidgetItem([node.name.rpartition('/')[2]],  0)
                parentitem.addChild(item)
                self.createTree(node, item)
                if self.showattrs:
                    for attrname, attrval in node.attrs.iteritems():
                        attritem=QTreeWidgetItem([self.attrstring(attrname, attrval)],  0)
                        item.addChild(attritem)

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
                    if c=='None':
                        c=None
                    else:
                        temp=c.lstrip('0')
                        if temp=='' and '0' in c:
                            c=0
                        else:
                            c=eval(temp)
                except:
                    pass
                a=a.strip("'").strip('"')
                print a, a in self.attrd.keys()
                self.attrd[a]=c
                print a,  ' updated to ',  `c`
                self.edited=True
                
class SegmentEditor(QDialog):
    def __init__(self, parent, SegmentData, cycledata, maxpts=9, markersize=2):
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
        self.firstderptsSpinBox.setValue(7)

        secderptsLabel=QLabel()
        secderptsLabel.setText('num pts/n2nd der')
        self.secderptsSpinBox=QSpinBox()
        self.secderptsSpinBox.setRange(5, 101)
        self.secderptsSpinBox.setValue(15)

        secdervalLabel=QLabel()
        secdervalLabel.setText('crit val\n2nd der')
        self.secdervalSpinBox=QDoubleSpinBox()
        self.secdervalSpinBox.setDecimals(4)
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
                    print '@@', minval
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
        
        self.arrComboBoxlist=[QComboBox() for i in range(3)]
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
        
        mainlayout=QGridLayout()
        mainlayout.addWidget(l0, 0, 0, 4, 1)
        mainlayout.addWidget(l1, 4, 0, 4, 1)
        for w, n in zip(self.arrComboBoxlist, [0, 4, 6]):
            mainlayout.addWidget(w, n, 1, 1, 1)
        mainlayout.addWidget(l2, 0, 2, 2, 1)
        mainlayout.addWidget(l3, 4, 2, 2, 1)
        for w, n in zip(self.dfltSpinBoxlist, [2, 3, 6]):
            mainlayout.addWidget(w, n, 2, 1, 1)
        for w, n in zip(self.useaveCheckBoxlist, [0, 4]):
            mainlayout.addWidget(w, n, 3, 4, 1)
            
            
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
            if 'CellResistance' in grp['analysis']:
                #self.imagenamelist+=[(grp.name).rpartition('/')[2]]
                for cb in self.arrComboBoxlist:
                    cb.insertItem(counter, (grp.name).rpartition('/')[2])
                self.arrpathlist+=[grp['analysis/CellResistance'].name]
                counter+=1
        self.h5file.close()
    def calcansave(self):
        self.h5file=h5py.File(self.h5path, mode='r+')
        sbvals=[sb.value() for sb in self.dfltSpinBoxlist]
        Rdfltbool=sbvals[0]>0. or len(self.arrpathlist)==0
        Adfltbool=sbvals[2]>0. or len(self.arrpathlist)==0
        
        pnts=[self.h5file[self.arrpathlist[cb.currentIndex()]] for cb in self.arrComboBoxlist if len(self.arrpathlist)>0]
        
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
            savearr[:, 0]=pnts[0][:]
            savearr[:, 1]=pnts[0].attrs['ambient_tempC']
            if self.useaveCheckBoxlist[0].isChecked():
                inds=numpy.where(savearr[:, 0]<=0.)
                notinds=numpy.where(savearr[:, 0]>0.)
                savearr[inds, 0]=numpy.mean(savearr[notinds, 0])
                savearr[inds, 1]=numpy.mean(savearr[notinds, 1])
        if Adfltbool:
            savearr[:, 2]=sbvals[2]
        else:
            r1=pnts[1][:]
            t1=pnts[1].attrs['ambient_tempC']
            r2=pnts[2][:]
            t2=pnts[2].attrs['ambient_tempC']
            
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
        self.h5file.close()

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
        
        plotw=plotwidget(self)
        mainlayout.addWidget(plotw, 0, 0, 1, 2)
        for count, arr in enumerate(data):
            if xdata is None:
                plotw.axes.plot(arr, style)
            else:
                plotw.axes.plot(xdata[count], arr, style)
        plotw.axes.set_xlabel('array index')
        plotw.axes.set_ylabel('array value')

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
            h5file[path].attrs[k]=v
        h5file.close()
