import time, copy
import os
import sys
import numpy
import h5py
#from PnSC_ui import *
#from PnSC_dataimport import *
from PnSC_SCui import *
#from PnSC_math import *
from PnSC_h5io import *
from PnSC_main import *
from matplotlib.ticker import FuncFormatter
import scipy.integrate
from scipy.interpolate import griddata
import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.ticker

#p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'

def myexpformat(x, pos):
    for ndigs in range(2):
        lab=(('%.'+'%d' %ndigs+'e') %x).replace('e+0','e').replace('e+','e').replace('e0','').replace('e-0','e-')
        if eval(lab)==x:
            return lab
    return lab
ExpTickLabels=FuncFormatter(myexpformat)

def make_ticklabels_invisible(ax, x=True, y=True, ticksalso=False):
    if x:
        for tl in ax.get_xticklabels():
            tl.set_visible(False)
        if ticksalso:
            ax.xaxis.set_major_locator(matplotlib.ticker.NullLocator())
    if y:
        for tl in ax.get_yticklabels():
            tl.set_visible(False)
        if ticksalso:
            ax.yaxis.set_major_locator(matplotlib.ticker.NullLocator())
savef='C:/Users/JohnnyG/Documents/HarvardWork/MG/PnSCplots/batchplotbycell_June2'

for i in range(1, 26):
    folder=os.path.join(savef, 'cell%02d' %i)
    if os.path.exists(folder):
        for fn in os.listdir(folder):
            if fn.startswith('Cpstack'):
                pylab.subplot(5, 5, i, frame_on=False)
                make_ticklabels_invisible(pylab.gca(), ticksalso=True)
                im=mpimg.imread(os.path.join(folder, fn))
                pylab.imshow(im)
pylab.show()
