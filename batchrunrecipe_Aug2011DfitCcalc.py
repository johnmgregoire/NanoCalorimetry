import time
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


#name, (cycles,daqHz)
exp_rec=[\
('heat1a',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat1b',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat1d',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat1c',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat2',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat3',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
#('heat4a',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat4b',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
#('heat6a',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat6b',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
('heat8',[(2, 'Dfit_t4_c2_existing'), (2, 'C_subDfit_interpars1')]),\
]

p=mm.h5path
#f=h5py.File(p, mode='r+')

h5expsource='heat1a'
#h5expsource='heat6c'
forcecopyrecipes=False
for exp, recipelist in exp_rec:
    for segind, rec in recipelist:
        print 'BATCH PROCESS OF ', exp, segind, rec
        id=SCanalysisDialog(None, p, exp, h5hpdflt=None)
        id.segComboBox.setCurrentIndex(segind)
        
        f=h5py.File(p, mode='r+')#r+ in case nSCrecipe needs to be created
        h5srcrec=getSCrecipegrp(f, exp)
        if rec in h5srcrec and not forcecopyrecipes:
            f.close()
        else:
            f.close()
            copySCrecipes(p, exp, h5expsource)
        id.recComboBox.clear()
        id.recComboBox.insertItem(0, rec)
        id.recComboBox.setCurrentIndex(0)
        print 'calcall starting'
        id.calcall()
print 'batch done'
