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
from PnSC_main import *


cell=25
p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/AuSiCu_pnsc_all.h5'

#name, (cycles,daqHz)
exp_rec=[\
('heat1a',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat1b',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat1d',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat1c',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat2',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat3',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat4a',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat4b',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat6a',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat6b',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat8',[(2, 'Dfit_t4_c2'), (2, 'C_subDfit_interpars1')]),\
]

if mm is None:
    mm=start()
#p=mm.h5path
#f=h5py.File(p, mode='r+')

h5expsource='heat2'
#h5expsource='heat6c'
forcecopyrecipes=True

for exp, recipelist in exp_rec:
    for segind, rec in recipelist:
        print 'BATCH PROCESS OF ', exp, segind, rec
        f=h5py.File(p, mode='r')
        h5hp=gethpgroup(f, exp)
        h5hpdflt=None
        for h5hpgrp in h5hp.itervalues():
            try:
                if h5hpgrp.attrs['CELLNUMBER']==cell:
                    h5hpdflt=h5hpgrp.name.rpartition('/')[2]
            except:
                pass
        f.close()
        if h5hpdflt is None:
            print 'no heat program found'
            continue
        print h5hpdflt
        id=SCanalysisDialog(None, p, exp, h5hpdflt=h5hpdflt)
        id.segComboBox.setCurrentIndex(segind)
        
        f=h5py.File(p, mode='r+')#r+ in case nSCrecipe needs to be created
        h5srcrec=getSCrecipegrp(f, exp)
        
        if not (rec in h5srcrec) or (forcecopyrecipes and exp!=h5expsource):
            f.close()
            copySCrecipes(p, exp, h5expsource)
        else:
            f.close()
        id.recComboBox.clear()
        id.recComboBox.insertItem(0, rec)
        id.recComboBox.setCurrentIndex(0)
        print 'calcall starting'
        id.calcall()
print 'batch done'
