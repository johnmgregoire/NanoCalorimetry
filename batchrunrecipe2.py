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


#name, (cycles,daqHz)
exp_rec=[\
#('heat1a',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat1b',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat1d',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat1c',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat1_cell21_other',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat2',[(2, 'CTpks_dflt')]),\
#('heat2',[(2, 'Dfit_t5_c1_existpt'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat2_cell21_other',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\

('heat3',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat4a',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat4b',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat6a',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat6b',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1')]),\
('heat8',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1')]),\

#('heat3',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat3_cell21_other',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat4a',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat4b',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat4_cell21_other',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat4c',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (3, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3')]),\
#('heat6a',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat6b',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\


#('heat6c',[(2, 'IVOL8_RTPDSG10o1_SdtSG15o1'), (3, 'IVOL8_RTPDSG10o1_SdtSG15o1')]),\
#('heat6_cell21_other',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
#('heat7',[(2, 'IVOL8_RTPDSG10o1_SdtSG15o1'), (3, 'IVOL8_RTPDSG10o1_SdtSG15o1')]),\
#('heat8',[(2, 'Dfit_T4_c2'), (2, 'C_subDfit_interpars1'), (2, 'CTpks_dflt')]),\
]

p=mm.h5path
#f=h5py.File(p, mode='r+')

h5expsource='heat2'
#h5expsource='heat6c'
forcecopyrecipes=True
for exp, recipelist in exp_rec:
    for segind, rec in recipelist:
        print 'BATCH PROCESS OF ', exp, segind, rec
        id=SCanalysisDialog(None, p, exp, h5hpdflt=None)
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
