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


rootfolder='E:/CHESS2010PnSC'

#name, (cycles,daqHz)
exp_rec=[\
#('heat1a',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1b',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1d',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1c',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat1_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat2',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat2_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat3',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat3_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4a',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4b',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat4c',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (3, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3')]),\
#('heat6a',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
#('heat6b',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
('heat6c',[(2, 'IVOL8_RTPDSG10o1_SdtSG15o1'), (3, 'IVOL8_RTPDSG10o1_SdtSG15o1')]),\
#('heat6_cell21_other',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
('heat7',[(2, 'IVOL8_RTPDSG10o1_SdtSG15o1'), (3, 'IVOL8_RTPDSG10o1_SdtSG15o1')]),\
#('heat8',[(2, 'OL10gap1_RTSG40o2_dtSG40o2_PDSG40o3'), (4, 'OL10gap1_RTPSG160o2_dtSG250o1_DSG160o2')]),\
]
p='E:/CHESS2010PnSC/AuSiCu_pnsc_all.h5'
#p=mm.h5path
#f=h5py.File(p, mode='r+')

#h5expsource='heat1a'
h5expsource='heat6c'
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
