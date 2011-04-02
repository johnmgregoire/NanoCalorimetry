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


rootfolder='E:/CHESS2010PnSC'

#name, (cycles,daqHz)
folders_pars=[\
('heat1a',(1,1e5, .00003)),\
('Rcal1a_before',(10,1e5, .0002)),\
('Rcal1a_after',(10,1e5, .0002)),\
('heat1b',(1,1e5, .00003)),\
('heat1d',(1,1e5, .00003)),\
('Rcal1d_after',(10,1e5, .0002)),\
('heat1c',(1,1e5, .00003)),\
('heat1_cell21_other',(1,1e5, .00003)),\
('heat2',(1,1e5, .00003)),\
('Rcal2_before',(1,1e5, .0002)),\
('Rcal2_after',(1,1e5, .0002)),\
('heat2_cell21_other',(1,1e5, .00003)),\
('Rcal2_cell21_other_before',(1,1e5, .0002)),\
('Rcal2_cell21_other_after',(1,1e5, .0002)),\
('heat3',(1,1e5, .00003)),\
('heat3_cell21_other',(1,1e5, .00003)),\
('Rcal3_cell21_other_before',(1,1e5, .0002)),\
('Rcal3_before',(1,1e5, .0002)),\
('Rcal3_after',(1,1e5, .0002)),\
('heat4a',(1,1e5, .00003)),\
('heat4b',(1,1e5, .00003)),\
('Rcal4a_before',(1,1e5, .0002)),\
('Rcal4bc_after',(1,1e5, .0002)),\
('heat4_cell21_other',(1,1e5, .00003)),\
('Rcal4_cell21_other_before',(1,1e5, .0002)),\
('heat4c',(1,1e5, .000001)),\
('heat6a',(1,1e5, .00003)),\
('heat6b',(1,1e5, .00003)),\
('heat6c',(1,1e3, .00003)),\
('Rcal6a_before',(1,1e5, .0002)),\
('Rcal6ab_after',(1,1e5, .0002)),\
('Rcal6c_before',(1,1e3, .0002)),\
('Rcal6c_after',(1,1e3, .0002)),\
('heat6_cell21_other',(1,1e5, .00003)),\
('Rcal6_cell21_other_before',(1,1e5, .0002)),\
('heat7',(1,1e3, .00003)),\
('Rcal7_before',(1,1e3, .0002)),\
('Rcal7_after',(1,1e3, .0002)),\
('heat8',(1,1e5, .00003)),\
('Rcal8_before',(1,1e5, .0002)),\
('Rcal8_after',(1,1e5, .0002))]

#p='E:/CHESS2010PnSC/2010Nov27_AuSiCu_pnsc.h5'
p=mm.h5path
#f=h5py.File(p, mode='r+')
#:2,2:8,8:20,20:
for folder, (ncycles, daqHz, sd) in [folders_pars[i] for i in [31, 32]]:#[20, 25]
    print 'BATCH PROCESS OF ', folder
    create_exp_grp(p,folder)
    for k,v in [('grpname',folder), ('ncycles', ncycles), ('daqHz', daqHz), ('secderptsSpinBox',30), ('secdervalSpinBox',sd), ('durSpinBox',None), ('nnoiseSpinBox',50), ('naboveSpinBox',20), ('protname','PatDAQ_SC'), ('firstderptsSpinBox',15),('nsigSpinBox',1.3)]:
        mm.batchattrdict[k]=v
    mm.batchrun_files(os.path.join(rootfolder,folder), skiperrors=False)

print 'batch done'
