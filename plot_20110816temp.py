import numpy, h5py, os
from PnSC_main import *
from PnSC_h5io import *

p='C:/Users/JohnnyG/Documents/PythonCode/Vlassak/NanoCalorimetry/20110816_Zr-Hf-B.h5'
h5f=h5py.File(p, mode='r')

ehl=[\
('quadlinheating2', 'pre_25mApluslinquad2_cell17_1_of_1', 'Zr-B, 1st'),\
('quadlinheating2', 'cell17_25malinquad2repeat_1_of_1', 'Zr-B, 2nd'),\
('quadlinheating2', 'cell17_25mAquadlincurrent_3rdtime_1_of_1', 'Zr-B, 3rd'), \
('quadlinheating2', 'pre_25mApluslinquad2_cell16_1_of_1', 'Hf-B, nth'), \
('quadlinheating2', 'cell11_25malinquad2_1_of_1', 'empty'), \
]
tarrs=[]
pprarrs=[]
for i, (e, h, l) in enumerate(ehl):
    hpsdl=CreateHeatProgSegDictList(p, e, h)
    T=hpsdl[2]['sampletemperature'][0, :]
    ppr=hpsdl[2]['samplepowerperrate'][0, :]
    pylab.plot(T, ppr*1.e6, label=l)

pylab.legend()
pylab.xlabel('Temperature (C)')
pylab.ylabel('Power per rate ($\mu$J/K)')
pylab.show()

