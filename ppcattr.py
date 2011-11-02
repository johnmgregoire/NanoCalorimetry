import h5py,numpy, os, operator, pylab, time

p='F:/CHESS2011_h5MAIN'
px='F:/CHESS2011XRD_asimported'

fnxrd_fnpnsc=[
##('2011Jun01b_AuSiCu.h5',['AuSiCuheat1','AuSiCuheat2','AuSiCuheat3','AuSiCuheat4']),\
##('2011Jun01b_NiTiHf.h5',['NiTiHfheat1','NiTiHfheat1_MA','NiTiHfheat1_fast','NiTiHfheat1_slow','NiTiHfheat2','NiTiHfheat2_MA']),\
##('2011Jun01a.h5', ['ZrCuAlheat1', 'ZrCuAlheat2', 'ZrCuAlheat3', 'ZrCuAlheat4', 'ZrCuAlheat5']), \
##('2011Oct02D_AuSiCu.h5', ['AuSiCuheats']), \
#('2011Oct02D_BiInSn.h5', ['Bi_ACheats', 'In_ACheats', 'Sn_ACheats']), \
##('2011Oct10B_NiTiHf.h5', ['NiTiHfheat1', 'NiTiHfheat2', 'NiTiHfheat3', 'NiTiHfheat1_MA', 'NiTiHfheat2_MA', 'NiTiHfheat3_MA']), \
#('2011Oct10B_FeNi.h5', ['ACheats']), \
#('2011Oct10C.h5', ['ACheats']), \
#('2011Oct10D.h5', ['ACheats']), \
#('BackgroundImages.dat.h5', serp), \
#('nosampleconfigurations.dat.h5', serp), \
]

savebool=0
getfromorigfilename=0
x=[]
y=[]
for fn, expgrplist in fnxrd_fnpnsc:
    if not fn.endswith('.h5'):
        continue
    f=h5py.File(os.path.join(p, fn), mode='r+')
    hppnt_epoch=[]
    #for node in f['Calorimetry'].values():
    for expgrp in expgrplist:
        node=f['Calorimetry'][expgrp]
#        for nn in node.values():
#            print nn, 'samplecurrent' in nn
#            if 'samplecurrent' in nn:
#                nam=nn.name.rpartition('/')[2]
#                f.copy(node[nam], node['measurement/HeatProgram'])
#                del node[nam]
        if not 'measurement/HeatProgram' in node:
            continue
        for node2 in node['measurement/HeatProgram'].itervalues():
            ppc=None
            nam=node2.name.rpartition('/')[2]
            getfromattr=True
            if getfromorigfilename:
                if 'ppc' in nam:
                    a, b, c=nam.partition('ppc')
                    a=a[::-1]
                    s=''
                    while a[0].isdigit() and len(a)>1:
                        s+=a[0]
                        a=a[1:]
                    try:
                        ppc=eval(s[::-1])
                        getfromattr=False
                    except:
                        pass
            getfromattr=getfromattr and 'ppc' in node2.attrs.keys()
            if getfromattr:
                ppc=node2.attrs['ppc']
                temp=''
            else:
                temp=', read from filename'
            if not ppc is None:
                print nam, ' : ', ppc, temp
                if savebool:
                    node2.attrs['pts_sincycle']=ppc
                    node2.attrs['ppc']=ppc
    f.close()

print 'done'

