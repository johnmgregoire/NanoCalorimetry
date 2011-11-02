import h5py,numpy, os, operator, pylab, time

p='F:/CHESS2011_h5MAIN'
px='F:/CHESS2011XRD_asimported'

fnxrd_fnpnsc=[
#('2011Jun01b_AuSiCu.h5',['AuSiCuheat1','AuSiCuheat2','AuSiCuheat3','AuSiCuheat4']),\
#('2011Jun01b_NiTiHf.h5',['NiTiHfheat1','NiTiHfheat1_MA','NiTiHfheat1_fast','NiTiHfheat1_slow','NiTiHfheat2','NiTiHfheat2_MA']),\
#('2011Jun01A_ZrCuAl_heat0.dat.h5', '2011Jun01a.h5', ['ZrCuAlheat1', 'ZrCuAlheat2', 'ZrCuAlheat3', 'ZrCuAlheat4', 'ZrCuAlheat5']), \
#('2011Oct02D_AuSiCu.h5', ['AuSiCuheats']), \
('2011Oct02D_BiInSn.h5', ['Bi_DCheats', 'Bi_ACheats', 'In_DCheats', 'In_ACheats', 'Sn_DCheats', 'Sn_ACheats']), \
#('2011Oct10B_NiTiHf.h5', ['NiTiHfheat1', 'NiTiHfheat2', 'NiTiHfheat3', 'NiTiHfheat1_MA', 'NiTiHfheat2_MA', 'NiTiHfheat3_MA']), \
#('2011Oct10B_FeNi.h5', ['DCheats', 'ACheats']), \
#('2011Oct10C.h5', ['borides']), \
#('2011Oct10D.h5', ['DCheats', 'ACheats']), \
#('BackgroundImages.dat.h5', serp), \
#('nosampleconfigurations.dat.h5', serp), \
]


tryattr=lambda pnt, s:(s in pnt.attrs.keys() and (pnt.attrs[s],) or (None,))[0]
x=[]
y=[]
for fn, expgrplist in fnxrd_fnpnsc:
    if not fn.endswith('.h5'):
        continue
    f=h5py.File(os.path.join(p, fn), mode='r')
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
            sptup=tuple()
            for k in ['specscan','prespecscan','postspecscan']:
                t=tryattr(node2,k)
                if not t is None:
                    pnt=f[`t`]['measurement/scalar_data/Seconds']
                    try:
                        sec=pnt[0]
                    except:
                        sec=pnt.value
                    t=(t, sec)
                sptup+=(t,)
            
            hppnt_epoch+=[(node2.attrs['epoch'], node2.attrs['CELLNUMBER'], node2.attrs['segment_ms'][-1], numpy.max(node2['samplecurrent'][:, :]))+sptup]
    
    
    f.close()
    

    eparr=numpy.array([tup[0] for tup in hppnt_epoch])
    hppnt_epoch=[hppnt_epoch[i] for i in numpy.argsort(eparr)]
    
    carr=numpy.array([tup[1] for tup in hppnt_epoch])
    allcells=sorted(list(set(carr)))
    maxnumexp=max([(carr==c).sum() for c in allcells])
    sc=[]
    se=[[] for i in range(maxnumexp)]
    for count, cell in enumerate(allcells):
        inds=numpy.where(carr==cell)[0]
        sc+=['cell%02d' %cell]
        for l, i in zip(se, inds):
            e, c, maxms, maxma, stup, prestup, poststup=hppnt_epoch[i]
            s='%.2fs PnSC up to %.1fmA' %(maxms/1000., maxma)
            if not stup is None:
                s+=' with insitu XRD'
            if not prestup is None:
                s+=', %ds preXRD' %prestup[1]
            if not poststup is None:
                s+=', %ds postXRD' %poststup[1]
            l+=[s]
    print fn
    print '\t'.join(sc)
    for l in se:
        print '\t'.join(l)
print 'done'

