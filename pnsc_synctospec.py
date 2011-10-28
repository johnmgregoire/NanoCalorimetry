import h5py,numpy, os, operator, pylab, time

p='F:/CHESS2011_h5MAIN'
px='F:/CHESS2011XRD_asimported'

fnxrd_fnpnsc=[
#('2011Jun01B.dat.h5','2011Jun01b_AuSiCu.h5',['AuSiCuheat1','AuSiCuheat2','AuSiCuheat3','AuSiCuheat4']),\
#('2011Jun01B.dat.h5','2011Jun01b_NiTiHf.h5',['NiTiHfheat1','NiTiHfheat1_MA','NiTiHfheat1_fast','NiTiHfheat1_slow','NiTiHfheat2','NiTiHfheat2_MA']),\
#('20101127AuSiCu_cell11.dat.h5', revstrip), \
#('2011Jun01A_ZrCuAl_heat0.dat.h5', '2011Jun01a.h5', ['ZrCuAlheat1', 'ZrCuAlheat2', 'ZrCuAlheat3', 'ZrCuAlheat4', 'ZrCuAlheat5']), \
#('2011Jun01B.dat.h5', '2011Jun01b.h5'), \
#('2011Oct02D_AuSiCu.dat.h5', '2011Oct02D.h5'), \
#('2011Oct02D_InSnBi.dat.h5', '2011Oct02D.h5'), \
#('2011Oct10B.dat.h5', '2011Oct10B_NiTiHf.h5'), \
#('2011Oct10B_FeNi.dat.h5', '2011Oct10B_FeNi.h5'), \
#('2011Oct10C.dat.h5', '2011Oct10C.h5', ['borides']), \
#('2011Oct10D_NiTiHf.dat.h5', '2011Oct10B_NiTiHf.h5'), \
#('BackgroundImages.dat.h5', serp), \
#('nosampleconfigurations.dat.h5', serp), \
]

savebool=False
predeleteattrs=False
x=[]
y=[]
for fnx, fn, expgrplist in fnxrd_fnpnsc:
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
            hppnt_epoch+=[(node2, node2.attrs['epoch'], node2.attrs['CELLNUMBER'])]
            if predeleteattrs:
                for k in ['specscan','prespecscan','postspecscan']:
                    if k in node2.attrs.keys():
                        node2.attrs['backup'+k]=node2.attrs[k]
                        del node2.attrs[k]
    
    specname_ep_cell_insitu=[]
    fx=h5py.File(os.path.join(px, fnx), mode='r')
    for node in fx.itervalues():
        node2=node['measurement/scalar_data/Epoch']
        try:
            node3=node['analysis/otherdata/cellnumber']
        except:
            continue
        if len(node2)>1 and len(node3)>1 and node3[0]!=node3[1]:
            for ep, c in zip(node2[:], node3[:]):
                specname_ep_cell_insitu+=[(node.name.rpartition('/')[2], ep, c, 0)]
        elif len(node3)>1 and node.attrs['acquisition_command'].startswith('tseries 1'):
            specname_ep_cell_insitu+=[(node.name.rpartition('/')[2], node2.value, node3[0], 1)]
        else:
            try:
                ep=node2[0]
            except:
                ep=node2.value
            specname_ep_cell_insitu+=[(node.name.rpartition('/')[2], ep, node3[0], 0)]
    #y+=specname_epoch
    fx.close()
    

    eparr=numpy.array([ep for pnt, ep, c in hppnt_epoch])
    carr=numpy.array([c for pnt, ep, c in hppnt_epoch])
    for sp, epsp, csp, insitu in specname_ep_cell_insitu:
        ind=numpy.argmin((eparr-epsp)**2)
        pnt=hppnt_epoch[ind][0]
        c=hppnt_epoch[ind][2]
        diff=eparr[ind]-epsp
        if insitu:
            if numpy.abs(diff)>20.:
                print 'cannot find match for insitu spec %s. The closest is %s which occured %ds later' %(sp, pnt.name.rpartition('/')[2], diff)
                continue
            elif csp!=c:
                print 'in situ spec %s match time with %s but spec says cell %d and onsc says cell %d' %(sp, pnt.name.rpartition('/')[2], csp, c)
                continue
            print 'insitu spec %s with %s' %(sp, pnt.name.rpartition('/')[2])
            if savebool:
                pnt.attrs['specscan']=eval(sp)
        else:
            presp=numpy.where((carr==csp)&(eparr>epsp))
            if len(presp[0])>0:
                prespind=presp[0][numpy.argmin(eparr[presp]-epsp)]
                pnt=hppnt_epoch[prespind][0]
                print 'prespec %s with %s' %(sp, pnt.name.rpartition('/')[2])
                if savebool:
                    pnt.attrs['prespecscan']=eval(sp)
            postp=numpy.where((carr==csp)&(eparr<epsp))
            if len(postp[0])>0:
                postspind=postp[0][numpy.argmin(epsp-eparr[postp])]
                pnt=hppnt_epoch[postspind][0]
                print 'postspec %s with %s' %(sp, pnt.name.rpartition('/')[2])
                if savebool:
                    pnt.attrs['postspecscan']=eval(sp)
#    for pnt, ep in hppnt_epoch:
#        if not 'specscan' in pnt.attrs.keys():
#            continue
#        epx=[t for n, t in specname_epoch if n==`pnt.attrs['specscan']`]
#        if len(epx)!=1:
#            print 'PROBLEM: ', epx, pnt.name, pnt.attrs['specscan']
#            continue
#        epx=epx[0]
#        print epx-ep, pnt.name.rpartition('/')[2]

    #x+=[(node.name.rpartition('/')[2], ep) for node, ep in hppnt_epoch]
    f.close()

#el=[ep for node, ep in x]
#el=numpy.sort(numpy.array(el))
#pylab.subplot(211)
#pylab.plot(el/60./60./24.)
#
#for elv in el:
#    print time.localtime(elv)
#    
#el=[numpy.mean(numpy.array(ep)) for node, ep in y]
#el=numpy.sort(numpy.array(el))
#pylab.subplot(212)
#pylab.plot(el/60./60./24.)


#for i in numpy.argsort(el):
#    print el[i]-min(el), x[i][0]
print 'done'

