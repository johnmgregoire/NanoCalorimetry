import h5py,numpy, os, operator, pylab, time

p='F:/CHESS2011_h5MAIN'
px='F:/CHESS2011XRD_asimported'

savebool=1
predeleteattrs=1
x=[]
y=[]

for fnx in os.path.listdir(px):
    if not fnx.endswith('.h5'):
        continue
    pxr=os.path.join(px, fnx)
    f=h5py.File(pxr, mode='r+')
    for node in f['Calorimetry'].values():
        if 'Ro' in node.name:
            continue
        
        node=[expgrp]
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

