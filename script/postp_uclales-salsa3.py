def openr(filepref,ff):
    import netCDF4 as nc
    fname = '%s.%04u%04u.nc' % (filepref, ff[0], ff[1])
    return nc.Dataset(fname,'r')

def openw(filepref):
    import netCDF4 as nc
    fname = '%s_pp.nc' % (filepref)
    return nc.Dataset(fname,'w',format='NETCDF3_CLASSIC')

def opena(filepref):
    import netCDF4 as nc
    fname = '%s_pp.nc' % (filepref)
    return nc.Dataset(fname,'a')

def getVariable(fid,vname):
    import numpy as np
    
    val = fid.variables[vname][:]
    dims = fid.variables[vname].dimensions
    coords = np.zeros(len(dims),dtype=np.object)
    c = 0
    for d in dims:
        coords[c] = fid.variables[d][:]
        c+=1

    #print coords
    return val,dims,coords

def getDimension(fid,dname):
    import numpy as np

    # Check that both dimension and axis variable exist
    if dname in fid.variables.keys():
        val = fid.variables[dname][:]
        varfound = True
    else:
        val = np.array([-999.])
        varfound = False

    dlen = len(val)

    return val,dlen,varfound

def getNames(fid,itype):

    # itype == 1 : Dimensions
    # itype == 2 : Variables

    if itype == 1:
        return fid.dimensions.keys()
    if itype == 2:
        return fid.variables.keys()
    

def makeVariables(fidout,filepref,imaax,jmax,dimnames,vars):
    import numpy as np
    import gc

    if vars == []:
        idin = openr(filepref,(0,0))
        vars = getNames(idin,2)
        idin.close()

    ii = range(imax)
    jj = range(jmax)
    ifn = []
    for i in ii:
        for j in jj:
            ifn.append((i,j))
            
    del ii,jj

    itfile = iter(ifn)
    # Iterate files
    try:
        while 1:
            fn = itfile.next()
            print 'file ', fn
            idin = openr(filepref,fn)

            # Iterate variables
            itvars = iter(vars)
            try:
                while 1:
                    vn = itvars.next()

                    print vn
                    print '--------'

                    vv = []

                    if vn in dimnames:
                        print 'axis variable, skip'
                        continue   
                    
                    val,dims,coords = getVariable(idin,vn)
           
                    if len(dims) > 4:
                        continue

                    if fn == (0,0):
                        vv = fidout.createVariable(vn,'d',dims)
                    else:
                        vv = fidout.variables[vn]

                    # Find out the axis coordinates for current block
                    mask1d = np.zeros(len(dims),dtype=np.object)
                    m = 0
                    for d in dims:
                        mask1d[m] = np.in1d(fidout.variables[d][:],coords[m])
                        m+=1    
               
                    # Different dimensionalities (shoudl be a nicer way to do this...)
                    if len(dims) == 1:
                        l1 = np.where( mask1d[0] )[0]
                        vv[l1[0]:l1[-1]+1] = val
                        del l1
                    elif len(dims) == 2:
                        l1 = np.where( mask1d[0] )[0]
                        l2 = np.where( mask1d[1] )[0]
                        vv[l1[0]:l1[-1]+1,l2[0]:l2[-1]+1] = val
                        del l1,l2
                    elif len(dims) == 3:
                        l1 = np.where( mask1d[0] )[0]
                        l2 = np.where( mask1d[1] )[0]
                        l3 = np.where( mask1d[2] )[0]
                        vv[l1[0]:l1[-1]+1,l2[0]:l2[-1]+1,l3[0]:l3[-1]+1] = val
                        del l1,l2,l3
                    elif len(dims) == 4:
                        l1 = np.where( mask1d[0] )[0]
                        l2 = np.where( mask1d[1] )[0]
                        l3 = np.where( mask1d[2] )[0]
                        l4 = np.where( mask1d[3] )[0]
                        vv[l1[0]:l1[-1]+1,l2[0]:l2[-1]+1,l3[0]:l3[-1]+1,l4[0]:l4[-1]+1] = val
                        del l1,l2,l3,l4

                    del val,dims,coords,mask1d

                    gc.collect()

            except StopIteration:
                pass

            idin.close()

            fidout.sync()
            del vv

    except StopIteration:
        pass



def makeStatVariables(fidout,filepref,imax,jmax,dimnames,vars):

    import numpy as np
    import gc
    import sys

    # These are converted to all-sky averages! Can then be renormalized by cfrac in
    # the post-processed data
    incloud = np.array(['Nc_ic','Na_int','SO4_ic','SO4_int',   \
                        'OC_ic','OC_int','BC_ic','BC_int',      \
                        'SS_ic','SS_int','DU_ic','DU_int',      \
                        'NO3_ic','NO3_int','NH3_ic','NH3_int'])
    outcloud = np.array(['Na_oc','SO4_oc','OC_oc','BC_oc','SS_oc','DU_oc','NO3_oc','NH3_oc'])

    if vars == []:
        idin = openr(filepref,(0,0))
        vars = getNames(idin,2)
        idin.close()

    ii = range(imax)
    jj = range(jmax)
    ifn = []
    for i in ii:
        for j in jj:
            ifn.append((i,j))
            
    del ii,jj

    # Weight factor for averaging. For now constant == naive averaging
    wght = 1./np.real(len(ifn))
    
    #print wght

    itfile = iter(ifn)
    # Iterate files
    try:
        while 1:
            fn = itfile.next()
            print 'file ', fn
            idin = openr(filepref,fn)

            # Iterate variables
            itvars = iter(vars)
            try:
                while 1:
                    vn = itvars.next()

                    print vn
                    print '--------'

                    vv = []

                    if vn in dimnames:
                        print 'axis variable, skip'
                        continue   
                    
                    val,dims,coords = getVariable(idin,vn)

                    if len(dims) > 4:
                        continue

                    if fn == (0,0):
                        vv = fidout.createVariable(vn,'d',dims)
                        vv[:] = 0.
                    else:
                        vv = fidout.variables[vn]

                    if np.any(incloud == vn):
                        frac,fdims,fcoords = getVariable(idin,'cfrac')
                        loc = np.where(np.isnan(val))[0]
                        val[loc] = 0.
                        vv[:] = vv[:] + val*frac*wght
                        del frac,fdims,fcoords,loc
                    #elif np.any(outcloud == vn):
                    #    frac,fdims,fcoords = getVariable(idin,'cfrac')
                    #    loc = np.where(np.isnan(val))[0]
                    #    val[loc] = 0.
                    #    vv[:] = vv[:] + val*(1.-frac)*wght
                    #    del frac,fdims,fcoords,loc
                    else:
                        vv[:] = vv[:] + wght*val


                    del val,dims,coords

                    gc.collect()

            except StopIteration:
                pass

            idin.close()

            fidout.sync()
            del vv

    except StopIteration:
        pass



def makeDimensions(fidout,filepref,imax,jmax):

    import numpy as np
    import gc

    idin = openr(filepref,(0,0))
    dimnames = getNames(idin,1)
    itdims = iter(dimnames)
    idin.close()

    ii = range(imax)
    jj = range(jmax)
    ifn = []
    for i in ii:
        for j in jj:
            ifn.append((i,j))

    del ii,jj

    try:
        while 1:
            dim = itdims.next()

            print dim
            print '------'

            dst_val = np.array([])
            dst_size = 0

            itfile = iter(ifn)
            try:
                while 1:
                    fn = itfile.next()
                    print 'file ',fn
                    
                    if dim == 'time' and fn != (0,0):
                        break

                    idin = openr(filepref,fn)
                    val,dlen,dvfound = getDimension(idin,dim)
                    idin.close()
                    
                    print '1'

                    if dvfound:
                        # If the values are already in the destination array, do nothing
                        # Otherwise append
                        if np.any(np.in1d(val,dst_val)):
                            continue

                        else:
                            dst_val = np.append(dst_val,val)
                            dst_size = dst_size + dlen

                    else:
                        # No corresponding axis variable for current dimensions
                        dst_size = dlen
                        continue

                    print '2'

                    del val,dlen
                    gc.collect()

            except StopIteration:
                pass
                
            fidout.createDimension(dim,dst_size)
            if dvfound:
                vv = fidout.createVariable(dim,'d',(dim,))
                vv[:] = dst_val

                fidout.sync()
                del vv
            
            del dst_val,dst_size,itfile

    except StopIteration:
        pass


    return dimnames

import os
import sys

filepref = sys.argv[1]
#filetype = sys.argv[1]


# Check that at least one file exists
if not os.path.isfile('%s.%04u%04u.nc' % (filepref, 0, 0)): raise Exception('Data not found!')

# Limits for the data file indices
imax,jmax=0,0
while imax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (filepref, imax, 0) ):
    imax+=1
    #
while jmax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (filepref, 0, jmax) ):
    jmax+=1


#varbls = ['u','v','w','t','p','q','l','r','S_RH','S_Na','S_Nb','S_Nc','S_Np']


# Open outputfile
fidout = openw(filepref)
dimnames = makeDimensions(fidout,filepref,imax,jmax)
fidout.close()
if filepref[-3:] == '.ts' or filepref[-3:] == '.ps':
    fidout = opena(filepref)
    makeStatVariables(fidout,filepref,imax,jmax,dimnames,[])
else:
    fidout = opena(filepref)
    makeVariables(fidout,filepref,imax,jmax,dimnames,[])


fidout.close()
