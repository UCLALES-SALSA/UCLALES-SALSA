
import netCDF4 as nc
import numpy as np
import gc
import sys
latdims = ['xt','xm','yt','ym']
skip = ['u0','v0','dn0']

class NCFile(object):

    def __init__(self,ifile,mode="r"):
        self.ncid = nc.Dataset(ifile,mode)

    def NcClose(self):
        self.ncid.close()


class NcInput(NCFile,object):

    def __init__(self,ifile,useBinned=True,dimNames=[],varNames=[]):
        super(NcInput,self).__init__(ifile,mode="r")
        
        # Inputs: ifile - input filename (full path)
        #         varnames - List of variable names to be read
        #         dimnames - List of dimension names to be read 
        
        self.dimensions = self._readData(dimNames,useBinned,dims=True)

        self.data = self._readData(varNames,useBinned,dims=False)

    
    def _readData(self,names,useBinned,dims=False):

        values = dict()

        if dims:
            if len(names) == 0:
                #If no names specified, read all dimensions
                names = self.ncid.dimensions.keys()
                
            for nn in names:
                values.update(self._fetchData(nn))

        elif not dims:
            
            if len(names) == 0:
                # If no names specified, read all variables
                names = self.ncid.variables.keys()
                # Remove dimension names from this list 
                rem = self.ncid.dimensions.keys()
                for r in rem: names.remove(r)
                
                # Remove binned variables if needed
                if not useBinned:
                    rem = []
                    for nn in names:
                        tmp = self.ncid.variables[nn]
                        if len(tmp.dimensions) > 4:
                            rem.append(nn)
                    
                    for r in rem: names.remove(r)
                    del rem, tmp
                    gc.collect()

            # Fetch the listed variables
            for nn in names:
                entry = self._fetchData(nn)
                values.update(self._fetchData(nn))
                
        return values

        
    def _fetchData(self,var):
                
        tmp = self.ncid.variables[var]
        return {var:Variable(tmp[:],tmp.dimensions)}


class NcNew(NCFile,object):

    def __init__(self,filename,dim,var,nFX,nFY):
        super(NcNew,self).__init__(filename,mode="w")
    
        self._makeDimensions(dim,nFX,nFY)
        self._makeVariables(var)
        
        self.ncid.close()

    def _makeDimensions(self,dim,nFX,nFY):
        #
        # Create dimensions that are the same for all raw data files,
        # i.e. everything except lateral dimensions
        #
        dimiter = dim.iteritems()
        while True:
            try:
                dimitem = dimiter.next()
                
                if dimitem[0] in latdims:
                    print("creating dimension",dimitem[0])
                    # In this case, the length of the dimension 
                    # has to be multiplied by the number of blocks
                    if dimitem[0] in latdims[0:2]:
                        length = dimitem[1].values.size*nFX
                    elif dimitem[0] in latdims[2:]:
                        length = dimitem[1].values.size*nFY
                    self.ncid.createDimension(dimitem[0],length)
                    self.ncid.createVariable(dimitem[0],'d',dimensions=dimitem[1].dims)

                    # Make the lateral axes by calculating from number of points and resolution
                    print("Writing lateral axes:",dimitem[0])
                    print ""
                    d0 = dimitem[1].values[0]
                    dd = np.diff(dimitem[1].values)[0]
                    axes = self._makeLateralAxes(d0,dd,length)
                    self._putInitialValues([dimitem[0],axes])

                else:
                    print("creating dimension",dimitem[0])
                    self.ncid.createDimension(dimitem[0],dimitem[1].values.size)
                    self.ncid.createVariable(dimitem[0],'d',dimensions=dimitem[1].dims)
                    
                    print("writing initial values for",dimitem[0])
                    print ""
                    self._putInitialValues(dimitem)
                    
            except StopIteration:
                break


    def _makeLateralAxes(self,x0,dx,nx):

        # x0: first x coordinate
        # dx: resolution
        # nx: number of points per file

        tmp = np.arange(x0,x0+nx*dx,dx)
        axes = Variable(tmp,())

        return axes


    def _makeVariables(self,var):

        variter = var.iteritems()
        while True:
            try:
                varitem = variter.next()
                if varitem[0] in skip: continue

                print "creating variable", varitem[0]
                self.ncid.createVariable(varitem[0],'f',dimensions=varitem[1].dims)

            except StopIteration:
                break
        print ""
                
    def _putInitialValues(self,var):

        # Used for putting initial values for dimensions
        self.ncid.variables[var[0]][:] = var[1].values
        
class NcExtend(NCFile,object):
    def __init__(self,filename):
        super(NcExtend,self).__init__(filename,mode="a")

    def processData(self,dim,var):

        borders = self._findTargetSpace(dim)
        print "------------------------"
        print "[xmin,xmax,ymin,ymax] (indices) for current block:",borders
        print "------------------------"
        print ""
        self._writeBlockData(var,borders)
        del borders
        del dim
        del var
        gc.collect()


    def _findTargetSpace(self,dim):

        # Find the indices for x and y ranges convered by the input block
        # xstart,xend,ystart,yend
        borders = [0,0,0,0]
        
        xtlo = dim['xt'].values[0]
        xthi = dim['xt'].values[-1]
        xtout = self.ncid.variables['xt'][:]

        print "Find location for current subdomain:"
        print "-----------------"
        print "Border locations in X-direction",xtlo,xthi
        print "within"
        print xtout
        print ""

        borders[0] = np.where(xtout == xtlo)[0][0]
        borders[1] = np.where(xtout == xthi)[0][0]
        
        ytlo = dim['yt'].values[0]
        ythi = dim['yt'].values[-1]
        ytout = self.ncid.variables['yt'][:]

        print "Border locations in Y-direction",ytlo,ythi
        print "within"
        print ytout
        print "-----------------"
        print ""

        borders[2] = np.where(ytout == ytlo)[0][0]
        borders[3] = np.where(ytout == ythi)[0][0]

        return borders

    def _writeBlockData(self,var,borders):
 
        variter = var.iteritems()
        while True:
            try:
                varitem = variter.next()
                if varitem[0] in skip: continue

                print "Writing block values for",varitem[0]
                print ""
                tmp = self.ncid.variables[varitem[0]]
                if len( tmp.dimensions ) == 4:
                    tmp[:,                       \
                        borders[2]:borders[3]+1, \
                        borders[0]:borders[1]+1, \
                        :] = varitem[1].values
                elif len( tmp.dimensions ) == 5:
                    tmp[:,                       \
                        :,                       \
                        borders[2]:borders[3]+1, \
                        borders[0]:borders[1]+1, \
                        :] = varitem[1].values
                
            except StopIteration:
                break

        self.ncid.sync()
        gc.collect()

class Variable(object):

    def __init__(self,vals,dims):
        self.values = vals
        self.dims = dims # This will be tuple of dimensions names
