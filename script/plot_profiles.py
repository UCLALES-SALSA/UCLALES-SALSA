# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 15:54:25 2016

@author: aholaj
"""
##########################################
###                                    ###
### PURPOSE OF THE SCRIPT              ###
### plot profiles from sound_in files  ###
###                                    ###
##########################################

import numpy as np

import matplotlib.pyplot as plt

class PlotProfiles:

    def __init__(self, file, subfolder, folder = '/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/' ):
        self.folder=folder
        self.file = file
        self.subfolder = subfolder
        #file = 'sound_in_DYCOMSIIRF02'
        #file = 'sound_in'
        filu = self.folder + self.subfolder + self.file
        f = open(filu, 'r')
    
        self.z = []
        self.t = []
        self.q = []
        self.u = []
        self.v = []
        
        
        
        
        for line in f:
            zA, tA, qA, uA, vA = line.split()
            self.z.append(float(zA))
            self.t.append(float(tA))
            self.q.append(float(qA))
            self.u.append(float(uA))
            self.v.append(float(vA))
        
            
        f.close()
        self.z[0]=0.
        
        self.zu = np.column_stack(( self.z, self.u ))
        self.zv = np.column_stack(( self.z, self.v ))
    
    def getU(self):
        return self.u
    
    def getV(self):
        return self.v

    def getZU(self):
        return self.zu

    def getZV(self):
        return self.zv

    def returnWindAppr(self, height, wind):
        found = False
        i = 0
        indexUpper=0
        while ( i < len(self.z) and (not found) ) :
            if self.z[i] > height:
                found = True
                indexUpper = i
            i += 1
        
        found = False
        i = len(self.z)-1
        indexLower=0
        while ( i >= 0 and (not found) ) :
            if self.z[i] < height:
                found = True
                indexLower = i
            i -= 1
#        print 'indexLower ' + str(indexLower)
#        print 'indexUpper ' + str(indexUpper)
        if ( indexUpper - indexLower == 2):
            WindAppr = wind[indexLower+1]
        else:
            WindAppr = ( wind[indexUpper]-wind[indexLower] )/( self.z[indexUpper] - self.z[indexLower] )*( height-self.z[indexLower] ) + wind[indexLower]
#        print 'WindAppr' + str(WindAppr)
        return WindAppr
    
    def returnUAppr(self, height):
        if (height <= 3000.):
            u = self.returnWindAppr( height,self.u)
        else:
            u = 0. #10.*np.random.random()
        return u

    def returnVAppr(self, height):
        if (height <= 3000.):
            v = self.returnWindAppr( height,self.v)
        else:
            v = 0. # -10.*np.random.random()
        return v
                    
            

    def plot(self):
    
        
        
        plt.figure()
        plt.plot( self.t, self.z )
        
        plt.figure()
        plt.plot( self.q, self.z )
        
        plt.figure()
        plt.plot( self.u, self.z )
        
        print 'u'
        plt.figure()
        plt.plot( self.u, self.z )
        
        plt.show()



#def giveU():
#    return
#
#
#if __name__ == "__main__":
#    main()
#print z