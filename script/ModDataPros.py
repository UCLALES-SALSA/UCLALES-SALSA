#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 13:50:26 2016

@author: aholaj
"""

###### read ice

#################################
### imports and global values ###
#################################


import numpy as np
import sys
import os

import matplotlib.pyplot as plt
from itertools import cycle
import matplotlib.patches as mpatches


#################################
### subroutines               ###
#################################

#################################
### time steps to strings     ###
#################################
def tstepH(time):
  return str(time*360/3600.)

#################################
### set values >= 10e-10      ###
#################################
def nollaa(data):
  for i in xrange(np.shape(data)[0]):
    if data[i]<10e-10:
      data[i]=10e-10  
  
  return data

#################################
### read variable from data   ###
#################################  
def read_Data(filename,var):
  from netCDF4 import Dataset
  fileH = Dataset(filename,mode='r')
  data = fileH.variables[var][:]
  fileH.close()
  return data

def read_NamelistValue( filename = 'NAMELIST', var = 'Tspinup' ):
    import re
    f = open( filename )
    value = None
    for line in f:
        
        aa = line.split('=')
        k = 0
        nimi = ''
        arvo = ''
        for sentence in aa:
            if   np.mod(k,2) == 0:
                nimi = re.sub(r"\s+", "", sentence, flags=re.UNICODE)
            elif np.mod(k,2) == 1:
                kakka = re.sub(r"\s+", "", sentence, flags=re.UNICODE)
                arvo  = kakka.split('!', 1)[0]
            k+=1
         
        #print 'nimi', nimi, 'arvo', arvo
        if nimi == var:
            value = arvo
            break
        


       
        
    f.close()
    return float(value)


######################################################
### count min and max values of 5 dimensional data ###
### time x y z bin                                 ###
######################################################
def laske_minimi_maksimi_bini(nimi,data):
  print ' '
  print nimi
  maxi =  data[ 0, 0, 0, 0, 0 ]
  mini =  data[ 0, 0, 0, 0, 0 ]
  maxikoord = [ 0, 0, 0, 0, 0 ]
  minikoord = [ 0, 0, 0, 0, 0 ]
  
  for z in  xrange( np.shape(data)[4] ):
    for y in xrange( np.shape(data)[3] ):
        for x in xrange( np.shape(data)[2] ):
            for bini in xrange( np.shape(data)[1] ):
                for time in xrange( np.shape(data)[0] ):
		  
		  if (data[ time, bini, x, y, z ] > maxi):
		    maxi =  data[ time, bini, x, y, z ] 
		    maxikoord = [ time, bini, x, y, z ]
		    
		  if (data[ time, bini, x, y, z ] < mini):
		    mini     = data[ time, bini, x, y, z ]
		    minikoord =    [ time, bini, x, y, z ]
		    
  print 'maxi ' + str(maxi)
  print 'mini ' + str(mini)
  print 'maxikoord' + str(maxikoord)
  print 'minikoord' + str(minikoord)
  print ' '

######################################################
### count min and max values of 4 dimensional data ###
### time x y z                                     ###
######################################################  
def laske_minimi_maksimi(nimi,data):
  print ' '
  print nimi
  maxi= data[0,0,0,0]
  mini= data[0,0,0,0]
  maxikoord=[0,0,0,0]
  minikoord=[0,0,0,0]
  for z in  xrange(np.shape(data)[3]):
    for y in xrange(np.shape(data)[2]):
        for x in xrange(np.shape(data)[1]):
	  for time in xrange(np.shape(data)[0]):
	    if (data[time,x,y,z] > maxi):
	      maxi = data[time,x,y,z]
	      maxikoord=[time,x,y,z]
	    if (data[time,x,y,z] < mini):
	      mini = data[time,x,y,z]
	      minikoord=[time,x,y,z]
  print 'maxi ' + str(maxi)
  print 'mini ' + str(mini)
  print 'maxikoord' + str(maxikoord)
  print 'minikoord' + str(minikoord)		  
  print ' '

##########################################################
### handle timeseries data by summing them columnwise  ###
### 4 dimensions: time, x, y, z                        ###
##########################################################
def vertical_sum_timeseries( nimi, data, aika, tulostus = False, piirra = False, uusikuva = True, label = None ):
    print ' '
    print nimi, 'aikasarja'
    kokoaika = np.shape(data)[0]
    aikasarja = np.zeros(kokoaika)

    for time in xrange(np.shape(data)[0]):
        for z in  xrange(np.shape(data)[3]):
            for y in xrange(np.shape(data)[2]):
                for x in xrange(np.shape(data)[1]):
                    aikasarja[time] = aikasarja[time] + data[time,x,y,z]
            
    if tulostus:
        for time in xrange(kokoaika):
            print 'ajanhetki: ' + str(aika[time]) + ' ' + ' arvo : ' + str(aikasarja[time])
    print ' '

    ## drawing ##

    uusikuva = ( piirra and uusikuva )
    plot_alustus() if uusikuva else False
    plottaa( aika, timeSeries, nimi, 'aika', 'path', label = label)  if piirra else False   

##########################################################
### handle timeseries data by summing them columnwise  ###
### 5 dimensions: time, x, y, z, bin                   ###
##########################################################    
def vertical_sum_timeseries_bini( nimi, data, aika, tulostus = False, piirra = False, uusikuva = True, label = None ):
    print ' '
    print nimi, 'aikasarja'
    kokoaika = np.shape(data)[0]
    aikasarja = np.zeros(kokoaika)

    for time in xrange(np.shape(data)[0]):
        for z in  xrange(np.shape(data)[4]):
            for y in xrange(np.shape(data)[3]):
                for x in xrange(np.shape(data)[2]):
                    for bini in xrange(np.shape(data)[1]):
                        aikasarja[time] = aikasarja[time] + data[time,bini,x,y,z]
	  
    if tulostus:
        for time in xrange(kokoaika):
            print 'ajanhetki: ' + str(aika[time]) + ' ' + ' arvo : ' + str(aikasarja[time])
    print ' '
  
    ## drawing ##

    uusikuva = ( piirra and uusikuva )
    plot_alustus() if uusikuva else False
    plottaa( aika, timeSeries, nimi, 'aika', 'path', label = label)  if piirra else False   


##########################################################
### print filename of data file                        ###
###                                                    ###
##########################################################
def print_filename(fname):
  head, tail = os.path.split(fname)
  print ' '
  print 'file: ' + tail  

##########################################################
### print shape of data file                           ###
###                                                    ###
##########################################################
def print_shape( var, data ):
  print ' '
  print 'variable: ' + var
  print 'shape var1: '+ str(np.asarray(np.shape(data)))



##########################################################
### return area of the domain                          ###
###                                                    ###
##########################################################    
def area( xm_data, ym_data ):
    x_dim = int(np.asarray(np.shape(xm_data)))
    y_dim = int(np.asarray(np.shape(ym_data)))
    x_size = 0.
    y_size = 0.;
    if ( x_dim == 1 ):
      x_size = xm_data[0]*2
    else:
      x_size = xm_data[ x_dim-1 ]*2
    if ( y_dim == 1 ):
      y_size = ym_data[0]*2
    else:
      y_size = ym_data[ y_dim-1 ]*2
      
    return x_size*y_size

##########################################################
### calculate vertical path of a variable according    ###
### to the air density (g/m^2)                         ###
###                                                    ###
### input variable: mix ratio ( kg / kg )              ###
### 4 dimensions: time, x, y z                         ###
##########################################################
def laske_path_aikasarjaXYZ( mixRatData, dn0_data, korkeus, aika = None, muunnosKerroin = 1000.,  onlyCloudy = False, tulostus = False, piirra = False, uusikuva = True, nimi = 'path aikasarja', xlabel = 'aika [s]', tightXAxis=False, label = None):
    
    #fig.laske_path_aikasarjaXYZ = None
    #ax.laske_path_aikasarjaXYZ  = None
    
    mixRatData = mixRatData*muunnosKerroin # kg/kg -> g/kg

    timeDim  = np.shape( mixRatData )[0]
    xDim     = np.shape( mixRatData )[1]
    yDim     = np.shape( mixRatData )[2]
    zDim     = np.shape( mixRatData )[3]
    nCol     = xDim * yDim

    dn0Kork  = dn0_data * korkeus

    onlyCloudyTXY = np.zeros( ( timeDim, xDim, yDim ) ) # 
    onesTXY       = np.ones(  ( timeDim, xDim, yDim ) )
    timeSeriesTXY = np.zeros( ( timeDim, xDim, yDim ) )
    timeSeries    = np.zeros(   timeDim             ) 

    #[a[i] < b[i] for i in range(5)]
    for t in  xrange(timeDim):
        for i in   xrange(xDim):
            for j in xrange(yDim):
                timeSeriesTXY[t, i, j] = np.dot( mixRatData[ t, i, j, : ], dn0Kork )
                if ( onlyCloudy and timeSeriesTXY[t, i, j] > 0.0 ):
                    onlyCloudyTXY[t, i, j] = 1.0

    timeSeries  = np.sum( np.sum( timeSeriesTXY, axis = 1), axis = 1 )
    onlyCloudyT = np.sum( np.sum(onlyCloudyTXY , axis = 1), axis = 1 )

    if onlyCloudy:
        timeSeries = np.where( onlyCloudyT > 0.0, timeSeries / onlyCloudyT , 0.0 )
    else:
        timeSeries = timeSeries / np.sum( np.sum( onesTXY,       axis = 1) , axis = 1 )




    if tulostus:
        print 'dimensiot aikasarjaXYZ'+ str(np.shape(aika))+ ' timeseries '+ str(np.shape(timeSeries))
        print ' '
        print nimi
        for t in xrange(timeDim):
            print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(timeSeries[t])
        print ' '    


        ## drawing ##

    uusikuva = ( piirra and uusikuva )
    
    if uusikuva:
        laske_path_aikasarjaXYZ.fig, laske_path_aikasarjaXYZ.ax = plot_alustus()
        
    plottaa( aika, timeSeries, nimi, xlabel , 'path [g/m^2]', tightXAxis = tightXAxis, label = label)  if piirra else False   
    
    
 #   if uusikuva:
    return laske_path_aikasarjaXYZ.fig, laske_path_aikasarjaXYZ.ax
    #else:
        #return None, None

##########################################################
### calculate vertical path of a variable according    ###
### to the air density (g/m^2)                         ###
###                                                    ###
### input variable: mix ratio ( kg / kg )              ###
### 2 dimensions: time, z                              ###
##########################################################
def laske_path_aikasarjaZ( mixRatData, dn0_data, korkeus, aika, tulostus = False, piirra = False, uusikuva = True, nimi = 'path aikasarja', label = None ):
  
    print ' '

    mixRatData = mixRatData * 1000.0 # kg/kg -> g/kg # ( timeDim, zdim ) 

    timeDim  = np.shape( mixRatData )[0]
    zDim     = np.shape( mixRatData )[1]

    dn0Kork  = dn0_data * korkeus

    timeSeriesTZ = np.zeros( ( timeDim, zDim ) )
    timeSeries   = np.zeros(   timeDim         ) 

    timeSeries   = np.dot(mixRatData, dn0Kork)

    print 'timeDim '+ str(timeDim)
    if ( np.shape( aika )[0] - np.shape( timeSeries )[0] == 1 ):
        timeSeries = np.insert( timeSeries, 0, 0 )
    elif ( np.shape( timeSeries )[0] - np.shape( aika )[0] == 1):
        aika       = np.insert( aika, 0, 0 )
    elif ( np.abs( np.shape( timeSeries )[0] - np.shape( aika )[0] ) > 1): 
        sys.exit( "something went really wrong with dimensions in laske_path_aikasarjaZ()" )
        

    if tulostus:
        print ' '
        print nimi
        for t in xrange(timeDim):
            print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(timeSeries[t])
        print ' '    

        print 'dimensiot aikasarjaZ'+ str(np.shape(aika))+ ' timeseries '+ str(np.shape(timeSeries))

        ## drawing ##

    uusikuva = ( piirra and uusikuva )
    plot_alustus() if uusikuva else False
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'path [g/m^2]', label = label)  if piirra else False   

##########################################################
### calculate content of a variable according          ###
### to the air density (g/m^3)                         ###
###                                                    ###
### input variable: mix ratio ( kg / kg )              ###
### 4 dimensions: time, x, y z                         ###
##########################################################
def laske_WC( mixRatData, dn0_data, aika, tulostus = False, piirra = False, uusikuva = True, nimi = 'water content aikasarja', label = None ):
  
    mixRatData = mixRatData*1000.0 # kg/kg -> g/kg

    timeDim  = np.shape( mixRatData )[0]
    xDim     = np.shape( mixRatData )[1]
    yDim     = np.shape( mixRatData )[2]
    zDim     = np.shape( mixRatData )[3]
    nCol     = xDim * yDim  

    onlyCloudyTXYZ = np.zeros( ( timeDim, xDim, yDim, zDim ) ) #

    timeSeriesTXYZ =  np.zeros( ( timeDim, xDim, yDim, zDim ) )
    timeSeries     =  np.zeros(   timeDim             ) 

    timeSeriesTXYZ = np.multiply( mixRatData, dn0_data )

    onlyCloudyTXYZ = np.where( mixRatData > 0.0, 1.0, 0.0)

    timeSeries  = np.sum( np.sum( np.sum( timeSeriesTXYZ, axis = 1), axis = 1), axis = 1 )
    onlyCloudyT = np.sum( np.sum( np.sum( onlyCloudyTXYZ, axis = 1), axis = 1), axis = 1 )

    #print np.shape(timeSeries)
    #print np.shape(onlyCloudyT)

    timeSeries = np.where( onlyCloudyT > 0.0, timeSeries / onlyCloudyT , 0.0 )


    if tulostus:
        print ' '
        print nimi
        for t in xrange(timeDim):
            print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(timeSeries[t])
    print ' '    

    ## drawing ##

    uusikuva = ( piirra and uusikuva )
    plot_alustus() if uusikuva else False
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'Content [g/m^3]', label = label)  if piirra else False    

#######################################################################
### calculate number concentration of a variable according          ###
### to the reference data                                           ###
###                                                                 ###
### input variable: mix ratio ( kg / kg )                           ###
### 4 dimensions: time, x, y z                                      ###
#######################################################################

############## NEEDS REVISION ############################
def laske_NumberConcentration( Ndata, refNdata, dn0_data, aika, tulostus = False, uusikuva = True, piirra =False, nimi = 'number concentration', label = None ):
    
    #dn0_data = dn0_data/1000.0 # kg/m^3 -> kg/l
    
    timeDim  = np.shape( Ndata )[0]
    xDim     = np.shape( Ndata )[1]
    yDim     = np.shape( Ndata )[2]
    zDim     = np.shape( Ndata )[3]
    nCol     = xDim * yDim  
  
    onlyCloudyTXYZ = np.zeros( ( timeDim, xDim, yDim, zDim ) ) #
  
    timeSeriesTXYZ =  np.zeros( ( timeDim, xDim, yDim, zDim ) )
    timeSeries     =  np.zeros(   timeDim             ) 
  
    #timeSeriesTXYZ = np.multiply( Ndata, dn0_data )
    timeSeriesTXYZ = np.where( Ndata/1000.0    > 1e-10, Ndata, 0.0 )
    onlyCloudyTXYZ = np.where( refNdata > 1e-10, 1.0  , 0.0 )
    
    
    timeSeries  = np.sum( np.sum( np.sum( timeSeriesTXYZ, axis = 1), axis = 1), axis = 1 )
    onlyCloudyT = np.sum( np.sum( np.sum( onlyCloudyTXYZ, axis = 1), axis = 1), axis = 1 )
  

  
    timeSeries = np.where( onlyCloudyT > 0.0, timeSeries / onlyCloudyT , 0.0 )
    
    
    
    if tulostus:
        
        print ' '
        print nimi
        for t in xrange(timeDim):
            print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(timeSeries[t])
        print ' '    

    ## drawing ##
      
    uusikuva = ( piirra and uusikuva )
    plot_alustus() if uusikuva else False
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'Number concentration [#/kg]', label = label)  if piirra else False
        

#######################################################################
### calculate mean diameter in a bin                                ###
### according to the reference data                                 ###
###                                                                 ###
### input variable: mix ratio ( kg / kg )                           ###
### 5 dimensions: time, x, y z, bin                                 ###
#######################################################################
def laske_MeanDiameterInBin( RadiusBinData, bini, refNdata, aika, tulostus = False, piirra =False, uusikuva = True, nimi = 'Mean diameter in bin ', label = None ):
    biniNimi = str(bini+1)
    nimi     = nimi + biniNimi
    
    timeDim  = np.shape( RadiusBinData )[0]
    binDim   = np.shape( RadiusBinData )[1]
    xDim     = np.shape( RadiusBinData )[2]
    yDim     = np.shape( RadiusBinData )[3]
    zDim     = np.shape( RadiusBinData )[4]
    
    nCol     = xDim * yDim  
  
    onlyCloudyTXYZ = np.zeros( ( timeDim, xDim, yDim, zDim ) ) #
  
    timeSeriesTXYZ = np.zeros( ( timeDim, xDim, yDim, zDim ) )
    timeSeries     = np.zeros(   timeDim             ) 
  
    timeSeriesTXYZ =  2.0 * RadiusBinData[ :, bini, :, : , : ]*1e6 # select only one bin and change to diameter in um
    
    onlyCloudyTXYZ = np.where( refNdata > 1e-10, 1.0, 0.0)
    
    timeSeries  = np.sum( np.sum( np.sum( timeSeriesTXYZ, axis = 1), axis = 1), axis = 1 )
    onlyCloudyT = np.sum( np.sum( np.sum( onlyCloudyTXYZ, axis = 1), axis = 1), axis = 1 )
  

  
    timeSeries = np.where( onlyCloudyT > 0.0, timeSeries / onlyCloudyT , 0.0 )

    
    if tulostus:
        
        print ' '
        print 'bini: ' + biniNimi
        print nimi
        for t in xrange(timeDim):
            print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(timeSeries[t])
        print ' '    

    ## drawing ##
      
    uusikuva = ( piirra and uusikuva )
    plot_alustus() if uusikuva else False
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'diameter [um]', label = label )  if piirra else False

#######################################################################
### calculate column mean PSD divided into bins                     ###
### at a spesific time                                              ###
###                                                                 ###
### input variable: mix ratio ( kg / kg )                           ###
### 5 dimensions: time, x, y z, bin                                 ###
#######################################################################

def laske_PSD_TimestepColumnAverage( data, aikaPisteet = 0, tulostus = False, piirra =False, uusikuva = True, nimi = 'PSD ', label = None ):
    #bini in xrange( np.shape(S_Rwiba_data)[1] ):
    biniNimi = str(bini+1)
    nimi     = nimi + biniNimi
    
    timeDim  = np.shape( data )[0]
    binDim   = np.shape( data )[1]
    xDim     = np.shape( data )[2]
    yDim     = np.shape( data )[3]
    zDim     = np.shape( data )[4]
    
    if not isinstance( aikaPisteet, np.ndarray):
        aikaPisteet = np.asarray( aikaPisteet )
    
    aikaPisteetDim = np.shape( aikaPisteet)[0]
    
    nXY      = xDim * yDim  
    
    nXYZ     = nXY*zDim
  
    aikaBinData  = np.zeros(( aikaPisteetDim, binDim ))
    
   
    for aika in aikaPisteet:
        for bini in binDim:
            aikaBinData = data( aika, bini)

    #timeSeriesTXYZ = np.zeros( ( timeDim, xDim, yDim, zDim ) )
    #timeSeries     = np.zeros(   timeDim             ) 
  
    timeSeriesTXYZ =  2.0 * data[ :, bini, :, : , : ]*1e6 # select only one bin and change to diameter in um
    if isinstance( aikaPisteet, np.ndarray):
      for t in aikaPisteet:
          print 'o'
    #onlyCloudyTXYZ = np.where( refNdata > 1e-10, 1.0, 0.0)
    
    timeSeries  = np.sum( np.sum( np.sum( timeSeriesTXYZ, axis = 1), axis = 1), axis = 1 )
    #onlyCloudyT = np.sum( np.sum( np.sum( onlyCloudyTXYZ, axis = 1), axis = 1), axis = 1 )
  

  
    #timeSeries = np.where( onlyCloudyT > 0.0, timeSeries / onlyCloudyT , 0.0 )
    #print ' '
    #print nimi, 'aikasarja'
    #kokoaika = np.shape(data)[0]
    #aikasarja = np.zeros(kokoaika)

    #for time in xrange(np.shape(data)[0]):
        #for z in  xrange(np.shape(data)[4]):
            #for y in xrange(np.shape(data)[3]):
                #for x in xrange(np.shape(data)[2]):
                    #for bini in xrange(np.shape(data)[1]):
                        #aikasarja[time] = aikasarja[time] + data[time,bini,x,y,z]
	  
    #if tulostus:
        #for time in xrange(kokoaika):
            #print 'ajanhetki: ' + str(aika[time]) + ' ' + ' arvo : ' + str(aikasarja[time])
    #print ' '
  
    ### drawing ##

    #uusikuva = ( piirra and uusikuva )
    #plot_alustus() if uusikuva else False
    #plottaa( aika, timeSeries, nimi, 'aika', 'path')  if piirra else False   
    
    if tulostus:
        
        print ' '
        print 'bini: ' + biniNimi
        print nimi
        for t in xrange(timeDim):
            print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(timeSeries[t])
        print ' '    

    ## drawing ##
      
    uusikuva = ( piirra and uusikuva )
    plot_alustus() if uusikuva else False
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'diameter [um]', label = label)  if piirra else False

########################################
### calculate root mean square error ###
###                                  ###
########################################
def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

###########################################
### print/draw timeseries of a variable ###
### most useful with .ts.nc files       ###
###                                     ###
###########################################
def aikasarjaTulostus( data, aika = 0, tulostus = False, piirra = False, uusikuva = True, nimi = 'aikasarja', xnimi = 'x-akseli', ynimi= 'y-akseli', changeColor=True, tightXAxis=False, LEGEND=True, omavari = False, label = None ):
  #fig.aikasarjaTulostus = None
  #ax.aikasarjaTulostus  = None
  
  if not isinstance(aika, np.ndarray):
    aika=np.zeros( (np.shape(data)[0]))
  
  if tulostus:
      print ' '
      print nimi
      for t in xrange( np.shape(data)[0] ):
         print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(data[t])
      print ' '    

## drawing ##
      
  uusikuva = ( piirra and uusikuva )
  
  #if uusikuva:
      #aikasarjaTulostus.fig, aikasarjaTulostus.ax = plot_alustus()
      
  aikasarjaTulostus.fig, aikasarjaTulostus.ax = plottaa( aika, data, nimi, xnimi, ynimi, changeColor = changeColor, tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = omavari, label = label, uusikuva = uusikuva )

  return aikasarjaTulostus.fig, aikasarjaTulostus.ax

###########################################
### print/draw timeseries of a variable ###
### most useful with .ts.nc files       ###
###                                     ###
###########################################
def profiiliTulostus( data, aikaPisteet = 0, korkeus=0,label = None, tulostus = False, piirra = False, uusikuva = True, nimi = 'profiili', xnimi = 'x-akseli', ynimi= 'y-akseli', changeColor=True, tightXAxis=False, tightYAxis=False, LEGEND=True, omavari = False, loc = 2):
  
  if not isinstance(korkeus, np.ndarray):
    korkeus=np.arange( (np.shape(data)[0]) )
  
  if tulostus:
      print ' '
      print nimi
      for t in xrange( np.shape(data)[0] ):
         print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(data[t,:])
      print ' '    

## drawing ##
  uusikuva = ( piirra and uusikuva )

  if uusikuva:
      profiiliTulostus.fig, profiiliTulostus.ax = plot_alustus()
  
  if isinstance( aikaPisteet, np.ndarray) or isinstance( aikaPisteet, list):
      for t in aikaPisteet:
        plottaa( data[t,:], korkeus, nimi, xnimi, ynimi, label = label, changeColor = changeColor, tightXAxis=tightXAxis, tightYAxis=tightYAxis, LEGEND=LEGEND, omavari = omavari, loc = loc )
  else:
      plottaa( data, korkeus, nimi, xnimi, ynimi, label = label, changeColor = changeColor, tightXAxis=tightXAxis, tightYAxis=tightYAxis, LEGEND=LEGEND, omavari = omavari, loc = loc )
          

  return profiiliTulostus.fig, profiiliTulostus.ax
####################################
###                              ###
### convert None string to empty ###
###                              ###
####################################
def xstr(s):
    if s is None:
        return ''
    return str(s)

########################################
### colorPool object class           ###
###                                  ###
########################################  
class colorPool:
    
    def __init__( self, colorNumber, shuffling = False ):
        colorMap = plt.cm.gist_ncar
        if colorNumber > 7:
            
            listing = [colorMap(i) for i in np.linspace(0, 0.95, colorNumber)]
            
            if shuffling:
                np.random.shuffle(listing)
                
            self.colors    = cycle( listing )
        else:
            self.colors = cycle( ['r','b','k','g','c','m','y'][:colorNumber] )

        self.currentColor = next(self.colors)

    def getCurrentColor(self):
        return self.currentColor
    
    def getNextColor(self):
        self.currentColor = next(self.colors)
        return self.currentColor

########################################
### colorPool object                 ###
### needs to be called if plots are  ###
### drawn                            ###
########################################  

def initializeColors(colorNRO=6, shuffling = False):
  global colorChoice
  colorChoice = colorPool(colorNRO, shuffling)

  return colorChoice

########################################
### add a vertical line              ###
###                                  ###
########################################   
def plot_vertical( x, vari = 'k', viivatyyli = '--' ):
    plt.axvline( x, color = vari , linestyle = viivatyyli )


########################################
### add a horizontal line            ###
###                                  ###
########################################   
def plot_horizontal( y, vari = 'k', viivatyyli = '--' ):
    plt.axhline( y, color = vari , linestyle = viivatyyli )



########################################
### initialize a new figure          ###
###                                  ###
########################################   
def plot_alustus(a=24,b=15, uusifig = True, uusisub = True, sub_i = 1, sub_j = 1, sub_k = 1, figuuri = None):
    
  if uusifig:
    plot_alustus.fig = plt.figure( figsize = (a,b) )
  if uusisub:
    if figuuri is None:
        kuva = plot_alustus.fig
    else:
        kuva = figuuri
    plot_alustus.ax  = kuva.add_subplot( sub_i, sub_j, sub_k )

  
  return kuva, plot_alustus.ax

  
########################################
### show figures                     ###
###                                  ###
########################################
def plot_lopetus():
  plt.show()
  
########################################
### close figure                     ###
###                                  ###
########################################
def plot_suljetus( joojoo = True ):
    if joojoo:
        plt.close()
        
########################################
### change y-limits of the plot      ###
###                                  ###
########################################
def plot_setYlim( minimiY, maksimiY, extendBelowZero = True, A = 0.05 ):
    from sys import float_info
    
    # A = extensionparametri
    #print 'minimiY '  + str(minimiY)
    #print 'maksimiY ' + str(maksimiY)
    if (abs(minimiY-maksimiY) >  float_info.epsilon*10 ):
        if extendBelowZero:
            ymin = minimiY - A*(maksimiY-minimiY) 
        else:
            ymin = 0.0
    
        ymax = maksimiY + A*(maksimiY-minimiY)
        #print 'y limit min ' + str(ymin)
        #print 'y limit max ' + str(ymax)
        plt.ylim( ymin, ymax )

########################################
### change x-limits of the plot      ###
###                                  ###
########################################
def plot_setXlim( minimiX, maksimiX, extendBelowZero = True, A = 0.05 ):
    from sys import float_info
    
    # A = extensionparametri
    #print 'minimiY '  + str(minimiY)
    #print 'maksimiY ' + str(maksimiY)
    if (abs(minimiX-maksimiX) >  float_info.epsilon*10 ):
        if extendBelowZero:
            xmin = minimiX - A*(maksimiX-minimiX) 
        else:
            xmin = 0.0
    
        xmax = maksimiX + A*(maksimiX-minimiX)
        #print 'y limit min ' + str(ymin)
        #print 'y limit max ' + str(ymax)
        plt.xlim( xmin, xmax )
    
########################################
### plot data                        ###
###                                  ###
########################################
def plottaa( x, y, tit = ' ', xl = ' ', yl = ' ', label=None, log=False, currentColor = 'b', changeColor=True, tightXAxis=False, tightYAxis = False, markers=False, LEGEND=True, omavari = False, scatter=False, uusikuva = False, gridi = True, loc = 3, markersize = 10, marker = 'o', a = 24, b=15, sub_i = 1, sub_j = 1, sub_k = 1, uusifig = True, uusisub = True, figuuri = None):
    
  
  if uusikuva:
      plottaa.fig, plottaa.ax = plot_alustus(a=a,b=b, uusifig = uusifig, uusisub = uusisub, sub_i = sub_i, sub_j = sub_j, sub_k = sub_k, figuuri = figuuri )
  
  global color
  if  label is None:
      label = tit
  if ( (omavari is False) and ('colorChoice' in globals() )):
    if changeColor:
        currentColor = colorChoice.getNextColor()
    else:
        currentColor = colorChoice.getCurrentColor()
  else:    
      currentColor = omavari

  


  if markers and not scatter:
      
      plt.plot( x, y, color = currentColor, label=label, linestyle='-', marker=marker, markersize = markersize )
  elif not markers and not scatter:
      plt.plot( x, y, color = currentColor, label=label) # default
  elif scatter:
      plt.scatter( x, y, color = currentColor, label=label, s=markersize**2, marker = marker)
  
  if LEGEND:
      if loc == 2: # right side
        plt.legend(bbox_to_anchor=(1.02, 1), loc=loc, borderaxespad=0., fancybox = True, shadow = True )
      elif loc == 3: # top
        plt.legend(bbox_to_anchor=(0., 1.06, 1., .102), loc=loc, ncol=6, fancybox = True, shadow = True , mode="expand" ) # upper center ,  prop={'size': 18}      bbox_to_anchor = ( 0., 1.1, 0.6, 20.102 ),loc=9,  ncol=2, mode="expand", borderaxespad=0., fancybox = True, shadow = True
  
  plt.xlabel( xl ) #r'$\#/m^3$'
  plt.ylabel( yl )
  
  #plt.xticks( xtikut )
  plt.grid( gridi )
  
  plt.autoscale( enable=True, axis='x', tight=tightXAxis )
  
  if tightYAxis:
    plt.autoscale( enable=True, axis='y', tight=tightYAxis )
  
  #patch = mpatches.Patch(color=c, label=legend)
  #plt.legend(handles=[patch])

#
  plt.title(tit)
  
      
  if (log):
    plt.xscale('log')

  return plottaa.fig, plottaa.ax
    
#y = zt_data[1:]
#c = next(colorpool)
#plt.plot(Nc1,y, color=c, linewidth=3, label='Number of cloud droplets' )
#c = next(colorpool)
#plt.plot(Nic1,y, color=c, linewidth=3, label='Number of ice particles')




#red_patch = mpatches.Patch(color='red', label='Number of cloud droplets')
#plt.legend(handles=[red_patch])


#blue_patch = mpatches.Patch(color='blue', label='Number of ice particles')
#plt.legend(handles=[blue_patch])
#plt.legend()
#plt.figure()



