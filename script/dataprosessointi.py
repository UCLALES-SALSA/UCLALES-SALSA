#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 13:50:26 2016

@author: aholaj
"""

###### read ice

from netCDF4 import Dataset
import numpy as np
import sys
import os

import matplotlib.pyplot as plt
from itertools import cycle
import matplotlib.patches as mpatches

def main( piirra = False ):
  for filebase in sys.argv[1:]:
    
    filenameNC = filebase + '.nc'
    filenamePS = filebase + '.ps.nc'
    filenameTS = filebase + '.ts.nc'
    
    global colors
    global colorpool    
    colors = ['r','b','g','c','m','y','k']
    colorpool = cycle(colors)

    os.path.isfile( filenameNC )
  
    S_Niba_data  = read_Data( filenameNC, 'S_Niba'  )
    S_Rwiba_data = read_Data( filenameNC, 'S_Rwiba' )

    S_Ncba_data  = read_Data( filenameNC, 'S_Ncba'  )
    S_Rwcba_data = read_Data( filenameNC, 'S_Rwcba' )

    S_Naba_data  = read_Data( filenameNC, 'S_Naba'  )
    S_Rwaba_data = read_Data( filenameNC, 'S_Rwaba' )


    S_Nic_data   = read_Data( filenameNC, 'S_Nic'   )
    
    print laske_minimi_maksimi('s_nic_data', S_Nic_data )
    print ' '
    S_Nc_data    = read_Data( filenameNC, 'S_Nc'    )

    zt_data      = read_Data( filenameNC, 'zt'      )
    zm_data      = read_Data( filenameNC, 'zm'      )
    
    yt_data      = read_Data( filenameNC, 'yt'      )
    xt_data      = read_Data( filenameNC, 'xt'      )

    ym_data      = read_Data( filenameNC, 'ym'      )
    xm_data      = read_Data( filenameNC, 'xm'      )

    t_data       = read_Data( filenameNC, 't'       )
    
    w_data       = read_Data( filenameNC, 'w'       )
    
    f_data       = read_Data( filenameNC, 'f'       ) # total ice mixing ratio (ice+snow)
    l_data       = read_Data( filenameNC, 'l'       ) # total liquid water mixing ratio
    
    time_data    = read_Data( filenameNC, 'time'    )
    
    S_RHI_data   = read_Data( filenameNC, 'S_RHI'   )
    
    dn0_data     = read_Data( filenameNC, 'dn0'     )
    
    #P_rv_data	 = read_Data( filenamePS, 'P_rv'    )
    
    liqWP_Tdata  = read_Data( filenameTS, 'lwp_bar'  )
    rainWP_Tdata  = read_Data( filenameTS, 'rwp_bar' )
    
    
    #top_data     = read_Data( filenameTS, 'zc'      )
    #base_data    = read_Data( filenameTS, 'zb'      )
    
    
    korkeus = ( zm_data - zt_data )*2.0
    
    #depth   = top_data - base_data
    
    #aikasarjaTulostus( depth, nimi= 'cloud depth' )

    #print_shape( 'S_Rwiba', S_Rwiba_data )
    #print_shape( 'S_Rwcba', S_Rwcba_data )
    #print_shape( 'S_Rwaba', S_Rwaba_data )
    #print_shape( 'S_Nic',   S_Nic_data   )
    #print_shape( 'S_Nc',    S_Nc_data    )
    #print_shape( 'xt',	    xt_data      )
    #print_shape( 'yt',	    yt_data      )
    #print_shape( 'f',	    f_data       )
    #print_shape( 'dn0',     dn0_data     )
    #print_shape( 'height',  korkeus      )
    #print_shape( 'aika',    time_data    )
    ##print_shape( 'depth',   depth         )
    #print_shape( 'water-vap mix rat', P_rv_data )
    
    print time_data
    
    #area( xm_data, ym_data )
    
    #laske_minimi_maksimi_bini('S_Rwiba',S_Rwiba_data)
    #laske_minimi_maksimi_bini('S_Rwcba',S_Rwcba_data)
    #laske_minimi_maksimi('S_Nc',S_Nc_data)
    #laske_minimi_maksimi('S_Nic',S_Nic_data)
    #laske_minimi_maksimi('t',t_data)
    #laske_minimi_maksimi('w',w_data)
    #laske_minimi_maksimi('f',f_data)
    #laske_minimi_maksimi('l',l_data)
    #laske_minimi_maksimi('S_RHI',S_RHI_data)
    #time=1
    #Zi=66
    #print 'Ajankohta: '+str(time)
    #print 'Korkeus: '+str(Zi)
    #print 
    
    #print 'S_RHI', S_RHI_data[time,0,0,Zi]
    #print 't', t_data[time,0,0,Zi]
    #print 'w', w_data[time,0,0,Zi]
    #print 'S_Nic', S_Nic_data[time,0,0,Zi]
    
    #vertical_sum_timeseries(      'S_Nic',   S_Nic_data,   time_data,  tulostus=True)
    #vertical_sum_timeseries(      'S_Nc',    S_Nc_data,    time_data,  tulostus=True)
    #vertical_sum_timeseries(      'f',       f_data,       time_data,  tulostus=True)
    #vertical_sum_timeseries_bini( 'S_Rwiba', S_Rwiba_data, time_data,  tulostus=True)
    #vertical_sum_timeseries_bini( 'S_Rwcba', S_Rwcba_data, time_data,  tulostus=True)
    
    #laske_path_aikasarjaXYZ( f_data,    dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, nimi = 'Ice water path')
    #laske_path_aikasarja( f_data, dn0_data, korkeus, time_data, onlyCloudy = True,  tulostus = True, piirra = False )
    laske_path_aikasarjaXYZ( l_data,    dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, nimi = 'Liquid water path' )
    #laske_path_aikasarjaZ( P_rv_data, dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, nimi = 'Water vapor path')
    
    laske_WC( f_data, dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'Ice water content')
    laske_WC( l_data, dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'Liquid water content')
    
    
    laske_NumberConcentration( S_Nic_data, S_Nc_data,  dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'number concentration, ice in liquid clouds (mixed-phase)' )
    
    laske_NumberConcentration( S_Nic_data, S_Nic_data, dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'number concentration, ice in ice clouds' )
    
    laske_NumberConcentration( S_Nc_data, S_Nc_data,   dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'number concentration, liquid in liquid clouds' )
    
    #laske_NumberConcentration( l_data, f_data, dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'number concentration, liquid in ice clouds' )

    for bini in xrange( np.shape(S_Rwiba_data)[1] ):
        laske_MeanDiameterInBin( S_Rwiba_data, bini, S_Nc_data, time_data, tulostus = True, piirra =piirra, nimi = 'Mean diameter of ice particles in liquid clouds in bin ' )
    
    for bini in xrange( np.shape(S_Rwcba_data)[1] ):
        laske_MeanDiameterInBin( S_Rwcba_data, bini, S_Nc_data, time_data, tulostus = True, piirra =piirra, nimi = 'Mean diameter of liquid particles in liquid clouds in bin ' )            
    
    for bini in xrange( np.shape(S_Rwiba_data)[1] ):
        vertical_sum_timeseries( 'Rwiba bini'+str(bini), S_Rwiba_data[ :, bini, :, :, : ], time_data, tulostus = True, piirra = piirra )
    
    
    
    for bini in xrange( np.shape(S_Rwcba_data)[1] ):
        vertical_sum_timeseries( 'Rwcba bini'+str(bini), S_Rwcba_data[ :, bini, :, :, : ], time_data, tulostus = True, piirra = piirra )
    
    vertical_sum_timeseries(      'f',       f_data,       time_data,  tulostus=True, piirra=piirra )
    vertical_sum_timeseries(      'l',       l_data,       time_data,  tulostus=True, piirra=piirra )
    
    #Nc = S_Nc_data[time,0,0,1:]
    #Nic = S_Nic_data[time,0,0,1:]

    #nollaa(Nc)  
    #nollaa(Nic)
    
    #titpohja="PSD aika: "+tstepH(time)
    #xl='sade'
    #yl='lkm'
    
    #plot_alustus()
    #tit=titpohja+" Aero"
    #y = S_Naba_data[time,:,0,0,korkeus]
    #x = S_Rwaba_data[time,:,0,0,korkeus]
    #plottaa(x,y,tit,xl,yl,True)
    
    #plot_alustus()
    #tit=titpohja+" Cloud"
    #y = S_Ncba_data[time,:,0,0,korkeus]
    #x = S_Rwcba_data[time,:,0,0,korkeus]
    #plottaa(x,y,tit,xl,yl,True)
    
    #plot_alustus()
    #tit=titpohja+" Ice"
    #y = S_Niba_data[time,:,0,0,korkeus]
    #x = S_Rwiba_data[time,:,0,0,korkeus]
    #plottaa(x,y,tit,xl,yl,True)
    
    #time=38
    #titpohja="Vertical profile at "+tstepH(time)
    
    #plot_alustus()
    #tit=titpohja+" Cloud"
    #y = S_Ncba_data[time,:,0,0,korkeus]
    #x = S_Rwcba_data[time,:,0,0,korkeus]
    #plottaa(x,y,tit,xl,yl,True)
    
    #plot_alustus()
    #tit=titpohja+" Ice"
    #y = S_Niba_data[time,:,0,0,korkeus]
    #x = S_Rwiba_data[time,:,0,0,korkeus]
    #plottaa(x,y,tit,xl,yl,True)
    
     
    if piirra:
      plot_lopetus()

def tstepH(time):
  return str(time*360/3600.)

def nollaa(data):
  for i in xrange(np.shape(data)[0]):
    if data[i]<10e-10:
      data[i]=10e-10  
  
  return data
  
def read_Data(filename,var):

  fileH = Dataset(filename,mode='r')
  data = fileH.variables[var][:]
  fileH.close()
  return data

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
  
def vertical_sum_timeseries( nimi, data, aika, tulostus = False, piirra = False ):
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

  if piirra:
    plot_alustus()
    plottaa( aika, aikasarja, nimi, 'aika', 'path')
    
def vertical_sum_timeseries_bini( nimi, data, aika, tulostus = False, piirra = False ):
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
  
  if piirra:
    plot_alustus()
    plottaa( aika, aikasarja, nimi, 'aika', 'path')

def print_filename(fname):
  head, tail = os.path.split(fname)
  print ' '
  print 'file: ' + tail  

def print_shape(var,data):
  print ' '
  print 'variable: ' + var
  print 'shape var1: '+ str(np.asarray(np.shape(data)))



    
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

def laske_path_aikasarjaXYZ( mixRatData, dn0_data, korkeus, aika, onlyCloudy = False, tulostus = False, piirra = False, nimi = 'path aikasarja' ):
  
  mixRatData = mixRatData*1000.0 # kg/kg -> g/kg
  
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
    
  
  print 'dimensiot aikasarjaXYZ'+ str(np.shape(aika))+ ' timeseries '+ str(np.shape(timeSeries))
  
  if tulostus:
    print ' '
    print nimi
    for t in xrange(timeDim):
      print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(timeSeries[t])
  print ' '    

  if piirra:
    plot_alustus()
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'path [g/m^2]')

def laske_path_aikasarjaZ( mixRatData, dn0_data, korkeus, aika, tulostus = False, piirra = False, nimi = 'path aikasarja' ):
  
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

  if piirra:
    plot_alustus()
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'path [g/m^2]')

def laske_WC( mixRatData, dn0_data, aika, tulostus = False, piirra = False, nimi = 'water content aikasarja' ):
  
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

  if piirra:
    plot_alustus()
    plottaa( aika, timeSeries, nimi, 'aika [s]', 'Content [g/m^3]')  

def laske_NumberConcentration( Ndata, refNdata, dn0_data, aika, tulostus = False, piirra =False, nimi = 'number concentration' ):
    
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

    if piirra:
        plot_alustus()
        plottaa( aika, timeSeries, nimi, 'aika [s]', 'Number concentration [#/kg]')  

def laske_MeanDiameterInBin( RadiusBinData, bini, refNdata, aika, tulostus = False, piirra =False, nimi = 'Mean diameter in bin ' ):
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

    if piirra:
        plot_alustus()
        plottaa( aika, timeSeries, nimi, 'aika [s]', 'diameter [um]')          

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def aikasarjaTulostus( data, aika = 0, nimi = 'aikasarja' ):
  if aika == 0:
    aika=np.zeros( (np.shape(data)[0]))
  print ' '
  print nimi
  for t in xrange( np.shape(data)[0] ):
    print 'ajanhetki: ' + str(aika[t]) + ' ' + ' arvo : ' + str(data[t])
  print ' '    

    
def plot_alustus():
  plt.figure()

def plot_lopetus():
  plt.show()

def plottaa(x,y,tit,xl,yl,log=False):


  c = next(colorpool)
  plt.plot(x,y, color=c)
  plt.xlabel(xl) #r'$\#/m^3$'
  plt.ylabel(yl)
  plt.title(tit)
  plt.xticks( [7200, 7201, 7210, 7300] )
  plt.grid(True)
  if (log):
    plt.xscale('log')  
  
main(True)   
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



