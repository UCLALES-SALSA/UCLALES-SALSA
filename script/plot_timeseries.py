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

#print matplotlib.matplotlib_fname()

filename1 = sys.argv[1]

S_Niba='S_Niba'
S_Rwiba='S_Rwiba'

S_Ncba='S_Ncba'
S_Rwcba='S_Rwcba'

S_Naba='S_Naba'
S_Rwaba='S_Rwaba'

S_Nc='S_Nc'
S_Nic='S_Nic'

zt='zt'

t='t'

fileH1 = Dataset(filename1,mode='r')
fileH2 = Dataset(filename2,mode='r')

#S_Niba_data = fileH.variables[S_Niba][:]
#S_Rwiba_data = fileH.variables[S_Rwiba][:]

#S_Ncba_data = fileH.variables[S_Ncba][:]
#S_Rwcba_data = fileH.variables[S_Rwcba][:]

#S_Naba_data = fileH.variables[S_Naba][:]
#S_Rwaba_data = fileH.variables[S_Rwaba][:]

S_Nc_data1 = fileH1.variables[S_Nc][:]
S_Nic_data1 = fileH1.variables[S_Nic][:]

S_Nc_data2 = fileH2.variables[S_Nc][:]


#S_Nic_data2 = fileH2.variables[S_Nic][:]

zt_data = fileH1.variables[zt][:]

t1_data = fileH1.variables[t][:]
t2_data = fileH2.variables[t][:]


fileH1.close()
fileH2.close()






time=30
korkeus=61
print 'Ajankohta: '+str(time)
print ' '



colors = ['r','b','g','c','m','y','k']
colorpool = cycle(colors)

def print_filename(fname):
  head, tail = os.path.split(fname)
  print ' '
  print 'file: ' + tail  

def print_shape(var,data):
  print ' '
  print 'variable: ' + var
  print 'shape var1: '+ str(np.shape(data))

print_filename(filename1)
print_filename(filename2)



##for bini in xrange(7):
  ##print 'bini: ' + str(bini)
  ##print 'maxi: '+ str(np.max(S_Ncba_data[time,bini,0,0,:]))
  ##print 'maxi indeksi: '+ str(np.argmax(S_Ncba_data[time,bini,0,0,:]))
  
#c = next(colorpool)
#y = S_Ncba_data[time,:,0,0,korkeus]
#x = S_Rwcba_data[time,:,0,0,korkeus]
#print y
#print x
#plt.plot(x,y, color=c)


#c = next(colorpool)
#y = S_Naba_data[time,:,0,0,korkeus]
#x = S_Rwaba_data[time,:,0,0,korkeus]
#plt.plot(x,y, color=c)

#c = next(colorpool)
#y = S_Niba_data[time,:,0,0,korkeus]
#x = S_Rwiba_data[time,:,0,0,korkeus]
#plt.plot(x,y, color=c)


#plt.xscale('log')


#print_shape(S_Nc, S_Nc_data1)

#print_shape(S_Nic, S_Nic_data1)

#print_shape(zt, zt_data)

#print zt_data

Nc1 = S_Nc_data1[20:,0,0,1:]
Nic1 = S_Nic_data1[20:,0,0,1:]   
#Nc2 = S_Nc_data2[time,0,0,1:]


print_shape('Nc',Nc1)
print_shape('Nic',Nic1)

for i in xrange(np.shape(Nic1)[0]):
  if Nic1[i]<10e-10:
    Nic1[i]=10e-10

for i in xrange(np.shape(Nc1)[0]):
  if Nc1[i]<10e-10:
    Nc1[i]=10e-10
    
#for i in xrange(np.shape(Nc2)[0]):
  #if Nc2[i]<10e-10:
    #Nc2[i]=10e-10
    
#print 'erotusvektorin summa'    
#print sum(Nc1 -Nc2)

#print_shape('t',t1_data)
#print 'jää_lkm_maks'
#print np.max(Nic1)

plt.figure()    
y = zt_data[1:]
c = next(colorpool)
plt.plot(Nc1,y, color=c, linewidth=3, label='Number of cloud droplets' )
c = next(colorpool)
plt.plot(Nic1,y, color=c, linewidth=3, label='Number of ice particles')

c = next(colorpool)
#plt.plot(Nc2,y, color=c)

plt.xlabel(r'$\#/m^3$')
plt.ylabel('height (m)')
plt.xscale('log')

#red_patch = mpatches.Patch(color='red', label='Number of cloud droplets')
#plt.legend(handles=[red_patch])


#blue_patch = mpatches.Patch(color='blue', label='Number of ice particles')
#plt.legend(handles=[blue_patch])
plt.legend()
#plt.figure()

#t1 = t1_data[time,0,0,1:]
#t2 = t2_data[time,0,0,1:]
#c = next(colorpool)
#plt.plot(t1,y, color=c)


#c = next(colorpool)
#plt.plot(t2,y, color=c)
plt.title("Vertical profile at 3 hours")
plt.show()
#print_shape(S_Ncba,S_Ncba_data)

#print_shape(S_Rwcba,S_Rwcba_data)

#print_shape(S_Naba,S_Naba_data)

#print_shape(S_Rwaba,S_Rwaba_data)