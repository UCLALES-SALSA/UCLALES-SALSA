# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 13:33:06 2017

@author: aholaj
"""
import os
import sys
#f = open(, 'r')

def Muunnos(aika):
   AIKA=float(aika)
   DAY=(24*60*60)
   HOUR=(60*60)
   MINUTE=60

   PAIVA=AIKA//DAY
   TUNTI=(AIKA%DAY)//HOUR
   MINUUTTI=round(((AIKA%DAY)%HOUR)/MINUTE,0)
   
   PAIVA=int(PAIVA)
   TUNTI=int(TUNTI)
   MINUUTTI=int(MINUUTTI)
   
   PAIVA=str(PAIVA)
   TUNTI=str(TUNTI)
   MINUUTTI=str(MINUUTTI)
   

  ###
  ### Muuta seuraavan rivin print komentoa, mikali kaytat python 2 -versiota
  #print PAIVA + "d" + ":" + TUNTI + "h" + ":" + MINUUTTI + "min"
   print( PAIVA + "d " + TUNTI+"h "+ MINUUTTI +"min")

if ( len(sys.argv) > 1):
    filu = sys.argv[1]
    
else:
    filu = '/home/aholaj/koodit_IL/bash/emulatortime/log'

folder = os.path.dirname(os.path.realpath( filu ))
f = open('/home/aholaj/koodit_IL/bash/emulatortime/log', 'r')

unixtimes=[]
for line in f:
    if (line[0] != 'L' and line[0] !=' '):
        unixtimes.append(int(line))
f.close()
aloitusaika  = min( unixtimes )
lopetusaika  = max( unixtimes )
suoritusaika = lopetusaika - aloitusaika

Muunnos(suoritusaika)
