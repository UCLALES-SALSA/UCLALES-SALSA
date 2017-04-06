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
   MINUUTTI=((AIKA%DAY)%HOUR)//MINUTE
   SEKUNTI=(((AIKA%DAY)%HOUR)%MINUTE)
   
   PAIVA    = str(int(PAIVA))
   TUNTI    = str(int(TUNTI))
   MINUUTTI = str(int(MINUUTTI))
   SEKUNTI  = str(int(SEKUNTI))
   

  ###
  ### Muuta seuraavan rivin print komentoa, mikali kaytat python 2 -versiota
   return PAIVA + "d" + " " + TUNTI + "h" + " " + MINUUTTI + "m"+ " " + SEKUNTI + "s"

if ( len(sys.argv) > 1):
    filu = sys.argv[1]
    
else:
    filu = '/home/aholaj/koodit_IL/bash/emulatortime/log'

folder = os.path.dirname(os.path.realpath( filu ))
f = open(filu, 'r+')

unixtimes=[]
for line in f:
    if (line[0] != 'L' and line[0] !=' ' and line[0] !='S'):
        unixtimes.append(int(line))

aloitusaika  = min( unixtimes )
lopetusaika  = max( unixtimes )
suoritusaika = lopetusaika - aloitusaika

f.write( 'Suoritusaika ' +  Muunnos(suoritusaika)+'\n' )




f.close()