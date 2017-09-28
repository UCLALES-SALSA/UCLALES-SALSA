# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 16:28:50 2017

@author: aholaj
"""

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

def waterBallMassToDiam(m):
    from math import pi
    from numpy import power
    diam = power( 6*m/(pi*1000.), (1./3.))*1e6
    
    return diam # um

def waterBallVolumeToDiam(V):
    from math import pi
    from numpy import power
    diam = power( 6*V/pi, (1./3.) )*1e6
    
    return diam # um

def myRound(x, base = 5):
    return int(base*round(float(x)/base))

def myRoundFloat(x, prec=2, base=.5):
  return round(base * round(float(x)/base),prec)