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