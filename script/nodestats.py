#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 13:21:16 2017

@author: aholaj
"""
import os
import sys
import time
import subprocess

def bool_convert(s):
    if s=="True":
        r = True
    elif s=="False":
        r = False
    else:
        r = True
    return r



if ( len(sys.argv) > 1):
    folder = sys.argv[1]

else:
    folder = '/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/script'


qlog = 'qstatlogi'
nodelog = 'nodelog'
filu = folder + qlog
nodefilu = folder + nodelog


command = 'qstat -u aholaj > ' + filu

os.system(  command )
k=0

while ( os.stat( filu ).st_size != 0 ):
#    while ( k < 2):

    unixtime = int( time.time() )

    nooditvapaat = subprocess.check_output( 'xtnodestat | grep "compute nodes" | cut -c49-60  | tr -dc "[:alnum:]"', shell=True)

    les_r_nodes   = 0
    les_q_nodes   = 0
    postp_r_nodes = 0
    postp_q_nodes = 0
    
    les_r_nro   = 0
    les_q_nro   = 0
    postp_r_nro = 0
    postp_q_nro = 0
    
    threads_active  = 0
    threads_passive = 0
    
    cases = []
    i=0
    f   = open(filu, 'r')
    for line in f:
        if i > 4:
            nodes  = line[54]
            job    = line[34:37]
            status = line[73]
            emul = ''
            if  ( job == 'LES' ):
                emul = line[42:44]        
                if status == 'R' :
                    #print 'LES ajossa'
                    les_r_nro   += 1
                    les_r_nodes += int(nodes)
                    if emul not in cases:
                        threads_active += 1
                        cases.append(emul)
                elif status == 'Q':
                    #print 'LES jonossa'
                    les_q_nro   += 1
                    les_q_nodes += int(nodes)
                    if emul not in cases:
                        threads_passive += 1
                        cases.append(emul)
            elif  ( job == 'nc_' or job == 'ps_' or job == 'ts_' ):
                emul = line[41:43]            
                if status == 'R' :
                    #print 'postprosessointi ajossa'
                    postp_r_nro   += 1
                    postp_r_nodes += int(nodes)
                    if emul not in cases:
                        threads_active += 1
                        cases.append(emul)
                elif status == 'Q':
                    #print 'postprosessointi jonossa'
                    postp_q_nro   += 1
                    postp_q_nodes += int(nodes)
                    if emul not in cases:
                        threads_passive += 1
                        cases.append(emul)                
                    
                    
        i += 1
    f.close()            
    data =  str(unixtime) +','+ str(les_r_nodes)  +','+ str(les_q_nodes)    +','+ str(postp_r_nodes)   +','+ str(postp_q_nodes) +',' \
                              + str(les_r_nro)    +','+ str(les_q_nro)      +','+ str(postp_r_nro)     +','+ str(postp_q_nro)   +',' \
                              + str(nooditvapaat) +','+ str(threads_active) +','+ str(threads_passive) +'\n'
    print data

    nod = open(nodefilu , 'w')
    nod.write(data)
    nod.close()

    time.sleep(15)               
    os.system(  command )
    k += 1