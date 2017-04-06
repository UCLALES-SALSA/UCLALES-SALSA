# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:35:32 2017

@author: aholaj
"""
import ModDataPros as mdp
import emulatorRunTime as ert
import sys
import os
import numpy as np




if ( len(sys.argv) > 1):
    filu = sys.argv[1]
    folder = os.path.dirname(os.path.realpath( filu ))
else:
    folder = '/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/script'
    filu   = 'nodelog'
    
f = open(filu, 'r')
#        data =  str(unixtime) +','+ str(les_r_nodes)  +','+ str(les_q_nodes)    +','+ str(postp_r_nodes)   +','+ str(postp_q_nodes) +',' \
#                                  + str(les_r_nro)    +','+ str(les_q_nro)      +','+ str(postp_r_nro)     +','+ str(postp_q_nro)   +',' \
#                                  + str(nooditvapaat) +','+ str(threads_active) +','+ str(threads_passive) +'\n'

unixtime = []

les_r_nodes   = []
les_q_nodes   = []
postp_r_nodes = []
postp_q_nodes = []

les_r_nro   = []
les_q_nro   = []
postp_r_nro = []
postp_q_nro = []

nooditvapaat = []

threads_active  = []
threads_passive = []

for line in f:
    A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11 = line.split(',')
    
    unixtime.append( int(A0) )
    
    les_r_nodes.append( int(A1) )
    les_q_nodes.append( int(A2) )
    postp_r_nodes.append( int(A3) )
    postp_q_nodes.append( int(A4) )
    
    les_r_nro.append( int(A5) )
    les_q_nro.append( int(A6) )
    postp_r_nro.append( int(A7) )
    postp_q_nro.append( int(A8) )
    
    nooditvapaat.append( int(A9) )
    
    threads_active.append( int(A10) )
    threads_passive.append( int(A11) )


time = np.asarray( unixtime) - unixtime[0]*np.ones( len( unixtime) )

suoritusaika = ert.Muunnos( time[-1] )

#############
### nodes ###
#############
mdp.initializeColors(6)

RUN = np.asarray(les_r_nodes) + np.asarray(postp_r_nodes)
QU  = np.asarray(les_q_nodes) + np.asarray(postp_q_nodes)

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(RUN), np.max(QU) )+1. )
mdp.plottaa( time, les_r_nodes,   ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='les_run',     changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, les_q_nodes,   ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='les_queue',   changeColor = True, markers=True, tightXAxis = True )

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(RUN), np.max(QU) )+1. )
mdp.plottaa( time, postp_r_nodes, ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='postp_run',   changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, postp_q_nodes, ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='postp_queue', changeColor = True, markers=True, tightXAxis = True )


mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(RUN), np.max(QU) )+1. )
mdp.plottaa( time, RUN, 'nodes ' + suoritusaika, xl = 'time '+ suoritusaika, yl = '# nodes', label='run',   changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, QU, 'nodes ' + suoritusaika, xl = 'time '+ suoritusaika, yl = '# nodes', label='queue', changeColor = True, markers=True, tightXAxis = True )

#############
### nro   ###
#############
mdp.initializeColors(6)

Rnro = np.asarray(les_r_nro) + np.asarray(postp_r_nro)
Qnro = np.asarray(les_q_nro) + np.asarray(postp_q_nro)

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(Rnro), np.max(Qnro) )+1. )
mdp.plottaa( time, les_r_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='les_run',   changeColor = True, markers=True,     tightXAxis = True )
mdp.plottaa( time, les_q_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='les_queue', changeColor = True, markers=True,     tightXAxis = True )

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(Rnro), np.max(Qnro) )+1. )
mdp.plottaa( time, postp_r_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='postp_run',   changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, postp_q_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='postp_queue', changeColor = True, markers=True, tightXAxis = True )

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(Rnro), np.max(Qnro) )+1. )
mdp.plottaa( time, Rnro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='run',   changeColor = True, markers=True, tightXAxis = True, )
mdp.plottaa( time, Qnro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='queue', changeColor = True, markers=True, tightXAxis = True, )

###############
### threads ###
###############
mdp.initializeColors(2)

TH = np.asarray(threads_active) + np.asarray(threads_passive)

mdp.plot_alustus()
mdp.plot_setYlim ( 0., np.max(TH)+1. )
mdp.plottaa( time, threads_active,   ' ', xl = 'time '+ suoritusaika, yl = '# jobs', label='threads_active',   changeColor = True, markers=True,     tightXAxis = True, tightYAxis = False   )
mdp.plottaa( time, threads_passive,  ' ', xl = 'time '+ suoritusaika, yl = '# jobs', label='threads_passive', changeColor = True, markers=True,     tightXAxis = True , tightYAxis = False  )



mdp.plot_lopetus()