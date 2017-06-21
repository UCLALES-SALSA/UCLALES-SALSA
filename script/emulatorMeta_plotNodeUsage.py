#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:35:32 2017

@author: aholaj
"""
#
# plottaa nodelogi -tiedoston perusteella noodien kayton
#
import ModDataPros as mdp
from PythonMethods import Muunnos
import sys
import os
import numpy as np
from datetime import datetime
import tzlocal
import matplotlib.pyplot as plt



if ( len(sys.argv) > 1):
    filu = sys.argv[1]
    folder = os.path.dirname(os.path.realpath( filu ))
else:
    folder = '/home/aholaj/mounttauskansiot/voimahomemount/UCLALES-SALSA/script/'
    folder = '/home/aholaj/mounttauskansiot/lustremount/UCLALES-SALSA/'
    folder = '/home/aholaj/mounttauskansiot/ibrixmount/'
   
    subfolder = ''#'case_emulator_DESIGN_v1.5.1_LES_cray.dev07042017_LVL3_BSPsmooth0.2/'
    subfolder = 'case_emulator_DESIGN_v1.5.1_LES_cray.dev07042017_LVL3_BSPsmooth0.1/'
    filu   = folder + subfolder + 'nodelog'
    
f = open(filu, 'r')
#        data =  str(unixtime) +','+ str(les_r_nodes)  +','+ str(les_q_nodes)    +','+ str(postp_r_nodes)   +','+ str(postp_q_nodes) +',' \
#                                  + str(les_r_nro)    +','+ str(les_q_nro)      +','+ str(postp_r_nro)     +','+ str(postp_q_nro)   +',' \
#                                  + str(nooditvapaat) +','+ str(threads_active) +','+ str(threads_passive) +'\n'

picturefolder= folder + '/nodestats_pictures/'

if not os.path.exists( picturefolder ):
    os.makedirs( picturefolder )

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

local_timezone = tzlocal.get_localzone() 
aloitusaika = datetime.fromtimestamp(  unixtime[0], local_timezone )
lopetusaika = datetime.fromtimestamp( unixtime[-1], local_timezone )

time = np.asarray( unixtime) - unixtime[0]*np.ones( len( unixtime) )
pvm = aloitusaika.strftime("%a %x %T") + ' -> ' + lopetusaika.strftime("%a %x %T")
kesto = Muunnos( time[-1] )

suoritusaika = pvm + ' == ' + kesto

#############
### nodes ###
#############
mdp.initializeColors(6)

RUN = np.asarray(les_r_nodes) + np.asarray(postp_r_nodes)
QU  = np.asarray(les_q_nodes) + np.asarray(postp_q_nodes)

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(RUN), np.max(QU) )+1. )
mdp.plottaa( time, les_r_nodes,   ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='# les run nodes',     changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, les_q_nodes,   ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='# les queue nodes',   changeColor = True, markers=True, tightXAxis = True )
plt.savefig( picturefolder + '/' + 'nodes_les' + '.png' )

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(RUN), np.max(QU) )+1. )
mdp.plottaa( time, postp_r_nodes, ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='# postp run nodes',   changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, postp_q_nodes, ' ' , xl = 'time '+ suoritusaika, yl = '# nodes', label='# postp queue nodes', changeColor = True, markers=True, tightXAxis = True )
plt.savefig( picturefolder + '/' + 'nodes_postp' + '.png' )

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max( RUN ), np.max( QU ), np.max( nooditvapaat ) ), +1. )
mdp.plottaa( time, RUN, ' ', xl = 'time '+ suoritusaika, yl = '# nodes', label='# run nodes',   changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, QU,  ' ', xl = 'time '+ suoritusaika, yl = '# nodes', label='# queue nodes', changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, nooditvapaat, ' ', xl = 'time '+ suoritusaika, yl = '# nodes', label='# free nodes', changeColor = True, markers=True, tightXAxis = True )
plt.savefig( picturefolder + '/' + 'nodes_total_free' + '.png' )

#############
### nro   ###
#############
mdp.initializeColors(6)

Rnro = np.asarray(les_r_nro) + np.asarray(postp_r_nro)
Qnro = np.asarray(les_q_nro) + np.asarray(postp_q_nro)

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(Rnro), np.max(Qnro) )+1. )
mdp.plottaa( time, les_r_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='# les run jobs',   changeColor = True, markers=True,     tightXAxis = True )
mdp.plottaa( time, les_q_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='# les queue jobs', changeColor = True, markers=True,     tightXAxis = True )
plt.savefig( picturefolder + '/' + 'jobs_les' + '.png' )

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(Rnro), np.max(Qnro) )+1. )
mdp.plottaa( time, postp_r_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='# postp run jobs',   changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, postp_q_nro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='# postp queue jobs', changeColor = True, markers=True, tightXAxis = True )
plt.savefig( picturefolder + '/' + 'jobs_postp' + '.png' )

mdp.plot_alustus()
mdp.plot_setYlim( 0.,  max(np.max(Rnro), np.max(Qnro) )+1. )
mdp.plottaa( time, Rnro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='# run jobs',   changeColor = True, markers=True, tightXAxis = True )
mdp.plottaa( time, Qnro, ' ' ,xl = 'time '+ suoritusaika, yl = '# jobs', label='# queue jobs', changeColor = True, markers=True, tightXAxis = True )
plt.savefig( picturefolder + '/' + 'jobs_total' + '.png' )

###############
### threads ###
###############
mdp.initializeColors(2)

TH = np.asarray(threads_active) + np.asarray(threads_passive)
print threads_active[10]+ threads_passive[10]
mdp.plot_alustus()
mdp.plot_setYlim ( 0., np.max(TH)+1. )
mdp.plottaa( time, threads_active,   ' ', xl = 'time '+ suoritusaika, yl = '# jobs', label='# threads active of '  + str(np.max(TH))+ " threads",   changeColor = True, markers=True,     tightXAxis = True   )
mdp.plottaa( time, threads_passive,  ' ', xl = 'time '+ suoritusaika, yl = '# jobs', label='# threads passive of  '+ str(np.max(TH))+ " threads", changeColor = True, markers=True,     tightXAxis = True  )
plt.savefig( picturefolder + '/' + 'threads' + '.png' )

mdp.plot_lopetus()

print suoritusaika

