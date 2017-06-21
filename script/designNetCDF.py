#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 13:54:50 2017

@author: aholaj
"""

import numpy as np
import matplotlib.pylab as plt
import ModDataPros as mdp

ncfolder = '/home/aholaj/mounttauskansiot/ibrixmount/DESIGNnetcdf/'

tag = 'v1.5.1'

def readDESIGN(tag, var):
    from ModDataPros import read_Data
    ncfile = ncfolder + 'design_'+tag + '.nc'
    return read_Data(ncfile, var)

def plottaaVariable(data, k, tit, xl, yl, tags):
    cases = np.arange(1,91)
    mdp.initializeColors(np.shape(data)[2])
    mdp.plot_alustus()
    for ver in xrange(np.shape(data)[2]):
        mdp.plottaa( cases, np.sort(data[:, k, ver ]), tit + ' '+tags[ver], xl, yl, markers=True )

def plottaaCross(data, i, xl, j, yl, tit, tags):
    mdp.initializeColors(np.shape(data)[2])
    mdp.plot_alustus()
    for ver in xrange(np.shape(data)[2]):
        mdp.plottaa( data[:, i, ver ], data[:, j, ver ], tit + ' '+tags[ver], xl, yl, markers=True, scatter = True )
#    mdp.plot_setYlim( 300, 3000, extendBelowZero = False)
def main():
    nroCases = 90

    tags = ['v1.4.0', 'v1.5.0', 'v1.5.1']
#    var = ['q_inv', 'tpot_inv', 'clw_max','tpot_pbl', 'pblh', 'num_pbl', 'q_pbl', 'cloudbase', 'thickness', 'lwp_max' ]
    var = ['q_inv', 'tpot_inv', 'clw_max','tpot_pbl', 'pblh', 'num_pbl', 'q_pbl', 'cloudbase', 'thickness']
    data = np.zeros( ( nroCases, len(var), len(tags) ) )
    for v in xrange(len(var)):
        for t in xrange(len(tags)):
            data[:, v, t ] = readDESIGN( tags[t], var[v] )
    
    plottaaVariable( data, 2,   'clw_max', 'sorted', 'g/kg', tags )
    plottaaVariable( data, 4,    '  pblh', 'sorted',    'm', tags )
    plottaaVariable( data, 8, 'thickness', 'sorted',    'm', tags )
#    plottaaVariable( data, 9,   'lwp_max', 'case', 'g/kg', tags )
    
#    plottaaCross(    data, 2,   'clw_max', 9, 'lwp_max',       'clw vs lwp', tags )
#    plottaaCross(    data, 8, 'thickness', 9, 'lwp_max', 'thickness vs lwp', tags )    
#def plottaa( x, y, tit, xl, yl, label=0, log=False, changeColor=True, tightXAxis=False, markers=False, LEGEND=True, omavari = False):


main()
mdp.plot_lopetus()
