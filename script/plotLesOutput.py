#!/usr/bin/python
import time
tic = time.clock()
import ModDataPros as mdp
import sys
from itertools import cycle
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
from netCDF4 import Dataset
from PythonMethods import myRound
global piirra, tulostus, tightXAxis, saveFig, LEGEND, tag

kylla = [ 'y', 'Y', 'yes', 'Yes', 'YES', 'True', 'true', '1' ]
emuloo = ['e', 'E', 'emul', 'y', 'Y', 'yes']
lvl = raw_input( 'Level: ')
EMUL = raw_input('Emulator/Regular pictures: e/r (y/n): ') in emuloo
tag  = raw_input('give tag: ')
icevari = 'k'
if len(tag) == 4 and tag[0:3]=='ice':
    if tag[-1] == '4':
        icevari = 'k'
    elif tag[-1] == '1':
        icevari = 'r'
else:
    icevari = None

if len(tag)>0:
    tag = tag + ' '

ICE = (not EMUL )
NORMIPUTTI = ( not EMUL ) and ( lvl == '4' )
importOnly = False
if lvl == '0':
    importOnly = True
#ICE  = raw_input('Ice plots: yes/no: ') in kylla
#NORMIPUTTI = raw_input('Regular LES drawings: yes/no: ') in kylla

    

#######################
##### setting up    ###    
##### DO NOT CHANGE ###
#######################
if not importOnly:
    ibrix = os.environ["IBRIXMOUNT"]
    global arguments
    global labelArray
    arguments = sys.argv
    
    labelArray = []

    filenameNC = []
    filenamePS = []
    filenameTS = []
    
    
    customLabels = raw_input("Do you want to give custom labels for plots (yes/no): ") in kylla
    LVLprint     = raw_input("Do you want to print LEVEL (yes/no): ") in kylla
    if LVLprint:
        LVLprintFig  = ' ' + 'LVL' + lvl
        LVLprintSave = '_' + 'LVL' + lvl
    else:
        LVLprintFig  = ''
        LVLprintSave = ''
        
    kkk = 1
    for filebase in arguments[1:]:

        filenameNC.append( filebase + '.nc'    )
        filenamePS.append( filebase + '.ps.nc' )
        filenameTS.append( filebase + '.ts.nc' )
        
        if not EMUL:
            if customLabels:
                teksti = "Give "+ str(kkk) +". label: "
                labelArray.append( raw_input(teksti) )
                kkk += 1
            else:
                labelArray.append( filebase )
        
        
            


##############################
#### YOU CAN CHANGE THESE: ###
##############################

##############################
#
# drawing GLOBAL  
#
# parameter settings 
#
##############################
cases = 1


piirra = True

tulostus = False

tightXAxis = True

saveFig=True


if not importOnly:
    if saveFig:
        picturefolder='./pictures/'
        if not os.path.exists( picturefolder ):
            os.makedirs( picturefolder )

    cases    = len(arguments)-1
    colorNRO = cases

    if ( colorNRO > 6 ) :
        LEGEND = False
    else:
        LEGEND = True

    mdp.initializeColors(colorNRO)

global jarjestys
jarjestys = np.arange( 1, cases + 1 )
ajanhetket = 3*3600. #7200.

spinup = mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'Tspinup' )
tmax   = mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'timmax'  )

ticksHours  = np.arange(0, tmax/3600. + 0.1, 0.5)
xLabelsHours = map(str, ticksHours )

xTicksSeconds = ticksHours*3600.


##########################
###                    ###
### GLOBAL SUBROUTINES ###
###                    ###
##########################

##############
def varibaari( variRefVektori, variKartta ):
    minRef = np.min( variRefVektori )
    maxRef = np.max( variRefVektori )
    step = 1    
    Z = [[ 0,0 ],[ 0,0 ]]
    levels = range( int(minRef), int(maxRef) + step, step )
    CS3 = plt.contourf( Z, levels, cmap = variKartta )
    plt.close()
    return CS3

##########################
def piirra_aikasarjasettii( muuttuja, variKartta = plt.cm.gist_rainbow, variRefVektori = jarjestys, colorBar = None, colorBarTickValues = [0,1], colorBarTickNames = ['< -1', '0', '> 1'], muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'time [h]', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, relative = False, savePrefix = None, omaVari = True, tit = ' ', nollaArvo = None, xlabels = None, ylabels = None, xticks = None, yticks = None, spinup = None, ylabelFont = False ):

    maksimi = None
    minimi  = None
    maxRef = np.max( variRefVektori )
    maksInd = None
    minInd  = None
    #vastausKaikkiin = False
    nono = None
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        time_data       = mdp.read_Data( filenameTS[i], 'time'  )
        muuttuja_Tdata  = mdp.read_Data( filenameTS[i], muuttuja )*muunnosKerroin
        

        #### relative    
        if relative:
            if ((nollaArvo is None) and (nono is None)):
                lahtoarvo = muuttuja_Tdata[0]
            elif ( isinstance( nollaArvo, np.ndarray) or isinstance(nollaArvo, list) ):
                lahtoarvo = nollaArvo[i]
            elif nono is not None:
                lahtoarvo = muuttuja_Tdata[nono]
            if lahtoarvo == 0.:
                    print 'VAROITUS: ei voida laskea relative koska lahtoarvo on nolla', muuttuja, 'case', i+1, 'nono', nono
                    if nono is None:
                        NEXT = raw_input('Haluatko etsia seuraavan ei-negatiivisen datasta kaikkiin kappyroihin: yes/no: ') in kylla
                    
                        if NEXT:
                            nono = 0
                            while( muuttuja_Tdata[nono] == 0. ):
                                nono += 1
                            lahtoarvo = muuttuja_Tdata[nono]
                        else:
                            print "Let's skip this", muuttuja, i+1
                            continue
                        
                    else:
                        print "Let's skip this", muuttuja, i+1
                        continue
                        
            muuttuja_Tdata = np.where( muuttuja_Tdata > -999.,  muuttuja_Tdata / lahtoarvo, 0. )
        #### end relative        
                
                
        if EMUL:
            case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
        else:
            case_indeksi = i
        
        
        colorMap = variKartta
        skal = variRefVektori[ case_indeksi ] / maxRef
        color = colorMap(skal)
        
        if omaVari:
            omaVari = color

        if maksInd is not None:
            if np.max(muuttuja_Tdata) > maksimi: # maksimia ei ole viela tassa kohtaa paivitetty
                # maksimi indeksi muuttuu koska uusi maksimi on suurempi kuin vanha maksimi
                maksInd = np.argmax( muuttuja_Tdata)
        else:
            maksInd = np.argmax( muuttuja_Tdata)
        
        if minInd is not None:
            if np.min(muuttuja_Tdata) < minimi: # maksimia ei ole viela tassa kohtaa paivitetty
                # maksimi indeksi muuttuu koska uusi maksimi on suurempi kuin vanha maksimi
                minInd = np.argmin( muuttuja_Tdata)
        else:
            minInd = np.argmin( muuttuja_Tdata)
        
            

        if maksimi is not None:
            maksimi = max( maksimi,   np.max( muuttuja_Tdata ) )
            
        else:
            maksimi = np.max( muuttuja_Tdata )
        
        if minimi is not None:
            minimi = min( minimi,  np.min( muuttuja_Tdata) ) 
        else:
            minimi = np.min( muuttuja_Tdata )
        


        if EMUL:
            nimi = longName + ' ' + tag + LVLprintFig
            label = filenameNC[i].split("/")[-2]
        else:
            nimi = longName + ' ' + tag
            label = labelArray[i]
        
        fig, ax = mdp.aikasarjaTulostus( muuttuja_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = omaVari, label = label )
        #######################
    print 'muuttuja', muuttuja, 'indeksi min', minInd, 'arvo min', minimi, 'indeksi max', maksInd, 'arvo max', maksimi
    
    xticksHours = np.arange(0, max(muuttuja_Tdata)/3600. + 0.1, 0.5)
    
    if xticks is None:
        xticks = xticksHours * 3600.
    
    if xlabels is None:
        xticks = map( str, xticksHours)
    
    if yticks is not None:
        ax.set_yticks( yticks)
    
    if ylabels is not None:
        ax.set_yticklabels( ylabels )
    
    
    ax.set_xticks( xticks )
    ax.set_xticklabels( xlabels )
    
    if ((ylabelFont is not None) and (yticks is not None)):
        for takka in ( ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() ):
            takka.label.set_fontsize(ylabelFont) 
    
    #if customYscale:
        #ax.set_yscale('symlog', basey=10, linthresy=[0,1])    
    
    if spinup is not None:
        mdp.plot_vertical( spinup )

    if colorBar is not None:
       cb = plt.colorbar( colorBar, shrink=.9, pad=.03, aspect=10, ticks = colorBarTickValues )#colorBar 
       cb.ax.set_yticklabels(colorBarTickNames)
       cb.ax.set_ylabel( tit )

    # jos ymin ja ymax arvoja ei ole ennalta annettu, niin kaytetaan kuvan raja-arvoina laskettuja arvoja
    if ymin is None:
        ymin = minimi
        
    if ymax is None:
        ymax = maksimi
    
    if ( asetaRajat ):
        mdp.plot_setYlim( ymin, ymax, extendBelowZero = extendBelowZero)
    
    if savePrefix is None:  
        savePrefix = muuttuja

    if saveFig:
        plt.savefig( picturefolder + savePrefix + '_' + tag + LVLprintSave + '.png')

#########################
def piirra_profiilisettii( muuttuja, variKartta = plt.cm.gist_rainbow, variRefVektori = jarjestys, colorBar = None, colorBarTickValues = [0,1], colorBarTickNames = ['< -1', '0', '> 1'], muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'xlabel', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, xmin = None, xmax = None, relative = False, savePrefix = None, ajanhetket = None, tit = ' ', nollaArvo = None, rajaKerros = None ):

    maksimi = None
    minimi  = None
    maxRef = np.max( variRefVektori )
    maksimifracZ = None
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        
        muuttuja_data      = mdp.read_Data( filenamePS[i], muuttuja )*muunnosKerroin   
        
        #### relative    
        if relative:
            if ((nollaArvo is None) and (nono is None)):
                lahtoarvo = muuttuja_data[0]
            elif ( isinstance( nollaArvo, np.ndarray) or isinstance(nollaArvo, list) ):
                lahtoarvo = nollaArvo[i]
            elif nono is not None:
                lahtoarvo = muuttuja_data[nono]
                
            if lahtoarvo == 0.:
                    print 'VAROITUS: ei voida laskea relative koska lahtoarvo on nolla', muuttuja, 'case', i+1, 'nono', nono
                    if nono is None:
                        NEXT = raw_input('Haluatko etsia seuraavan ei-negatiivisen datasta kaikkiin kappyroihin: yes/no: ') in kylla
                    
                        if NEXT:
                            nono = 0
                            while( muuttuja_data[nono] == 0. ):
                                nono += 1
                            lahtoarvo = muuttuja_data[nono]
                        else:
                            print "Let's skip this", muuttuja, i+1
                            continue
                        
                    else:
                        print "Let's skip this", muuttuja, i+1
                        continue
                        
            muuttuja_Tdata = muuttuja_data / lahtoarvo
        #### end relative 
        
        if EMUL:
            case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
        else:
            case_indeksi = i  

        time_data   = mdp.read_Data( filenameTS[i], 'time'    )

        height_data = mdp.read_Data( filenameNC[i], 'zt' )
        
        if ajanhetket is not None:
            aikaP = np.argmin( np.abs(ajanhetket - time_data) )
        
        if aikaP > np.shape(muuttuja_data)[0] -1:
            print 'VAROITUS: indeksi liian iso, skipataan - ', 'aikaP indeksi:', aikaP, 'muuttuja_data viimeinen indeksi:',  np.shape(muuttuja_data)[0] -1, 'muuttuja', muuttuja
            continue


        #### pbl height
        if rajaKerros is None:
            pbl_height = float(raw_input('Anna rajakerroksen korkeus lahtotilanteessa (float): ') )
        elif ( isinstance(rajaKerros, np.ndarray) or isinstance(rajaKerros, list) ):
            pbl_height = rajaKerros[case_indeksi]
        #### end pbl height

        dens = 20
        
        fracZ   = np.zeros( dens*(len(height_data)-1) +1 )
        normZ   = np.zeros( len(fracZ) )
        
        pSpline = np.zeros( ( 2, len(fracZ) ) )
        
        
        for k in xrange(len(fracZ)-1):
            h_indeksi = int(np.floor(k/dens) )
            normZ[k] =  ( height_data[ h_indeksi ] + np.mod(k,dens)*(height_data[ h_indeksi + 1] - height_data[ h_indeksi ])/float(dens) )
            fracZ[k] =  normZ[k] / pbl_height
    
        normZ[-1] = height_data[-1]
        fracZ[-1] = normZ[-1] / pbl_height
        
        tck0 = interpolate.splrep( height_data, muuttuja_data[ 0,: ]     )
        tckT = interpolate.splrep(height_data, muuttuja_data[ aikaP,:]  )
        for k in xrange( np.shape(pSpline)[1] ):
            pSpline[ 0,k ]  = interpolate.splev( normZ[k], tck0 )
            pSpline[ 1,k ]  = interpolate.splev( normZ[k], tckT )
        
        p_difference = pSpline[ 1,:] / pSpline[ 0,:] -1
            

        if maksimifracZ is not None:
            maksimifracZ = max( maksimifracZ,  np.max(fracZ) )
        else:
            maksimifracZ = np.max(fracZ)
        
        if maksimi is not None:
            maksimi = max( maksimi,  np.max(p_difference) )
        else:
            maksimi = np.max(p_difference)

        if minimi is not None:
            minimi = min( minimi,  np.min(p_difference) )
        else:
            minimi = np.min(p_difference)
        
        colorMap = variKartta
        skal = variRefVektori[ case_indeksi ] / maxRef
        color = colorMap(skal)

        
        if EMUL:
            nimi = longName + ' ' + tag + LVLprintFig
            label = filenameNC[i].split("/")[-2]
        else:
            nimi = longName + ' ' + tag
            label = labelArray[i]
        fig, ax = mdp.profiiliTulostus( p_difference, aikaPisteet = 0, korkeus = fracZ, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color, label = label )

        ####################
    # tikkien fonttikoko
    for takka in ( ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() ):
                takka.label.set_fontsize(18) 


    font = 26
    ax.xaxis.get_label().set_fontsize(font)
    ax.yaxis.get_label().set_fontsize(font)
    

                
    if colorBar is not None:
       cb = plt.colorbar( colorBar, shrink=.9, pad=.03, aspect=10, ticks = colorBarTickValues )#colorBar 
       cb.ax.set_yticklabels(colorBarTickNames)
       cb.ax.set_ylabel( tit )


        # jos ymin ja ymax arvoja ei ole ennalta annettu, niin kaytetaan kuvan raja-arvoina laskettuja arvoja
    if ymin is None:
        ymin = minimi

    if ymax is None:
        ymax = maksimi

    if ( asetaRajat ):
        mdp.plot_setYlim( ymin, ymax, extendBelowZero = extendBelowZero)

    if (xmin is not None) or (xmax is not None):
        plt.xlim( xmin, xmax)
    if savePrefix is None:
        savePrefix = muuttuja

    if saveFig:
        plt.savefig( picturefolder + savePrefix + '_' + tag + LVLprintSave + '.png')




############################
def piirra_aikasarjaPathXYZ( muuttuja, longName = None, savePrefix = None, xaxislabel = 'time [h]', xlabels = None, ylabels = None, xticks = None, yticks = None, spinup = None ):
    
        
    if longName is None:
        longName = muuttuja
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False

        muuttuja_data = mdp.read_Data( filenameNC[i], muuttuja )
        time_data     = mdp.read_Data( filenameNC[i], 'time'   )
        dn0_data      = mdp.read_Data( filenameNC[i], 'dn0'    )
        zt_data       = mdp.read_Data( filenameNC[i], 'zt'     )
        zm_data       = mdp.read_Data( filenameNC[i], 'zm'     )
        korkeus = ( zm_data - zt_data )*2.0

        nimi = longName +' ' + filenameNC[i]

        fig,ax = mdp.laske_path_aikasarjaXYZ( muuttuja_data, dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xlabel = xaxislabel, tightXAxis=tightXAxis)
    
    
    ajat  =  np.round( mdp.read_Data( filenameNC[0], 'time' ), 0 )
    
    xticksHours = np.arange(0, max(ajat)/3600. + 0.1, 0.5)
    
    if xticks is None:
        xticks = xticksHours * 3600.
    
    if xlabels is None:
        xticks = map( str, xticksHours)
    
    
    ax.set_xticks( xticks )
    ax.set_xticklabels( xlabels )
    
    if spinup is not None:
        mdp.plot_vertical( spinup )
    
    if savePrefix is None:
        savePrefix = muuttuja
    if saveFig:
        plt.savefig( picturefolder + savePrefix + '_' + tag + LVLprintSave + '.png')


##########################
def piirra_maksimiKeissit( muuttuja, muunnosKerroin = 1.0, longName = 'pitka nimi', xlabel = 'case', ylabel = 'ylabel', savePrefix = None):
    lista = np.zeros( cases )
    xTikit = np.zeros( cases ) 

    
    for i in xrange(len(arguments)-1):
        muuttuja_Tdata  = mdp.read_Data( filenameTS[i], muuttuja )*muunnosKerroin
        time_data       = mdp.read_Data( filenameNC[i], 'time'   )
        if EMUL:
            tikki = int(filenameNC[i].split("/")[-2][-2:])
        else:
            
            tikki = int(i+1)

        xTikit[i] = tikki
        
        if ajanhetket is not None:
            aikaP = np.argmin( np.abs( ajanhetket - time_data) )
        else:
            aikaP = 0
        
        if aikaP > np.shape(muuttuja_Tdata)[0] -1:
            print 'VAROITUS: indeksi liian iso skipataan', 'aikaP', aikaP, 'np.shape(muuttuja_Tdata)[0]', np.shape(muuttuja_Tdata)[0], 'muuttuja', muuttuja
            continue
        #print 'np.shape(lista)', np.shape(lista), 'case_indeksi',i, 'aikaP', aikaP, 'np.shape(muuttuja_Tdata)[0]', np.shape(muuttuja_Tdata)[0], 'muuttuja', muuttuja
        
        lista[i] = np.max( muuttuja_Tdata[aikaP:] )

    fig, ax = mdp.plottaa( jarjestys, lista, tit = longName+ ' ' + tag + LVLprintFig, xl = xlabel, yl = ylabel, changeColor = True, markers = True, uusikuva = True, gridi = False, tightXAxis = False, tightYAxis = True, LEGEND = LEGEND )
    
    mdp.plot_setYlim( 0., max(lista), extendBelowZero = False, A = 0.05 )
    #print map(str, xTikit)
    ax.set_xticks( map(int, jarjestys) )
    for takka in ax.xaxis.get_major_ticks():
                takka.label.set_fontsize(9)
    #ax.set_xticklabels( map(str, xTikit) , fontsize = 8 )

    if savePrefix is None:
                savePrefix = longName.replace( " ", "_" )
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + tag + LVLprintSave + '.png')

#######################
def piirra_binijakauma( muuttuja, muunnosKerroin = 1.0, longName = None, savePrefix = None ):

    if longName is None:
        longName = muuttuja

    #for bini in xrange( np.shape(S_Rwiba_data)[1] ):
        #mdp.laske_MeanDiameterInBin( S_Rwiba_data, bini, S_Nc_data, time_data, tulostus = True, piirra =piirra, nimi = 'Mean diameter of ice particles in liquid clouds in bin ' )
        

#######################
def piirra_domainProfiili( muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = [0,1], colorBarTickNames = ['< -1', '0', '> 1'], xlabels = None, ylabels = None, xticks = None, yticks = None, variKartta = plt.cm.Blues, profiili = False, spinup = None ):
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib import colors
    
    if spinup is None:
        sys.exit("EXIT anna spinup @ piirra_domainProfiili")
    if profiili:
        tiedostonimi = filenamePS
    else:
        tiedostonimi = filenameNC
    
    
    if isinstance( variKartta, list):
        variKartta   = colors.ListedColormap(variKartta)
        bounds = colorBarTickValues
        norm   = colors.BoundaryNorm(bounds, variKartta.N)
    
    
    
    # Make plot with vertical (default) colorbar
    asp = 25./45.
    w, h = mpl.figure.figaspect( asp )
    fig, ax = plt.subplots(figsize = (24,15)) #figsize=(w, h)
    for i in xrange(len(arguments)-1):
        muuttuja_data  = np.multiply( mdp.read_Data( tiedostonimi[i], muuttuja ), muunnosKerroin)

        if useDN:
            dn00            = np.power( mdp.read_Data( filenameNC[i], 'dn0'    ), -1 ) 
            muuttuja_data = np.multiply ( muuttuja_data, dn00 )
    
        if profiili:
            muuttuja_meanProfile = muuttuja_data
        else:
            muuttuja_meanProfile = np.mean( np.mean( muuttuja_data, axis = 1), axis = 1)
            
        
        if transpose:
            muuttuja_meanProfile = np.rot90(muuttuja_meanProfile)  #zip(*muuttuja_meanProfile[::-1])#muuttuja_meanProfile.T

        cax = ax.imshow( muuttuja_meanProfile, interpolation='nearest', cmap=variKartta, aspect=asp, norm = norm ) 
        
        
    if longName is None:
        longName = muuttuja
        
    ax.set_title(longName)
    
    if xticks is None:
         sys.exit("anna xticks")
         
    #if yticks is None:
         #yticks  = np.arange( np.shape(muuttuja_data)[-1])
         
    #if xlabels is None:
         #xlabels = map(str,  np.multiply( xticks,  1./3600. ) ) 
    #if ylabels is None:
         #ylabels = map(str, mdp.read_Data( filenameNC[0], 'zt' )  )
    
    
    
    ajat  =  np.round( mdp.read_Data( filenameNC[0], 'time' ), 0 )
    oikeatXtikit = np.zeros( ( np.shape( xticks ) ))
    
    for i in xrange(len(xticks)):
        oikeatXtikit[i] = np.argmin( np.abs(ajat - 3600.*xticks[i]) )
    
    oikeatXleimat= map( int, np.multiply( ajat[[ map( int, oikeatXtikit)]], 1./3600) )

    ax.set_xticks( oikeatXtikit )
    ax.set_xticklabels( oikeatXleimat )
    
    
    
    k = 0
    for label in ax.xaxis.get_ticklabels():
        if np.mod(k,4) != 0:
            label.set_visible(False)
        k+=1
    
    
    
    
    korkeudet =  np.round( mdp.read_Data( filenameNC[0], 'zt' ), 0 ) 
    Ytikit = np.arange(0, 1500, 200 )  #np.zeros( ( np.shape( yticks ) ))
    
    #oikeatYtikit = np.zeros( np.shape(Ytikit))
    
    #yind = np.arange( np.shape( muuttuja_meanProfile)[0] )
    
    
    oikeatYleimat = map( str, Ytikit  )
    #print Ytikit
    
    #for i in xrange(len(Ytikit)):
        #oikeatYtikit[i] = np.interp( Ytikit[i], korkeudet, yind)
    #print oikeatYtikit
    
    oikeatYtikit = np.arange( np.shape( muuttuja_meanProfile)[0], 0, -20  )
    
    ax.set_yticks( oikeatYtikit )
    ax.set_yticklabels( oikeatYleimat )
    

    #j = 0
    #for label in ax.yaxis.get_ticklabels():
        #if np.mod(j,4) != 0:
            #label.set_visible(False)
        #j+=1

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixe  at 0.05 inch.
    #divider = make_axes_locatable(ax)
    #cax1 = divider.append_axes("right",  "5%", pad="3%")

    #Add colorbar, make sure to specify tick locations to match desired ticklabels
    #cbar = plt.colorbar(cax,fraction=0.046, pad=0.04, ticks = colorBarTickValues)
    #cbar = fig.colorbar(cax, ticks = colorBarTickValues) # , cax=cax1
    cbar = fig.colorbar(cax, ax=ax, shrink=.9, pad=.03, aspect=10, ticks = colorBarTickValues, norm = norm)
    cbar.ax.set_yticklabels(colorBarTickNames)  # vertically oriented colorbar
    

    

    spinupIND = np.argmin( np.abs( ajat - spinup ) )
    mdp.plot_vertical( spinupIND )
    plt.ylabel('z(m)')
    plt.xlabel('Time(h)')
    plt.tight_layout()
    
    if savePrefix is None:
        savePrefix = 'domainTimeseriesProfiili'    
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + tag + LVLprintSave + '.png')
    
###########################
#### NEEDS REVISION BIG TIME (INTERPOLATION)
######################
def piirra_hiukkausjakauma( tyyppi = 'ice', bini='a', ajanhetket = [0], korkeus = [0], useDN = True, color = 'b', savePrefix = None   ):
    from scipy.stats.mstats import gmean
    if tyyppi not in ['aerosol', 'cloud', 'precipitation', 'ice', 'snow' ]:
        sys.exit("wrong type")
        
    tyyppi1 = tyyppi[0:1]
    tyyppi2 = tyyppi[0:2]
    
    NumBin_name     = 'S_N'    + tyyppi1 + 'b' + bini
    RadBin_name     = 'S_Rw'   + tyyppi1 + 'b' + bini
    #SizBin_name =  tyyppi2 + bini
    
    
    
    for i in xrange(len(arguments)-1):
        NumBin  = mdp.read_Data( filenameNC[i], NumBin_name)
        RadBin  = mdp.read_Data( filenameNC[i], RadBin_name)
        #SizBin = mdp.read_Data( filenameNC[i], SizBin_name)
        
        if useDN:
            dn00   = np.power( mdp.read_Data( filenameNC[i], 'dn0'    ), -1 ) 
            NumBin = np.multiply ( NumBin, dn00 )
        
        ######### aika slaissaus
        time = mdp.read_Data( filenameNC[i], 'time')
        aikaindeksit = []
        for t in ajanhetket:
            aikaindeksit.append( np.argmin( np.abs(time - t*3600.) ))
        
        Tslize = map( int, np.arange( min(aikaindeksit), max(aikaindeksit)+0.5 ) )
        
        TslizeSTR = r'$t_0$' + ' = ' + str( time[ min(aikaindeksit) ]/3600. ) 
        if len(aikaindeksit)>0:
            TslizeSTR += ' ' + r'$t_1$' + ' = ' + str( time[ max(aikaindeksit) ]/3600. )
        print TslizeSTR
        ###############################
        
        ######### korkeus slaissaus
        zt = mdp.read_Data( filenameNC[i], 'zt')
        korkeusindeksit = []
        for h in korkeus:
            korkeusindeksit.append( np.argmin( np.abs(zt - h) ))
        
        Hslize = map( int, np.arange( min(korkeusindeksit), max(korkeusindeksit)+0.5 ) )
        
        HslizeSTR = r'$h_0$' + ' = ' + str( zt[ min(korkeusindeksit) ] ) 
        if len(korkeusindeksit)>0:
            HslizeSTR += + ' ' + r'$h_1$' + ' = ' + str( zt[ max(korkeusindeksit) ] )
        #print HslizeSTR
        #####################
        
        NumBinSlize  = NumBin[ Tslize, :, :, :, Hslize ]
        RadBinSlize  = RadBin[ Tslize, :, :, :, Hslize ]
        
        if len(Tslize)>1:
            print 'Tslize pituus', len(Tslize)
            NumBinSlize = np.mean( NumBinSlize, axis = 0 )
            RadBinSlize = np.mean( RadBinSlize, axis = 0 )
        if len(Hslize)>1:
            print 'Hslize pituus', len(Hslize)
            NumBinSlize = np.mean( NumBinSlize, axis = -1 )
            RadBinSlize = np.mean( RadBinSlize, axis = -1 )
        
        
        NumBinSlizeMean  = gmean( gmean(NumBinSlize, axis=1), axis=1) # #/L
        DiaBinSlizeMean  = gmean( gmean(RadBinSlize, axis=1), axis=1)*1000.*2. #  m -> mm, radius -> diameter
        

        tit = tyyppi + ' ' + TslizeSTR + ' ' + HslizeSTR
        
        mdp.plot_alustus()
        mdp.plottaa( SizBin, NumBinSlizeMean/DiaBinSlizeMean, tit , 'Diameter [mm]', r'$dN/dP(L^{-1}mm^{-1})$', log = True, changeColor=False, tightXAxis=True, markers=True, LEGEND=True, omavari = color, scatter=False )
    
    if savePrefix is None:
        savePrefix = 'Hiukkasjakauma'        
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + tag + LVLprintSave + '.png')
        
###########################
def piirra_MeanSize( tyyppi = 'ice', bini='a', ajanhetket = [0], korkeus = [0], useDN = True, color = 'b', savePrefix = None   ):
    from scipy.stats.mstats import gmean
    if tyyppi not in ['aerosol', 'cloud', 'precipitation', 'ice', 'snow' ]:
        sys.exit("wrong type")
        
    tyyppi1 = tyyppi[0:1]
    tyyppi2 = tyyppi[0:2]
    
    Rad_name = 'S_Rw'   + tyyppi1  + bini
    
    if tyyppi == 'ice':
        Num_name = 'S_Nic'
    else:
        Num_name = 'S_N'   + tyyppi1
    
    print tyyppi, Rad_name, Num_name
    
    
    for i in xrange(len(arguments)-1):
        Num  = np.multiply( mdp.read_Data( filenameNC[i], Num_name), 1./1000. )
        Rad  = mdp.read_Data( filenameNC[i], Rad_name)

        
        if useDN:
            dn00   = np.power( mdp.read_Data( filenameNC[i], 'dn0'    ), -1 ) 
            Num = np.multiply ( Num, dn00 )
            
        ######### aika slaissaus
        time = mdp.read_Data( filenameNC[i], 'time')
        aikaindeksit = []
        for t in ajanhetket:
            aikaindeksit.append( np.argmin( np.abs(time - t*3600.) ))
        
        Tslize = map( int, np.arange( min(aikaindeksit), max(aikaindeksit)+0.5 ) )
        
        TslizeSTR = r'$t_0$' + ' = ' + str( time[ min(aikaindeksit) ]/3600. ) 
        if len(aikaindeksit)>0:
            TslizeSTR += ' ' + r'$t_1$' + ' = ' + str( time[ max(aikaindeksit) ]/3600. )
        print TslizeSTR
        ###############################
        
        ######### korkeus slaissaus
        zt = mdp.read_Data( filenameNC[i], 'zt')
        korkeusindeksit = []
        for h in korkeus:
            korkeusindeksit.append( np.argmin( np.abs(zt - h) ))
        
        Hslize = map( int, np.arange( min(korkeusindeksit), max(korkeusindeksit)+0.5 ) )
        
        HslizeSTR = r'$h_0$' + ' = ' + str( zt[ min(korkeusindeksit) ] ) 
        if len(korkeusindeksit)>0:
            HslizeSTR += + ' ' + r'$h_1$' + ' = ' + str( zt[ max(korkeusindeksit) ] )
        print HslizeSTR
        #####################
        
        NumSlize  = Num[ Tslize,  :, :, Hslize ]
        RadSlize  = Rad[ Tslize,  :, :, Hslize ]
        
        if len(Tslize)>1:
            NumSlize = np.mean( NumSlize, axis = 0 )
            RadSlize = np.mean( RadSlize, axis = 0 )
        if len(Hslize)>1:
            NumSlize = np.mean( NumSlize, axis = -1 )
            RadSlize = np.mean( RadSlize, axis = -1 )

        
        #NumSlizeMean     = np.mean( np.mean(NumSlize, axis=1), axis=1) # #/L
        #DiaBinSlizeMean  = np.mean( np.mean(RadSlize, axis=1), axis=1)*1000.*2. #  m -> mm, radius -> diameter
        DiaBinSlize = RadSlize*1000.*2. #  m -> mm, radius -> diameter
        
        DwiAvg = np.mean( sum( DiaBinSlize * NumSlize ) / sum( NumSlize ) )
        NumAvg = np.mean( NumSlize )
        
        print 'DwiAvg', DwiAvg, 'NumAvg/DwiAwg', NumAvg/DwiAvg
        tit = tyyppi + ' ' + TslizeSTR + ' ' + HslizeSTR
        
        #mdp.plot_alustus()
        #mdp.plottaa( SizBin, NumSlizeMean/DiaBinSlizeMean, tit , 'Diameter [mm]', r'$dN/dP(L^{-1}mm^{-1})$', log = True, changeColor=False, tightXAxis=True, markers=True, LEGEND=True, omavari = color, scatter=False )        
    #if savePrefix is None:
        #savePrefix = 'MeanSize'
    #if saveFig:
        #plt.savefig( picturefolder + savePrefix+ '_' + tag + LVLprintSave + '.png')

#############################        
def piirra_domainMeanProfiili( muuttuja, muunnosKerroin = 1.0, ajanhetket = [0], useDN = True, profiili = False, xAxisL = '', color = 'k', savePrefix = None    ):
        
    minimi  = None
    maksimi = None
    
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        data  = np.multiply( mdp.read_Data( filenameNC[i], muuttuja), muunnosKerroin )
        zt = mdp.read_Data( filenameNC[i], 'zt')
        if useDN:
            dn00   = np.power( mdp.read_Data( filenameNC[i], 'dn0'    ), -1 ) 
            data = np.multiply ( data, dn00 )
        
        
        ######### aika slaissaus
        time = mdp.read_Data( filenameNC[i], 'time')
        aikaindeksit = []
        for t in ajanhetket:
            aikaindeksit.append( np.argmin( np.abs(time - t*3600.) ))
        
        Tslize = map( int, np.arange( min(aikaindeksit), max(aikaindeksit)+0.5 ) )
        
        TslizeSTR = r'$t_0$' + ' = ' + str( time[ min(aikaindeksit) ]/3600. ) 
        if len(aikaindeksit)>0:
            TslizeSTR += ' ' + r'$t_1$' + ' = ' + str( time[ max(aikaindeksit) ]/3600. )
        print TslizeSTR
        ###############################
        
        
        dataSlize  = data[ Tslize,  :, :, : ]
        
        if len(Tslize)>1:
            dataSlize = np.mean( dataSlize, axis = 0 )
        
        dataSlizeMean = np.mean( np.mean( dataSlize, axis = 0), axis = 0)
        
        tit = muuttuja + ' ' + TslizeSTR
        
        if maksimi is not None:
            maksimi = max( maksimi,  np.max(dataSlizeMean) )
        else:
            maksimi = np.max(dataSlizeMean)

        if minimi is not None:
            minimi = min( minimi,  np.min(dataSlizeMean) )
        else:
            minimi = np.min(dataSlizeMean)
        
        mdp.plottaa( dataSlizeMean, zt, tit , xl = xAxisL, yl='height [m]', changeColor=True, tightXAxis=True, tightYAxis = True, markers=False, LEGEND=True, omavari = color, scatter=False, uusikuva=uusikuva )       
        mdp.plot_setXlim( minimi, maksimi, extendBelowZero = False, A = 0.05 )
    
    if savePrefix is None:
        savePrefix = 'domainMeanProfiili'
    
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + tag + LVLprintSave + '.png')

def animoi_path(muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = [0,1], colorBarTickNames = ['< -1', '0', '> 1'], xlabels = None, ylabels = None, xticks = None, yticks = None, variKartta = plt.cm.Blues, profiili = False, spinup = None ):        
    
    import matplotlib.animation as animation
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib import colors
    from itertools import cycle

    if profiili:
        tiedostonimi = filenamePS
    else:
        tiedostonimi = filenameNC
    
    
    if isinstance( variKartta, list):
        variKartta   = colors.ListedColormap(variKartta)
        bounds = colorBarTickValues
        norm   = colors.BoundaryNorm(bounds, variKartta.N)
    else:
        norm = None
    
    
    # Make plot with vertical (default) colorbar
    asp = 1.
    w, h = mpl.figure.figaspect( asp )
    fig, ax = plt.subplots(figsize = (24,15)) #figsize=(w, h)

    def onClick(event):
        global pause
        pause ^= True

    #for i in xrange(len(arguments)-1):
    global muuttuja_meanDomain
    if (len(arguments)-1)>1:
        sys.exit("animation doesn't work at the moment with multiple arguments")
    
    
    muuttuja_data  = np.multiply( mdp.read_Data( tiedostonimi[0], muuttuja ), muunnosKerroin)
    dn0_data = mdp.read_Data( filenameNC[0], 'dn0'    )
    zt_data       = mdp.read_Data( filenameNC[0], 'zt'     )
    zm_data       = mdp.read_Data( filenameNC[0], 'zm'     )
    korkeus = ( zm_data - zt_data )*2.0
    
    dn0Kork  = dn0_data * korkeus
    
    muuttuja_meanDomain = np.sum( muuttuja_data * dn0Kork, axis = -1 )
    
    print 'muoto txyz', np.shape(muuttuja_data)
    print 'muoto txy', np.shape(muuttuja_meanDomain)
    
    
    
    if transpose:
        muuttuja_meanDomain = np.rot90(muuttuja_meanDomain)  #zip(*muuttuja_meanDomain[::-1])#muuttuja_meanDomain.T

    data = muuttuja_meanDomain[0,:,:]




    def simData():
        ajat  =  np.round( mdp.read_Data( filenameNC[0], 'time' ), 0 )
        indeksijoukko = cycle( np.arange(np.shape(muuttuja_meanDomain)[0]) )
        t_max = ajat[-1]
        dt = ajat[2]-ajat[1]
        t = 0.0

        while t < t_max:
            if not pause:
                k = next(indeksijoukko)
                data = muuttuja_meanDomain[k,:,:]
                t = t + dt
                
            yield data, t

    def simPoints(simData):
        data, t = simData[0], simData[1]
        time_text.set_text(time_template%(t))
        cax.set_array( data )
        return cax, time_text

    global pause
    pause = False

    cax = ax.imshow( [[],[]] , interpolation='nearest', cmap=variKartta, aspect=asp, norm = norm, animated=True ) 
    
    if longName is None:
        longName = muuttuja
        
    ax.set_title(longName)
    
    xt = mdp.read_Data( filenameNC[0], 'xt'     )
    oikeatXtikit = np.arange( np.shape( data)[1], 0, -20  )
    ax.set_yticks( oikeatXtikit  )
    ax.set_xticklabels( map(str, xt) )

    h = 0
    for label in ax.xaxis.get_ticklabels():
        if np.mod(h,20) != 0:
            label.set_visible(False)
        h+=1
    

    yt = mdp.read_Data( filenameNC[0], 'yt'     )
    
    oikeatYtikit = np.arange( np.shape(data)[0], 0, -20  )
    ax.set_yticks( oikeatYtikit  )
    ax.set_yticklabels( map(str, yt) )

    j = 0
    for label in ax.yaxis.get_ticklabels():
        if np.mod(j,6) != 0:
            label.set_visible(False)
        j+=1


    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    cbar = fig.colorbar(cax, ax=ax, shrink=.45, pad=.03, aspect=10, ticks = colorBarTickValues, norm = norm)
    cbar.ax.set_yticklabels(colorBarTickNames)

    time_template = 'Time = %.1f s'
    time_text = ax.text(-1.00, 0.9, '', transform=ax.transAxes )
    fig.canvas.mpl_connect('button_press_event', onClick)
    ani = animation.FuncAnimation(fig, simPoints, simData, blit=False, interval=10, repeat=True)
      
        
##############################
#
# drawing STYLES  
#
# parameter settings 
#
##############################

if EMUL and not importOnly:
    CFRAC = True
    CLOUD = True
    PRCP  = True
    LWP   = True
    SHF   = False
    LHF   = False
    QDIFF = True
    TDIFF = True
    CDNC  = True
    WMAX  = True
    
    
    #CFRAC = False
    #CLOUD = False
    #PRCP  = False
    #LWP   = False
    #SHF   = False
    #LHF   = False
    #QDIFF = False
    #TDIFF = False
    #CDNC  = True

    NCFwrite = False

    emulatorname = os.getcwd().split("/")[-1]

    #print 'emulatorname', emulatorname
    vers = emulatorname[21:27]
    design_input = 'Do you want to give an other design version other than: ' + vers + ' (yes/no)? '
    designLogical = raw_input( design_input ) not in kylla
    
    if not designLogical:
        vers= str( raw_input( 'Give a design (esim: v1.5.1 ): '))
    
    if int(lvl) < 4:
        CDNC = False
        prcp = 'prcp'
        sadekerroin = 1.
    else:
        prcp = 'rmH2Opr'
        sadekerroin = 2.5e+06 # latent heat of vaporization J/kg



    maksimisateet   = np.zeros(cases)

    maksimiLatent   = np.zeros(cases)

    maksimiSensible = np.zeros(cases)

    maksimiLWPlist  = np.zeros(cases)

    designroot = ibrix + '/DESIGN/'

    cwd = os.getcwd()

    os.chdir(designroot)


    thickness_color = plt.cm.gist_rainbow #Paired Blues
    pblh_color      = plt.cm.cool
    num_pbl_color   = plt.cm.Wistia

    for file in glob.glob("*.csv"):
        designbasename=file

    filu = designroot + designbasename
    os.chdir(cwd)



    #####################
    ncfolder = ibrix+'/DESIGNnetcdf/'
    print "DESIGN VERSION", vers
    ncfilename = ncfolder + 'design_'+vers + '.nc'
    print 'ncfilename', ncfilename
    ncfile = Dataset( ncfilename, 'r+' )
    
    q_inv_design     = ncfile.variables['q_inv'][:]
    tpot_inv_design  = ncfile.variables['tpot_inv'][:]
    clw_max_design   = ncfile.variables['clw_max'][:]
    tpot_pbl_design  = ncfile.variables['tpot_pbl'][:]
    pblh_design      = ncfile.variables['pblh'][:]
    num_pbl_design   = ncfile.variables['num_pbl'][:]
    q_pbl_design     = ncfile.variables['q_pbl'][:]
    cloudbase_design = ncfile.variables['cloudbase'][:]
    thickness_design = ncfile.variables['thickness'][:]
    
    ##########
    thickBAR = varibaari( thickness_design, thickness_color )
    pblhBAR  = varibaari( pblh_design,      pblh_color      )
    cdncBAR  = varibaari( num_pbl_design,    num_pbl_color   )
    thickTIT = 'cloud thickness'
    pblhTIT  = 'PBL height'
    cdncTIT  = 'CDNC #/kg'
    
    cbvalThick =  map( int, np.arange( myRound(min(thickness_design), 100),  myRound(max(thickness_design), 100), 100.) )  #map( int, np.arange( 0, myRound( max(thickness_design), 200 ), 200.) )
    cbvalThickStr = map( str, cbvalThick )
    
    cbvalPblh = map( int, np.arange( myRound(min(pblh_design), 100),  myRound(max(pblh_design), 100), 200.) )
    cbvalPblhStr = map( str, cbvalPblh )
    
    cbvalQ = map(int, np.arange( 0, max(num_pbl_design), 50.) )
    cbvalQStr = map( str, cbvalQ )
    
    zcTicks = [ 0, 0.5, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.20, 1.25, 1.3, 1.35, 1.4]
    zcTicksStr = map( str, zcTicks)

    
    ##########

    #mdp.plot_alustus()
    
    # thickness histogram
    #plt.hist(thickness_design, 20, normed=1, facecolor='green', alpha=0.75)
    # thickBAR
    if CFRAC:
        piirra_aikasarjasettii( muuttuja = 'cfrac',  variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName= 'Cloud fraction',               ylabel = 'Cloud fraction',      variKartta = thickness_color, ymin = 0.0, ymax = 1.0, tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds )
    
    if CLOUD:
        piirra_aikasarjasettii( muuttuja = 'zc',     variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName= 'Relative change of cloud top', ylabel = 'relative change',     variKartta = thickness_color, relative = True, savePrefix = 'cloud_top_rel_change', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, yticks = zcTicks, ylabelFont = 18  )     
        
    if PRCP:
        piirra_aikasarjasettii( muuttuja = prcp,     muunnosKerroin = sadekerroin, variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Surface precipitation',       ylabel = 'Precipitation W/m^2', variKartta = thickness_color, ymin = 0.0, savePrefix = 'prcp', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds )
        
        piirra_maksimiKeissit( muuttuja = prcp, muunnosKerroin = sadekerroin, longName = 'Maximum precipitation after '+str(int(ajanhetket/3600.))+'h', ylabel = 'precipitation W/m^2',      savePrefix = 'prcp_max' )
    
    if SHF:
        piirra_aikasarjasettii( muuttuja = 'shf_bar', variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Sensible heat flux',        ylabel = 'Sensible heat flux W/m^2', variKartta = thickness_color, savePrefix = 'heat_flux_sensible', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds )
        
        piirra_maksimiKeissit( maksimiSensible, longName = "Maximum sensible heat",                                       ylabel = 'Sensible heat flux W/m^2', savePrefix = 'heat_flx_sensible_max' ) 
    
    if LHF:
        piirra_aikasarjasettii( muuttuja = 'lhf_bar', variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Latent heat flux',        ylabel = 'Latent heat flux W/m^2', variKartta = thickness_color, savePrefix = 'heat_flux_latent', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds )
        
        piirra_maksimiKeissit( maksimiLatent,   longName = "Maximum latent heat",                                         ylabel = 'Latent heat flux W/m^2',   savePrefix = 'heat_flx_latent_max' ) 

    if LWP:
        piirra_aikasarjasettii( muuttuja = 'lwp_bar', muunnosKerroin = 1000.0, variKartta = thickness_color,  variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'LWP', ylabel = 'LWP g/m^2', ymin = 0.0,  savePrefix = 'lwp', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds )

    if CDNC:
        piirra_aikasarjasettii( muuttuja = 'Nc_ic', variKartta = num_pbl_color, variRefVektori = num_pbl_design, colorBar = cdncBAR, colorBarTickValues = cbvalQ, colorBarTickNames = cbvalQStr, longName = 'Relative change of in-cloud CDNC', ylabel = 'Relative change of In-cloud CDNC #/kg', tit = cdncTIT, relative = True, nollaArvo = 1.e6*num_pbl_design,  spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds) #

    if QDIFF:
        piirra_profiilisettii( 'q', variKartta = pblh_color, variRefVektori = pblh_design, colorBar = pblhBAR, colorBarTickValues = cbvalPblh, colorBarTickNames = cbvalPblhStr, longName =  'Tot. ' + r'$H_2{O}$' + ' mix. rat. ' + ' relative change ' +'0h - '+str((ajanhetket/3600.)) + 'h', xlabel = r'$\frac{q_{t=2h}}{q_{t=0h}}-1$', ylabel = 'z/pblh', ymin = 0.0, ymax = 1.01, xmin = -0.35, xmax = 0.025, savePrefix = 'q_diff', ajanhetket = ajanhetket, tit = pblhTIT, rajaKerros = pblh_design, relative = True, nollaArvo = q_pbl_design )

    if TDIFF:
        piirra_profiilisettii( 't', variKartta = pblh_color, variRefVektori = pblh_design, colorBar = pblhBAR, colorBarTickValues = cbvalPblh, colorBarTickNames = cbvalPblhStr, longName =  r'$\theta$' + ' relative change 0h - '+str((ajanhetket/3600.)) + 'h', xlabel = r'$\frac{\theta_{t=2h}}{\theta_{t=0h}}-1$', ylabel = 'z/pblh', ymin = 0.0, ymax = 1.01, xmin = -0.05, xmax = 0.025, savePrefix = 't_diff', ajanhetket = ajanhetket, tit = pblhTIT, rajaKerros = pblh_design, relative = True, nollaArvo = tpot_pbl_design )
    if WMAX:
        piirra_aikasarjasettii( muuttuja = 'wmax', variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Maximum vertical velocity',        ylabel = 'Maximum vertical velocity m/s', variKartta = thickness_color, savePrefix = 'w_max', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds )
    

    if NCFwrite:
        maksimiLWPlist_ncf    = ncfile.createVariable( 'lwp_max',     np.dtype('float32').char, ('case') )
        maksimiLWPlist_ncf[:] = maksimiLWPlist
        print 'eniten sadetta keissi', np.argmax(maksimiLWPlist)
    ncfile.close()
################

if ICE and not importOnly:
    

    
    korkeustikit =np.arange(0, 1250, 50)
    ylabels    = map(str, korkeustikit )
    
    profiiliVariLIQ = [ '#000099', '#00ccff', '#00e600', '#f9f906', '#ff9900', '#ff0000' ]
    profiiliVariICE = [ '#000099', '#00ccff', '#29a385', '#00e600', '#f9f906', '#f9f906', '#ff9900', '#ff0000', '#990000', '#660000' ]
    cbvalLIQ    = np.arange(0, 0.241, 0.04)
    cbvalLIQStr = map(str, cbvalLIQ)

    if tag[:-1]=='ice1':
        cbvalICE    = np.arange( 0, 0.51, 0.05) # np.arange(0, 1.4, 0.1)
    else:
        cbvalICE    =  np.arange(0, 1.41, 0.1)
    cbvalICEStr = map(str, cbvalICE)
    
    cbvalLIQPATH = np.arange(0, 61, 5)
    cbvalLIQPATHStr = map(str, cbvalLIQPATH)
    
    if int(lvl)>=4:
    
        piirra_aikasarjaPathXYZ( 'l', longName = 'Liquid Water Path', savePrefix = 'lwp', xaxislabel = 'time [h]' )
        #mdp.plot_vertical( spinup )
        #plt.xticks( ticksHours, xLabesHours )
        
        ##piirra_aikasarjasettii( muuttuja = 'lwp_bar', muunnosKerroin = 1000.0, longName = 'LWP', ylabel = 'LWP g/m^2', ymin = 0.0,  savePrefix = 'lwpTS', omaVari = False, xlabel = 'time [h]', spinup = spinup )
        #mdp.plot_vertical( spinup )
        #plt.xticks( ticksHours, xLabesHours )
        
        piirra_domainProfiili( 'l', muunnosKerroin = 1000., longName = tag + "Liquid water mixing ratio  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQ, colorBarTickNames = cbvalLIQStr, xlabels = xLabesHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = profiiliVariLIQ, spinup = spinup )
        
        #animoi_path( 'l', muunnosKerroin = 1000., longName = tag + "Liquid water path  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQPATH, colorBarTickNames = cbvalLIQPATHStr, xlabels = xLabesHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.Reds, spinup = spinup )
        
        #piirra_MeanSize(tyyppi = 'cloud', ajanhetket = [6], korkeus = [700], color = 'r')
        
        #piirra_domainMeanProfiili( 'S_Nc',  muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        
        #piirra_domainProfiili( 'w_2', longName = tag + "vertical velocity squared " + r'$m^{2}/s^{-2}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xLabesHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.RdPu, profiili = True, spinup = spinup )
        
    if int(lvl)>= 5:
        piirra_aikasarjaPathXYZ( 'f', longName = 'Ice Water Path', savePrefix = 'iwp', xaxislabel = 'time [h]',xlabels = xLabesHours, ylabels = ylabels, xticks = ticksHours, spinup = spinup )
        mdp.plot_vertical( spinup )
        #plt.xticks( ticksHours, xLabesHours )
    

        piirra_domainProfiili( 'f', muunnosKerroin = 1000.*np.power(10.,2), longName = tag + 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xLabesHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup  )
        
        #animoi_path( 'f', muunnosKerroin = 1000., longName = tag + "Ice water path  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQPATH, colorBarTickNames = cbvalLIQPATHStr, xlabels = xLabesHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.Blues, spinup = spinup )
        ## muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = [0,1], colorBarTickNames = ['< -1', '0', '> 1'], xlabels = None, ylabels = None, xticks = None, yticks = None
        
        piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [700] )
        piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [400] )
        piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [200] )
        piirra_domainMeanProfiili( 'S_Nic',  muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        piirra_domainMeanProfiili( 'S_Rwia', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$D [{\mu}m]$', color = icevari )   


toc = time.clock()
print toc - tic
########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra and not importOnly:
    mdp.plot_lopetus()