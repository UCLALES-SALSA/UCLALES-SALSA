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
    arguments = sys.argv


    filenameNC = []
    filenamePS = []
    filenameTS = []

    for filebase in arguments[1:]:

        filenameNC.append( filebase + '.nc'    )
        filenamePS.append( filebase + '.ps.nc' )
        filenameTS.append( filebase + '.ts.nc' )


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
    plt.clf()
    return CS3

##########################
def piirra_aikasarjasettii( muuttuja, variKartta = plt.cm.gist_rainbow, variRefVektori = jarjestys, colorBar = None, muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'time [s]', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, relative = False, savePrefix = None, ajanhetket = None, omaVari = True, tit = None ):

    maksimi = None
    minimi  = None
    maxRef = np.max( variRefVektori )
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        time_data    = mdp.read_Data( filenameTS[i], 'time'  )
        muuttuja_Tdata  = mdp.read_Data( filenameTS[i], muuttuja )*muunnosKerroin

        if relative:
            muuttuja_Tdata = muuttuja_Tdata / muuttuja_Tdata[0]

        if EMUL:
            case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
        else:
            case_indeksi = i
        
        
        colorMap = variKartta
        skal = variRefVektori[ case_indeksi ] / maxRef
        color = colorMap(skal)
        
        if omaVari:
            omaVari = color


        if maksimi is not None:
            maksimi = max( maksimi,   np.max( muuttuja_Tdata ) )
        else:
            maksimi = np.max( muuttuja_Tdata )
        
        if minimi is not None:
            minimi = min( minimi,  np.min( muuttuja_Tdata) ) 
        else:
            minimi = np.min( muuttuja_Tdata )

        if LEGEND and EMUL:
            nimi = longName + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
        else:
            nimi = longName + ' ' + tag + ' ' + 'LVL' + lvl

        mdp.aikasarjaTulostus( muuttuja_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = omaVari )
        #######################

    if colorBar is not None:
       cb = plt.colorbar( colorBar )
    if tit is not None:
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
        plt.savefig( picturefolder + savePrefix + '_' + tag + '_' + 'LVL'+lvl + '.png')

#########################
def piirra_profiilisettii( muuttuja, variKartta = plt.cm.gist_rainbow, variRefVektori = jarjestys, colorBar = None, muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'xlabel', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, xmin = None, xmax = None, relative = False, savePrefix = None, ajanhetket = None, tit = None):

    maksimi = None
    minimi  = None
    maxRef = np.max( variRefVektori )
    maksimifracZ = None
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        
        muuttuja_data      = mdp.read_Data( filenamePS[i], muuttuja )*muunnosKerroin   
        if relative:
            muuttuja_Tdata = muuttuja_Tdata / muuttuja_Tdata[0]

        time_data   = mdp.read_Data( filenameTS[i], 'time'    )

        height_data = mdp.read_Data( filenameNC[i], 'zt' )
        
        if ajanhetket is not None:
            aikaP = np.argmin( np.abs(ajanhetket - time_data) )
        
        
        if EMUL:
            case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
        else:
            case_indeksi = i  # tama ei toimi

        dens = 20
        
        fracZ   = np.zeros( dens*(len(height_data)-1) +1 )
        normZ   = np.zeros( len(fracZ) )
        
        pSpline = np.zeros( ( 2, len(fracZ) ) )
        
        
        for k in xrange(len(fracZ)-1):
            h_indeksi = int(np.floor(k/dens) )
            normZ[k] =  ( height_data[ h_indeksi ] + np.mod(k,dens)*(height_data[ h_indeksi + 1] - height_data[ h_indeksi ])/float(dens) )
            fracZ[k] =  normZ[k] / pblh_design[case_indeksi]
    
        normZ[-1] = height_data[-1]
        fracZ[-1] = normZ[-1] / pblh_design[case_indeksi]
        
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

        
        if LEGEND and EMUL:
            nimi = longName + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
        else:
            nimi = longName + ' ' + tag + ' ' + 'LVL' + lvl

        mdp.profiiliTulostus( p_difference, aikaPisteet = 0, korkeus = fracZ, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )
        ####################
        
    if colorBar is not None:
        cb = plt.colorbar( colorBar )
    if tit is not None:
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
        plt.savefig( picturefolder + savePrefix + '_' + tag + '_' + 'LVL'+lvl + '.png')




############################
def piirra_aikasarjaPathXYZ( muuttuja, longName = None, savePrefix = None, xlabel = 'time [s]' ):
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

        mdp.laske_path_aikasarjaXYZ( muuttuja_data, dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xlabel = xlabel, tightXAxis=tightXAxis)
    if savePrefix is None:
        savePrefix = muuttuja
    if saveFig:
        plt.savefig( picturefolder + savePrefix + '_' + tag + '_' + 'LVL'+lvl + '.png')


##########################
def piirra_maksimiKeissit( muuttuja, muunnosKerroin = 1.0, longName = 'pitka nimi', xlabel = 'case', ylabel = 'ylabel', savePrefix = None):
    lista = np.zeros( cases )
    mdp.plot_alustus()
    for i in xrange(len(arguments)-1):
        muuttuja_Tdata  = mdp.read_Data( filenameTS[i], muuttuja )*muunnosKerroin
        time_data       = mdp.read_Data( filenameNC[i], 'time'   )
        if EMUL:
            case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
        else:
            case_indeksi = i

        if ajanhetket is not None:
            aikaP = np.argmin( np.abs( ajanhetket - time_data) )
        else:
            aikaP = 0

        lista[case_indeksi] = np.max( muuttuja_Tdata[aikaP:] )

    mdp.plottaa( jarjestys, lista, tit = longName+ ' ' + tag + ' ' + 'LVL' + lvl, xl = xlabel, yl = ylabel, changeColor = True, markers = True )

    if savePrefix is None:
                savePrefix = longName.replace( " ", "_" )
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + tag + '_' + 'LVL'+lvl + '.png')

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
    cbar = fig.colorbar(cax, ax=ax, shrink=.45, pad=.03, aspect=10, ticks = colorBarTickValues, norm = norm)
    cbar.ax.set_yticklabels(colorBarTickNames)  # vertically oriented colorbar
    

    

    spinupIND = np.argmin( np.abs( ajat - spinup ) )
    mdp.plot_vertical( spinupIND )
    plt.ylabel('z(m)')
    plt.xlabel('Time(h)')
    plt.tight_layout()
    
    if savePrefix is None:
        savePrefix = 'domainTimeseriesProfiili'    
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + tag + '_' + 'LVL'+lvl + '.png')
    
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
        print HslizeSTR
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
        plt.savefig( picturefolder + savePrefix+ '_' + tag + '_' + 'LVL'+lvl + '.png')
        
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
        #plt.savefig( picturefolder + savePrefix+ '_' + tag + '_' + 'LVL'+lvl + '.png')
        
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
        #print TslizeSTR
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
        plt.savefig( picturefolder + savePrefix+ '_' + tag + '_' + 'LVL'+lvl + '.png')
        
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

    NCFwrite = False

    emulatorname = os.getcwd().split("/")[-1]


    tag = emulatorname[21:27]

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
    pblh_color      = plt.cm.Reds
    num_pbl_color   = plt.cm.Oranges

    for file in glob.glob("*.csv"):
        designbasename=file

    filu = designroot + designbasename
    os.chdir(cwd)



    #####################
    ncfolder = ibrix+'/DESIGNnetcdf/'
    print "DESIGN VERSION", tag
    ncfile = Dataset( ncfolder + 'design_'+tag + '.nc', 'r+' )
    
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
    pblhTIT  = 'pblh'
    cdncTIT  = 'num_pbl'
    ##########

    mdp.plot_alustus()
    plt.hist(thickness_design, 20, normed=1, facecolor='green', alpha=0.75)
    
    if CFRAC:
        piirra_aikasarjasettii( muuttuja = 'cfrac',  variRefVektori = thickness_design, colorBar = thickBAR, longName= 'Cloud fraction',               ylabel = 'Cloud fraction',      variKartta = thickness_color, ymin = 0.0, ymax = 1.0, tit = thickTIT )
    
    if CLOUD:
        piirra_aikasarjasettii( muuttuja = 'zc',     variRefVektori = thickness_design, colorBar = thickBAR, longName= 'Relative change of cloud top', ylabel = 'relative change',     variKartta = thickness_color, relative = True, savePrefix = 'cloud_top_rel_change', tit = thickTIT )     

    if PRCP:
        piirra_aikasarjasettii( muuttuja = prcp,     muunnosKerroin = sadekerroin, variRefVektori = thickness_design, colorBar = thickBAR, longName = 'Surface precipitation',       ylabel = 'Precipitation W/m^2', variKartta = thickness_color, ymin = 0.0, ajanhetket = ajanhetket, savePrefix = 'prcp', tit = thickTIT )
        piirra_maksimiKeissit( muuttuja = prcp, muunnosKerroin = sadekerroin, longName = 'Maximum precipitation after '+str(int(ajanhetket/3600.))+'h', ylabel = 'precipitation W/m^2',      savePrefix = 'prcp_max' )
    
    if SHF:
        piirra_aikasarjasettii( muuttuja = 'shf_bar', variRefVektori = thickness_design, colorBar = thickBAR, longName = 'Sensible heat flux',        ylabel = 'Sensible heat flux W/m^2', variKartta = thickness_color, savePrefix = 'heat_flux_sensible', tit = thickTIT )
        piirra_maksimiKeissit( maksimiSensible, longName = "Maximum sensible heat",                                       ylabel = 'Sensible heat flux W/m^2', savePrefix = 'heat_flx_sensible_max' ) 
    
    if LHF:
        piirra_aikasarjasettii( muuttuja = 'lhf_bar', variRefVektori = thickness_design, colorBar = thickBAR, longName = 'Latent heat flux',        ylabel = 'Latent heat flux W/m^2', variKartta = thickness_color, savePrefix = 'heat_flux_latent', tit = thickTIT )
        piirra_maksimiKeissit( maksimiLatent,   longName = "Maximum latent heat",                                         ylabel = 'Latent heat flux W/m^2',   savePrefix = 'heat_flx_latent_max' ) 

    if LWP:
        piirra_aikasarjasettii( muuttuja = 'lwp_bar', muunnosKerroin = 1000.0, variKartta = thickness_color,  variRefVektori = thickness_design, colorBar = thickBAR, longName = 'LWP', ylabel = 'LWP g/m^2', ymin = 0.0,  savePrefix = 'lwp', tit = thickTIT )

    if CDNC:
        piirra_aikasarjasettii( muuttuja = 'Nc_ic', relative = True, variKartta = num_pbl_color, variRefVektori = num_pbl_design, longName = 'Relative change of In-cloud CDNC #/kg', ylabel = 'Relative change of In-cloud CDNC #/kg', tit = cdncTIT)

    if QDIFF:
        piirra_profiilisettii( 'q', variKartta = pblh_color, variRefVektori = pblh_design, colorBar = pblhBAR, longName =  'Total water mix. rat g/m^2 difference between 0h & '+str((ajanhetket/3600.)) + 'h', xlabel = r'$\frac{q_{t=2h}}{q_{t=0h}}-1$', ylabel = 'z/pblh', ymin = 0.0, ymax = 1.01, xmin = -0.3, xmax = 0.025, savePrefix = 'q_diff', ajanhetket = ajanhetket, tit = pblhTIT )

    if TDIFF:
        piirra_profiilisettii( 't', variKartta = pblh_color, variRefVektori = pblh_design, colorBar = pblhBAR, longName =  'Temperature [K] difference between 0h & '+str((ajanhetket/3600.)) + 'h', xlabel = r'$\frac{t_{t=2h}}{t_{t=0h}}-1$', ylabel = 'z/pblh', ymin = 0.0, ymax = 1.01, xmin = -0.05, xmax = 0.025, savePrefix = 't_diff', ajanhetket = ajanhetket, tit = pblhTIT )

    

    if NCFwrite:
        maksimiLWPlist_ncf    = ncfile.createVariable( 'lwp_max',     np.dtype('float32').char, ('case') )
        maksimiLWPlist_ncf[:] = maksimiLWPlist
        print 'eniten sadetta keissi', np.argmax(maksimiLWPlist)
    ncfile.close()
################

if ICE and not importOnly:
    
    spinup = mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'Tspinup' )
    

    tuntileimat = np.arange(0, 8.5, 0.5)
    xlabels    = map(str, tuntileimat )
    
    korkeusleimat =np.arange(0, 1250, 50)
    ylabels    = map(str, korkeusleimat )
    
    profiiliVariLIQ = [ '#000099', '#00ccff', '#00e600', '#f9f906', '#ff9900', '#ff0000' ]
    profiiliVariICE = [ '#000099', '#00ccff', '#29a385', '#00e600', '#f9f906', '#f9f906', '#ff9900', '#ff0000', '#990000', '#660000' ]
    cbvalLIQ    = np.arange(0, 0.241, 0.04)
    cbvalLIQStr = map(str, cbvalLIQ)

    if tag[:-1]=='ice1':
        cbvalICE    = np.arange( 0, 0.51, 0.05) # np.arange(0, 1.4, 0.1)
    else:
        cbvalICE    =  np.arange(0, 1.41, 0.1)
    cbvalICEStr = map(str, cbvalICE)
    
    if int(lvl)>=4:
    
        #piirra_aikasarjaPathXYZ( 'l', longName = 'Liquid Water Path', savePrefix = 'lwp', xlabel = 'time [h]' )
        #mdp.plot_vertical( spinup )
        #plt.xticks( xtikit, xlabels )
        
        ##piirra_aikasarjasettii( muuttuja = 'lwp_bar', muunnosKerroin = 1000.0, longName = 'LWP', ylabel = 'LWP g/m^2', ymin = 0.0,  savePrefix = 'lwpTS', omaVari = False, xlabel = 'time [h]' )
        #mdp.plot_vertical( spinup )
        #plt.xticks( xtikit, xlabels )
        
        piirra_domainProfiili( 'l', muunnosKerroin = 1000., longName = "ice1 Liquid water mixing ratio  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQ, colorBarTickNames = cbvalLIQStr, xlabels = xlabels, ylabels = ylabels, xticks = tuntileimat, yticks = korkeusleimat,  variKartta = profiiliVariLIQ, spinup = spinup )
        
        #piirra_MeanSize(tyyppi = 'cloud', ajanhetket = [6], korkeus = [700], color = 'r')
        
        piirra_domainMeanProfiili( 'S_Nc',  muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        
        #piirra_domainProfiili( 'w_2', longName = tag + "vertical velocity squared " + r'$m^{2}/s^{-2}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xlabels, ylabels = ylabels, xticks = tuntileimat, yticks = korkeusleimat,  variKartta = plt.cm.RdPu, profiili = True, spinup = spinup )
        
    if int(lvl)>= 5:
        #piirra_aikasarjaPathXYZ( 'f', longName = 'Ice Water Path', savePrefix = 'iwp', xlabel = 'time [h]' )
        #mdp.plot_vertical( spinup )
        #plt.xticks( xtikit, xlabels )
    

        piirra_domainProfiili( 'f', muunnosKerroin = 1000.*np.power(10.,2), longName = tag + 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xlabels, ylabels = ylabels, xticks = tuntileimat, yticks = korkeusleimat, variKartta = profiiliVariICE, spinup = spinup  )
        

        ## muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = [0,1], colorBarTickNames = ['< -1', '0', '> 1'], xlabels = None, ylabels = None, xticks = None, yticks = None
        
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [700] )
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [400] )
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [200] )
        #piirra_domainMeanProfiili( 'S_Nic',  muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        #piirra_domainMeanProfiili( 'S_Rwia', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$D [{\mu}m]$', color = icevari )   


toc = time.clock()
print toc - tic
########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra and not importOnly:
    mdp.plot_lopetus()


























































































###################################
#if CLOUD:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#
#        time_data    = mdp.read_Data( filenameTS[i], 'time'    )
#        top_Tdata    = mdp.read_Data( filenameTS[i], 'zc' )
#        top_Tdata = top_Tdata/top_Tdata[0]
#
#        case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
#        
#        colorMap = thickness_color
#        skal = thickness_design[ case_indeksi ] / np.max( thickness_design)
#        color = colorMap(skal)
#
#        if 'maksimiCLOUD' in locals():
#            maksimiCLOUD = max( maksimiCLOUD,   np.max( top_Tdata ) )
#        else:
#            maksimiCLOUD = np.max( top_Tdata )
#        
#        if 'minimiCLOUD' in locals():
#            minimiCLOUD = min( minimiCLOUD,  np.min(top_Tdata) ) 
#        else:
#            minimiCLOUD = np.min( top_Tdata )
#            
#
#        if LEGEND:
#            nimi = 'Relative change of cloud top' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi = 'Relative change of cloud top' + ' ' + tag + ' ' + 'LVL' + lvl
#            
#        mdp.aikasarjaTulostus(  top_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva,    nimi = nimi, xnimi = 'time [s]', ynimi = 'relative change', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )
#
#
#    mdp.plot_setYlim( minimiCLOUD, maksimiCLOUD, extendBelowZero = True)
#
#if saveFig:
#    plt.savefig( picturefolder + 'change_of_cloud_top'+ '_' + tag + '_' + 'LVL'+lvl + '.png')
###################################
#if CFRAC:
    #for i in xrange(len(arguments)-1):
        #uusikuva = True if i == 0 else  False
        #time_data    = mdp.read_Data( filenameTS[i], 'time'  )
        #cfrac_Tdata  = mdp.read_Data( filenameTS[i], 'cfrac' )
        
        #case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
        
        #colorMap = thickness_color
        #skal = thickness_design[ case_indeksi ] / np.max( thickness_design)
        #color = colorMap(skal)

        #if LEGEND:
            #nimi = 'Cloud fraction' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
        #else:
            #nimi = 'Cloud fraction' + ' ' + tag + ' ' + 'LVL' + lvl

        #if uusikuva:
            #fig, ax = mdp.aikasarjaTulostus( cfrac_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'cloud fraction', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )
        #else:
            #mdp.aikasarjaTulostus( cfrac_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'cloud fraction', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )

    ##cax.imshow(
    #mdp.plot_setYlim( 0.0, 1.0, extendBelowZero = True)

    
    #if saveFig:
        #plt.savefig( picturefolder + 'cfrac'+ '_' + tag + '_' + 'LVL'+lvl + '.png')

###################################
#if PRCP:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#        time_data    = mdp.read_Data( filenameTS[i], 'time'    )
#        prcp_Tdata  = mdp.read_Data( filenameTS[i], prcp )*sadekerroin
#        
#        aikaP = np.argmin( np.abs(ajanhetket - time_data) )
#        
#        case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
#        
#        colorMap = thickness_color
#        skal = thickness_design[ case_indeksi ] / np.max( thickness_design)
#        color = colorMap(skal)        
#    
#        maksimisateet[case_indeksi] = np.max(prcp_Tdata[aikaP:])
#        if 'maksimiPRCP' in locals():
#            maksimiPRCP = max( maksimiPRCP,   np.max(prcp_Tdata) ) 
#        else:
#            maksimiPRCP = np.max(prcp_Tdata)
#
#        if LEGEND:
#            nimi = 'Surface precipitation' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi = 'Surface precipitation' + ' ' + tag + ' ' + 'LVL' + lvl
#
#        mdp.aikasarjaTulostus( prcp_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'precipitation W/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color   )
#
#
#    mdp.plot_setYlim( 0.0, maksimiPRCP, extendBelowZero = True)
#
#    if saveFig:
#        plt.savefig( picturefolder + 'prcp'+ '_' + tag + '_' + 'LVL'+lvl + '.png')
#print ' '
#print 'SENSIBLE'
###################################
#if SHF:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#
#        time_data = mdp.read_Data( filenameTS[i], 'time'    )
#        shf_Tdata = mdp.read_Data( filenameTS[i], 'shf_bar' )
#
#        case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
#        maksimiSensible[ case_indeksi ] = np.max(shf_Tdata)
#        
#        if 'maksimiSHF' in locals():
#            maksimiSHF = max( maksimiSHF,  np.max(shf_Tdata) )
#        else:
#            maksimiSHF = np.max(shf_Tdata)
#        
#        if 'minimiSHF' in locals():
#            minimiSHF = min( minimiSHF,  np.min(shf_Tdata) )
#        else:
#            minimiSHF = np.min(shf_Tdata)
#
#        if LEGEND:
#            nimi = 'Sensible heat flux' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi = 'Sensible heat flux' + ' ' + tag + ' ' + 'LVL' + lvl
#            
#        mdp.aikasarjaTulostus( shf_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Sensible heat flux W/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND )
#        
#    
#    mdp.plot_setYlim( minimiSHF, maksimiSHF )
#    if saveFig:
#        plt.savefig( picturefolder + 'heat_flux_sensible'+ '_' + tag + '_' + 'LVL'+lvl + '.png')
#
#print ' '
#print 'LATENT'
###################################
#if LHF:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#
#        time_data    = mdp.read_Data( filenameTS[i], 'time'    )
#        lhf_Tdata = mdp.read_Data( filenameTS[i], 'lhf_bar' )
#
#        case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
#        maksimiLatent[ case_indeksi ] = np.max( lhf_Tdata )
#        if 'maksimiLHF' in locals():
#            maksimiLHF = max( maksimiLHF,  np.max(lhf_Tdata) )
#        else:
#            maksimiLHF = np.max(lhf_Tdata)
#        
#        if 'minimiLHF' in locals():
#            minimiLHF = min( minimiLHF,  np.min(lhf_Tdata) )
#        else:
#            minimiLHF = np.min(lhf_Tdata)
#
#        if LEGEND:
#            nimi = 'Latent heat flux' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi ='Latent heat flux' + ' ' + tag + ' ' + 'LVL' + lvl
#            
#        mdp.aikasarjaTulostus( lhf_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Latent heat flux W/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND )    
#
#    mdp.plot_setYlim( minimiLHF, maksimiLHF )
#    if saveFig:
#        plt.savefig( picturefolder + 'heat_flux_latent'+ '_' + tag + '_' + 'LVL'+lvl + '.png')
#
#############################################
#if LWP:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#        
#        time_data    = mdp.read_Data( filenameTS[i], 'time'    )
#        liqWP_Tdata  = np.multiply( mdp.read_Data( filenameTS[i], 'lwp_bar' ), 1000.0)
#        case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
#        maksimiLWPlist[ case_indeksi ] = np.max(liqWP_Tdata)
#        
#        colorMap = thickness_color
#        skal = thickness_design[ case_indeksi ] / np.max( thickness_design)
#        color = colorMap(skal)
#        
#        if 'maksimiLWP' in locals():
#            if maksimiLWP <  np.max(liqWP_Tdata):
#                maksimiLWP    =    np.max(liqWP_Tdata)
#        else:
#            maksimiLWP    =    np.max(liqWP_Tdata)
#        
#        
#        
#        
#        #if 'minimiLWP' in locals():
#            #minimiLWP = min( minimiLWP,  np.min(liqWP_Tdata) )
#        #else:
#            #minimiLWP = np.min(liqWP_Tdata)
#            
#        if LEGEND:
#            nimi = 'LWP' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi = 'LWP' + ' ' + tag + ' ' + 'LVL' + lvl
#
#        mdp.aikasarjaTulostus( liqWP_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'LWP g/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )
#
#
#    mdp.plot_setYlim( 0.0, maksimiLWP )
#    if saveFig:
#        plt.savefig( picturefolder + 'lwp'+ '_' + tag + '_' + 'LVL'+lvl + '.png')

#print ' '
#print 'q profile changes'
#if QDIFF:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#        
#        q_data      = mdp.read_Data( filenamePS[i], 'q' )   
#        time_data   = mdp.read_Data( filenameTS[i], 'time'    )
#
#        height_data = mdp.read_Data( filenameNC[i], 'zt' )
#        
#
#        aikaP = np.argmin( np.abs(ajanhetket - time_data) )
#        
#        #q_difference = q_data[ aikaP, : ] / q_data[ 0, : ] -1
#
#
#        case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
#        dens = 20
#        
#        fracZ   = np.zeros( dens*(len(height_data)-1) +1 )
#        normZ   = np.zeros( len(fracZ) )
#        
#        qSpline = np.zeros( ( 2, len(fracZ) ) )
#        
#        
#        for k in xrange(len(fracZ)-1):
#            h_indeksi = int(np.floor(k/dens) )
#            normZ[k] =  ( height_data[ h_indeksi ] + np.mod(k,dens)*(height_data[ h_indeksi + 1] - height_data[ h_indeksi ])/float(dens) )
#            fracZ[k] =  normZ[k] / pblh_design[case_indeksi]
#    
#        normZ[-1] = height_data[-1]
#        fracZ[-1] = normZ[-1] / pblh_design[case_indeksi]
#        #np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})
#        tck0 = interpolate.splrep( height_data, q_data[ 0,: ]     )
#        tckT = interpolate.splrep(height_data, q_data[ aikaP,:]  )
#        for k in xrange( np.shape(qSpline)[1] ):
#            qSpline[ 0,k ]  = interpolate.splev( normZ[k], tck0 )
#            qSpline[ 1,k ]  = interpolate.splev( normZ[k], tckT )
#        
#        q_difference = qSpline[ 1,:] / qSpline[ 0,:] -1
#            
#
#        #print 'minimi '  + str(np.min(liqWP_Tdata))
#        #print 'maksimi ' + str(np.max(liqWP_Tdata))
#        if 'maksimifracZ' in locals():
#            maksimifracZ = max( maksimifracZ,  np.max(fracZ) )
#        else:
#            maksimifracZ = np.max(fracZ)
#        
#        if 'maksimiQ' in locals():
#            maksimiQ = max( maksimiQ,  np.max(q_difference) )
#        else:
#           maksimiQ = np.max(q_difference)
#
#        if 'minimiQ' in locals():
#            minimiQ = max( minimiQ,  np.max(q_difference) )
#        else:
#            minimiQ = np.max(q_difference)
#        
#        colorMap = pblh_color #plt.cm.Reds
#        skal = pblh_design[ case_indeksi ] / np.max( pblh_design)
#        color = colorMap(skal)
#        #if 'minimiLWP' in locals():
#            #minimiLWP = min( minimiLWP,  np.min(liqWP_Tdata) )
#        #else:
#            #minimiLWP = np.min(liqWP_Tdata)
#        
#        if LEGEND:
#            nimi = 'Total water mix. rat g/m^2 difference between 0h & '+str((ajanhetket/3600.)) + 'h' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi = 'Total water mix. rat g/m^2 difference between 0h & '+str((ajanhetket/3600.)) + 'h' + ' ' + tag + ' ' + 'LVL' + lvl
#
#        mdp.profiiliTulostus( q_difference, aikaPisteet = 0, korkeus = fracZ, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = r'$\frac{q_{t=2h}}{q_{t=0h}}-1$', ynimi = 'z/pblh', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )
#
#    mdp.plot_setYlim( 0.0, 1.01 )
#    plt.xlim( -0.3, 0.025)
#    if saveFig:
#        plt.savefig( picturefolder + 'q_diff'+ '_' + tag + '_' + 'LVL'+lvl + '.png')

#############################################
##print ' '
##print 't profile changes'
#if TDIFF:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#        
#        t_data      = mdp.read_Data( filenamePS[i], 't' )   
#        time_data   = mdp.read_Data( filenameTS[i], 'time'    )
#
#        height_data = mdp.read_Data( filenameNC[i], 'zt' )
#        
#        aikaP = np.argmin( np.abs(ajanhetket - time_data) )
#        
#        #t_difference = t_data[ aikaP, : ] / t_data[ 0, : ] -1
#
#
#        case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
#        dens = 20
#        
#        fracZ   = np.zeros( dens*(len(height_data)-1) +1 )
#        normZ   = np.zeros( len(fracZ) )
#        
#        tSpline = np.zeros( ( 2, len(fracZ) ) )
#        
#        
#        for k in xrange(len(fracZ)-1):
#            h_indeksi = int(np.floor(k/dens) )
#            normZ[k] =  ( height_data[ h_indeksi ] + np.mod(k,dens)*(height_data[ h_indeksi + 1] - height_data[ h_indeksi ])/float(dens) )
#            fracZ[k] =  normZ[k] / pblh_design[case_indeksi]
#    
#        normZ[-1] = height_data[-1]
#        fracZ[-1] = normZ[-1] / pblh_design[case_indeksi]
#        #np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})
#        tck0 = interpolate.splrep( height_data, t_data[ 0,: ]     )
#        tckT = interpolate.splrep(height_data, t_data[ aikaP,:]  )
#        for k in xrange( np.shape(tSpline)[1] ):
#            tSpline[ 0,k ]  = interpolate.splev( normZ[k], tck0 )
#            tSpline[ 1,k ]  = interpolate.splev( normZ[k], tckT )
#        
#        t_difference = tSpline[ 1,:] / tSpline[ 0,:] -1
#            
#
#        #print 'minimi '  + str(np.min(liqWP_Tdata))
#        #print 'maksimi ' + str(np.max(liqWP_Tdata))
#        if 'maksimifracZ' in locals():
#            maksimifracZ = max( maksimifracZ,  np.max(fracZ) )
#        else:
#            maksimifracZ = np.max(fracZ)
#        
#        if 'maksimiT' in locals():
#            maksimiT = max( maksimiT,  np.max(t_difference) )
#        else:
#            maksimiT = np.max(t_difference)
#
#        if 'minimiT' in locals():
#            minimiT = max( minimiT,  np.max(t_difference) )
#        else:
#            minimiT = np.max(t_difference)
#        
#        colorMap = pblh_color # plt.cm.Reds
#        skal = pblh_design[ case_indeksi ] / np.max( pblh_design)
#        color = colorMap(skal)
#        #if 'minimiLWP' in locals():
#            #minimiLWP = min( minimiLWP,  np.min(liqWP_Tdata) )
#        #else:
#            #minimiLWP = np.min(liqWP_Tdata)
#        
#        if LEGEND:
#            nimi = 'Temperature [K] difference between 0h & '+str((ajanhetket/3600.)) + 'h' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi = 'Temperature [K] difference between 0h & '+str((ajanhetket/3600.)) + 'h' + ' ' + tag + ' ' + 'LVL' + lvl
#
#        mdp.profiiliTulostus( t_difference, aikaPisteet = 0, korkeus = fracZ, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = r'$\frac{t_{t=2h}}{t_{t=0h}}-1$', ynimi = 'z/pblh', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )
#
#    mdp.plot_setYlim( 0.0, 1.01 )
#    plt.xlim( -0.05, 0.025)
#    if saveFig:
#        plt.savefig( picturefolder + 't_diff'+ '_' + tag + '_' + 'LVL'+lvl + '.png')
#
###################################
#if CDNC:
#    for i in xrange(len(arguments)-1):
#        uusikuva = True if i == 0 else  False
#        time_data    = mdp.read_Data( filenameTS[i], 'time'    )
#        Nc_ic_Tdata  = mdp.read_Data( filenameTS[i], 'Nc_ic' )
#        case_indeksi = int(filenameTS[i].split("/")[-2][-2:])-1
#        Nc_ic_Tdata = Nc_ic_Tdata / (np.ones(len(Nc_ic_Tdata))*num_pbl_design[case_indeksi]*1e6)
#        
#        colorMap = num_pbl_color # plt.cm.Oranges
#        skal = num_pbl_design[ case_indeksi ] / np.max( num_pbl_design )
#        
#        color = colorMap(skal)
#
#        #maksimicdnc[i] = np.max(Nc_ic_Tdata)
#        if 'maksimiCDNC' in locals():
#            maksimiCDNC = max( maksimiCDNC,   np.max(Nc_ic_Tdata) ) 
#        else:
#            maksimiCDNC = np.max(Nc_ic_Tdata)
#
#        if LEGEND:
#            nimi = 'Relative change of In-cloud CDNC #/kg' + ' ' + tag + ' ' + 'LVL' + lvl +  ' ' + filenameNC[i].split("/")[-2]
#        else:
#            nimi = 'Relative change of In-cloud CDNC #/kg' + ' ' + tag + ' ' + 'LVL' + lvl
#
#        mdp.aikasarjaTulostus( Nc_ic_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Relative change of In-cloud CDNC #/kg', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color  )
#
#
#    mdp.plot_setYlim( 0.0, maksimiCDNC, extendBelowZero = True)
#
#    if saveFig:
#        plt.savefig( picturefolder + 'cdnc'+ '_' + tag + '_' + 'LVL'+lvl + '.png')
#
