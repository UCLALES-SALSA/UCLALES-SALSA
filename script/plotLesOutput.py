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
import matplotlib.colors as colors
import os
import glob
from netCDF4 import Dataset
from PythonMethods import myRound
global piirra, tulostus, tightXAxis, saveFig, LEGEND, tag, kylla


if __name__ == "__main__":
    
    kylla = [ 'y', 'Y', 'yes', 'Yes', 'YES', 'True', 'true', '1' ]
    emuloo = ['e', 'E', 'emul', 'y', 'Y', 'yes']
    lvl = raw_input( 'Level: ')

    
    EMUL             = False
    ICE              = False
    
    if lvl == '0':
        print "Level comparisons, WARNING make sure that the files are given in order of ice0, ice1, ice4"
    else:
        EMUL = raw_input('Emulator/Regular pictures: e/r (y/n): ') in emuloo
        
    tag  = raw_input('give tag: ')
    LVLprint     = raw_input("Do you want to print LEVEL (yes/no): ") in kylla
    
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
        saveTag = tag.replace(" ","_")
        if not LVLprint:
            saveTag=saveTag[:-1]
    else:
        saveTag = ''

    ICE = (not EMUL )
    
    customLabels = raw_input("Do you want to give custom labels for plots (yes/no): ") in kylla

    
    if EMUL:
        EMULCASES = raw_input("Do you want to give custom emul cases (yes/no): ") in kylla
    
        
    #ICE  = raw_input('Ice plots: yes/no: ') in kylla


    

#######################
##### setting up    ###    
##### DO NOT CHANGE ###
#######################
    ibrix = os.environ["IBRIXMOUNT"]
    global arguments
    global labelArray
    arguments = sys.argv
    
    labelArray = []
    emulCaseIndexes  = []

    filenameNC = []
    filenamePS = []
    filenameTS = []
    
    

    
    if LVLprint:
        LVLprintFig  = ' ' + 'LVL' + lvl
        LVLprintSave = '_' + 'LVL' + lvl
    else:
        LVLprintFig  = ''
        LVLprintSave = ''
        
    kkk = 1
    uuu = 1
    applyForAll = False
    for filebase in arguments[1:]:

        filenameNC.append( filebase + '.nc'    )
        filenamePS.append( filebase + '.ps.nc' )
        filenameTS.append( filebase + '.ts.nc' )
        
        
        if customLabels:
            teksti = "Give "+ str(kkk) +". label: (" + str(filebase)+") : "
            labelArray.append( raw_input(teksti) )
            kkk += 1
        elif ( (not customLabels) and (not EMUL) ):
            labelArray.append( filebase )
    
        if EMUL and EMULCASES :
                emulteksti = "Give "+ str(uuu) +". emul case (" + str(filebase)+") : "
                if not applyForAll:
                    emulcase = int(raw_input(emulteksti))-1
                    
                    applyteksti = "Do you want to apply this case " + str(emulcase+1) + " for all the rest? (yes/no): "
                    applyForAll = raw_input(applyteksti) in kylla
                
                emulCaseIndexes.append( emulcase )
                uuu += 1            
            


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

naytaPlotit = False

piirra = True

tulostus = False

tightXAxis = True

saveFig=True

global jarjestys


if __name__ == "__main__":
    if saveFig:
        
        kuvavaliliite = raw_input("Give custom picturefolder tag: ")
        if len(kuvavaliliite) > 0:
            kuvavaliliite = '_' + kuvavaliliite
            
        picturefolder='./pictures' + kuvavaliliite + '/'
        if not os.path.exists( picturefolder ):
            os.makedirs( picturefolder )

    cases    = len(arguments)-1
    colorNRO = cases
    
    jarjestys = np.arange( 1, cases + 1 )

    if ( colorNRO > 6 ) :
        LEGEND = False
    else:
        LEGEND = True

    colorChoice = mdp.initializeColors(colorNRO)



    ajanhetket = 2*3600. #7200.

    spinup = mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'Tspinup' )
    tmax   = mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'timmax'  )
    
    slaissaaXY   = max( mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'nxp'  ), mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'nyp'  ))


    ticksHours  = np.arange(0., tmax/3600. + 0.1, 0.5)
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

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

##########################
def piirra_aikasarjasettii( muuttuja, variKartta = plt.cm.gist_rainbow, variRefVektori = jarjestys, colorBar = None, colorBarTickValues = [0,1], colorBarTickNames = ['0','1'], muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'time [h]', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, relative = False, savePrefix = None, omaVari = True, tit = ' ', nollaArvo = None, xlabels = None, ylabels = None, xticks = None, yticks = None, spinup = None, ylabelFont = False, askChangeOfVariable = False, piilotaOsaXlabel = False, legenda = True ):
    origmuuttuja = muuttuja
    maksimi = None
    minimi  = None
    maxRef = np.max( variRefVektori )
    maksInd = None
    minInd  = None
    #vastausKaikkiin = False
    nono = None
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        muuttuja = origmuuttuja
        if EMUL:
            if not EMULCASES:
                case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
            else:
                case_indeksi = emulCaseIndexes[i]
        else:
            case_indeksi = i
            
        if askChangeOfVariable:
            printti = "Case is "+ arguments[i+1] + ", the variable is currently: " + muuttuja + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. prcp/rmH2Opr Nc_ic/CCN): '
            apumuuttuja = raw_input(printti)
            
            if apumuuttuja != '':
                muuttuja = apumuuttuja
            
            if muuttuja == 'prcp':
                muunnosKerroin = 1.
            elif muuttuja == 'rmH2Opr':
                muunnosKerroin = 2.5e+06
            
            print 'variable is', muuttuja, 'muunnosKerroin', muunnosKerroin
            print ' '
            

            
        time_data       = mdp.read_Data( filenameTS[i], 'time'  )
        muuttuja_Tdata  = mdp.read_Data( filenameTS[i], muuttuja )*muunnosKerroin
        

        #### relative    
        if relative:
            
           
            if ((nollaArvo is None) and (nono is None)):
                lahtoarvo = muuttuja_Tdata[0]
                if lahtoarvo == -999.:
                    print 'WARNING lahtoarvo = -999. lasketaan seuraavasta arvosta', muuttuja
                    lahtoarvo = muuttuja_Tdata[1]
                    
            elif ( isinstance( nollaArvo, np.ndarray) or isinstance(nollaArvo, list) ):
                lahtoarvo = nollaArvo[case_indeksi]
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
                
                

        
        
        colorMap = variKartta
        skal = variRefVektori[ case_indeksi ] / maxRef
        color = colorMap(skal)
        
        if omaVari:
            omavari = color
        else:
            omavari = omaVari

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
            if customLabels:
                label = labelArray[i]
            else:    
                label = str(case_indeksi+1)
        else:
            nimi = longName + ' ' + tag
            label = labelArray[i]
        
        fig, ax = mdp.aikasarjaTulostus( muuttuja_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=tightXAxis, LEGEND=legenda, omavari = omavari, label = label )
        #######################
    #print 'muuttuja', muuttuja, 'indeksi min', minInd, 'arvo min', minimi, 'indeksi max', maksInd, 'arvo max', maksimi
    
    xticksHours = np.arange(0., max(time_data)/3600. + 0.1, 0.5)

    xticks = xticksHours * 3600.
    
    if xlabels is None:
        xlabels = map( str, xticksHours)
    
    if yticks is not None:
        ax.set_yticks( yticks)
    
    if ylabels is not None:
        ax.set_yticklabels( ylabels )
    

    ax.set_xticks( xticks )
    ax.set_xticklabels( xlabels )
    
    if piilotaOsaXlabel:
        k = 0
        for label in ax.xaxis.get_ticklabels():
            if np.mod(k,4) != 0:
                label.set_visible(False)
            k+=1
    
    if ((ylabelFont is not None) and (yticks is not None)):
        for takka in ( ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() ):
            takka.label.set_fontsize(ylabelFont) 
    
    #if customYscale:
        #ax.set_yscale('symlog', basey=10, linthresy=[0,1])    
    
    if spinup is not None:
        mdp.plot_vertical( spinup )

    if colorBar is not None and ( omaVari == True):
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
        plt.savefig( picturefolder + savePrefix + '_' + saveTag + LVLprintSave + '.png')

def piirra_aikasarjavertailusettii( muuttuja, variKartta = plt.cm.gist_rainbow, variRefVektori = jarjestys, colorBar = None, colorBarTickValues = [0,1], colorBarTickNames = ['0','1'], muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'time [h]', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, relative = False, savePrefix = None, omaVari = True, tit = ' ', nollaArvo = None, xlabels = None, ylabels = None, xticks = None, yticks = None, spinup = None, ylabelFont = False, askChangeOfVariable = False, piilotaOsaXlabel = False ):
    
    origmuuttuja = muuttuja
    maksimi = None
    minimi  = None
    maxRef = np.max( variRefVektori )
    maksInd = None
    minInd  = None
    #vastausKaikkiin = False
    nono = None
    
    time_data       = mdp.read_Data( filenameTS[0], 'time'  )
    
    for i in xrange(1,len(arguments)-1):
    
        uusikuva = True 
        case_indeksi = i
            
       
        if i == 1:
            longName = "ice1-ice0"
            stri = "ice1"
        elif i== 2:
            longName = "ice4-ice0"
            stri = "ice4"
        
        if muuttuja == "lwp":
            
            muuttuja_Tdata_ice0 = mdp.read_Data( filenameTS[0], 'lwp_bar' )*muunnosKerroin
            muuttuja_Tdata      = mdp.read_Data( filenameTS[i], 'lwp_bar' )*muunnosKerroin
            ylabel = r'$\Delta$LWP [$g/m^2$]'
            
            if i == 1:
                omavari = 'red'
            elif i == 2:
                omavari = 'darkviolet'
                
            
            
        elif muuttuja == "lwpiwp":
            
            muuttuja_Tdata_ice0 = ( mdp.read_Data( filenameTS[0], 'lwp_bar' )  ) *muunnosKerroin
            muuttuja_Tdata      = ( mdp.read_Data( filenameTS[i], 'lwp_bar' ) +  mdp.read_Data( filenameTS[i], 'iwp_bar' ) ) *muunnosKerroin
            ylabel = r'$\Delta$(LWP+IWP) [$g/m^2$]'
            
            if i == 1:
                omavari = 'blue'
            elif i == 2:
                omavari = 'darkblue'
            
        else:
            sys.exit("anna oikea muuttuja")
        
        muuttuja_Tdata = muuttuja_Tdata - muuttuja_Tdata_ice0



        maksimi = np.max( muuttuja_Tdata )
        
        minimi  = np.min( muuttuja_Tdata )
        




        nimi = longName + ' ' + tag

        
        fig, ax = mdp.aikasarjaTulostus( muuttuja_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=tightXAxis, LEGEND=False, omavari = omavari)
        #######################
    #print 'muuttuja', muuttuja, 'indeksi min', minInd, 'arvo min', minimi, 'indeksi max', maksInd, 'arvo max', maksimi
    
        xticksHours = np.arange(0., max(time_data)/3600. + 0.1, 0.5)

        xticks = xticksHours * 3600.
        
        if xlabels is None:
            xlabels = map( str, xticksHours)
        
        if yticks is not None:
            ax.set_yticks( yticks)
        
        if ylabels is not None:
            ax.set_yticklabels( ylabels )
        

        ax.set_xticks( xticks )
        ax.set_xticklabels( xlabels )
        
        if piilotaOsaXlabel:
            k = 0
            for label in ax.xaxis.get_ticklabels():
                if np.mod(k,4) != 0:
                    label.set_visible(False)
                k+=1
        
        if ((ylabelFont is not None) and (yticks is not None)):
            for takka in ( ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() ):
                takka.label.set_fontsize(ylabelFont) 
        
        #if customYscale:
            #ax.set_yscale('symlog', basey=10, linthresy=[0,1])    
        
        if spinup is not None:
            mdp.plot_vertical( spinup )
        
        mdp.plot_horizontal( 0. )

        if colorBar is not None and ( omaVari == True):
            cb = plt.colorbar( colorBar, shrink=.9, pad=.03, aspect=10, ticks = colorBarTickValues )#colorBar 
            cb.ax.set_yticklabels(colorBarTickNames)
            cb.ax.set_ylabel( tit )


        
        if ( asetaRajat ):
            mdp.plot_setYlim( minimi, maksimi, extendBelowZero = extendBelowZero)
        

            
            
        if savePrefix is None:  
            savePrefix = muuttuja+longName+stri
        else:
            savePrefix = savePrefix + stri

        if saveFig:
            plt.savefig( picturefolder + savePrefix + '_' + saveTag + LVLprintSave + '.png')




#########################
def piirra_profiilisettii( muuttuja, variKartta = plt.cm.gist_rainbow, variRefVektori = jarjestys, colorBar = None, colorBarTickValues = [0,1], colorBarTickNames = ['0','1'], muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'xlabel', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, xmin = None, xmax = None, relative = False, savePrefix = None, ajanhetket = None, tit = ' ', nollaArvo = None, rajaKerros = None, omaVari = True ):

    maksimi = None
    minimi  = None
    maxRef = np.max( variRefVektori )
    maksimifracZ = None
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        
        muuttuja_data      = mdp.read_Data( filenamePS[i], muuttuja )*muunnosKerroin   
        
        if EMUL:
            if not EMULCASES:
                case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
            else:
                case_indeksi = emulCaseIndexes[i]
        else:
            case_indeksi = i  
        
        #### relative    
        if relative:
            if ((nollaArvo is None) and (nono is None)):
                lahtoarvo = muuttuja_data[0]
                apu = lahtoarvo[0,:] == -999.
                if apu.all():
                    print 'WARNING lahtoarvo = -999. lasketaan seuraavasta arvosta', muuttuja
                    lahtoarvo = muuttuja_data[1]
            elif ( isinstance( nollaArvo, np.ndarray) or isinstance(nollaArvo, list) ):
                lahtoarvo = nollaArvo[case_indeksi]
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

        time_data   = mdp.read_Data( filenameTS[i], 'time'    )

        height_data = mdp.read_Data( filenameNC[i], 'zt' )
        
        if ajanhetket is not None:
            aikaP = np.argmin( np.abs(ajanhetket - time_data) )
        
        if aikaP > np.shape(muuttuja_data)[0] -1:
            #print 'VAROITUS: indeksi liian iso, skipataan - ', 'aikaP indeksi:', aikaP, 'muuttuja_data viimeinen indeksi:',  np.shape(muuttuja_data)[0] -1, 'muuttuja', muuttuja
            if aikaP - (np.shape(muuttuja_data)[0] -1) > 1:
                print 'VAROITUS: indeksi liian iso, skipataanko, koska ero viimeiseen indeksiin on >1 : ', aikaP - (np.shape(muuttuja_data)[0] -1) , 'aikaP indeksi:', aikaP, 'muuttuja_data viimeinen indeksi:',  np.shape(muuttuja_data)[0] -1, 'muuttuja', muuttuja
                ippu = 'Do you want use the last index ' + str(aikaP) + ' (yes/no): '
                lastIND = raw_input( ippu ) in kylla
                if not lastIND:
                    continue
            # asetetaan aikaP viimeiseksi indeksiksi
            aikaP = np.shape(muuttuja_data)[0] -1
            
            


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
        
        if omaVari:
            omavari = color
        else:
            omavari = omaVari
        
        if EMUL:
            nimi = longName + ' ' + tag + LVLprintFig
            if customLabels:
                label = labelArray[i]
            else:    
                label = str(case_indeksi+1)
        else:
            nimi = longName + ' ' + tag
            label = labelArray[i]
        fig, ax = mdp.profiiliTulostus( p_difference, aikaPisteet = 0, korkeus = fracZ, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = omavari, label = label, loc = 3 )

        ####################
    # tikkien fonttikoko
    for takka in ( ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() ):
                takka.label.set_fontsize(18) 


    font = 26
    ax.xaxis.get_label().set_fontsize(font)
    ax.yaxis.get_label().set_fontsize(font)
    

                
    if colorBar is not None and ( omaVari == True):
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
        plt.savefig( picturefolder + savePrefix + '_' + saveTag + LVLprintSave + '.png')


def piirra_profiiliKehitys(  muuttuja, variKartta = plt.cm.gist_rainbow, colorBar = None, colorBarTickValues = [0,1], colorBarTickNames = ['0','1'], muunnosKerroin = 1.0, longName = 'titteli', xlabel = 'xlabel', ylabel='ylabel', extendBelowZero = True, asetaRajat = True, ymin = None, ymax = None, xmin = None, xmax = None,  savePrefix = None, aikaPisteet = None, tit = ' ', nollaArvo = None, rajaKerros = None, paksuus = None, xlabels = None, ylabels = None, xticks = None, yticks = None, tempConversion = False, askChangeOfVariable = False):
    #import matplotlib.patches as mpatches
    from PythonMethods import myRoundFloat
    from emulator_inputs import absT
    if aikaPisteet is None:
        sys.exit("You did not give any time points for profile evolution, exiting")
    
    time_data   = mdp.read_Data( filenameTS[0], 'time'    )
    maxRef = np.max( time_data )
    variRefVektori = time_data
    
    lastIND = None
    for i in xrange(len(arguments)-1):
        maksimi = None
        minimi  = None
        maksimifracZ = None
        if EMUL:
            if not EMULCASES:
                case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
            else:
                case_indeksi = emulCaseIndexes[i]
        else:
            case_indeksi = i

        if askChangeOfVariable:
            printti = "Case is "+ arguments[i+1] + ", the variable is currently: " + muuttuja + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. [prcp/rmH2Opr] [Nc_ic/CCN] [P_rl/l]): '
            apumuuttuja = raw_input(printti)
        
            if apumuuttuja != '':
                muuttuja = apumuuttuja
        
            print 'variable is', muuttuja
            print ' '
    
    

        for t in xrange(len(aikaPisteet)):
            uusikuva = True if t == 0 else  False

            muuttuja_data      = mdp.read_Data( filenamePS[i], muuttuja )*muunnosKerroin
            
            if tempConversion:
                paine_data    =  mdp.read_Data( filenamePS[i], 'p' )
                for aa in xrange( np.shape(paine_data)[0] ):
                    for bb in xrange( np.shape(paine_data)[1] ):
                        muuttuja_data[aa, bb] = absT( muuttuja_data[aa,bb], paine_data[aa,bb], 1.)
            
            height_data = mdp.read_Data( filenameNC[i], 'zt' )
            
            base = mdp.read_Data( filenameTS[i], 'zb' )
            top  = mdp.read_Data( filenameTS[i], 'zc' )
            
            
            aikaP = np.argmin( np.abs(aikaPisteet[t] - time_data) )
            
                
            if aikaP > np.shape(muuttuja_data)[0] -1:
                #print 'VAROITUS: indeksi liian iso, skipataan - ', 'aikaP indeksi:', aikaP, 'muuttuja_data viimeinen indeksi:',  np.shape(muuttuja_data)[0] -1, 'muuttuja', muuttuja
                if aikaP - (np.shape(muuttuja_data)[0] -1) > 1:
                    if i == 0:
                        print ' '
                    print 'VAROITUS: indeksi liian iso, skipataanko, koska ero viimeiseen indeksiin on >1 : ', aikaP - (np.shape(muuttuja_data)[0] -1) , 'aikaP indeksi:', aikaP, 'muuttuja_data viimeinen indeksi:',  np.shape(muuttuja_data)[0] -1, 'muuttuja', muuttuja
                    ippu = 'Do you want use the last index ' + str(aikaP) + ' (yes/no): '
                    if lastIND is None:
                        lastIND = raw_input( ippu ) in kylla
                        if not lastIND:
                            continue
                # asetetaan aikaP viimeiseksi indeksiksi
                aikaP = np.shape(muuttuja_data)[0] -1
            
            
            
            ##### pbl height
            #if rajaKerros is None:
                #pbl_height = float(raw_input('Anna rajakerroksen korkeus lahtotilanteessa (float): ') )
            #elif ( isinstance(rajaKerros, np.ndarray) or isinstance(rajaKerros, list) ):
                #pbl_height = rajaKerros[case_indeksi]
            ##### end pbl height

            #dens = 20
            
            #fracZ   = np.zeros( dens*(len(height_data)-1) +1 )
            #normZ   = np.zeros( len(fracZ) )
            
            #pSpline = np.zeros( ( 2, len(fracZ) ) )
            
            
            #for k in xrange(len(fracZ)-1):
                #h_indeksi = int(np.floor(k/dens) )
                #normZ[k] =  ( height_data[ h_indeksi ] + np.mod(k,dens)*(height_data[ h_indeksi + 1] - height_data[ h_indeksi ])/float(dens) )
                #fracZ[k] =  normZ[k] / pbl_height
        
            #normZ[-1] = height_data[-1]
            #fracZ[-1] = normZ[-1] / pbl_height
            
            #tck0 = interpolate.splrep( height_data, muuttuja_data[ 0,: ]     )
            #tckT = interpolate.splrep(height_data, muuttuja_data[ aikaP,:]  )
            #for k in xrange( np.shape(pSpline)[1] ):
                #pSpline[ 0,k ]  = interpolate.splev( normZ[k], tck0 )
                #pSpline[ 1,k ]  = interpolate.splev( normZ[k], tckT )
            
            #p_difference = pSpline[ 1,:] / pSpline[ 0,:] -1
                

            #if maksimifracZ is not None:
                #maksimifracZ = max( maksimifracZ,  np.max(fracZ) )
            #else:
                #maksimifracZ = np.max(fracZ)

            if maksimi is not None:
                maksimi = max( maksimi,  np.max(muuttuja_data) )
            else:
                maksimi = np.max(muuttuja_data)

            if minimi is not None:
                minimi = min( minimi,  np.min(muuttuja_data) )
            else:
                minimi = np.min(muuttuja_data)
            
            colorMap = variKartta
            skal = variRefVektori[ aikaP ] / maxRef
            color = colorMap(skal)

            
            
            if EMUL and not customLabels:
                aputagi = str(case_indeksi+1)
                
            else:
                aputagi = labelArray[i]
            
            legend = str(myRoundFloat( time_data[aikaP]/3600.) ) + ' [h]'
            nimi = longName + ' ' + aputagi + ' ' + tag
                
                

            fig, ax = mdp.profiiliTulostus( muuttuja_data[aikaP,:], aikaPisteet = 0, korkeus = height_data, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = xlabel, ynimi = ylabel, tightXAxis=False, tightYAxis = True, LEGEND=True, omavari = color, label = legend, loc = 2 )
            
            if aikaP >0:
                plt.axhline( base[aikaP], color = color, linestyle= 'dashed' )
                plt.axhline(  top[aikaP], color = color, linestyle= 'dashed' )
            elif aikaP == 0 and EMUL:
                plt.axhline( rajaKerros[case_indeksi], color = color, linestyle= 'dashed' )
                plt.axhline( rajaKerros[case_indeksi] - paksuus[case_indeksi], color = color, linestyle= 'dashed' )
            
        ####################
        ## tikkien fonttikoko
        for takka in ( ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() ):
            takka.label.set_fontsize(18) 

        
        font = 26
        ax.xaxis.get_label().set_fontsize(font)
        ax.yaxis.get_label().set_fontsize(font)
        #### end tikkien fonttikoko
        
        korkeustikit = map( int, np.arange( 0, max(height_data), 100.) )
        #if yticks is not None:
        ax.set_yticks( korkeustikit )
        #if ylabes
        ax.set_yticklabels( map(str, korkeustikit))
        
        k = 0
        for label in ax.yaxis.get_ticklabels():
            if np.mod(k,2) != 0:
                label.set_visible(False)
            k+=1
        
        #patch = mpatches.Patch(color=color, label=label)
        #plt.legend(handles=[patch])
                    
        #if colorBar is not None:
            #cb = plt.colorbar( colorBar, shrink=.9, pad=.03, aspect=10, ticks = colorBarTickValues )#colorBar 
            #cb.ax.set_yticklabels(colorBarTickNames)
            #cb.ax.set_ylabel( tit )


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
            valitagi = '_'
            if customLabels:
                valitagi = labelArray[i]
            elif applyForAll and not customLabels:
                valitagi = str(case_indeksi+1) + str(i)
            else:
                valitagi = str(case_indeksi+1)
            
            subfolder = picturefolder + '/' + savePrefix + '/'
            if not os.path.exists( subfolder ):
                os.makedirs( subfolder )
            
            plt.savefig( subfolder + savePrefix + '_' + valitagi + '_' + saveTag + LVLprintSave + '.png')
        


############################
def piirra_aikasarjaPathXYZ( muuttuja, muunnosKerroin = 1.0, longName = None, savePrefix = None, xaxislabel = 'time [h]', xlabels = None, ylabels = None, xticks = None, yticks = None, spinup = None, piilotaOsaXlabel = False ):
    
        
    if longName is None:
        longName = muuttuja
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False

        if EMUL:
            if not EMULCASES:
                case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
            else:
                case_indeksi = emulCaseIndexes[i]
        else:
            case_indeksi = i
            
        muuttuja_data = mdp.read_Data( filenameNC[i], muuttuja )*muunnosKerroin
        time_data     = mdp.read_Data( filenameNC[i], 'time'   )
        dn0_data      = mdp.read_Data( filenameNC[i], 'dn0'    )
        zt_data       = mdp.read_Data( filenameNC[i], 'zt'     )
        zm_data       = mdp.read_Data( filenameNC[i], 'zm'     )
        korkeus = ( zm_data - zt_data )*2.0
        
        
        if customLabels:
            aputagi = labelArray[i]
        else: 
            aputagi = str(case_indeksi+1)
        
        nimi = longName + ' grid data'
        label = aputagi

        fig,ax = mdp.laske_path_aikasarjaXYZ( muuttuja_data, dn0_data, korkeus, time_data, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xlabel = xaxislabel, tightXAxis=tightXAxis, label = label)
    
    
    ajat  =  np.round( mdp.read_Data( filenameNC[0], 'time' ), 0 )
    
    xticksHours = np.arange(0, max(ajat)/3600. + 0.1, 0.5)
    
    if xticks is None:
        xticks = xticksHours * 3600.
    
    if xlabels is None:
        xlabels = map( str, xticksHours)
    
    print 'xticks', xticks
    print 'xlabels', xlabels
    ax.set_xticks( xticks )
    ax.set_xticklabels( xlabels )
    
    if piilotaOsaXlabel:
        k = 0
        for label in ax.xaxis.get_ticklabels():
            if np.mod(k,4) != 0:
                label.set_visible(False)
            k+=1
    
    if spinup is not None:
        mdp.plot_vertical( spinup )
    
    if savePrefix is None:
        savePrefix = muuttuja
    if saveFig:
        plt.savefig( picturefolder + savePrefix + '_' + saveTag + LVLprintSave + '.png')


##########################
def piirra_maksimiKeissit( muuttuja, muunnosKerroin = 1.0, longName = 'pitka nimi', xlabel = 'case', ylabel = 'ylabel', savePrefix = None, askChangeOfVariable = False, ajanhetket = None ):
    lista = np.zeros( cases )
    xTikit = np.zeros( cases ) 

    
    for i in xrange(len(arguments)-1):
        if askChangeOfVariable:
            printti = "Case is "+ arguments[i+1] + ", the variable is currently: " + muuttuja + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. [prcp/rmH2Opr] [Nc_ic/CCN] [P_rl/l]): '
            apumuuttuja = raw_input(printti)
        
            if apumuuttuja != '':
                muuttuja = apumuuttuja
            
            if muuttuja == 'prcp':
                muunnosKerroin = 1.
            elif muuttuja == 'rmH2Opr':
                muunnosKerroin = 2.5e+06
            
            print 'variable is', muuttuja, 'muunnosKerroin', muunnosKerroin
            print ' '
                
        muuttuja_Tdata  = mdp.read_Data( filenameTS[i], muuttuja )*muunnosKerroin
        time_data       = mdp.read_Data( filenameNC[i], 'time'   )
        
        if EMUL:
            if not EMULCASES:
                tikki = int(filenameNC[i].split("/")[-2][-2:])
            else:
                tikki = emulCaseIndexes[i]+1    
        else:
            
            tikki = int(i+1)

        xTikit[i] = tikki
        
        if ajanhetket is not None:
            aikaP = np.argmin( np.abs( ajanhetket - time_data) )
        else:
            aikaP = 0
        
        if aikaP > np.shape(muuttuja_Tdata)[0] -1:
            #print 'VAROITUS: indeksi liian iso, skipataan - ', 'aikaP indeksi:', aikaP, 'muuttuja_Tdata viimeinen indeksi:',  np.shape(muuttuja_Tdata)[0] -1, 'muuttuja', muuttuja
            if aikaP - (np.shape(muuttuja_Tdata)[0] -1) > 1:
                print 'VAROITUS: indeksi liian iso, skipataanko, koska ero viimeiseen indeksiin on >1 : ', aikaP - (np.shape(muuttuja_Tdata)[0] -1) , 'aikaP indeksi:', aikaP, 'muuttuja_Tdata viimeinen indeksi:',  np.shape(muuttuja_Tdata)[0] -1, 'muuttuja', muuttuja
                ippu = 'Do you want use the last index ' + str(aikaP) + ' (yes/no): '
                lastIND = raw_input( ippu ) in kylla
                if not lastIND:
                    continue
            # asetetaan aikaP viimeiseksi indeksiksi
            aikaP = np.shape(muuttuja_Tdata)[0] -1
            
            
            
        #print 'np.shape(lista)', np.shape(lista), 'case_indeksi',i, 'aikaP', aikaP, 'np.shape(muuttuja_Tdata)[0]', np.shape(muuttuja_Tdata)[0], 'muuttuja', muuttuja
        
        lista[i] = np.max( muuttuja_Tdata[aikaP:] )

    fig, ax = mdp.plottaa( jarjestys, lista, tit = longName+ ' ' + tag + LVLprintFig, xl = xlabel, yl = ylabel, changeColor = True, markers = True, uusikuva = True, gridi = False, tightXAxis = False, tightYAxis = True, LEGEND = False )
    
    mdp.plot_setYlim( 0., max(lista), extendBelowZero = False, A = 0.05 )
    #print map(str, xTikit)
    ax.set_xticks( map(int, jarjestys) )
    for takka in ax.xaxis.get_major_ticks():
                takka.label.set_fontsize(9)
    #ax.set_xticklabels( map(str, xTikit) , fontsize = 8 )

    if savePrefix is None:
                savePrefix = longName.replace( " ", "_" )
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + saveTag + LVLprintSave + '.png')

#######################
def piirra_binijakauma( muuttuja, muunnosKerroin = 1.0, longName = None, savePrefix = None ):

    if longName is None:
        longName = muuttuja

    #for bini in xrange( np.shape(S_Rwiba_data)[1] ):
        #mdp.laske_MeanDiameterInBin( S_Rwiba_data, bini, S_Nc_data, time_data, tulostus = True, piirra =piirra, nimi = 'Mean diameter of ice particles in liquid clouds in bin ' )
        

#######################
def piirra_domainProfiili( muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = None, colorBarTickNames = None, xlabels = None, ylabels = None, xticks = None, yticks = None, variKartta = plt.cm.Blues, profiili = False, radCool = False, cloudBaseTop = False, spinup = None, testi = False, sliceXbeg = None, sliceXend = None, sliceYbeg = None, sliceYend = None, askChangeOfVariable = False ):
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib import colors
    
    if spinup is None:
        sys.exit("EXIT anna spinup @ piirra_domainProfiili")
    if profiili:
        tiedostonimi = filenamePS
    else:
        tiedostonimi = filenameNC
    
    


    norm = None
    
    maxi = None
    mini = None
    
    slices = ((sliceXbeg is not None) or (sliceXend is not None) or (sliceYbeg is not None) or (sliceYend is not None)) # slices
    
    if (sliceXbeg is not None) and (sliceXend is None):
        sliceXend = sliceXbeg
    if (sliceXbeg is None) and (sliceXend is not None):
        sliceXbeg = sliceXend    
    if (sliceXbeg is not None) and (sliceXend is None):
        sliceXend = sliceXbeg
    if (sliceXbeg is None) and (sliceXend is not None):
        sliceXbeg = sliceXend

    for i in xrange(len(arguments)-1):
        if askChangeOfVariable:
            printti = "Case is "+ arguments[i+1] + ", the variable is currently: " + muuttuja + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. prcp/rmH2Opr Nc_ic/CCN): '
            apumuuttuja = raw_input(printti)
            
            if apumuuttuja != '':
                muuttuja = apumuuttuja
            print 'variable is', muuttuja
            print ' '
            
        if EMUL:
            if not EMULCASES:
                case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
            else:
                case_indeksi = emulCaseIndexes[i]
        else:
            case_indeksi = i
        # Make plot with vertical (default) colorbar
        asp = 25./45.
        
        #w, h = mpl.figure.figaspect( asp )
        #fig, ax = plt.subplots(figsize = (24,15)) #figsize=(w, h)
        fig, ax = mdp.plot_alustus()
        muuttuja_data  = np.multiply( mdp.read_Data( tiedostonimi[i], muuttuja ), muunnosKerroin)

        if useDN:
            dn0           = np.power( mdp.read_Data( filenamePS[i], 'dn0' ))
            muuttuja_data = np.multiply( muuttuja_data, dn0 )
    
        if profiili:
            muuttuja_meanProfile = muuttuja_data
        elif not profiili:
            muuttuja_meanProfile = np.mean( np.mean( muuttuja_data, axis = 1), axis = 1)
        elif not profiili and slices:
            muuttuja_meanProfile = np.mean( muuttuja_data[:,sliceYbeg:sliceYend+1,sliceXbeg:sliceXend+1,:], axis = 1)
            
        
        if testi:
            muuttuja_meanProfile = np.where(mdp.read_Data( tiedostonimi[i], 'P_rl' )*1000. > 0.001, muuttuja_meanProfile, 0.)
            muuttuja_meanProfile = np.where(mdp.read_Data( tiedostonimi[i], 'P_RHi' ) > 1.05, muuttuja_meanProfile, 0.)
            #muuttuja_meanProfile = np.where(muuttuja_meanProfile < 1.05, muuttuja_meanProfile, 0.)
        
        if transpose:
            muuttuja_meanProfile = np.rot90(muuttuja_meanProfile)  #zip(*muuttuja_meanProfile[::-1])#muuttuja_meanProfile.T
        
        if maxi is not None:
            maxi = max( maxi,  np.max(muuttuja_meanProfile) )
        else:
            maxi = np.max(muuttuja_meanProfile)
        
        if mini is not None:
            mini = min( mini,  np.min(muuttuja_meanProfile) )
        else:
            mini = np.min(muuttuja_meanProfile)
        
        
        #print 'muuttuja', muuttuja, 'mini', mini, 'maxi', maxi
        if radCool:
            for pp in xrange(np.shape(muuttuja_meanProfile)[0]):
                for rr in xrange(1, np.shape(muuttuja_meanProfile)[1]):
                    muuttuja_meanProfile[pp,rr-1] = muuttuja_meanProfile[pp,rr] - muuttuja_meanProfile[pp,rr-1]
                muuttuja_meanProfile[pp,-1] = 0.
            
            mini = np.min(muuttuja_meanProfile)
            maxi = np.max(muuttuja_meanProfile)

            midpoint = np.abs( mini ) / np.abs( maxi-mini )
            print 'midpoint', midpoint
            
            apu = (1-midpoint) / midpoint
            start = 0.
            stop  = 1.
            if apu<1.:
                stop  = min( midpoint + apu*(1-midpoint), 1.)
            else:
                start = max( 0, midpoint - (1./apu*midpoint) ) 
            print 'start', start, 'stop', stop
            
            variKartta = shiftedColorMap( variKartta, start=start, stop=stop, midpoint = midpoint,   name='shrunk') #
        
        if colorBarTickValues is None:
            baseN = 5
            if ( mini < 0 and maxi > 0 ):
                posi = np.arange( 0., min(10., maxi)+.5,  1. )
                nega = np.arange( max(-10, mini), 0., 1. )
                kokonaisluvut = map(int, np.arange( int(mini/baseN)*baseN, int(maxi/baseN)*baseN, baseN  ))
                colorBarTickValues = map(int, np.concatenate((np.asarray([int(mini)]), nega, kokonaisluvut, posi, np.asarray([int(maxi)]) )) )
            else:
                colorBarTickValues =map(int, np.arange( int(mini/baseN)*baseN, int(maxi/baseN)*baseN+baseN, baseN  ))# map(int, np.arange( int(mini), int(maxi), 1  ))

         
        #print 'colorbartickvalues', colorBarTickValues 
        if isinstance( variKartta, list):
            variKartta   = colors.ListedColormap(variKartta)
            bounds = colorBarTickValues
            norm   = colors.BoundaryNorm(bounds, variKartta.N)
        
        cax = ax.imshow( muuttuja_meanProfile, interpolation='nearest', cmap=variKartta, vmin= min(colorBarTickValues), vmax=max(colorBarTickValues), aspect=asp, norm = norm ) 
        
        
        
        
        if longName is None:
            longName = muuttuja
        
        
        ax.set_title(longName)
        ax.title.set_fontsize(35)
        
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
        
        for tt in xrange(len(xticks)):
            oikeatXtikit[tt] = np.argmin( np.abs(ajat - 3600.*xticks[tt]) )
        
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
        cbar = fig.colorbar(cax, ax=ax, shrink=.383, pad=.03, aspect=10, ticks = colorBarTickValues, norm = norm)#,fraction=0.046, pad=0.04
        
        if colorBarTickNames is None:
            colorBarTickNames = map( str, colorBarTickValues)
            
        cbar.ax.set_yticklabels(colorBarTickNames)  # vertically oriented colorbar
        cbar.ax.tick_params(labelsize=18)

        

        spinupIND = np.argmin( np.abs( ajat - spinup ) )
        mdp.plot_vertical( spinupIND )
        plt.ylabel('z [m]')
        plt.xlabel('time [h]')
        plt.tight_layout()
        
        
        if savePrefix is None:
            savePrefix = 'domainTimeseriesProfiili'    
        if saveFig:
            valitagi = '_'
            if customLabels:
                valitagi = labelArray[i]
            elif applyForAll and not customLabels:
                valitagi = str(case_indeksi+1) + str(i)
            else:
                valitagi = str(case_indeksi+1)
            plt.plot(0,0, label = valitagi)
            plt.legend(bbox_to_anchor=(0., 1.13, 1., .102), loc=3, ncol=6, fancybox = True, shadow = True , mode="expand" )
            
            if slices:
                kuvakansio = picturefolder + savePrefix + '_slices/'
                
                if not os.path.exists( picturefolder ):
                    os.makedirs( picturefolder )
            
                savePrefix = savePrefix + mdp.xstr(sliceXbeg) + mdp.xstr(sliceXend) + mdp.xstr(sliceYbeg) + mdp.xstr(sliceYend)
                
            else:
            	kuvakansio = picturefolder
            
            plt.savefig( kuvakansio + savePrefix + '_' + muuttuja + '_' + valitagi + '_' + saveTag + LVLprintSave + '.png')

            
        if cloudBaseTop:
            top  = mdp.read_Data( filenameTS[i], 'zc'   )
            base = mdp.read_Data( filenameTS[i], 'zb'   )
            aika = mdp.read_Data( filenameTS[i], 'time' )
            zt   = mdp.read_Data( filenamePS[i], 'zt'   )
            
            pilvi = np.zeros(np.shape(muuttuja_meanProfile))

            for mm in xrange(np.shape(top)[0]):
                indT = np.argmin( np.abs(  top[ mm ] - zt ))
                indB = np.argmin( np.abs( base[ mm ] - zt ))
                #print 'mm', mm,'base', base[mm], zt[indB], indB, 'top', top[mm], zt[indT], indT
                
                pilvi[-indT, mm] = 100.
                pilvi[-indB, mm] = 100.
            #asp = 25./45.
            #w, h = mpl.figure.figaspect( asp )
            #fig, ax = plt.subplots(figsize = (24,15)) #figsize=(w, h)    
            plt.imshow( pilvi, cmap='Greys',  interpolation='none', alpha= 0.1 )
            #mdp.plottaa( aika, top, uusikuva = True  )
            #mdp.plottaa( aika, base )
            
###################################################################
def piirra_kokojakauma( muuttujaR = 'S_Rwiba', muuttujaN = 'S_Niba', typename = 'Ice', muunnosKerroinR = 2.e3, korkeusH = None, aikaT = None, xlabel = 'Diameter [mm]', ylabel='#', asetaRajat = True, xmax = None, ymax = None, savePrefix = None, legenda = True, interplo = True, sekoita = True, xCustSize = None, yCustSize = None, askChangeOfVariable = False ):
    from scipy.interpolate import interp1d
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.ensemble import AdaBoostRegressor
    from random import shuffle
    from math import ceil
    
    if korkeusH is None:
        korkeusH = int(raw_input("Anna z [m]: "))
    if aikaT is None:
        aikaT = int(raw_input("Anna t [h]: "))
    
    
    
    LKM = len(arguments)-1
    #fig , ax = mdp.plot_alustus(a=24,b=15*LKM, sub_i = 1, sub_j = LKM, )
    #uusifig = True
    figuuri = None
    for i in xrange(len(arguments)-1):
        
        uusikuva = True
        
        if askChangeOfVariable:
            printti = "Case is "+ arguments[i+1] + ", the variable is currently: " + muuttujaR + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. prcp/rmH2Opr Nc_ic/CCN): '
            apumuuttujaR = raw_input(printti)
            
            if apumuuttujaR != '':
                muuttujaR = apumuuttujaR
            print 'variable is', muuttujaR
            print ' '
            
            printti = "Case is "+ arguments[i+1] + ", the variable is currently: " + muuttujaN + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. prcp/rmH2Opr Nc_ic/CCN): '
            apumuuttujaN = raw_input(printti)
            
            if apumuuttujaN != '':
                muuttujaN = apumuuttujaN
            print 'variable is', muuttujaN
            print ' '
            
            
            
        if EMUL:
            if not EMULCASES:
                case_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
            else:
                case_indeksi = emulCaseIndexes[i]
        else:
            case_indeksi = i
            
        dn0_data     = mdp.read_Data( filenameNC[i], 'dn0'     )
        time_data    = mdp.read_Data( filenameNC[i], 'time'    )
        zt_data      = mdp.read_Data( filenameNC[i], 'zt'      )
        
        Radius_data = mdp.read_Data( filenameNC[i], muuttujaR )
        Number_data = mdp.read_Data( filenameNC[i], muuttujaN  )
        
        t =  np.argmin(np.abs(aikaT*3600-time_data))
        z =  np.argmin(np.abs(korkeusH-zt_data))
                
        
        maksimi = np.max( Number_data )

        ysize = np.shape(Radius_data)[2]
        xsize = np.shape(Radius_data)[3]
        
        Diameter_data = Radius_data[t,:,:,:,z]* muunnosKerroinR
        Number_data  = Number_data[t,:,:,:,z]* dn0_data[z]
        
        
        Diameter_data_flatten = Diameter_data.flatten()
        Number_data_flatten  = Number_data.flatten()
        
        jarj = np.argsort(Diameter_data_flatten)
        
        Diameter_data_flatten = np.reshape( np.sort(Diameter_data_flatten) , (-1,1) )
        koko =  np.shape(Number_data_flatten)
        apu = np.zeros(koko)
        
        for kkk in xrange(koko[0]):
            apu[kkk] = Number_data_flatten[ jarj[kkk] ]
        
        Number_data_flatten = np.reshape(apu, (-1,1))
        
        ### decision tree regression
        xnew = np.linspace(np.min(Diameter_data_flatten),np.max(Diameter_data_flatten), 1000, endpoint = True)[:,np.newaxis]
        
        rng = np.random.RandomState(1)
        regr_2 = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4),n_estimators=500, random_state=rng)
        regr_2.fit(Diameter_data_flatten, Number_data_flatten)
        y_2 = regr_2.predict(xnew)
        ### end decision tree regression
        
        
            
        #if customLabels:
                #label = labelArray[i]
         #else:    
                #label = str(case_indeksi+1)
        #else:
        
        label = labelArray[i]
        
        #tit = typename + ' size distribution, slize t = ' + str(round(time_data[t]/3600.,1)) + ' [h], h = ' + str(zt_data[z]) + ' [m] ' + label
        tit ='slize t = ' + str(round(time_data[t]/3600.,1)) + ' [h], h = ' + str(zt_data[z]) + ' [m] ' + label
        if xCustSize is None:
            xlist = range(xsize)
        else:
            a = xCustSize/2
            b = xCustSize - a
            xlist = map(int, np.arange( ceil(int(xsize/2)-a+.1), ceil(int(xsize/2)+b+.1),1))
        
        if yCustSize is None:
            ylist = range(ysize)
        else:
            a = yCustSize/2
            b = yCustSize - a
            ylist = map(int, np.arange( ceil(int(ysize/2)-a+.1), ceil(int(ysize/2)+b+.1),1))
            
        print 'kokojakauma x ja y indeksit', xlist, ylist    
        mdp.initializeColors(len(xlist)*len(ylist), shuffling = sekoita)
        
        oikeatXtikit = np.arange( 0, xmax+0.01,  0.2)
        for x in xlist:
            for y in ylist:
                if figuuri is None:
                    uusifig = True
                else:
                    uusifig = False
                    ylabel = ' '
                
                fig, ax =  mdp.plottaa(Diameter_data[:,y,x], Number_data[:,y,x],xl=xlabel, yl=ylabel, tit = ' ', uusikuva = uusikuva, scatter = True, LEGEND = False, tightXAxis = True, tightYAxis = True, gridi = False, sub_i = 1, sub_j = LKM, sub_k = i+1, figuuri = None, uusisub = True, uusifig = uusifig )
                
                figuuri = fig
                ax.title.set_text(tit)
                ax.title.set_fontsize(22)
                
                
                ax.set_xticklabels( map(str, oikeatXtikit ))
                ax.set_xticks( oikeatXtikit )
                
                j = 0
                for xxlabel in ax.xaxis.get_ticklabels():
                    if np.mod(j-1,4) != 0 or j == 0:
                    #if j==0:
                        xxlabel.set_visible(False)
                    j+=1
                xxlabel.set_visible(True)
                if i > 0:
                    for yylabel in ax.yaxis.get_ticklabels():
                        yylabel.set_visible(False)
                    
                
        if interplo:
            mdp.plottaa(xnew, y_2,xl=xlabel, yl=ylabel, tit = ' ', uusikuva = uusikuva, omavari='k', LEGEND = False, sub_i = 1, sub_j = LKM, sub_k = i+1)
        
        figuuri.suptitle( typename + ' size distribution') #, fontsize=45

        
        # jos ymin ja ymax arvoja ei ole ennalta annettu, niin kaytetaan kuvan raja-arvoina laskettuja arvoja
        if ymax is None:
            ymax = maksimi
        
        if xmax is None:
            xmax = np.max(Diameter_data_flatten)
        
        if ( asetaRajat ):
            plt.ylim( 0., 1.1*ymax )
            plt.xlim( 0., xmax )
        
        
        

        
        
            
        
        
            
    if savePrefix is None:  
        savePrefix = typename+ '_size_distribution'

    if saveFig:
        plt.savefig( picturefolder + savePrefix  + '_' + saveTag+ '_'+ str(round(time_data[t]/3600.,1)) + 'h_' + str(zt_data[z]) + 'm' + LVLprintSave + '.png')


   
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
            dn0    = mdp.read_Data( filenameNC[i], 'dn0' )
            NumBin = np.multiply( NumBin, dn0 )
        
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
        plt.savefig( picturefolder + savePrefix+ '_' + saveTag + LVLprintSave + '.png')
        
###########################
def piirra_MeanSize( tyyppi = 'ice', bini='a', ajanhetket = [0], korkeus = [0], useDN = True, color = 'b', savePrefix = None   ):
    from scipy.stats.mstats import gmean
    if tyyppi not in ['aerosol', 'cloud', 'precipitation', 'ice', 'snow' ]:
        sys.exit("wrong type")
        
    tyyppi1 = tyyppi[0:1]
    tyyppi2 = tyyppi[0:2]
    
    Rad_name = 'S_Rw'   + tyyppi1  + bini
    
    Num_name = 'S_N'   + tyyppi1
    
    print tyyppi, Rad_name, Num_name
    
    
    for i in xrange(len(arguments)-1):
        Num  = np.multiply( mdp.read_Data( filenameNC[i], Num_name), 1./1000. )
        Rad  = mdp.read_Data( filenameNC[i], Rad_name)

        
        if useDN:
            dn0 = mdp.read_Data( filenameNC[i], 'dn0' )
            Num = np.multiply( Num, dn0 )
            
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
            HslizeSTR +=  ' ' + r'$h_1$' + ' = ' + str( zt[ max(korkeusindeksit) ] )
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
        #plt.savefig( picturefolder + savePrefix+ '_' + saveTag + LVLprintSave + '.png')

#############################        
def piirra_domainMeanProfiili( muuttuja, muuttujaPainotus =  'S_Niba', muuttujaPainotusPotenssi =  'S_Rwiba', nimi = None, muunnosKerroin = 1.0, ajanhetket = [0], useDN = True, profiili = False, binidata = False, xAxisL = '', color = 'k', savePrefix = None, akselit = False, xmax = None, askChangeOfVariable = False   ):
        
    minimi  = None
    maksimi = None
    
    if nimi is None:
        nimi = muuttuja
        
    if profiili:
        tiedostonimi = filenamePS
    elif (not profiili)  or binidata:
        tiedostonimi = filenameNC        
    
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        if askChangeOfVariable:
            printti = "Case is "+ arguments[i+1] + ", the variable is currently: " + muuttuja + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. prcp/rmH2Opr Nc_ic/CCN): '
            apumuuttuja = raw_input(printti)
            if apumuuttuja != '':
                muuttuja = apumuuttuja
            print 'variable is', muuttuja
            print ' '
            
            printti = "Case is "+ arguments[i+1] + ", the muuttujaPainotus is currently: " + muuttujaPainotus + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. prcp/rmH2Opr Nc_ic/CCN): '
            muuttujaPainotusApu = raw_input(printti)
            if muuttujaPainotusApu != '':
                muuttujaPainotus = muuttujaPainotusApu
            print 'muuttujaPainotus is', muuttujaPainotus
            print ' '
            
            printti = "Case is "+ arguments[i+1] + ", the muuttujaPainotus is currently: " + muuttujaPainotusPotenssi + ' if you want to keep it press [Enter], otherwise give the name of the variable (e.g. prcp/rmH2Opr Nc_ic/CCN): '
            muuttujaPainotusPotenssiApu = raw_input(printti)
            if muuttujaPainotusPotenssiApu != '':
                muuttujaPainotusPotenssi = muuttujaPainotusPotenssiApu
            print 'muuttujaPainotus is', muuttujaPainotusPotenssi
            print ' '
                
        data  = np.multiply( mdp.read_Data( tiedostonimi[i], muuttuja), muunnosKerroin )
        zt    = np.asmatrix( mdp.read_Data( tiedostonimi[i], 'zt') )
        
        dn0 = mdp.read_Data( tiedostonimi[i], 'dn0' )
        if useDN:
            data = np.multiply( data, dn0 )
        
        
        ######### aika slaissaus
        time = mdp.read_Data( filenameNC[i], 'time')
        aikaindeksit = []
        for t in ajanhetket:
            aikaindeksit.append( np.argmin( np.abs(time - t*3600.) ))
        
        Tslize = map( int, np.arange( min(aikaindeksit), max(aikaindeksit)+0.5 ) )
        
        TslizeSTR = r'$t_{0} = $' + str( round( time[ min(aikaindeksit) ]/3600.,1) ) + ' h'
        if len(aikaindeksit)>0:
            TslizeSTR += ' to ' + r'$t_{1} = $' + str( round( time[ max(aikaindeksit) ]/3600.,1) ) + ' h'
            
            TslizeSTR = 'from ' + TslizeSTR
        else:
            TslizeSTR = 'at ' + TslizeSTR
        #print TslizeSTR
        ###############################
        
        if profiili:
            dataSlize  = data[ Tslize,   : ]
        elif binidata:
            dataSlize  = data[ Tslize, :, :, :, : ]
        elif not profiili:
            dataSlize  = data[ Tslize, :, :, : ]
        
        if len(Tslize)>1:
            dataSlize = np.mean( dataSlize, axis = 0 )
        
        if profiili:
            dataSlizeMean = dataSlize
        elif binidata:
            from math import pi
            rho = 44.2*6./pi
            painotus =  np.multiply( mdp.read_Data( tiedostonimi[i], muuttujaPainotus), np.multiply( np.power(mdp.read_Data( tiedostonimi[i], muuttujaPainotusPotenssi ),3), dn0*44.2) )  
            painotusSlize = painotus[ Tslize, :, :, :, : ]
            if len(Tslize)>1:
                painotusSlize = np.mean( painotusSlize, axis = 0 )
                
            nimittaja = np.sum( np.sum( np.sum(painotusSlize, axis = 0), axis= 0), axis=0 )
            nimittaja = np.where( nimittaja > 0.,  np.power(nimittaja, -1), 0.)
            
            osoittaja = np.sum( np.sum( np.sum(np.multiply(dataSlize, painotusSlize), axis = 0), axis= 0), axis=0 )
            
            dataSlizeMean = np.multiply( osoittaja, nimittaja)
        elif not profiili:
            dataSlizeMean = np.mean( np.mean( dataSlize, axis = 0), axis = 0)
        
        tit = nimi + ' ' + TslizeSTR
        
        if maksimi is not None:
            maksimi = max( maksimi,  np.max(dataSlizeMean) )
        else:
            maksimi = np.max(dataSlizeMean)

        if minimi is not None:
            minimi = min( minimi,  np.min(dataSlizeMean) )
        else:
            minimi = np.min(dataSlizeMean)
        
        dataSlizeMean = np.asmatrix(dataSlizeMean)
        
        if EMUL:
            if customLabels:
                label = labelArray[i]
            else:    
                label = str(case_indeksi+1)
        else:
            label = labelArray[i]       

        fig, ax = mdp.plottaa( dataSlizeMean.A1, zt.A1, tit , xl = xAxisL, yl='height [m]', changeColor=True, tightXAxis=True, tightYAxis = True, markers=False, LEGEND=True, label = label, omavari = color, scatter=False, uusikuva=uusikuva )       
    ax.title.set_fontsize(35)
    mdp.plot_setXlim( minimi, maksimi, extendBelowZero = False, A = 0.05 )
    
    if akselit or (xmax is not None):
        if xmax is None:
            xmax = np.max(dataSlizeMean)
        
        oikeatXtikit = map(int, np.arange( 0, xmax+0.01,  50))
        ax.set_xticklabels( map(str, oikeatXtikit ))
        ax.set_xticks( oikeatXtikit )
        
    if savePrefix is None:
        savePrefix = 'domainMeanProfiili'
    
    if saveFig:
        plt.savefig( picturefolder + savePrefix+ '_' + muuttuja + '_' + saveTag + LVLprintSave + '.png')
        
def piirra_domainMeanScatter( muuttujaX, muuttujaY, muuttujaYY = None, muunnosKerroinX = 1.0, muunnosKerroinY = 1.0, muunnosKerroinYY = 1.0, nimi = None, ajanhetket = [0], useDNX = False,useDNY = False, profiiliX = False, profiiliY = False, xaxislabel = None, yaxislabel = None, color = None, savePrefix = None, minimiX =  None, maksimiX = None, minimiY =  None, maksimiY = None    ):
        
    minimi  = None
    maksimi = None
    
    if nimi is None:
        nimi = muuttujaX + muuttujaY + mdp.xstr(muuttujaYY)
        
    #if profiiliX:
        #tiedostonimiX = filenamePS
    #else:
        #tiedostonimiX = filenameNC        
        
    #if profiiliY:
        #tiedostonimiY = filenamePS
    #else:
        #tiedostonimiY = filenameNC
        
    varit = cycle( ['b','r','k','g','c','m','y'][:(len(arguments)-1)] )
    
    if (len(arguments)-1)>7:
        changeColor = True
        color = False
    else:
        changeColor = False
        color = True
    
    for i in xrange(len(arguments)-1):
        uusikuva = True if i == 0 else  False
        dataX  = np.multiply( mdp.read_Data( filenameTS[i], muuttujaX), muunnosKerroinX )
        dataY  = np.multiply( mdp.read_Data( filenameTS[i], muuttujaY), muunnosKerroinY )
        
        
        
        if muuttujaYY is not None:
            dataY = dataY + np.multiply( mdp.read_Data( filenameTS[i], muuttujaYY), muunnosKerroinYY )
        
        #zt    = np.asmatrix( mdp.read_Data( tiedostonimi[i], 'zt') )
        if useDNX:
            dn0  = mdp.read_Data( filenameTS[i], 'dn0' )
            dataX = np.multiply( dataX, dn0 )
        
        if useDNY:
            dn0   = mdp.read_Data( filenameTS[i], 'dn0' )
            dataY = np.multiply( dataY, dn0 )
        
        
        ######### aika slaissaus
        time = mdp.read_Data( filenameTS[i], 'time')
        aikaindeksit = []
        for t in ajanhetket:
            aikaindeksit.append( np.argmin( np.abs(time - t*3600.) ))
        
        Tslize = map( int, np.arange( min(aikaindeksit), max(aikaindeksit)+0.5 ) )
        
        TslizeSTR = r'$t_{0} = $' + str( round( time[ min(aikaindeksit) ]/3600.,1) ) + ' h'
        if len(aikaindeksit)>0:
            TslizeSTR += ' to ' + r'$t_{1} = $' + str( round( time[ max(aikaindeksit) ]/3600.,1) ) + ' h'
            
            TslizeSTR = 'from ' + TslizeSTR
        else:
            TslizeSTR = 'at ' + TslizeSTR
        print ' '
        print TslizeSTR
        ###############################
        dataX = np.where( dataX > -999.,  dataX, 0. )
        dataY = np.where( dataY > -999.,  dataY, 0. )
        
        dataSlizeX  = dataX[ Tslize ]
        dataSlizeY  = dataY[ Tslize ]
        
        if len(Tslize)>1:
            dataSlizeMeanX = np.mean( dataSlizeX, axis = 0 )
        
        if len(Tslize)>1:
            dataSlizeMeanY = np.mean( dataSlizeY, axis = 0 )
        
        
        
        tit = nimi + ' ' + TslizeSTR
        
#        if maksimi is not None:
#            maksimi = max( maksimi,  np.max(dataSlizeMean) )
#        else:
#            maksimi = np.max(dataSlizeMean)
#
#       if minimi is not None:
#            minimi = min( minimi,  np.min(dataSlizeMean) )
#        else:
#            minimi = np.min(dataSlizeMean)
#        
#        dataSlizeMean = np.asmatrix(dataSlizeMean)
        
#        if EMUL:
#            if customLabels:
#                label = labelArray[i]
#            else:    
#                label = str(case_indeksi+1)
#        else:
#            label = labelArray[i]       
        if color is not False:
            color = next(varit)
            
        print 'case', filenameTS[i], 'muuttujaX', muuttujaX, 'dataSlizeMeanX', dataSlizeMeanX, 'dataSlizeMeanY', dataSlizeMeanY, 'color', color
        
        fig, ax = mdp.plottaa( dataSlizeMeanX, dataSlizeMeanY, a = 20, b = 20, tit = tit , xl = xaxislabel, yl = yaxislabel, changeColor=changeColor, tightXAxis=True, tightYAxis = True, LEGEND=False, omavari = color, scatter=True, uusikuva=uusikuva, marker= '+', markersize = 90 )       

        ax.title.set_fontsize(35)
        
        if (minimiX is not None) and (maksimiX is not None):
            plt.xlim( minimiX, maksimiX )
        
        if (minimiY is not None) and (maksimiY is not None):    
            plt.ylim( minimiY, maksimiY )

        if savePrefix is None:
            savePrefix = 'scatterSlize'

        if saveFig:
            plt.savefig( picturefolder + savePrefix+ '_' + muuttujaX + muuttujaY + '_' + saveTag + LVLprintSave + '.png')        

def animoi_path(muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = [0,1], colorBarTickNames = ['0','1'], xlabels = None, ylabels = None, xticks = None, yticks = None, variKartta = plt.cm.Blues, profiili = False, spinup = None ):        
    
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
if __name__ == "__main__":

    aika_color      = truncate_colormap(  plt.cm.Blues, minval = 0.3)
    aikaPisteet     = xTicksSeconds
    aikaBAR         = varibaari( aikaPisteet,  aika_color)
    aikaTIT         = 'time [h]'
    cbvalT          =  xTicksSeconds
    cbvalTStr       = map(str, ticksHours)

if EMUL:
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
    RAD   = True
    
    refColorbarSwitch = True
    askChangeOfVariable = False
    
    
    #CFRAC = False
    #CLOUD = True
    #PRCP  = False
    #LWP   = False
    #SHF   = False
    #LHF   = False
    #QDIFF = False
    #TDIFF = False
    #CDNC  = False
    #WMAX  = False

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
        liqW='l'
    else:
        prcp = 'rmH2Opr'
        sadekerroin = 2.5e+06 # latent heat of vaporization J/kg
        liqW = 'P_rl'



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
        piirra_aikasarjasettii( muuttuja = 'cfrac',  variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName= 'Cloud fraction',               ylabel = 'Cloud fraction',      variKartta = thickness_color, ymin = 0.0, ymax = 1.0, tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, omaVari = refColorbarSwitch, legenda = LEGEND )
        
        mdp.plot_suljetus(naytaPlotit)
        
    if CLOUD:
        piirra_aikasarjasettii( muuttuja = 'zc',     variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName= 'Relative change of cloud top', ylabel = 'relative change',     variKartta = thickness_color, relative = True, savePrefix = 'cloud_top_rel_change', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, yticks = zcTicks, ylabelFont = 18, omaVari = refColorbarSwitch, legenda = LEGEND  )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_profiiliKehitys( liqW, muunnosKerroin = 1000.0, variKartta = aika_color, colorBar = aikaBAR, colorBarTickValues = cbvalT, colorBarTickNames = cbvalTStr, longName =  'Cloud ' + r'$H_2{O}$' + ' mix. rat. evolution', xlabel = 'q_cloud [g/kg]', ylabel = 'z [m]', savePrefix = 'q_cloud_evol', aikaPisteet = aikaPisteet, tit = aikaTIT, rajaKerros = pblh_design, asetaRajat = False, paksuus = thickness_design, askChangeOfVariable = askChangeOfVariable )
        
        mdp.plot_suljetus(naytaPlotit)
        
    if PRCP:
        piirra_aikasarjasettii( muuttuja = prcp,     muunnosKerroin = sadekerroin, variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Surface precipitation',       ylabel = 'Precipitation [' + r'$W/m^2$' + ']', variKartta = thickness_color, ymin = 0.0, savePrefix = 'prcp', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, omaVari = refColorbarSwitch, askChangeOfVariable = askChangeOfVariable, legenda = LEGEND )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_maksimiKeissit( muuttuja = prcp, muunnosKerroin = sadekerroin, longName = 'Maximum precipitation after '+str(int(ajanhetket/3600.))+'h', ylabel = 'precipitation W/m^2',      savePrefix = 'prcp_max', askChangeOfVariable = askChangeOfVariable, ajanhetket = ajanhetket )
        
        mdp.plot_suljetus(naytaPlotit)
    
    if SHF:
        piirra_aikasarjasettii( muuttuja = 'shf_bar', variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Sensible heat flux',        ylabel = 'Sensible heat flux W/m^2', variKartta = thickness_color, savePrefix = 'heat_flux_sensible', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, omaVari = refColorbarSwitch, legenda = LEGEND )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_maksimiKeissit( maksimiSensible, longName = "Maximum sensible heat",                                       ylabel = 'Sensible heat flux W/m^2', savePrefix = 'heat_flx_sensible_max', ajanhetket = ajanhetket ) 
        
        mdp.plot_suljetus(naytaPlotit)
        
    if LHF:
        piirra_aikasarjasettii( muuttuja = 'lhf_bar', variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Latent heat flux',        ylabel = 'Latent heat flux W/m^2', variKartta = thickness_color, savePrefix = 'heat_flux_latent', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, omaVari = refColorbarSwitch, legenda = LEGEND )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_maksimiKeissit( maksimiLatent,   longName = "Maximum latent heat",                                         ylabel = 'Latent heat flux W/m^2',   savePrefix = 'heat_flx_latent_max', ajanhetket = ajanhetket ) 
        
        mdp.plot_suljetus(naytaPlotit)
        
    if LWP:
        piirra_aikasarjasettii( muuttuja = 'lwp_bar', muunnosKerroin = 1000.0, variKartta = thickness_color,  variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'LWP', ylabel = 'LWP g/m^2', ymin = 0.0,  savePrefix = 'lwp', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, omaVari = refColorbarSwitch, legenda = LEGEND )
        
        mdp.plot_suljetus(naytaPlotit)
        
    if CDNC:
        piirra_aikasarjasettii( muuttuja = 'Nc_ic', variKartta = num_pbl_color, variRefVektori = num_pbl_design, colorBar = cdncBAR, colorBarTickValues = cbvalQ, colorBarTickNames = cbvalQStr, longName = 'Relative change of in-cloud CDNC', ylabel = 'Relative change of In-cloud CDNC #/kg', tit = cdncTIT, relative = True, nollaArvo = 1.e6*num_pbl_design,  spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, omaVari = refColorbarSwitch, askChangeOfVariable = askChangeOfVariable, legenda = LEGEND) #
        
        mdp.plot_suljetus(naytaPlotit)
        
    if QDIFF:
        piirra_profiilisettii( 'q', variKartta = pblh_color, variRefVektori = pblh_design, colorBar = pblhBAR, colorBarTickValues = cbvalPblh, colorBarTickNames = cbvalPblhStr, longName =  'Tot. ' + r'$H_2{O}$' + ' mix. rat. ' + ' relative change ' +'0h - '+str((ajanhetket/3600.)) + 'h', xlabel = r'$\frac{q_{t=2h}}{q_{t=0h}}-1$', ylabel = 'z/pblh', ymin = 0.0, ymax = 1.01, xmin = -0.35, xmax = 0.025, savePrefix = 'q_diff', ajanhetket = ajanhetket, tit = pblhTIT, rajaKerros = pblh_design, relative = True, nollaArvo = q_pbl_design, omaVari = refColorbarSwitch )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_profiiliKehitys( 'q', muunnosKerroin = 1000.0, variKartta = aika_color, colorBar = aikaBAR, colorBarTickValues = cbvalT, colorBarTickNames = cbvalTStr, longName =  'Tot. ' + r'$H_2{O}$' + ' mix. rat. evolution', xlabel = 'q [g/kg]', ylabel = 'z [m]', savePrefix = 'q_evol', aikaPisteet = aikaPisteet, tit = aikaTIT, rajaKerros = pblh_design, asetaRajat = False, paksuus = thickness_design )
        
        mdp.plot_suljetus(naytaPlotit)

    if TDIFF:
        piirra_profiilisettii( 't', variKartta = pblh_color, variRefVektori = pblh_design, colorBar = pblhBAR, colorBarTickValues = cbvalPblh, colorBarTickNames = cbvalPblhStr, longName =  r'$\theta$' + ' relative change 0h - '+str((ajanhetket/3600.)) + 'h', xlabel = r'$\frac{\theta_{t=2h}}{\theta_{t=0h}}-1$', ylabel = 'z/pblh', ymin = 0.0, ymax = 1.01, xmin = -0.05, xmax = 0.025, savePrefix = 't_diff', ajanhetket = ajanhetket, tit = pblhTIT, rajaKerros = pblh_design, relative = True, nollaArvo = tpot_pbl_design, omaVari = refColorbarSwitch )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_profiiliKehitys( 't',  variKartta = aika_color, colorBar = aikaBAR, colorBarTickValues = cbvalT, colorBarTickNames = cbvalTStr, longName =  'Potential temperature', xlabel = r'$\theta$' + ' [K]', ylabel = 'z [m]', savePrefix = 'theta_evol', aikaPisteet = aikaPisteet, tit = aikaTIT, rajaKerros = pblh_design, asetaRajat = False, paksuus = thickness_design )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_profiiliKehitys( 't',  variKartta = aika_color, colorBar = aikaBAR, colorBarTickValues = cbvalT, colorBarTickNames = cbvalTStr, longName =  'Absolute temperature', xlabel = 'theta' + ' [K]', ylabel = 'z [m]', savePrefix = 'temp_evol', aikaPisteet = aikaPisteet, tit = aikaTIT, rajaKerros = pblh_design, asetaRajat = False, paksuus = thickness_design, tempConversion = True )
        
        mdp.plot_suljetus(naytaPlotit)
        
    if WMAX:
        piirra_aikasarjasettii( muuttuja = 'wmax', variRefVektori = thickness_design, colorBar =  thickBAR, colorBarTickValues = cbvalThick, colorBarTickNames = cbvalThickStr, longName = 'Maximum vertical velocity',        ylabel = 'Maximum vertical velocity m/s', variKartta = thickness_color, savePrefix = 'w_max', tit = thickTIT, spinup = spinup, xlabels = xLabelsHours, xticks = xTicksSeconds, omaVari = refColorbarSwitch, legenda = LEGEND )
        
        mdp.plot_suljetus(naytaPlotit)
    
    #if RAD:
        

    if NCFwrite:
        maksimiLWPlist_ncf    = ncfile.createVariable( 'lwp_max',     np.dtype('float32').char, ('case') )
        maksimiLWPlist_ncf[:] = maksimiLWPlist
        print 'eniten sadetta keissi', np.argmax(maksimiLWPlist)
    ncfile.close()
################

if ICE:
    

    naytaPlotit = False
    korkeustikit =np.arange(0, 1250, 50)
    ylabels    = map(str, korkeustikit )
    
    profiiliVariLIQ = [ '#000099', '#00ccff', '#00e600', '#f9f906', '#ff9900', '#ff0000' ]
    profiiliVariICE = [ '#000099', '#00ccff', '#29a385', '#00e600', '#c8e600', '#f9f906', '#ff9900', '#ff0000' ] # , '#990000', '#660000'
    cbvalLIQ    = np.arange(0, 0.241, 0.04)
    cbvalLIQStr = map(str, cbvalLIQ)

    if tag[:-1] == 'ice1':
        cbvalICE    = np.arange( 0, 0.41, 0.05) # np.arange(0, 1.4, 0.1)
        nicLKM = 1.41
    else:
        cbvalICE    =  np.arange(0, 1.41, 0.1)
        nicLKM = 4.41
        
    cbvalICEStr = map(str, cbvalICE)
    
    cbvalLIQPATH = np.arange(0, 61, 5)
    cbvalLIQPATHStr = map(str, cbvalLIQPATH)
    
    piilotaOsaXlabel = True
    
    cbvalDiam = np.arange(0,600, 75)
    cbvalDiamStr = map(str,cbvalDiam)
    
    diamVari = profiiliVariICE #[plt.cm.tab20(3),plt.cm.tab20(2),plt.cm.tab20c(1), plt.cm.tab20c(0),plt.cm.tab20c(7),plt.cm.tab20c(6),plt.cm.tab20c(5),plt.cm.tab20c(4),plt.cm.tab20c(11),plt.cm.tab20c(10),plt.cm.tab20c(9),plt.cm.tab20c(8) ]
    askChangeOfVariable = raw_input( "Snow stuff included? (yes/no) Then will ask askChangeOfVariable " ) in kylla
    if len(arguments)-1>1:
        icevari = False
    
    
    if int(lvl)>=4:
    
        #piirra_aikasarjaPathXYZ( 'l', longName = 'Liquid Water Path', savePrefix = 'lwp', xaxislabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel )
        #mdp.plot_vertical( spinup )
        #plt.xticks( ticksHours, xLabelsHours )
        
        piirra_aikasarjasettii( muuttuja = 'lwp_bar', muunnosKerroin = 1000.0, longName = 'Liquid water path', ylabel = 'path ' + r'[$g/m^2$]', ymin = 0.0,  savePrefix = 'lwpTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True  )
        
        mdp.plot_suljetus(naytaPlotit)
        
        #mdp.plot_vertical( spinup )
        #plt.xticks( ticksHours, xLabelsHours )
        
        #piirra_aikasarjasettii( muuttuja = 'shf_bar',  longName = 'Sensible heat flux', ylabel = 'W/m^2', ymin = 0.0,  savePrefix = 'heat_flux_sensible', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel  )
        
        #piirra_aikasarjasettii( muuttuja = 'lhf_bar',  longName = 'Latent heat flux', ylabel = 'W/m^2', ymin = 0.0,  savePrefix = 'heat_flux_latent', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel  )
        
        
        
        #piirra_domainProfiili( 'l', muunnosKerroin = 1000., longName = "Liquid water mixing ratio  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQ, colorBarTickNames = cbvalLIQStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = profiiliVariLIQ, spinup = spinup )
        
        piirra_domainProfiili( 'P_rl', muunnosKerroin = 1000., longName = "Liquid water mixing ratio  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQ, colorBarTickNames = cbvalLIQStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = profiiliVariLIQ, spinup = spinup, profiili = True ) # variKartta = profiiliVariLIQ
        
        mdp.plot_suljetus(naytaPlotit)
        
        
        #animoi_path( 'l', muunnosKerroin = 1000., longName = "Liquid water path  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQPATH, colorBarTickNames = cbvalLIQPATHStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.Reds, spinup = spinup )
        
        #piirra_MeanSize(tyyppi = 'cloud', ajanhetket = [6], korkeus = [700], color = 'r')
        
        
        
        
        #piirra_domainMeanProfiili( 'S_Nc',nimi = 'Cloud number concentration averaged',  muunnosKerroin=1.e-6, ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        
        piirra_domainMeanProfiili( 'P_Nca',nimi = 'Cloud number concentration averaged',  muunnosKerroin=1.e-6, ajanhetket = [6,8], useDN = True, profiili = True, xAxisL = r'[#$/cm^{3}$]', color = icevari )
        
        mdp.plot_suljetus(naytaPlotit)
        
        
        
        #piirra_domainProfiili( 'w_2', longName = "vertical velocity squared " + r'$m^{2}/s^{-2}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.RdPu, profiili = True, spinup = spinup )
        
        #piirra_domainProfiili( 'rflx', longName = "Radiative cooling " + r'$W/m^{2}$', useDN = False, transpose = True, profiili = True, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.bwr, spinup = spinup, radCool = True, cloudBaseTop = True )
        
        #piirra_profiiliKehitys( 'theta',  variKartta = aika_color, colorBar = aikaBAR, colorBarTickValues = cbvalT, colorBarTickNames = cbvalTStr, longName =  'Potential temperature', xlabel = r'$\theta$' + ' [K]', ylabel = 'z [m]', savePrefix = 'theta_evol', aikaPisteet = aikaPisteet, tit = aikaTIT, asetaRajat = False )
        
        #piirra_aikasarjasettii( muuttuja = 'rmH2Opr',     muunnosKerroin = 2.5e+06,  longName = "Removal of water by sedimentation of rain",       ylabel = 'flux ' + r'[$W/m^2$]', ymin = 0.0,  savePrefix = 'depRain', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True  )
        
        piirra_aikasarjasettii( muuttuja = 'Nc_ic', muunnosKerroin = 1.e-6, longName = 'in-cloud CDNC', ylabel = '[#/mg]', ymin = 0.0,  savePrefix = 'cdncTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True  )
        
        mdp.plot_suljetus(naytaPlotit)
        
        #piirra_domainProfiili( 'theta', longName = 'Potential temperature', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Reds, spinup = spinup, profiili = True  )

        #piirra_domainProfiili( 'thl', longName = 'Liquid water potential temperature', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Reds, spinup = spinup, profiili = True  )
        
        #piirra_aikasarjasettii( muuttuja = 'vtke', muunnosKerroin = 1.0, longName = 'Vertical integral of total TKE', ylabel = '[kg/s]', ymin = 0.0,  savePrefix = 'vtkeTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True  )
        
        mdp.plot_suljetus(naytaPlotit)
        
    if int(lvl)== 5:
        #piirra_aikasarjaPathXYZ( 'i', longName = 'Ice Water Path', savePrefix = 'iwp', xaxislabel = 'time [h]',xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel )
        #mdp.plot_vertical( spinup )
        #plt.xticks( ticksHours, xLabelsHours )
        
        piirra_aikasarjasettii( muuttuja = 'iwp_bar', muunnosKerroin = 1000.0, longName = 'Ice water path', ylabel = 'path ' + r'[$g/m^2$]', ymin = 0.0,  savePrefix = 'iwpTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True, askChangeOfVariable = askChangeOfVariable  )
        
        mdp.plot_suljetus(naytaPlotit)
        
        #piirra_aikasarjasettii( muuttuja = 'swp_bar', muunnosKerroin = 1000.0, longName = 'Ice water path', ylabel = 'path ' + r'[$g/m^2$]', ymin = 0.0,  savePrefix = 'swpTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True  )

        #piirra_domainProfiili( 'i', muunnosKerroin = 1000.*np.power(10.,2), longName = 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup  )
        
        piirra_domainProfiili( 'P_ri', muunnosKerroin = 1000.*np.power(10.,2), longName = 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup, profiili = True, askChangeOfVariable = askChangeOfVariable  ) #variKartta = profiiliVariICE, colorBarTickValues = cbvalICE plt.cm.Blues
        
        mdp.plot_suljetus(naytaPlotit)
        
        #piirra_domainProfiili( 'P_ri', muunnosKerroin = 1000.*np.power(10.,4), longName = 'Ice mixing ratio ' + r'$10^{4}g/kg^{-1}$', useDN = False, transpose = True, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Blues, spinup = spinup, profiili = True  ) #variKartta = profiiliVariICE, colorBarTickValues = cbvalICE plt.cm.Blues
        
        nic = np.arange(0, nicLKM, 0.2)
        piirra_domainProfiili( 'P_Nia', muunnosKerroin = 1.e-3, longName = 'Ice number concentration ', useDN = False, transpose = True, colorBarTickValues = nic, colorBarTickNames = map(str,nic), xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup, profiili = True, askChangeOfVariable = askChangeOfVariable  ) #variKartta = profiiliVariICE, colorBarTickValues = cbvalICE 
        
        mdp.plot_suljetus(naytaPlotit)
        #piirra_domainProfiili( 'P_ri', muunnosKerroin = 1000.*np.power(10.,2), longName = 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Blues, spinup = spinup, profiili = True  )
        
        #animoi_path( 'f', muunnosKerroin = 1000., longName = tag + "Ice water path  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQPATH, colorBarTickNames = cbvalLIQPATHStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.Blues, spinup = spinup )
        ## muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = [0,1], colorBarTickNames = ['0','1'], xlabels = None, ylabels = None, xticks = None, yticks = None
        
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [700] )
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [400] )
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [200] )
        #piirra_domainMeanProfiili( 'S_Ni',   nimi = 'Ice number concentration averaged',   muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        ##piirra_domainMeanProfiili( 'S_Ni',   nimi = 'Ice number concentration averaged',   muunnosKerroin=1./1000., ajanhetket = [3,6], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        #piirra_domainMeanProfiili( 'S_Rwia', nimi = 'Ice particle mean diameter averaged', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$D [{\mu}m]$', color = icevari )   

        piirra_domainMeanProfiili( 'P_Nia',  nimi = 'Ice number concentration averaged',   muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = True, xAxisL = r'[#$/L$]', color = icevari, askChangeOfVariable = askChangeOfVariable )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_domainMeanProfiili( 'P_Rwia', nimi = 'Ice particle mean diameter averaged', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, profiili = True, xAxisL = r'[${\mu}m$]', color = icevari, askChangeOfVariable = askChangeOfVariable )   
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_domainMeanProfiili( 'S_Rwiba', nimi = 'Ice particle mass mean diameter averaged testi', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, binidata = True, xAxisL = r'[${\mu}m$]', color = icevari, savePrefix = 'domainMassMeanProfiili_salsa', askChangeOfVariable = askChangeOfVariable )
        
        mdp.plot_suljetus(naytaPlotit)
        #piirra_domainMeanProfiili( 'S_Rwiba', nimi = 'Ice particle mean diameter averaged testi', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, binidata = True, xAxisL = r'[${\mu}m$]', color = icevari, savePrefix = 'domainMeanProfiili_salsa' )

        piirra_domainMeanProfiili( 'P_ri',  nimi = 'Ice mixing ratio averaged',   muunnosKerroin=1000., ajanhetket = [6,8], useDN = False, profiili = True, xAxisL = r'[g/kg]', color = icevari, askChangeOfVariable = askChangeOfVariable )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_domainProfiili( 'P_Rwia', muunnosKerroin = 2.e6, longName = 'Ice particle mean diameter '+r'[${\mu}m$]', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, spinup = spinup, profiili = True, colorBarTickValues = cbvalDiam, colorBarTickNames = cbvalDiamStr, variKartta = diamVari, askChangeOfVariable = askChangeOfVariable  )
        
        mdp.plot_suljetus(naytaPlotit)
        
        #piirra_domainProfiili( 'thil', longName = 'Ice-Liquid water potential temperature', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Reds, spinup = spinup, profiili = True  )
        
        piirra_aikasarjasettii( muuttuja = 'rmH2Oic',     muunnosKerroin = 2.5e+06,  longName = 'Removal of water by sedimentation of ice',       ylabel = 'flux ' + r'[$W/m^2$]', ymin = 0.0,  savePrefix = 'depIce', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, askChangeOfVariable = askChangeOfVariable  )
        
        mdp.plot_suljetus(naytaPlotit)
        
    if int(lvl)== 6:  # plain snow stuff
        print 'plain snow stuff'
        #piirra_aikasarjaPathXYZ( 'i', longName = 'Ice Water Path', savePrefix = 'iwp', xaxislabel = 'time [h]',xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel )
        #mdp.plot_vertical( spinup )
        #plt.xticks( ticksHours, xLabelsHours )
        
        piirra_aikasarjasettii( muuttuja = 'swp_bar', muunnosKerroin = 1000.0, longName = 'Ice water path (snow)', ylabel = 'path ' + r'[$g/m^2$]', ymin = 0.0,  savePrefix = 'swpTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True, askChangeOfVariable = askChangeOfVariable  )
        
        mdp.plot_suljetus(naytaPlotit)
        
        #piirra_aikasarjasettii( muuttuja = 'swp_bar', muunnosKerroin = 1000.0, longName = 'Ice water path', ylabel = 'path ' + r'[$g/m^2$]', ymin = 0.0,  savePrefix = 'swpTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True  )

        #piirra_domainProfiili( 'i', muunnosKerroin = 1000.*np.power(10.,2), longName = 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup  )
        
        piirra_domainProfiili( 'P_rs', muunnosKerroin = 1000.*np.power(10.,2), longName = 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalICE, colorBarTickNames = cbvalICEStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup, profiili = True, askChangeOfVariable = askChangeOfVariable  ) #variKartta = profiiliVariICE, colorBarTickValues = cbvalICE plt.cm.Blues
        
        mdp.plot_suljetus(naytaPlotit)
        
        #piirra_domainProfiili( 'P_ri', muunnosKerroin = 1000.*np.power(10.,4), longName = 'Ice mixing ratio ' + r'$10^{4}g/kg^{-1}$', useDN = False, transpose = True, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Blues, spinup = spinup, profiili = True  ) #variKartta = profiiliVariICE, colorBarTickValues = cbvalICE plt.cm.Blues
        
        nic = np.arange(0, nicLKM, 0.2)
        piirra_domainProfiili( 'P_Ns', muunnosKerroin = 1.e-3, longName = 'Ice number concentration ', useDN = False, transpose = True, colorBarTickValues = nic, colorBarTickNames = map(str,nic), xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup, profiili = True, askChangeOfVariable = askChangeOfVariable  ) #variKartta = profiiliVariICE, colorBarTickValues = cbvalICE 
        
        mdp.plot_suljetus(naytaPlotit)
        #piirra_domainProfiili( 'P_ri', muunnosKerroin = 1000.*np.power(10.,2), longName = 'Ice mixing ratio ' + r'$10^{2}g/kg^{-1}$', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Blues, spinup = spinup, profiili = True  )
        
        #animoi_path( 'f', muunnosKerroin = 1000., longName = tag + "Ice water path  " + r'$g/kg^{-1}$', useDN = False, transpose = True, colorBarTickValues = cbvalLIQPATH, colorBarTickNames = cbvalLIQPATHStr, xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit,  variKartta = plt.cm.Blues, spinup = spinup )
        ## muuttuja, muunnosKerroin = 1.0, transpose = False, longName = None , savePrefix = None, useDN = False, colorBarTickValues = [0,1], colorBarTickNames = ['0','1'], xlabels = None, ylabels = None, xticks = None, yticks = None
        
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [700] )
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [400] )
        #piirra_MeanSize(tyyppi = 'ice', ajanhetket = [6], korkeus = [200] )
        #piirra_domainMeanProfiili( 'S_Ni',   nimi = 'Ice number concentration averaged',   muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        ##piirra_domainMeanProfiili( 'S_Ni',   nimi = 'Ice number concentration averaged',   muunnosKerroin=1./1000., ajanhetket = [3,6], useDN = True, profiili = False, xAxisL = r'$N [L^{-1}]$', color = icevari )
        #piirra_domainMeanProfiili( 'S_Rwia', nimi = 'Ice particle mean diameter averaged', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = True, profiili = False, xAxisL = r'$D [{\mu}m]$', color = icevari )   

        piirra_domainMeanProfiili( 'P_Ns',  nimi = 'Ice number concentration averaged',   muunnosKerroin=1./1000., ajanhetket = [6,8], useDN = True, profiili = True, xAxisL = r'[#$/L$]', color = icevari, askChangeOfVariable = askChangeOfVariable )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_domainMeanProfiili( 'P_Rws', nimi = 'Ice particle mean diameter averaged', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, profiili = True, xAxisL = r'[${\mu}m$]', color = icevari, askChangeOfVariable = askChangeOfVariable )   
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_domainMeanProfiili( 'S_Rwsba', muuttujaPainotus = 'S_Nsba', muuttujaPainotusPotenssi =  'S_Rwsba',  nimi = 'Ice particle mass mean diameter averaged testi', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, binidata = True, xAxisL = r'[${\mu}m$]', color = icevari, savePrefix = 'domainMassMeanProfiili_salsa', askChangeOfVariable = askChangeOfVariable )
        
        mdp.plot_suljetus(naytaPlotit)
        #piirra_domainMeanProfiili( 'S_Rwiba', nimi = 'Ice particle mean diameter averaged testi', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, binidata = True, xAxisL = r'[${\mu}m$]', color = icevari, savePrefix = 'domainMeanProfiili_salsa' )

        piirra_domainMeanProfiili( 'P_rs',  nimi = 'Ice mixing ratio averaged',   muunnosKerroin=1000., ajanhetket = [6,8], useDN = False, profiili = True, xAxisL = r'[g/kg]', color = icevari, askChangeOfVariable = askChangeOfVariable )
        
        mdp.plot_suljetus(naytaPlotit)
        
        piirra_domainProfiili( 'P_Rws', muunnosKerroin = 2.e6, longName = 'Ice particle mean diameter '+r'[${\mu}m$]', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, spinup = spinup, profiili = True, colorBarTickValues = cbvalDiam, colorBarTickNames = cbvalDiamStr, variKartta = diamVari, askChangeOfVariable = askChangeOfVariable  )
        
        mdp.plot_suljetus(naytaPlotit)
        
        #piirra_domainProfiili( 'thil', longName = 'Ice-Liquid water potential temperature', useDN = False, transpose = True,  xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = plt.cm.Reds, spinup = spinup, profiili = True  )
        
        piirra_aikasarjasettii( muuttuja = 'rmH2Osn',     muunnosKerroin = 2.5e+06,  longName = 'Removal of water by sedimentation of ice',       ylabel = 'flux ' + r'[$W/m^2$]', ymin = 0.0,  savePrefix = 'depIce', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, askChangeOfVariable = askChangeOfVariable  )
        
        mdp.plot_suljetus(naytaPlotit)
        
    if int(lvl) == 0:
        piirra_aikasarjavertailusettii( muuttuja = 'lwp',    muunnosKerroin = 1000.0,  savePrefix = 'changeLWP',    omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel  )
        piirra_aikasarjavertailusettii( muuttuja = 'lwpiwp', muunnosKerroin = 1000.0,  savePrefix = 'changeLWPIWP', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel  )
        
        piirra_domainMeanScatter(muuttujaX = 'lwp_bar', muuttujaY = 'zc', muuttujaYY = 'zb', muunnosKerroinX = 1000., muunnosKerroinY = 1.0, muunnosKerroinYY = -1.0, nimi = None, ajanhetket = [6,8], xaxislabel = r'LWP [$g/m^2$]', yaxislabel = 'CLD depth [m]', savePrefix = None, minimiX =  0, maksimiX = 70, minimiY  = 0., maksimiY = 500 )
        
        piirra_domainMeanScatter(muuttujaX = 'zb', muuttujaY = 'zc', muuttujaYY = 'zb', muunnosKerroinX = 1., muunnosKerroinY = 1.0, muunnosKerroinYY = -1.0, nimi = None, ajanhetket = [6,8], xaxislabel = r'Z cld base [m]', yaxislabel = 'CLD depth [m]', savePrefix = None, minimiX =  400, maksimiX = 750, minimiY  = 0., maksimiY = 500 )
        
        piirra_domainMeanScatter(muuttujaX = 'zc', muuttujaY = 'zc', muuttujaYY = 'zb', muunnosKerroinX = 1., muunnosKerroinY = 1.0, muunnosKerroinYY = -1.0, nimi = None, ajanhetket = [6,8], xaxislabel = r'Z cld top [m]', yaxislabel = 'CLD depth [m]', savePrefix = None, minimiX =  700, maksimiX = 900, minimiY  = 0., maksimiY = 500 )
    
    if int(lvl) == -1:
        nic = np.arange(0, 1.41, 0.2)
        
        piirra_domainProfiili( 'P_Nia', muunnosKerroin = 1.e-3, longName = 'Ice number concentration ', useDN = True, transpose = True, colorBarTickValues = nic, colorBarTickNames = map(str,nic), xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup, profiili = True  ) #variKartta = profiiliVariICE, colorBarTickValues = cbvalICE
    

        nyp = mdp.read_NamelistValue( os.path.dirname(os.path.realpath(arguments[1]))+"/NAMELIST" ,var = 'nyp'  )-4
        
        for slaissi in xrange(0,nyp):
            piirra_domainProfiili( 'S_Ni', muunnosKerroin = 1.e-3, longName = 'Ice number concentration ', useDN = True, transpose = True, colorBarTickValues = nic, colorBarTickNames = map(str,nic), xlabels = xLabelsHours, ylabels = ylabels, xticks = ticksHours, yticks = korkeustikit, variKartta = profiiliVariICE, spinup = spinup, profiili = False, sliceYbeg = slaissi, sliceXend = slaissi+1, savePrefix = 'NicSlices'  ) 
    
    if int(lvl) == -2:
        
        piirra_kokojakauma( muuttujaR = 'S_Rwiba', muuttujaN = 'S_Niba', typename = 'Ice', korkeusH = 700, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False, askChangeOfVariable = askChangeOfVariable ) # , xCustSize = 2
        piirra_kokojakauma( muuttujaR = 'S_Rwiba', muuttujaN = 'S_Niba', typename = 'Ice', korkeusH = 400, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False, askChangeOfVariable = askChangeOfVariable ) # , xCustSize = 2
        piirra_kokojakauma( muuttujaR = 'S_Rwiba', muuttujaN = 'S_Niba', typename = 'Ice', korkeusH = 200, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False, askChangeOfVariable = askChangeOfVariable ) # , xCustSize = 2

    if int(lvl) == -3:
        piirra_domainMeanProfiili( 'P_Rwia', nimi = 'Ice particle mean diameter averaged', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, profiili = True, xAxisL = r'[${\mu}m$]', color = icevari, xmax=450 )           
        piirra_domainMeanProfiili( 'S_Rwiba', nimi = 'Ice particle mass mean diameter averaged testi', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, binidata = True, xAxisL = r'[${\mu}m$]', color = icevari, savePrefix = 'domainMassMeanProfiili_salsa' )
        
    if int(lvl) == -4:
        
        piirra_kokojakauma( muuttujaR = 'S_Rwsba', muuttujaN = 'S_Nsba', typename = 'Snow', korkeusH = 700, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False, askChangeOfVariable = askChangeOfVariable ) # , yCustSize = 2
        piirra_kokojakauma( muuttujaR = 'S_Rwsba', muuttujaN = 'S_Nsba', typename = 'Snow', korkeusH = 400, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False, askChangeOfVariable = askChangeOfVariable ) # , yCustSize = 2
        piirra_kokojakauma( muuttujaR = 'S_Rwsba', muuttujaN = 'S_Nsba', typename = 'Snow', korkeusH = 200, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False, askChangeOfVariable = askChangeOfVariable ) # , yCustSize = 2        
    if int(lvl) == -5:
        
        piirra_kokojakauma( muuttujaR = 'S_Rwaba', muuttujaN = 'S_Naba', typename = 'Aero', korkeusH = 1100, aikaT = 6,  yCustSize = 2, sekoita = False, interplo = False ) # , yCustSize = 2
        #piirra_kokojakauma( muuttujaR = 'S_Rwiba', muuttujaN = 'S_Niba', typename = 'Aero', korkeusH = 400, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False ) # , yCustSize = 2
        #piirra_kokojakauma( muuttujaR = 'S_Rwiba', muuttujaN = 'S_Niba', typename = 'Aero', korkeusH = 200, aikaT = 6, ymax = 11000, xmax = 1.6, yCustSize = 2, sekoita = False, interplo = False ) # , yCustSize = 2
    if int(lvl) == -6:
        piirra_aikasarjasettii( muuttuja = 'iwp_bar', muunnosKerroin = 1000.0, longName = 'Ice water path', ylabel = 'path ' + r'[$g/m^2$]', ymin = 0.0,  savePrefix = 'iwpTS', omaVari = False, xlabel = 'time [h]', spinup = spinup, piilotaOsaXlabel = piilotaOsaXlabel, legenda = True  )
    
    if int(lvl) == -7:
        piirra_domainMeanProfiili( 'S_Rwiba', nimi = 'Ice particle mass mean diameter averaged testi', muunnosKerroin=2.e6  ,   ajanhetket = [6,8], useDN = False, binidata = True, xAxisL = r'[${\mu}m$]', color = icevari, savePrefix = 'domainMassMeanProfiili_salsa', askChangeOfVariable = askChangeOfVariable )
toc = time.clock()
print toc - tic
########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if __name__ == "__main__":
    if piirra and naytaPlotit:
        mdp.plot_lopetus()