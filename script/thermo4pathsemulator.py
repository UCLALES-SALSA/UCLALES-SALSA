import ModDataPros as mdp
import sys
from itertools import cycle
import numpy as np


#######################
##### setting up    ###    
##### DO NOT CHANGE ###
#######################
filenameNC = []
filenamePS = []
filenameTS = []

for filebase in sys.argv[1:]:

    filenameNC.append( filebase + '.nc'    )
    filenamePS.append( filebase + '.ps.nc' )
    filenameTS.append( filebase + '.ts.nc' )





##############################
#### YOU CAN CHANGE THESE: ###
##############################

##############################
#
# drawing  
#
# parameter settings 
#
##############################

piirra = True

tulostus = False

tightXAxis = True

colorNRO = len(sys.argv)-1

if ( colorNRO > 6 ) :
    LEGEND = False
else:
    LEGEND = True

mdp.initializeColors(colorNRO)

maksimisateet   = np.zeros(colorNRO)

maksimiLatent   = np.zeros(colorNRO)

maksimiSensible = np.zeros(colorNRO)

###################################
#
#
# cfrac - loop over given files 
#
#
###################################
#print ' '
#print 'LWP'

###################################
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False
    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    cfrac_Tdata  = mdp.read_Data( filenameTS[i], 'cfrac' )

    if LEGEND:
        nimi = 'Cloud fraction ' + "".join(filenameTS[i].split(".")[0])[-2:]
    else:
        nimi = 'Cloud fraction'

    mdp.aikasarjaTulostus( cfrac_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'cloud fraction', tightXAxis=tightXAxis, LEGEND=LEGEND )


mdp.plot_setYlim( 0.0, 1.0, extendBelowZero = True)

###################################
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    #base_Tdata   = mdp.read_Data( filenameTS[i], 'zb' )
    top_Tdata    = mdp.read_Data( filenameTS[i], 'zc' )
    top_Tdata = top_Tdata - np.ones(len(top_Tdata))*top_Tdata[0]

    #print 'minimi BASE '  + str(np.min(base_Tdata))
    #print 'maksimi BASE ' + str(np.max(base_Tdata))
    #print 'minimi TOP '  + str(np.min(top_Tdata))
    #print 'maksimi TOP ' + str(np.max(top_Tdata))
    #print np.min(base
    if 'maksimiCLOUD' in locals():
        maksimiCLOUD = max( maksimiCLOUD,   np.max( top_Tdata ) ) 
    else:
        maksimiCLOUD = np.max( top_Tdata )
    
    if 'minimiCLOUD' in locals():
        minimiCLOUD = min( minimiCLOUD,  np.min(top_Tdata) ) 
    else:
        minimiCLOUD = np.min( top_Tdata )
        

    if LEGEND:
        nimi = 'Change of cloud top ' + "".join(filenameTS[i].split(".")[0])[-2:]
    else:
        nimi = 'Change of cloud top'
        
    #mdp.aikasarjaTulostus( base_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'height [m]' )
    mdp.aikasarjaTulostus(  top_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva,    nimi = nimi, xnimi = 'time [s]', ynimi = 'height [m]', tightXAxis=tightXAxis, LEGEND=LEGEND )


mdp.plot_setYlim( minimiCLOUD, maksimiCLOUD, extendBelowZero = True)


###################################
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False
    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    prcp_Tdata  = mdp.read_Data( filenameTS[i], 'prcp' )
    

    maksimisateet[i] = np.max(prcp_Tdata[30:])
    if 'maksimiPRCP' in locals():
        maksimiPRCP = max( maksimiPRCP,   np.max(prcp_Tdata) ) 
    else:
        maksimiPRCP = np.max(top_Tdata)

    if LEGEND:
        nimi = 'Surface precipitation ' + "".join(filenameTS[i].split(".")[0])[-2:]
    else:
        nimi = 'Surface precipitation'

    mdp.aikasarjaTulostus( prcp_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'precipitation W/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND  )


mdp.plot_setYlim( 0.0, 200., extendBelowZero = True)


#print ' '
#print 'SENSIBLE'
###################################
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data = mdp.read_Data( filenameTS[i], 'time'    )
    shf_Tdata = mdp.read_Data( filenameTS[i], 'shf_bar' )
    #print 'minimi '  + str(np.min(shf_Tdata))
    #print 'maksimi ' + str(np.max(shf_Tdata))
    maksimiSensible[i] = np.max(shf_Tdata)
    
    if 'maksimiSHF' in locals():
        maksimiSHF = max( maksimiSHF,  np.max(shf_Tdata) )
    else:
        maksimiSHF = np.max(shf_Tdata)
    
    if 'minimiSHF' in locals():
        minimiSHF = min( minimiSHF,  np.min(shf_Tdata) )
    else:
        minimiSHF = np.min(shf_Tdata)

    if LEGEND:
        nimi = 'Sensible heat flux ' + "".join(filenameTS[i].split(".")[0])[-2:]
    else:
        nimi = 'Sensible heat flux'
        
    mdp.aikasarjaTulostus( shf_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Sensible heat flux W/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND )
    
   
mdp.plot_setYlim( minimiSHF, maksimiSHF )    

#print ' '
#print 'LATENT'
###################################
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    lhf_Tdata = mdp.read_Data( filenameTS[i], 'lhf_bar' )
    #print 'minimi '  + str(np.min(lhf_Tdata))
    #print 'maksimi ' + str(np.max(lhf_Tdata))
    maksimiLatent[i] = np.max(lhf_Tdata)
    if 'maksimiLHF' in locals():
        maksimiLHF = max( maksimiLHF,  np.max(lhf_Tdata) )
    else:
        maksimiLHF = np.max(lhf_Tdata)
    
    if 'minimiLHF' in locals():
        minimiLHF = min( minimiLHF,  np.min(lhf_Tdata) )
    else:
        minimiLHF = np.min(lhf_Tdata)

    if LEGEND:
        nimi = 'Latent heat flux ' + "".join(filenameTS[i].split(".")[0])[-2:]
    else:
        nimi ='Latent heat flux'
        
    mdp.aikasarjaTulostus( lhf_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Latent heat flux W/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND )    

mdp.plot_setYlim( minimiLHF, maksimiLHF )


############################

mdp.plot_alustus()
mdp.plottaa(range( 1,len(sys.argv) ), maksimisateet, 'maximum precipitation after 2h', 'case', 'precipitation W/m^2', changeColor = True, markers=True)
mdp.plot_alustus()
mdp.plottaa(range( 1,len(sys.argv) ), maksimiSensible, 'maximum sensible heat', 'case', 'Sensible heat flux W/m^2', changeColor = True, markers=True)
mdp.plot_alustus()
mdp.plottaa(range( 1,len(sys.argv) ), maksimiLatent, 'maximum latent heat', 'case', 'Latent heat flux W/m^2', changeColor = True, markers=True)

########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra:
    mdp.plot_lopetus()
