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

colorNRO = len(sys.argv)-1

mdp.initializeColors(colorNRO)

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


    nimi = 'CFRAC ' + filenameTS[i][4:6]

    mdp.aikasarjaTulostus( cfrac_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'cloud fraction' )


mdp.plot_setYlim( 0.0, 1.0, extendBelowZero = False)

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
        


    nimi = 'Cloud top ' + filenameTS[i][4:6]
    #mdp.aikasarjaTulostus( base_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'height [m]' )
    mdp.aikasarjaTulostus(  top_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva,    nimi = nimi, xnimi = 'time [s]', ynimi = 'height [m]')


mdp.plot_setYlim( minimiCLOUD, maksimiCLOUD, extendBelowZero = True)

###################################
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False
    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    prcp_Tdata  = mdp.read_Data( filenameTS[i], 'prcp' )

    if 'maksimiPRCP' in locals():
        maksimiPRCP = max( maksimiPRCP,   np.max(prcp_Tdata) ) 
    else:
        maksimiPRCP = np.max(top_Tdata)

    nimi = 'Surface precipitation ' + filenameTS[i][4:6]

    mdp.aikasarjaTulostus( prcp_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'precipitation' )


mdp.plot_setYlim( 0.0, maksimiPRCP, extendBelowZero = False)

########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra:
    mdp.plot_lopetus()
