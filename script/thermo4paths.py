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
# LWP - loop over given files 
#
#
###################################
#print ' '
#print 'LWP'
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False
    
    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    liqWP_Tdata  = np.multiply( mdp.read_Data( filenameTS[i], 'lwp_bar' ), 1000.0)
    #print 'minimi '  + str(np.min(liqWP_Tdata))
    #print 'maksimi ' + str(np.max(liqWP_Tdata))
    if 'maksimiLWP' in locals():
        maksimiLWP = max( maksimiLWP,  np.max(liqWP_Tdata) )
    else:
        maksimiLWP = np.max(liqWP_Tdata)
    
    #if 'minimiLWP' in locals():
        #minimiLWP = min( minimiLWP,  np.min(liqWP_Tdata) )
    #else:
        #minimiLWP = np.min(liqWP_Tdata)
        
    
    nimi = 'LWP ' + "_".join(filenameTS[i][14:-6].split("/")[::2]) # + filenameTS[i]

    mdp.aikasarjaTulostus( liqWP_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'LWP g/m^2')


mdp.plot_setYlim( 0.0, maksimiLWP, False)

#print ' '
#print 'RWP'
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    rainWP_Tdata = np.multiply( mdp.read_Data( filenameTS[i], 'rwp_bar' ), 1000.0)
    #print 'minimi '  + str(np.min(rainWP_Tdata))
    #print 'maksimi ' + str(np.max(rainWP_Tdata))
    if 'maksimiRWP' in locals():
        maksimiRWP = max( maksimiRWP,  np.max(rainWP_Tdata) )
    else:
        maksimiRWP = np.max(rainWP_Tdata)
    
    #if 'minimiRWP' in locals():
        #minimiRWP = min( minimiRWP,  np.min(rainWP_Tdata) )
    #else:
        #minimiRWP = np.min(rainWP_Tdata)
        
    nimi = 'RWP ' + "_".join(filenameTS[i][14:-6].split("/")[::2]) # + filenameTS[i]
    mdp.aikasarjaTulostus( rainWP_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'RWP g/m^2' )

mdp.plot_setYlim( 0.0, maksimiRWP, False)

#print ' '
#print 'CLOUD'
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    base_Tdata   = mdp.read_Data( filenameTS[i], 'zb' )
    top_Tdata    = mdp.read_Data( filenameTS[i], 'zc' )
    #print 'minimi BASE '  + str(np.min(base_Tdata))
    #print 'maksimi BASE ' + str(np.max(base_Tdata))
    #print 'minimi TOP '  + str(np.min(top_Tdata))
    #print 'maksimi TOP ' + str(np.max(top_Tdata))
    #print np.min(base
    if 'maksimiCLOUD' in locals():
        maksimiCLOUD = max( maksimiCLOUD,  max( np.max(base_Tdata), np.max(top_Tdata) ) )
    else:
        maksimiCLOUD = max( np.max(base_Tdata), np.max(top_Tdata) )
    
    #if 'minimiCLOUD' in locals():
        #minimiCLOUD = min( minimiCLOUD,  min( np.min(base_Tdata), np.min(top_Tdata) ) )
    #else:
        #minimiCLOUD = min( np.min(base_Tdata), np.min(top_Tdata) )
        


    nimi = 'Cloud base & top ' + "_".join(filenameTS[i][14:-6].split("/")[::2]) # + filenameTS[i]
    mdp.aikasarjaTulostus( base_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'height [m]' )
    mdp.aikasarjaTulostus(  top_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = False,    nimi = nimi, xnimi = 'time [s]', ynimi = 'height [m]', changeColor = False)


mdp.plot_setYlim( 0.0, maksimiCLOUD, False )

#print ' '
#print 'SENSIBLE'

for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data = mdp.read_Data( filenameTS[i], 'time'    )
    shf_Tdata = mdp.read_Data( filenameTS[i], 'shf_bar' )
    #print 'minimi '  + str(np.min(shf_Tdata))
    #print 'maksimi ' + str(np.max(shf_Tdata))
    
    if 'maksimiSHF' in locals():
        maksimiSHF = max( maksimiSHF,  np.max(shf_Tdata) )
    else:
        maksimiSHF = np.max(shf_Tdata)
    
    if 'minimiSHF' in locals():
        minimiSHF = min( minimiSHF,  np.min(shf_Tdata) )
    else:
        minimiSHF = np.min(shf_Tdata)

    nimi = 'SHflux ' + "_".join(filenameTS[i][14:-6].split("/")[::2]) # + filenameTS[i]
    mdp.aikasarjaTulostus( shf_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Sensible heat flux W/m^2' )
    
   
mdp.plot_setYlim( minimiSHF, maksimiSHF )    

#print ' '
#print 'LATENT'

for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    lhf_Tdata = mdp.read_Data( filenameTS[i], 'lhf_bar' )
    #print 'minimi '  + str(np.min(lhf_Tdata))
    #print 'maksimi ' + str(np.max(lhf_Tdata))
    if 'maksimiLHF' in locals():
        maksimiLHF = max( maksimiLHF,  np.max(lhf_Tdata) )
    else:
        maksimiLHF = np.max(lhf_Tdata)
    
    if 'minimiLHF' in locals():
        minimiLHF = min( minimiLHF,  np.min(lhf_Tdata) )
    else:
        minimiLHF = np.min(lhf_Tdata)

    nimi = 'LHflux ' + "_".join(filenameTS[i][14:-6].split("/")[::2]) # + filenameTS[i]
    mdp.aikasarjaTulostus( lhf_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Latent heat flux W/m^2' )    

mdp.plot_setYlim( minimiLHF, maksimiLHF )

########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra:
    mdp.plot_lopetus()
