import ModDataPros as mdp
import sys
from itertools import cycle



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

###################################
#
#
# LWP - loop over given files 
#
#
###################################

for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    liqWP_Tdata  = mdp.read_Data( filenameTS[i], 'lwp_bar' )
    

    #mdp.print_shape( 'lwp', liqWP_Tdata )

    #mdp.print_shape( 'time', time_data )
    
    nimi = 'LWP ' + filenameTS[i]

    mdp.aikasarjaTulostus( liqWP_Tdata*1000.0, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'aika [s]', ynimi = 'LWP g/m^2' )

for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    rainWP_Tdata = mdp.read_Data( filenameTS[i], 'rwp_bar' )

    nimi = 'RWP ' + filenameTS[i]
    mdp.aikasarjaTulostus( rainWP_Tdata*1000.0, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'aika [s]', ynimi = 'RWP g/m^2' )




########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra:
    mdp.plot_lopetus()
