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
    timeNC_data    = mdp.read_Data( filenameNC[i], 'time'    )
    
        
    dn0_data     = mdp.read_Data( filenameNC, 'dn0'     )


    #mdp.print_shape( 'lwp', liqWP_Tdata )

    #mdp.print_shape( 'time', time_data )
    
    nimi = 'IWP ' + filenameTS[i]

    mdp.laske_path_aikasarjaXYZ( f_data,    dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, uusikuva = uusikuva, nimi = 'Ice water path')



for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameNC[i], 'time'    )
    
    dn0_data     = mdp.read_Data( filenameNC[i], 'dn0'     )

    
    nimi = 'LWP ' + filenameTS[i]

    mpd.laske_path_aikasarjaXYZ( l_data,    dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, uusikuva = uusikuva, nimi = 'Liquid water path' )


for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False

    time_data    = mdp.read_Data( filenameTS[i], 'time'    )
    liqWP_Tdata  = mdp.read_Data( filenameTS[i], 'lwp_bar' )
    
    nimi = 'LWP ts' + filenameTS[i]

    mdp.aikasarjaTulostus( liqWP_Tdata*1000.0, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'aika [s]', ynimi = 'LWP g/m^2' )


########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra:
    mdp.plot_lopetus()

















  
    
    
    #P_rv_data	 = read_Data( filenamePS, 'P_rv'    )
    
    liqWP_Tdata  = mdp.read_Data( filenameTS, 'lwp_bar'  )
    rainWP_Tdata = mdp.read_Data( filenameTS, 'rwp_bar' )
    
    
    #top_data     = mdp.read_Data( filenameTS, 'zc'      )
    #base_data    = mdp.read_Data( filenameTS, 'zb'      )
    
    
    korkeus = ( zm_data - zt_data )*2.0
    
    #depth   = top_data - base_data
    
    #mdp.aikasarjaTulostus( depth, nimi= 'cloud depth' )

    
    print time_data
    
    #area( xm_data, ym_data )
    
    
    #mdp.laske_path_aikasarjaZ( P_rv_data, dn0_data, korkeus, time_data, tulostus = True, piirra = piirra, nimi = 'Water vapor path')
    
    mdp.laske_WC( f_data, dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'Ice water content')
    mdp.laske_WC( l_data, dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'Liquid water content')
    
    
    mdp.laske_NumberConcentration( S_Nic_data, S_Nc_data,  dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'number concentration, ice in liquid clouds (mixed-phase)' )
    
    mdp.laske_NumberConcentration( S_Nic_data, S_Nic_data, dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'number concentration, ice in ice clouds' )
    
    mdp.laske_NumberConcentration( S_Nc_data, S_Nc_data,   dn0_data, time_data, tulostus = True, piirra = piirra, nimi = 'number concentration, liquid in liquid clouds' )
    

    for bini in xrange( np.shape(S_Rwiba_data)[1] ):
        mdp.laske_MeanDiameterInBin( S_Rwiba_data, bini, S_Nc_data, time_data, tulostus = True, piirra =piirra, nimi = 'Mean diameter of ice particles in liquid clouds in bin ' )
    
    for bini in xrange( np.shape(S_Rwcba_data)[1] ):
        mdp.laske_MeanDiameterInBin( S_Rwcba_data, bini, S_Nc_data, time_data, tulostus = True, piirra =piirra, nimi = 'Mean diameter of liquid particles in liquid clouds in bin ' )            
    
    for bini in xrange( np.shape(S_Rwiba_data)[1] ):
        mdp.vertical_sum_timeseries( 'Rwiba bini'+str(bini), S_Rwiba_data[ :, bini, :, :, : ], time_data, tulostus = True, piirra = piirra )
    
    
    
    for bini in xrange( np.shape(S_Rwcba_data)[1] ):
        mdp.vertical_sum_timeseries( 'Rwcba bini'+str(bini), S_Rwcba_data[ :, bini, :, :, : ], time_data, tulostus = True, piirra = piirra )
    
