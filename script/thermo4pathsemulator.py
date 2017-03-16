import ModDataPros as mdp
import sys
from itertools import cycle
import numpy as np
from scipy import interpolate
from emulator_inputs import read_design
import matplotlib.pyplot as plt

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
        nimi = 'Cloud fraction ' + filenameNC[i].split("/")[-2]
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
        nimi = 'Change of cloud top ' + filenameNC[i].split("/")[-2]
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
        nimi = 'Surface precipitation ' + filenameNC[i].split("/")[-2]
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
        nimi = 'Sensible heat flux ' + filenameNC[i].split("/")[-2]
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
        nimi = 'Latent heat flux ' + filenameNC[i].split("/")[-2]
    else:
        nimi ='Latent heat flux'
        
    mdp.aikasarjaTulostus( lhf_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'Latent heat flux W/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND )    

mdp.plot_setYlim( minimiLHF, maksimiLHF )

#############################################


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
        
    
    nimi = 'LWP ' + filenameNC[i].split("/")[-2]

    mdp.aikasarjaTulostus( liqWP_Tdata, time_data,  tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'time [s]', ynimi = 'LWP g/m^2', tightXAxis=tightXAxis, LEGEND=LEGEND )


mdp.plot_setYlim( 0.0, maksimiLWP )

#############################################

#print ' '
#print 'q profile changes'
for i in xrange(len(sys.argv)-1):
    uusikuva = True if i == 0 else  False
    
    q_data      = mdp.read_Data( filenamePS[i], 'q' )   
    time_data   = mdp.read_Data( filenameTS[i], 'time'    )

    height_data = mdp.read_Data( filenameNC[i], 'zt' )
    
    ajanhetket = 2*3600. #7200.
    aikaP = np.argmin( np.abs(ajanhetket - time_data) )
    
    #q_difference = q_data[ aikaP, : ] / q_data[ 0, : ] -1
    filu ='/home/aholaj/mounttauskansiot/ibrixmount/DESIGN/corr_design.csv'
    case, q_inv, tpot_inv, lwc_max, tpot_pbl, pblh, num_pbl = read_design( filu )

    pbl_indeksi = int(filenameNC[i].split("/")[-2][-2:])-1
    dens = 20
    
    fracZ   = np.zeros( dens*(len(height_data)-1) +1 )
    normZ   = np.zeros( len(fracZ) )
    
    qSpline = np.zeros( ( 2, len(fracZ) ) )
    
    
    for k in xrange(len(fracZ)-1):
        h_indeksi = int(np.floor(k/dens) )
        normZ[k] =  ( height_data[ h_indeksi ] + np.mod(k,dens)*(height_data[ h_indeksi + 1] - height_data[ h_indeksi ])/float(dens) )
        fracZ[k] =  normZ[k] / pblh[pbl_indeksi]
   
    normZ[-1] = height_data[-1]
    fracZ[-1] = normZ[-1] / pblh[pbl_indeksi]
    #np.set_printoptions(formatter={'float': lambda x: "{0:0.1f}".format(x)})
    tck0 = interpolate.splrep( height_data, q_data[ 0,: ]     )
    tckT = interpolate.splrep(height_data, q_data[ aikaP,:]  )
    for k in xrange( np.shape(qSpline)[1] ):
        qSpline[ 0,k ]  = interpolate.splev( normZ[k], tck0 )
        qSpline[ 1,k ]  = interpolate.splev( normZ[k], tckT )
    
    q_difference = qSpline[ 1,:] / qSpline[ 0,:] -1
        

    #print 'minimi '  + str(np.min(liqWP_Tdata))
    #print 'maksimi ' + str(np.max(liqWP_Tdata))
    if 'maksimifracZ' in locals():
        maksimifracZ = max( maksimifracZ,  np.max(fracZ) )
    else:
        maksimifracZ = np.max(fracZ)
    
    if 'maksimiQ' in locals():
        maksimiQ = max( maksimiQ,  np.max(q_difference) )
    else:
        maksimiQ = np.max(q_difference)

    if 'minimiQ' in locals():
        minimiQ = max( minimiQ,  np.max(q_difference) )
    else:
        minimiQ = np.max(q_difference)
    colorMap = plt.cm.rainbow
    skal = pblh[ pbl_indeksi ] / np.max( pblh)
    color = colorMap(skal)
    #if 'minimiLWP' in locals():
        #minimiLWP = min( minimiLWP,  np.min(liqWP_Tdata) )
    #else:
        #minimiLWP = np.min(liqWP_Tdata)
        
    
    nimi = 'Total water mix. rat g/m^2 difference between 0h & '+str((ajanhetket/3600.)) + 'h ' + filenameNC[i].split("/")[-2]

    mdp.profiiliTulostus( q_difference, aikaPisteet = 0, korkeus = fracZ, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = r'$\frac{q_{t=2h}}{q_{t=0h}}-1$', ynimi = 'z/pblh', tightXAxis=tightXAxis, LEGEND=LEGEND, omavari = color )
    #mdp.profiiliTulostus( q_data[0,:],  aikaPisteet = 0, korkeus = height_data, tulostus = tulostus, piirra = piirra, uusikuva = uusikuva, nimi = nimi, xnimi = 'q [g/kg]', ynimi = 'height [m]', tightXAxis=tightXAxis, LEGEND=LEGEND )
    #mdp.profiiliTulostus( q_data[60,:], aikaPisteet = 0, korkeus = height_data, tulostus = tulostus, piirra = piirra, uusikuva = False,    nimi = nimi, xnimi = 'q [g/kg]', ynimi = 'height [m]', tightXAxis=tightXAxis, LEGEND=LEGEND )

mdp.plot_setYlim( 0.0, 1.01 )
plt.xlim( -0.25, 0.025)

#############################################

mdp.plot_alustus()
mdp.plottaa(range( 1,len(sys.argv) ), maksimisateet, 'maximum precipitation after 2h', 'case', 'precipitation W/m^2', changeColor = True, markers=True)
#mdp.plot_alustus()
#mdp.plottaa(range( 1,len(sys.argv) ), maksimiSensible, 'maximum sensible heat', 'case', 'Sensible heat flux W/m^2', changeColor = True, markers=True)
#mdp.plot_alustus()
#mdp.plottaa(range( 1,len(sys.argv) ), maksimiLatent, 'maximum latent heat', 'case', 'Latent heat flux W/m^2', changeColor = True, markers=True)

########################
### finishing up     ###
### DO NOT CHANGE    ###
########################
if piirra:
    mdp.plot_lopetus()
