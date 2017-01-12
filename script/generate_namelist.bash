#!/bin/bash

# Exit on error
set -e

nxp=${nxp:-204}
nyp=${nyp:-204}
nzp=${nzp:-200}
deltax=${deltax:-50.}
deltay=${deltay:-50.}
deltaz=${deltaz:-20.}
nxpart=${nxpart:-.true.}
dzmax=${dzmax:-3500.}
dzrat=${dzrat:-1.0}
dtlong=${dtlong:-1.}
distim=${distim:-100.}
timmax=${timmax:-10800.}
Tspinup=${Tspinup:-3600.}
runtype=${runtype:-'"INITIAL"'}
level=${level:-4}
CCN=${CCN:-600.e6}
prndtl=${prndtl:--0.3333333}
filprf=${filprf:-"'emul'"}
hfilin=${hfilin:-"'emul.rst'"}
ssam_intvl=${ssam_intvl:-120.}
savg_intvl=${savg_intvl:-120.}
mcflg=${mcflg:-.FALSE.}
lbinanl=${lbinanl:-.true.}
frqanl=${frqanl:-3600.}
corflg=${corflg:-.false.}
itsflg=${itsflg:-1}
strtim=${strtim:-180.0}
sed_aero=${sed_aero:-.FALSE.}
sed_cloud=${sed_cloud:-.TRUE.}
sed_precp=${sed_precp:-.TRUE.}
sed_ice=${sed_ice:-.FALSE.}
sed_snow=${sed_snow:-.FALSE.}
iradtyp=${iradtyp:-3}                   # ! 1 = no radiation, only large-scale forcing, 3 = radiation + large-scale forcing 
case_name=${case_name:-"'default'"}     # ! Case-specific large-scale forcing: none = not used, 
                              # ! default = simple divergence forcing with specified div 
div=${div:-1.5e-6}              # ! Divergence for e.g. case_name = 'default'
sfc_albedo=${sfc_albedo:-0.05}
radsounding=${radsounding:-"'datafiles/dsrt.lay'"}
dthcon=${dthcon:-20.}
drtcon=${drtcon:-248.296548667015}
ubmin=${ubmin:--0.25}
zrough=${zrough:-0.01}
th00=${th00:-289.}
umean=${umean:-0.1}
vmean=${vmean:--0.1}


# salsa	
nlcoag=${nlcoag:-.TRUE.}       # Master coagulation switch
nlcgcc=${nlcgcc:-.TRUE.}    #   ! Self-collection of cloud droplets
nlcgpp=${nlcgpp:-.TRUE.}    #   ! Self-collection of rain drops
nlcgpc=${nlcgpc:-.TRUE.}    #   ! Rain collection of cloud droplets
nlcgaa=${nlcgaa:-.FALSE.}   #   ! Aerosol coagulation
nlcgca=${nlcgca:-.TRUE.}   #   ! Cloud collection of aerosols
nlcgpa=${nlcgpa:-.TRUE.}    #   ! Rain collection of aerosols
nlcgia=${nlcgia:-.TRUE.}    #   ! Ice collection of aerosols
nlcgic=${nlcgic:-.TRUE.}    #   ! Ice collection of cloud droplets
nlcgii=${nlcgii:-.TRUE.}    #   ! Self-collection of ice
nlcgip=${nlcgip:-.TRUE.}    #   ! Ice collection of rain drops
nlcgsa=${nlcgsa:-.TRUE.}    #   ! Snow collection of aerosols
nlcgsc=${nlcgsc:-.TRUE.}    #   ! Snow collection of cloud droplets
nlcgsi=${nlcgsi:-.TRUE.}    #   ! Snow collection of ice particles
nlcgsp=${nlcgsp:-.TRUE.}    #   ! Snow collection of rain drops
nlcgss=${nlcgss:-.TRUE.}    #   ! Self-collection of snow

nlcnd=${nlcnd:-.TRUE.} #  ! Master condensation switch
nlcndgas=${nlcndgas:-.FALSE.} # ! --Aerosol precursor gas codensation
nlcndh2oae=${nlcndh2oae:-.TRUE.} # ! --Condensation of water on aerosols (if FALSE, equilibrium assumed)
nlcndh2ocl=${nlcndh2ocl:-.TRUE.} # ! --Condensation of water on cloud droplets (and drizzle)
nlcndh2oic=${nlcndh2oic:-.TRUE.} # ! --Condensation of water on ice particles
nlauto=${nlauto:-.TRUE.} # ! Master autoconversion switch
nlautosnow=${nlautosnow:-.FALSE.} # ! Master snow autoconversion switch
nlactiv=${nlactiv:-.TRUE.}  #! Master cloud activation switch
nlactbase=${nlactbase:-.FALSE.} # ! --Switch for parameterized cloud base activation
nlactintst=${nlactintst:-.TRUE.}  #! --Switch for interstitial activation based on host model Smax

nlichom=${nlichom:-.FALSE.}  #   ! Switch for homogeneous ice nucleation
nlichet=${nlichet:-.FALSE.}  #   ! Switch for heterogeneous ice nucleation
nlicimmers=${nlicimmers:-.FALSE.}  # ! Switch for ice nucleation by immersion
nlicmelt=${nlicmelt:-.FALSE.}  #  ! Switch for ice'n' snow melting

rhlim=${rhlim:-1.2}       #  ! RH limit for SALSA during initialization and spinup

isdtyp=${isdtyp:-0}
nspec=${nspec:-1} # 3
listspec=${listspec:-"'SO4','','','','','',''"} # "'SO4','DU','OC','','','',''"
volDistA=${volDistA:-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}   
volDistB=${volDistB:-0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
nf2a=${nf2a:-1.0}

sigmag=${sigmag:-1.2, 1.7, 2.0, 2.0, 2.0, 2.0, 2.0}  #! Stdev for initial aerosol size distribution for isdtyp == 0 (uniform)  
dpg=${dpg:-0.022, 0.12, 0.2, 0.2, 0.2, 0.2, 0.2}   #! Mode mean diameters in micrometers
n=${n:-125., 46.2502241367474 , 0., 0., 0., 0., 0.} #        ! Mode number concentrations in #/cm^3

##########################
###                    ###
### write the namelist ###
###                    ###
##########################

cat > ${dir}NAMELIST <<EOF
 &version
  ver="v1.0.2"
 /

 &model
  nxp =   ${nxp}
  nyp =   ${nyp}
  nzp =   ${nzp}
  deltax = ${deltax}
  deltay = ${deltay}
  deltaz = ${deltaz}
  nxpart = ${nxpart}
  dzmax  = ${dzmax}
  dzrat  = ${dzrat}
  dtlong = ${dtlong}
  distim = ${distim}
  timmax = ${timmax}
  Tspinup = ${Tspinup}
  runtype = ${runtype}
  level = ${level}
  CCN = ${CCN}
  prndtl = ${prndtl}
  filprf = ${filprf}
  hfilin = ${hfilin}
  ssam_intvl = ${ssam_intvl}
  savg_intvl = ${savg_intvl}
  mcflg = ${mcflg}
  lbinanl = ${lbinanl}
  frqanl = ${frqanl}
  corflg = ${corflg}
  itsflg = ${itsflg}
  strtim = ${strtim}
  sed_aero = ${sed_aero}
  sed_cloud = ${sed_cloud}
  sed_precp = ${sed_precp}
  sed_ice = ${sed_ice}
  sed_snow = ${sed_snow}
  iradtyp = ${iradtyp}                ! 1 = no radiation, only large-scale forcing, 3 = radiation + large-scale forcing 
  case_name = ${case_name}            ! Case-specific large-scale forcing: none = not used, 
                                      ! default = simple divergence forcing with specified div 
  div = ${div}              ! Divergence for e.g. case_name = 'default'
  sfc_albedo = ${sfc_albedo}
  radsounding = ${radsounding}
  dthcon = ${dthcon}
  drtcon = ${drtcon}
  ubmin  = ${ubmin}
  zrough = ${zrough}
  th00 = ${th00}
  umean =  ${umean}
  vmean = ${vmean}
 /

 &salsa	
   nlcoag = ${nlcoag}       ! Master coagulation switch
   nlcgcc = ${nlcgcc}       ! Self-collection of cloud droplets
   nlcgpp = ${nlcgpp}       ! Self-collection of rain drops
   nlcgpc = ${nlcgpc}       ! Rain collection of cloud droplets
   nlcgaa = ${nlcgaa}      ! Aerosol coagulation
   nlcgca = ${nlcgca}       ! Cloud collection of aerosols
   nlcgpa = ${nlcgpa}       ! Rain collection of aerosols
   nlcgia = ${nlcgia}       ! Ice collection of aerosols
   nlcgic = ${nlcgic}       ! Ice collection of cloud droplets
   nlcgii = ${nlcgii}       ! Self-collection of ice
   nlcgip = ${nlcgip}       ! Ice collection of rain drops
   nlcgsa = ${nlcgsa}       ! Snow collection of aerosols
   nlcgsc = ${nlcgsc}       ! Snow collection of cloud droplets
   nlcgsi = ${nlcgsi}       ! Snow collection of ice particles
   nlcgsp = ${nlcgsp}       ! Snow collection of rain drops
   nlcgss = ${nlcgss}       ! Self-collection of snow

   nlcnd       = ${nlcnd}  ! Master condensation switch
   nlcndgas    = ${nlcndgas}  ! --Aerosol precursor gas codensation
   nlcndh2oae  = ${nlcndh2oae}  ! --Condensation of water on aerosols (if FALSE, equilibrium assumed)
   nlcndh2ocl  = ${nlcndh2ocl}  ! --Condensation of water on cloud droplets (and drizzle)
   nlcndh2oic  = ${nlcndh2oic}  ! --Condensation of water on ice particles
   nlauto      = ${nlauto}  ! Master autoconversion switch
   nlautosnow  = ${nlautosnow} ! Master snow autoconversion switch
   nlactiv     = ${nlactiv}  ! Master cloud activation switch
   nlactbase   = ${nlactbase}  ! --Switch for parameterized cloud base activation
   nlactintst  = ${nlactintst}  ! --Switch for interstitial activation based on host model Smax

   nlichom     = ${nlichom}     ! Switch for homogeneous ice nucleation
   nlichet     = ${nlichet}     ! Switch for heterogeneous ice nucleation
   nlicimmers  = ${nlicimmers}   ! Switch for ice nucleation by immersion
   nlicmelt    = ${nlicmelt}    ! Switch for ice'n' snow melting

   rhlim = ${rhlim}          ! RH limit for SALSA during initialization and spinup

   isdtyp = ${isdtyp}
   nspec = ${nspec}
   listspec = ${listspec}
   volDistA = ${volDistA}
   volDistB = ${volDistB}
   nf2a = ${nf2a}

   sigmag = ${sigmag}  ! Stdev for initial aerosol size distribution for isdtyp == 0 (uniform)  
   dpg = ${dpg}   ! Mode mean diameters in micrometers
   n = ${n}         ! Mode number concentrations in #/cm^3
 /

EOF
 
exit