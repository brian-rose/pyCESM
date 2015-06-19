'''
Python translation of CESM 1.2.1 code
$CCSMROOT/models/csm_share/shr/shr_const_mod.F90

for the pyCESM package
Brian Rose
June 2015
'''
SHR_CONST_PI      = 3.14159265358979323846  #! pi
SHR_CONST_CDAY    = 86400.0      #! sec in calendar day ~ sec
SHR_CONST_SDAY    = 86164.0      #! sec in siderial day ~ sec
SHR_CONST_OMEGA   = 2.0*SHR_CONST_PI/SHR_CONST_SDAY #! earth rot ~ rad/sec
SHR_CONST_REARTH  = 6.37122e6    #! radius of earth ~ m
SHR_CONST_G       = 9.80616      #! acceleration of gravity ~ m/s^2

SHR_CONST_STEBOL  = 5.67e-8      #! Stefan-Boltzmann constant ~ W/m^2/K^4
SHR_CONST_BOLTZ   = 1.38065e-23  #! Boltzmann's constant ~ J/K/molecule
SHR_CONST_AVOGAD  = 6.02214e26   #! Avogadro's number ~ molecules/kmole
SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ       #! Universal gas constant ~ J/K/kmole
SHR_CONST_MWDAIR  = 28.966       #! molecular weight dry air ~ kg/kmole
SHR_CONST_MWWV    = 18.016       #! molecular weight water vapor
SHR_CONST_RDAIR   = SHR_CONST_RGAS/SHR_CONST_MWDAIR        #! Dry air gas constant     ~ J/K/kg
SHR_CONST_RWV     = SHR_CONST_RGAS/SHR_CONST_MWWV          #! Water vapor gas constant ~ J/K/kg
SHR_CONST_ZVIR    = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0    #! RWV/RDAIR - 1.0
SHR_CONST_KARMAN  = 0.4          #! Von Karman constant
SHR_CONST_PSTD    = 101325.0     #! standard pressure ~ pascals
SHR_CONST_PDB     = 0.0112372    #! ratio of 13C/12C in Pee Dee Belemnite (C isotope standard)
 
SHR_CONST_TKTRIP  = 273.16       #! triple point of fresh water        ~ K
SHR_CONST_TKFRZ   = 273.15       #! freezing T of fresh water          ~ K 
SHR_CONST_TKFRZSW = SHR_CONST_TKFRZ - 1.8 #! freezing T of salt water  ~ K
#! density of dry air at STP  ~ kg/m^3
SHR_CONST_RHODAIR = SHR_CONST_PSTD/(SHR_CONST_RDAIR*SHR_CONST_TKFRZ)  

SHR_CONST_RHOFW   = 1.000e3      #! density of fresh water     ~ kg/m^3
SHR_CONST_RHOSW   = 1.026e3      #! density of sea water       ~ kg/m^3
SHR_CONST_RHOICE  = 0.917e3      #! density of ice             ~ kg/m^3
SHR_CONST_CPDAIR  = 1.00464e3    #! specific heat of dry air   ~ J/kg/K
SHR_CONST_CPWV    = 1.810e3      #! specific heat of water vap ~ J/kg/K
SHR_CONST_CPVIR   = (SHR_CONST_CPWV/SHR_CONST_CPDAIR)-1.0 #! CPWV/CPDAIR - 1.0
SHR_CONST_CPFW    = 4.188e3      #! specific heat of fresh h2o ~ J/kg/K
SHR_CONST_CPSW    = 3.996e3      #! specific heat of sea h2o   ~ J/kg/K
SHR_CONST_CPICE   = 2.11727e3    #! specific heat of fresh ice ~ J/kg/K
SHR_CONST_LATICE  = 3.337e5      #! latent heat of fusion      ~ J/kg
SHR_CONST_LATVAP  = 2.501e6      #! latent heat of evaporation ~ J/kg
#! latent heat of sublimation ~ J/kg
SHR_CONST_LATSUB  = SHR_CONST_LATICE + SHR_CONST_LATVAP  

SHR_CONST_OCN_REF_SAL = 34.7     #! ocn ref salinity (psu)
SHR_CONST_ICE_REF_SAL =  4.0     #! ice ref salinity (psu)

SHR_CONST_SPVAL   = 1.0e30       #! special missing value
