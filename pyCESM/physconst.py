'''
Python translation of CESM 1.2.1 code
$CCSMROOT/models/atm/cam/src/control/physconst.F90

for the pyCESM package
Brian Rose
June 2015
'''
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from pyCESM.shr_const_mod import *

avogad      = SHR_CONST_AVOGAD     #! Avogadro's number (molecules/kmole)
boltz       = SHR_CONST_BOLTZ      #! Boltzman's constant (J/K/molecule)
cday        = SHR_CONST_CDAY       #! sec in calendar day ~ sec
cpair       = SHR_CONST_CPDAIR     #! specific heat of dry air (J/K/kg)
cpliq       = SHR_CONST_CPFW       #! specific heat of fresh h2o (J/K/kg)
karman      = SHR_CONST_KARMAN     #! Von Karman constant
latice      = SHR_CONST_LATICE     #  ! Latent heat of fusion (J/kg)
latvap      = SHR_CONST_LATVAP     #  ! Latent heat of vaporization (J/kg)
pi          = SHR_CONST_PI         #! 3.14...
pstd        = SHR_CONST_PSTD       #! Standard pressure (Pascals)
r_universal = SHR_CONST_RGAS       #! Universal gas constant (J/K/kmol)
rhoh2o      = SHR_CONST_RHOFW      #! Density of liquid water (STP)
spval       = SHR_CONST_SPVAL      #!special value
stebol      = SHR_CONST_STEBOL     #! Stefan-Boltzmann's constant (W/m^2/K^4)
h2otrip     = SHR_CONST_TKTRIP     #! Triple point temperature of water (K)
c0          = 2.99792458E8         #! Speed of light in a vacuum (m/s)
planck      = 6.6260755E-34        #! Planck's constant (J.s)
#   ! Molecular weights
mwco2       =  44.             #! molecular weight co2
mwn2o       =  44.             #! molecular weight n2o
mwch4       =  16.             #! molecular weight ch4
mwf11       = 136.             #! molecular weight cfc11
mwf12       = 120.             #! molecular weight cfc12
mwo3        =  48.             #! molecular weight O3
mwso2       =  64.
mwso4       =  96.
mwh2o2      =  34.
mwdms       =  62.
#   ! modifiable physical constants for aquaplanet
gravit       = SHR_CONST_G     # ! gravitational acceleration (m/s**2)
sday         = SHR_CONST_SDAY  # ! sec in siderial day ~ sec
mwh2o        = SHR_CONST_MWWV  # ! molecular weight h2o
cpwv         = SHR_CONST_CPWV  # ! specific heat of water vapor (J/K/kg)
mwdry        = SHR_CONST_MWDAIR# ! molecular weight dry air
rearth       = SHR_CONST_REARTH# ! radius of earth (m)
tmelt        = SHR_CONST_TKFRZ # ! Freezing point of water (K)
# !---------------  Variables below here are derived from those above -----------------------
rga          = 1./SHR_CONST_G                 #! reciprocal of gravit
ra           = 1./SHR_CONST_REARTH            #! reciprocal of earth radius
omega        = SHR_CONST_OMEGA                   #! earth rot ~ rad/sec
rh2o         = SHR_CONST_RWV                     #! Water vapor gas constant ~ J/K/kg
rair         = SHR_CONST_RDAIR   #! Dry air gas constant     ~ J/K/kg
epsilo       = SHR_CONST_MWWV/SHR_CONST_MWDAIR   #! ratio of h2o to dry air molecular weights
zvir         = SHR_CONST_ZVIR                    #! (rh2o/rair) - 1
cpvir        = SHR_CONST_CPVIR                   #! CPWV/CPDAIR - 1.0
rhodair      = SHR_CONST_RHODAIR                 #! density of dry air at STP  ~ kg/m^3
cappa        = (SHR_CONST_RGAS/SHR_CONST_MWDAIR)/SHR_CONST_CPDAIR  #! R/Cp
#ez           ! Coriolis expansion coeff -> omega/sqrt(0.375)
Cpd_on_Cpv   = SHR_CONST_CPDAIR/SHR_CONST_CPWV
