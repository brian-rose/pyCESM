'''
pyCESM:  Python translation of CESM 1.2.1 code
$CCSMROOT/models/atm/cam/src/physics/cam/wv_sat_methods.F90

Brian Rose
June 2015

Documentation from CESM code:
--------------------------------

module wv_sat_methods

! This portable module contains all CAM methods for estimating
! the saturation vapor pressure of water.
!
! wv_saturation provides CAM-specific interfaces and utilities
! based on these formulae.
!
! Typical usage of this module:
!
! Init:
! call wv_sat_methods_init(r8, tmelt, h2otrip, tboil, errstring)
!
! Get scheme index from a name string:
! scheme_idx = wv_sat_get_scheme_idx(scheme_name)
! if (.not. wv_sat_valid_idx(scheme_idx)) <throw some error>
!
! Get pressures:
! es = wv_sat_svp_water(t, scheme_idx)
! es = wv_sat_svp_ice(t, scheme_idx)
!
! Use ice/water transition range:
! es = wv_sat_svp_trice(t, ttrice, scheme_idx)
!
! Note that elemental functions cannot be pointed to, nor passed
! as arguments. If you need to do either, it is recommended to
! wrap the function so that it can be given an explicit (non-
! elemental) interface.
'''
import numpy as np
from wv_saturation import tmelt, h2otrip, tboil, ttrice, epsilo

omeps = 1. - epsilo

# !---------------------------------------------------------------------
# ! UTILITIES
# !---------------------------------------------------------------------

#! Get saturation specific humidity given pressure and SVP.
#! Specific humidity is limited to range 0-1.
def wv_sat_svp_to_qsat(es, p):
    #  ! If pressure is less than SVP, set qs to maximum of 1.
    qsat = np.where(p <= es, 1., epsilo*es / (p - omeps*es))
    return qsat


#! Get saturation "mass mixing ratio".
#! It is almost always preferable to use saturation
#! specific humidity rather than this mixing ratio,
#! which blows up at low pressure.
def wv_sat_svp_to_qmmr(es, p):
    # ! When this function is used in regions of very low
    # ! pressure, we set it to "huge" rather than using
    #  ! the small denominator.
    #if ( (p - es) < epsilon(1.0_r8)**2 ) then
    # qmmr = huge(1.0_r8)
    # else
     #qmmr = epsilo*es / (p - es)
    #end if
    #  BRIAN: didn't implement that switch
    return epsilo*es / (p - es)


def wv_sat_qsat_water(t, p, method='GoffGratch'):
    '''
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
    '''
    es = wv_sat_svp_water(t, method=method)
    qs = wv_sat_svp_to_qsat(es, p)
    # ! Ensures returned es is consistent with limiters on qs.
    es = np.where(es > p, p, es)
    return es, qs


def wv_sat_qsat_ice(t, p, method='GoffGratch'):
    '''
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
    '''
    es = wv_sat_svp_ice(t, method=method)
    qs = wv_sat_svp_to_qsat(es, p)
    # ! Ensures returned es is consistent with limiters on qs.
    es = np.where(es > p, p, es)
    return es, qs


def wv_sat_qsat_trans(t, p, method='GoffGratch'):
    '''
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index
    '''
    es = wv_sat_svp_trans(t, method=method)
    qs = wv_sat_svp_to_qsat(es, p)
    # ! Ensures returned es is consistent with limiters on qs.
    es = np.where(es > p, p, es)
    return es, qs

# !---------------------------------------------------------------------
# ! SVP INTERFACE FUNCTIONS
# !---------------------------------------------------------------------

def wv_sat_svp_water(t, method='GoffGratch'):
    if method is 'GoffGratch':
        es = GoffGratch_svp_water(t)
    elif method is 'MurphyKoop':
        es = MurphyKoop_svp_water(t)
    elif method is 'OldGoffGratch':
        es = OldGoffGratch_svp_water(t)
    elif method is 'Bolton':
        es = Bolton_svp_water(t)
    else:
        raise ValueError('method not recognized.')
    return es


def wv_sat_svp_ice(t, method='GottGratch'):
    if method is 'GoffGratch':
        es = GoffGratch_svp_ice(t)
    elif method is 'MurphyKoop':
        es = MurphyKoop_svp_ice(t)
    elif method is 'OldGoffGratch':
        es = OldGoffGratch_svp_ice(t)
    elif method is 'Bolton':
        es = Bolton_svp_water(t)
    else:
        raise ValueError('method not recognized.')
    return es


def wv_sat_svp_trans(t, method='GoffGratch'):
    eswater = np.where(t >= tmelt-ttrice, wv_sat_svp_water(t, method), 0.)
    esice = wv_sat_svp_ice(t, method)
    weight = np.where((tmelt - t) > ttrice, 1., (tmelt - t)/ttrice)
    estrans = weight*esice + (1.0 - weight)*eswater
    es = np.where(t < tmelt, estrans, eswater)
    return es  # in Pa

# ! ---------------------------------------------------------------------
# ! SVP METHODS
# !---------------------------------------------------------------------

#! Goff & Gratch (1946)
def GoffGratch_svp_water(t):
    '''Goff and Gratch (1946) formula for saturation vapor pressure
    OVER WATER.
    Input is temperature in K
    Output is saturation vapor pressure in Pa.

    Comments in CESM code say  ! uncertain below -70 C
    '''
    es = (10.**(-7.90298*(tboil/t - 1.) +
                5.02808 * np.log10(tboil/t) -
                1.3816E-7 * (10.**(11.344 * (1.-t/tboil)) - 1.) +
                8.1328E-3 * (10.**(-3.49149 * (tboil/t - 1.)) - 1.) +
                np.log10(1013.246)))*100.
    return es  # in Pa


def GoffGratch_svp_ice(t):
    '''Goff and Gratch (1946) formula for saturation vapor pressure
    OVER ICE.
    Input is temperature in K
    Output is saturation vapor pressure in Pa.

    Comments in CESM code say  ! good down to -100 C
    '''
    es = (10.**(-9.09718 * (h2otrip/t - 1.) -
               3.56654 * np.log10(h2otrip/t) +
               0.876793 * (1. - t/h2otrip) +
               np.log10(6.1071)))*100.
    return es

#
#! Murphy & Koop (2005)
#
def MurphyKoop_svp_water(t):
#  ! (good for 123 < T < 332 K)
    es = np.exp(54.842763 - (6763.22 / t) - (4.210 * np.log(t)) +
                (0.000367 * t) + (np.tanh(0.0415 * (t - 218.8)) *
                (53.878 - (1331.22 / t) - (9.44523 * np.log(t)) +
                 0.014025 * t)))
    return es  # in Pa


def MurphyKoop_svp_ice(t):
#  ! (good down to 110 K)
    es = np.exp(9.550426 - (5723.265 / t) + (3.53068 * np.log(t))
             - (0.00728332 * t))
    return es  # in Pa

#
#! Old CAM implementation, also labelled Goff & Gratch (1946)
#
#! The water formula differs only due to compiler-dependent order of
#! operations, so differences are roundoff level, usually 0.
#
#! The ice formula gives fairly close answers to the current
#! implementation, but has been rearranged, and uses the
#! 1 atm melting point of water as the triple point.
#! Differences are thus small but above roundoff.
#
#! A curious fact: although using the melting point of water was
#! probably a mistake, it mildly improves accuracy for ice svp,
#! since it compensates for a systematic error in Goff & Gratch.

def OldGoffGratch_svp_water(t):
    ps = 1013.246
    e1 = 11.344*(1.0 - t/tboil)
    e2 = -3.49149*(tboil/t - 1.0)
    f1 = -7.90298*(tboil/t - 1.0)
    f2 = 5.02808*np.log10(tboil/t)
    f3 = -1.3816*(10.0**e1 - 1.0)/10000000.0
    f4 = 8.1328*(10.0**e2 - 1.0)/1000.0
    f5 = np.log10(ps)
    f = f1 + f2 + f3 + f4 + f5
    es = (10.0**f)*100.
    return es  # in Pa

def OldGoffGratch_svp_ice(t):
    term1 = 2.01889049 / (tmelt/t)
    term2 = 3.56654 * np.log(tmelt/t)
    term3 = 20.947031 * (tmelt/t)
    es = 575.185606E10 * np.exp(-(term1 + term2 + term3))
    return es  # in Pa

#! Bolton (1980)
#! zm_conv deep convection scheme contained this SVP calculation.
#! It appears to be from D. Bolton, 1980, Monthly Weather Review.
#! Unlike the other schemes, no distinct ice formula is associated
#! with it. (However, a Bolton ice formula exists in CLUBB.)
#
#! The original formula used degrees C, but this function
#! takes Kelvin and internally converts.
#

def Bolton_svp_water(t):
    c1 = 611.2
    c2 = 17.67
    c3 = 243.5
    es = c1*np.exp((c2*(t - tmelt))/((t - tmelt)+c3))
    return es  # in Pa
