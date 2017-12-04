'''
pyCESM:  Python translation of CESM 1.2.1 code
$CCSMROOT/models/atm/cam/src/physics/cam/wv_saturation.F90

Brian Rose
June 2015

Documentation from CESM code:
--------------------------------

module wv_saturation

!--------------------------------------------------------------------!
! Module Overview:                                                   !
!                                                                    !
! This module provides an interface to wv_sat_methods, providing     !
! saturation vapor pressure and related calculations to CAM.         !
!                                                                    !
! The original wv_saturation codes were introduced by J. J. Hack,    !
! February 1990. The code has been extensively rewritten since then, !
! including a total refactoring in Summer 2012.                      !
!                                                                    !
!--------------------------------------------------------------------!
! Methods:                                                           !
!                                                                    !
! Pure water/ice saturation vapor pressures are calculated on the    !
! fly, with the specific method determined by a runtime option.      !
! Mixed phase SVP is interpolated from the internal table, estbl,    !
! which is created during initialization.                            !
!                                                                    !
! The default method for calculating SVP is determined by a namelist !
! option, and used whenever svp_water/ice or qsat are called.        !
!                                                                    !
!--------------------------------------------------------------------!
'''
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import numpy as np
from pyCESM.physconst import epsilo, latvap, latice, rh2o, cpair, tmelt, h2otrip

# ! transition range from es over H2O to es over ice
ttrice = 20.00
# ! Precalculated because so frequently used.
omeps = 1.0 - epsilo

from pyCESM.wv_sat_methods import (wv_sat_svp_water, wv_sat_svp_ice,
                            wv_sat_svp_trans,
                            wv_sat_qsat_water, wv_sat_qsat_ice)
from pyCESM.wv_sat_methods import wv_sat_svp_to_qmmr as svp_to_qmmr
from pyCESM.wv_sat_methods import wv_sat_svp_to_qsat as svp_to_qsat

#! Table of saturation vapor pressure values (estbl) from tmin to
#! tmax+1 Kelvin, in one degree increments.  ttrice defines the
#! transition region, estbl contains a combination of ice & water
#! values.
tmin = 127.16
tmax = 375.16

#  ! Set coefficients for polynomial approximation of difference
#  ! between saturation vapor press over water and saturation pressure
#  ! over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial is
#  ! valid in the range -40 < t < 0 (degrees C).
pcf = np.array([ 5.04469588506e-01,
                -5.47288442819e+00,
                -3.67471858735e-01,
                -8.95963532403e-03,
                -7.78053686625e-05])

def svp_water(t, method='GoffGratch'):
    '''! Compute saturation vapor pressure over water
       real(r8), intent(in) :: t ! Temperature (K)
       real(r8) :: es            ! SVP (Pa)
    '''
    es = wv_sat_svp_water(t, method=method)
    return es


def svp_ice(t, method='GoffGratch'):
    '''! Compute saturation vapor pressure over ice
       real(r8), intent(in) :: t ! Temperature (K)
       real(r8) :: es            ! SVP (Pa)
    '''
    es = wv_sat_svp_ice(t, method=method)
    return es


def svp_trans(t, method='GoffGratch'):
    '''! Compute saturation vapor pressure with an ice-water transition
       real(r8), intent(in) :: t ! Temperature (K)
       real(r8) :: es            ! SVP (Pa)
    '''
    es = wv_sat_svp_trans(t, method=method)
    return es


def tq_enthalpy(t, q, hltalt):
    '''! Get enthalpy based only on temperature
       ! and specific humidity.
       real(r8), intent(in) :: t      ! Temperature
       real(r8), intent(in) :: q      ! Specific humidity
       real(r8), intent(in) :: hltalt ! Modified hlat for T derivatives
       real(r8) :: enthalpy
    '''
    enthalpy = cpair * t + hltalt * q
    return enthalpy

#!---------------------------------------------------------------------
#! LATENT HEAT OF VAPORIZATION CORRECTIONS
#!---------------------------------------------------------------------

def no_ip_hltalt(t):
    '''!------------------------------------------------------------------!
       ! Purpose:                                                         !
       !   Calculate latent heat of vaporization of pure liquid water at  !
       !   a given temperature.                                           !
       !------------------------------------------------------------------!'''
    hltalt = latvap * np.ones_like(t)
    #! Account for change of latvap with t above freezing where
    #! constant slope is given by -2369 j/(kg c) = cpv - cw
    hltalt = np.where(t >= tmelt, latvap - 2369.0*(t-tmelt), latvap)
    return hltalt


def calc_hltalt(t):
    '''!------------------------------------------------------------------!
       ! Purpose:                                                         !
       !   Calculate latent heat of vaporization of water at a given      !
       !   temperature, taking into account the ice phase if temperature  !
       !   is below freezing.                                             !
       !   Optional argument also calculates a term used to calculate     !
       !   d(es)/dT within the water-ice transition range.                !
       !------------------------------------------------------------------!'''
    # t       ! Temperature in Kelvin
    hltalt = no_ip_hltalt(t)  # ! Appropriately modified hlat
    tc = t - tmelt  # ! Temperature in degrees C
    # ! Weighting of hlat accounts for transition from water to ice.
    weight = np.where(tc >= -ttrice, -tc/ttrice, 1.)
    hltalt += np.where(t < tmelt, weight*latice, 0.)
    #   ! polynomial expression approximates difference between es
    #   ! over water and es over ice from 0 to -ttrice (C) (max of
    #   ! ttrice is 40): required for accurate estimate of es
    #   ! derivative in transition range from ice to water
    tterm = np.zeros_like(t)
    for term in np.flipud(pcf):
        tterm = term + tc*tterm
    tterm = tterm/ttrice
    tterm = np.where(tc >= -ttrice, tterm, 0.)
    tterm = np.where(t < tmelt, tterm, 0.)

    return hltalt, tterm


#! Temperature derivative outputs, for qsat_*
def deriv_outputs(t, p, es, qs, hltalt, tterm):
    '''
  ! Inputs
  real(r8), intent(in) :: t      ! Temperature
  real(r8), intent(in) :: p      ! Pressure
  real(r8), intent(in) :: es     ! Saturation vapor pressure
  real(r8), intent(in) :: qs     ! Saturation specific humidity
  real(r8), intent(in) :: hltalt ! Modified latent heat
  real(r8), intent(in) :: tterm  ! Extra term for d(es)/dT in
                                 ! transition region.

  ! Outputs
  real(r8), intent(out), optional :: gam      ! (hltalt/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt    ! (d(qs)/dt)
  '''
    desdt = hltalt*es/(rh2o*t*t) + tterm
    dqsdt = np.where(qs == 1., 0., qs*p*desdt/(es*(p-omeps*es)))
    gam   = dqsdt * (hltalt/cpair)
    return gam, dqsdt


#!---------------------------------------------------------------------
#! QSAT (SPECIFIC HUMIDITY) PROCEDURES
#!---------------------------------------------------------------------

def qmmr(t, p):
    '''!------------------------------------------------------------------!
       ! Purpose:                                                         !
       !     Provide saturation mass mixing ratio.                        !
       !                                                                  !
       ! Note that qmmr is a ratio over dry air, whereas qsat is a ratio  !
       ! over total mass. These differ mainly at low pressures, where     !
       ! SVP may exceed the actual pressure, and therefore qmmr can blow  !
       ! up.                                                              !
       !------------------------------------------------------------------!
  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass)
    '''
    es = wv_sat_svp_water(t)
    qm = svp_to_qmmr(es, p)
    #! Ensures returned es is consistent with limiters on qmmr.
    es = np.where(es>p, p, es)
    return qm, es


def qsat(t, p, method='GoffGratch'):
    '''
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Look up and return saturation vapor pressure from precomputed  !
  !   table, then calculate and return saturation specific humidity. !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q
    '''
    #  es = estblf(t)
    #  Brian: not using lookup table. Just calculate directly
    es = svp_trans(t, method=method)
    qs = svp_to_qsat(es, p)
    #! Ensures returned es is consistent with limiters on qs.
    es = np.where(es > p, p, es)
    #es = np.min(es, p)
    #! "generalized" analytic expression for t derivative of es
    #! accurate to within 1 percent for 173.16 < t < 373.16
    hltalt, tterm = calc_hltalt(t)
    enthalpy = tq_enthalpy(t, qs, hltalt)
    gam, dqsdt = deriv_outputs(t, p, es, qs, hltalt, tterm)
    return es, qs, gam, dqsdt, enthalpy


def qsat_water(t, p, method='GoffGratch'):
    '''
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!
  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q
    '''
    es, qs = wv_sat_qsat_water(t, p, method=method)
    # ! "generalized" analytic expression for t derivative of es
    # ! accurate to within 1 percent for 173.16 < t < 373.16
    hltalt = no_ip_hltalt(t)  # ! Appropriately modified hlat
    enthalpy = tq_enthalpy(t, qs, hltalt)
    # ! For pure water/ice transition term is 0.
    gam, dqsdt = deriv_outputs(t, p, es, qs, hltalt, 0.)
    return es, qs, gam, dqsdt, enthalpy


def qsat_ice(t, p, method='GoffGratch'):
    '''
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !   Optionally return various temperature derivatives or enthalpy  !
  !   at saturation.                                                 !
  !------------------------------------------------------------------!
  ! Inputs
  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: p    ! Pressure
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure
  real(r8), intent(out) :: qs  ! Saturation specific humidity

  real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
  real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
  real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q
    '''
    es, qs = wv_sat_qsat_ice(t, p, method=method)
    # ! For pure ice, just add latent heats.
    hltalt = latvap + latice
    enthalpy = tq_enthalpy(t, qs, hltalt)
    # ! For pure water/ice transition term is 0.
    gam, dqsdt = deriv_outputs(t, p, es, qs, hltalt, 0.)
    return es, qs, gam, dqsdt, enthalpy


# !---------------------------------------------------------------------
# ! FINDSP (WET BULB TEMPERATURE) PROCEDURES
# !---------------------------------------------------------------------

#  BRIAN: haven't implemented these yet

#subroutine findsp_vc(q, t, p, use_ice, tsp, qsp)
#
#  use cam_logfile,  only: iulog
#  use abortutils,   only: endrun
#  ! Wrapper for findsp which is 1D and handles the output status.
#  ! Changing findsp to elemental restricted debugging output.
#  ! If that output is needed again, it's preferable *not* to copy findsp,
#  ! but to change the existing version.
#
#  ! input arguments
#  real(r8), intent(in) :: q(:)        ! water vapor (kg/kg)
#  real(r8), intent(in) :: t(:)        ! temperature (K)
#  real(r8), intent(in) :: p(:)        ! pressure    (Pa)
#  logical,  intent(in) :: use_ice     ! flag to include ice phase in calculations
#
#  ! output arguments
#  real(r8), intent(out) :: tsp(:)     ! saturation temp (K)
#  real(r8), intent(out) :: qsp(:)     ! saturation mixing ratio (kg/kg)
#
#  integer :: status(size(q))   ! flag representing state of output
#                               ! 0 => Successful convergence
#                               ! 1 => No calculation done: pressure or specific
#                               !      humidity not within usable range
#                               ! 2 => Run failed to converge
#                               ! 4 => Temperature fell below minimum
#                               ! 8 => Enthalpy not conserved
#
#  integer :: n, i
#
#  n = size(q)
#
#  call findsp(q, t, p, use_ice, tsp, qsp, status)
#
#  ! Currently, only 2 and 8 seem to be treated as fatal errors.
#  do i = 1,n
#     if (status(i) == 2) then
#        write(iulog,*) ' findsp not converging at i = ', i
#        write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
#        write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
#        call endrun ('wv_saturation::FINDSP -- not converging')
#     else if (status(i) == 8) then
#        write(iulog,*) ' the enthalpy is not conserved at i = ', i
#        write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
#        write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
#        call endrun ('wv_saturation::FINDSP -- enthalpy is not conserved')
#     endif
#  end do
#
#end subroutine findsp_vc
#
#elemental subroutine findsp (q, t, p, use_ice, tsp, qsp, status)
#!-----------------------------------------------------------------------
#!
#! Purpose:
#!     find the wet bulb temperature for a given t and q
#!     in a longitude height section
#!     wet bulb temp is the temperature and spec humidity that is
#!     just saturated and has the same enthalpy
#!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
#!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
#!
#! Method:
#! a Newton method is used
#! first guess uses an algorithm provided by John Petch from the UKMO
#! we exclude points where the physical situation is unrealistic
#! e.g. where the temperature is outside the range of validity for the
#!      saturation vapor pressure, or where the water vapor pressure
#!      exceeds the ambient pressure, or the saturation specific humidity is
#!      unrealistic
#!
#! Author: P. Rasch
#!
#!-----------------------------------------------------------------------
#!
#!     input arguments
#!
#
#  real(r8), intent(in) :: q        ! water vapor (kg/kg)
#  real(r8), intent(in) :: t        ! temperature (K)
#  real(r8), intent(in) :: p        ! pressure    (Pa)
#  logical,  intent(in) :: use_ice  ! flag to include ice phase in calculations
#!
#! output arguments
#!
#  real(r8), intent(out) :: tsp      ! saturation temp (K)
#  real(r8), intent(out) :: qsp      ! saturation mixing ratio (kg/kg)
#  integer,  intent(out) :: status   ! flag representing state of output
#                                    ! 0 => Successful convergence
#                                    ! 1 => No calculation done: pressure or specific
#                                    !      humidity not within usable range
#                                    ! 2 => Run failed to converge
#                                    ! 4 => Temperature fell below minimum
#                                    ! 8 => Enthalpy not conserved
#!
#! local variables
#!
#  integer, parameter :: iter = 8    ! max number of times to iterate the calculation
#  integer :: l                      ! iterator
#
#  real(r8) es                   ! sat. vapor pressure
#  real(r8) gam                  ! change in sat spec. hum. wrt temperature (times hltalt/cpair)
#  real(r8) dgdt                 ! work variable
#  real(r8) g                    ! work variable
#  real(r8) hltalt               ! lat. heat. of vap.
#  real(r8) qs                   ! spec. hum. of water vapor
#
#! work variables
#  real(r8) t1, q1, dt, dq
#  real(r8) qvd
#  real(r8) r1b, c1, c2
#  real(r8), parameter :: dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
#  real(r8), parameter :: dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
#  real(r8) enin, enout
#
#  ! Saturation specific humidity at this temperature
#  if (use_ice) then
#     call qsat(t, p, es, qs)
#  else
#     call qsat_water(t, p, es, qs)
#  end if
#
#  ! make sure a meaningful calculation is possible
#  if (p <= 5._r8*es .or. qs <= 0._r8 .or. qs >= 0.5_r8 &
#       .or. t < tmin .or. t > tmax) then
#     status = 1
#     ! Keep initial parameters when conditions aren't suitable
#     tsp = t
#     qsp = q
#     enin = 1._r8
#     enout = 1._r8
#
#     return
#  end if
#
#  ! Prepare to iterate
#  status = 2
#
#  ! Get initial enthalpy
#  if (use_ice) then
#     call calc_hltalt(t,hltalt)
#  else
#     call no_ip_hltalt(t,hltalt)
#  end if
#  enin = tq_enthalpy(t, q, hltalt)
#
#  ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
#  c1 = hltalt*c3
#  c2 = (t + 36._r8)**2
#  r1b = c2/(c2 + c1*qs)
#  qvd = r1b * (q - qs)
#  tsp = t + ((hltalt/cpair)*qvd)
#
#  ! Generate qsp, gam, and enout from tsp.
#  if (use_ice) then
#     call qsat(tsp, p, es, qsp, gam=gam, enthalpy=enout)
#  else
#     call qsat_water(tsp, p, es, qsp, gam=gam, enthalpy=enout)
#  end if
#
#  ! iterate on first guess
#  do l = 1, iter
#
#     g = enin - enout
#     dgdt = -cpair * (1 + gam)
#
#     ! New tsp
#     t1 = tsp - g/dgdt
#     dt = abs(t1 - tsp)/t1
#     tsp = t1
#
#     ! bail out if past end of temperature range
#     if ( tsp < tmin ) then
#        tsp = tmin
#        ! Get latent heat and set qsp to a value
#        ! that preserves enthalpy.
#        if (use_ice) then
#           call calc_hltalt(tsp,hltalt)
#        else
#           call no_ip_hltalt(tsp,hltalt)
#        end if
#        qsp = (enin - cpair*tsp)/hltalt
#        enout = tq_enthalpy(tsp, qsp, hltalt)
#        status = 4
#        exit
#     end if
#
#     ! Re-generate qsp, gam, and enout from new tsp.
#     if (use_ice) then
#        call qsat(tsp, p, es, q1, gam=gam, enthalpy=enout)
#     else
#        call qsat_water(tsp, p, es, q1, gam=gam, enthalpy=enout)
#     end if
#     dq = abs(q1 - qsp)/max(q1,1.e-12_r8)
#     qsp = q1
#
#     ! if converged at this point, exclude it from more iterations
#     if (dt < dttol .and. dq < dqtol) then
#        status = 0
#        exit
#     endif
#  end do
#
#  ! Test for enthalpy conservation
#  if (abs((enin-enout)/(enin+enout)) > 1.e-4_r8) status = 8
#
#end subroutine findsp
