"""
kernel.py

Code to carry out radiative feedback calculations using the
radiative kernel method.

This module is designed to work with xarray datasets

Call air temperature 'TA' instead of 'T', consistent with pyCESM
"""
from __future__ import division
from __future__ import print_function
import numpy as np
import xarray as xr
#  Use climlab thermodynamic routines to calculate saturation specific humidity
from climlab.utils.thermo import qsat

#  Here is a general purpose interpolation routine
#   for transforming between two different latitude - pressure level grids
#   WOULD BE GREAT TO MAKE THIS MORE ROBUST USING XARRAY LABELLED DIMENSIONS
def regrid(lat1, lev1, data, lat2, lev2):
    '''A general purpose interpolation routine for transforming
    between two different latitude - pressure level grids.'''
    from scipy.interpolate import griddata
    nlat = lat2.size
    nlev = lev2.size
    X1, Y1 = np.meshgrid(lat1, lev1)
    X2, Y2 = np.meshgrid(lat2, lev2)
    points_in = np.transpose(np.array([X1.flatten(), Y1.flatten()]))
    points_out = np.transpose(np.array([X2.flatten(), Y2.flatten()]))
    data_interp = griddata(points_in, data.flatten(),
                           points_out, method='linear').reshape(nlev, nlat)
    return data_interp


def toa_flux(ctrl, pert, kernels, cloud_mask_correction=0.):
    '''Compute changes in TOA flux using radiative kernels.
    Inputs are xarray dataset handles for model climatology files,
    and another dataset for the kernels.

    The sign convention is that the computed flux is POSITIVE DOWN
    (source of energy)
    thus the flux has the same sign as the corresponding feedback
    once normalized by surface temperature change.

    Currently implemented for CAM4 aquaplanet output.'''
    #   The kernels are not on the same grid as the CAM4 output
    lev_kernel = kernels.pfull
    dp_kernel = np.tile(np.transpose(np.atleast_2d(np.diff(kernels.phalf.data))),(1,90))

    # Flux is calculated by convolving the kernel with
    #  temperature / humidity response in the TROPOSPHERE

    #   define a mask from surface to 100 mb at equator,
    #    decreasing linearly to 300 mb at the poles
    #      (following Feldl and Roe 2013)
    maxlev = np.tile(100. + (200. * np.abs(kernels.lat)/90.), (lev_kernel.size, 1))
    p_kernel = np.tile(np.expand_dims(lev_kernel, axis=1), (1, kernels.lat.size))
    dp_masked = np.ma.masked_array(dp_kernel, np.where(p_kernel > maxlev, False, True))

    #  this custom function will compute vertical integrals
    #  taking advantage of xarray named coordinates
    def integral(field):
        return (field * dp_masked).sum(dim='pfull')

    #  Temperature anomalies at every point in the latitude - pressure plane
    DeltaT = pert.data_vars['TA'] - ctrl.data_vars['TA']
    #  Surface temperature anomalies
    DeltaTS = pert.data_vars['TS'] - ctrl.data_vars['TS']

    #  2D interpolation for the tropospheric temperature anomalies
    #   here I think I need to abandon the xarray dataset and just use numpy
    field = regrid(ctrl.lat.data, ctrl.lev.data, DeltaT.data,
                   kernels.lat.data, lev_kernel.data)
    DeltaT_interp = np.ma.masked_array(field, np.isnan(field))
    #  These use the "masked array" module (part of numpy)
    #   to make it easy to mask out missing data
    #  1D interpolation for the surface temperature anomalies
    field = np.interp(kernels.lat, DeltaTS.lat, DeltaTS)
    DeltaTS_interp = np.ma.masked_array(field, np.isnan(field))

    #  finite difference approximation to the slope d/dT (q_saturation)
    small = 0.01
    dqsatdT = (qsat(ctrl.data_vars['TA']+small, ctrl.lev) -
               qsat(ctrl.data_vars['TA']-small, ctrl.lev)) / (2*small)

    #  actual specific humidity anomalies
    DeltaQ = pert.Q - ctrl.Q
    #  relative humidity in control run (convert from percent to fraction)
    RH_ctrl = ctrl.RELHUM / 100.

    #  Equivalent temperature change
    #  (actual humidity change expressed as temperature change at fixed RH)
    DeltaTequiv = DeltaQ / (RH_ctrl * dqsatdT )
    #  Scaled by local surface temp. anomaly
    #DeltaTequiv_scaled = DeltaTequiv / DeltaTS
     #  But actually we are supposed to be using log(q)
    DeltaLogQ = np.log(pert.Q) - np.log(ctrl.Q)
    dlogqsatdT = (np.log(qsat(ctrl.data_vars['TA']+small, ctrl.lev)) -
               np.log(qsat(ctrl.data_vars['TA']-small, ctrl.lev))) / (2*small)
    DeltaTequiv_log = DeltaLogQ / (dlogqsatdT)
    #  Interpolated to kernel grid:
    field = regrid(ctrl.lat.data, ctrl.lev.data, DeltaTequiv_log.data,
                   kernels.lat.data, lev_kernel.data)
    DeltaTequiv_interp = np.ma.masked_array(field, np.isnan(field))

    #  create a new dictionary to hold all the feedbacks
    flux = {}
    #   Compute the feedbacks!
    flux['lw_q'] = integral(kernels['lw_q'] * DeltaTequiv_interp)
    #  shortwave water vapor
    flux['sw_q'] = integral(kernels['sw_q'] * DeltaTequiv_interp)
    #  longwave temperature  (Planck and lapse rate)
    flux['Planck'] =  (kernels['lw_ts'] +
                       integral(kernels['lw_t'])) * DeltaTS_interp
    flux['lapse'] = integral(kernels['lw_t'] * (DeltaT_interp-DeltaTS_interp))
    flux['lw_t'] = flux['Planck'] + flux['lapse']
    #  Add up the longwave feedbacks
    flux['lw_net'] = flux['lw_t'] + flux['lw_q']

    #  Compute cloud radiative forcing
        #   first, the all-sky TOA radiation fields
    DeltaOLR = pert.FLNT - ctrl.FLNT
    DeltaASR = pert.FSNT - ctrl.FSNT
        #  and the clear-sky diagnostics
    DeltaOLRclr = pert.FLNTC - ctrl.FLNTC
    DeltaASRclr = pert.FSNTC - ctrl.FSNTC
    #  Cloud radiative forcing is just the difference between these:
    DeltaOLRcld = DeltaOLR - DeltaOLRclr
    DeltaASRcld = DeltaASR - DeltaASRclr

    #  Soden et al. 2008 show how to compute the cloud feedback accurately
    #  accounting for cloud masking --  see their equation (25)

    # The cloud radiative forcing itself, interpolated
    flux['CRF_sw'] = np.interp( kernels.lat, ctrl.lat, DeltaASRcld )
    flux['CRF_lw'] = np.interp( kernels.lat, ctrl.lat, -DeltaOLRcld )
    flux['CRF_net'] = flux['CRF_sw'] + flux['CRF_lw']
    #  The corrections...
    #    water vapor
    flux['cld_q_lw'] = integral(kernels['cld_q_lw'] * DeltaTequiv_interp)
    flux['cld_q_sw'] = integral(kernels['cld_q_sw'] * DeltaTequiv_interp)
    flux['cld_q'] = flux['cld_q_lw'] + flux['cld_q_sw']
    #    temperature
    flux['cld_t'] = ( kernels['cld_ts']*DeltaTS_interp +
                    integral(kernels['cld_t'] * DeltaT_interp) )
    #  Corrected cloud feedbacks
    flux['cloud_lw'] = flux['CRF_lw'] + flux['cld_t'] + flux['cld_q_lw']
    flux['cloud_sw'] = flux['CRF_sw'] + flux['cld_q_sw']
    flux['cloud'] = flux['cloud_lw'] + flux['cloud_sw']
    #    Soden et al. (2008) argue that there should be an additional
    #     correction due to cloud masking of radiative forcing...
    #    This is only an issue for 2xCO2 forcing.
    #    Soden et al. suggest 0.69 W/m2/K globally
    flux['cloud'] += cloud_mask_correction

    #  Total net feedback
    flux['total'] = flux['cloud'] + flux['lw_net'] + flux['sw_q']

    # Alternate description using fixed relative humidity as state
    # variable following Held and Shell J. Clim. 2012

    # temperature kernel for fixed relative humidity is the sum of
    #  traditional temperature and traditional water vapor kernels
    kernels['lw_talt'] = kernels['lw_t'] + kernels['lw_q']
    flux['Planck_alt'] =  (kernels['lw_ts'] +
                           integral(kernels['lw_talt'])) * DeltaTS_interp
    flux['lapse_alt'] = integral(kernels['lw_talt'] *
                                 (DeltaT_interp-DeltaTS_interp))
    # Get RH feedback by subtracting original water vapor kernel times
    # atmospheric temperature response from traditional water vapor feedback.
    flux['RH'] = flux['lw_q'] - integral(kernels['lw_q'] * DeltaT_interp)

    #  package output into xarray datasets
    flux_dataset = xr.Dataset(data_vars=flux)
    return flux_dataset


def local_feedback(ctrl, pert, kernels, cloud_mask_correction=0.):
    '''Compute climate feedbacks with respect to local surface temperature
          using radiative kernels.
    Inputs are netCDF4 dataset handles for model climatology files.

    Currently implemented for CAM4 aquaplanet output.'''

    flux = toa_flux(ctrl, pert, kernels, cloud_mask_correction)
    #  Surface temperature anomalies
    DeltaTS = pert.data_vars['TS'] - ctrl.data_vars['TS']
    #  1D interpolation for the surface temperature anomalies
    field = np.interp(kernels.lat, DeltaTS.lat, DeltaTS)
    DeltaTS_interp = np.ma.masked_array(field, np.isnan(field))

    return flux / DeltaTS_interp


def global_feedback(ctrl, pert, kernels, cloud_mask_correction=0.):
    '''Compute climate feedbacks with respect to global mean surface temperature
         using radiative kernels.
    Inputs are netCDF4 dataset handles for model climatology files.

    Currently implemented for CAM4 aquaplanet output.'''

    flux = toa_flux(ctrl, pert, kernels, cloud_mask_correction)
    #  Surface temperature anomalies
    DeltaTS = pert.data_vars['TS'] - ctrl.data_vars['TS']
    DeltaTSglobal = np.average(DeltaTS, weights=np.cos(np.deg2rad(ctrl.lat)))

    return flux / DeltaTSglobal


def global_mean_feedback(ctrl, pert, kernels, cloud_mask_correction=0.):
    fb_global = global_feedback(ctrl, pert, kernels, cloud_mask_correction)
    coslat = np.cos(np.deg2rad(fb_global.lat))
    #fb_global_mean = {}
    #for name in fb_global.data_vars:
    #    fb_global_mean[name] = np.average(fb_global[name], weights=coslat)
    #return fb_global_mean
    return fb_global.apply(np.average, weights=coslat)
