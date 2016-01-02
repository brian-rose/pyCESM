'''
Custom diagnostics for CESM / CAM model output
This package is built on top of `xray` which provides the underlying
grid-aware data structures.

The method `open_dataset()` wraps the `xray.open_dataset()` method
and attemps to compute a bunch of useful diagnostics in addition to returning
a handle to the raw model output.
'''

import numpy as np
import xray
from scipy import integrate
#from climlab import thermo

import physconst
mb_to_Pa = 100.  # conversion factor from mb to Pa


def open_dataset(filename_or_ob, verbose=True, **kwargs):
    '''Convenience method to open and return an xray dataset handle,
    with precomputed CAM-specific diagnostic quantities.
    '''
    if verbose:
        print 'Opening dataset ', filename_or_ob
    dataset = xray.open_dataset(filename_or_ob, **kwargs)
    #  xray.Dataset has a property called .T which returns transpose()
    #  but this creates a conflict with data named 'T'
    #  Attempt to rename the T variable to TA
    if ('T' in dataset and 'TA' not in dataset):
        dataset = dataset.rename({'T':'TA'})
    if verbose:
        print 'Variable T renamed to TA.'
    fulldataset = compute_diagnostics(dataset)
    if verbose:
        print 'Gridpoint diagnostics have been computed.'
    try:
        fulldataset = compute_heat_transport_xray(fulldataset)
        if verbose:
            print 'Heat transport diagnostics have been computed.'
    except:
        if verbose:
            print 'Heat transport computation failed.'
        else:
            pass
    try:
        # overturning mass streamfunction (in 10^9 kg/s or "mass Sverdrup")
        V = fulldataset.V
        if 'lon' in V.dims:
            #  Take zonal average first
            V = V.mean(dim='lon')
        Psi = overturning(V)
        fulldataset = fulldataset.assign(Psi = Psi)
        if verbose:
            print 'Overturning streamfunction has been computed.'
    except:
        if verbose:
            print 'Overturning computation failed.'
        else:
            pass
    #  string data prevents us from doing arithmetic between two datasets
    #  So move any string data to attributes (where it belongs)
    #  Would be better to dynamically search for string data,
    #  but these are the two fields that usually appear in my CAM output
    string_fields = ['time_written', 'date_written']
    for item in string_fields:
        try:
            fulldataset.attrs[item] = fulldataset[item].values
            fulldataset = fulldataset.drop(item)
        except:
            pass

    return fulldataset


def inferred_heat_transport(energy_in, lat=None, latax=None):
    '''Compute heat transport as integral of local energy imbalance.
    Required input:
        energy_in: energy imbalance in W/m2, positive in to domain
    As either numpy array or xray.DataArray
    If using plain numpy, need to supply these arguments:
        lat: latitude in degrees
        latax: axis number corresponding to latitude in the data
            (axis over which to integrate)
    returns the heat transport in PW.
    Will attempt to return data in xray.DataArray if possible.
    '''
    if lat is None:
        try: lat = energy_in.lat
        except:
            raise InputError('Need to supply latitude array if input data is not self-describing.')
    lat_rad = np.deg2rad(lat)
    coslat = np.cos(lat_rad)
    field = coslat*energy_in
    if latax is None:
        try: latax = field.get_axis_num('lat')
        except:
            raise ValueError('Need to supply axis number for integral over latitude.')
    #  result as plain numpy array
    integral = integrate.cumtrapz(field, x=lat_rad, initial=0., axis=latax)
    result = (1E-15 * 2 * np.math.pi * physconst.rearth**2 * integral)
    if isinstance(field, xray.DataArray):
        result_xray = field.copy()
        result_xray.values = result
        return result_xray
    else:
        return result


def overturning(V, lat=None, lev=None, levax=0):
    '''compute overturning mass streamfunction (SV)
    Required input:
        V: meridional velocity (m/s) (zonal average)
    As either numpy array or xray.DataArray
    If using plain numpy, need to supply these arguments:
        lat: latitude in degrees
        lev: pressure levels in mb or hPa
        levax: axis number corresponding to pressure in the data
            (axis over which to integrate)
        levax argument is ignored if V is self-describing xray.DataArray
    Returns the overturning streamfunction in SV or 10^9 kg/s.
    Will attempt to return data in xray.DataArray if possible.
    '''
    if lat is None:
        try: lat = V.lat
        except:
            raise ValueError('Need to supply latitude array if input data is not self-describing.')
    if lev is None:
        try:
            lev = V.lev
        except:
            raise ValueError('Need to supply pressure array if input data is not self-describing.')
    lat_rad = np.deg2rad(lat)
    coslat = np.cos(lat_rad)
    field = coslat*V
    try: levax = field.get_axis_num('lev')
    except: pass
    #  result as plain numpy array
    result = (2*np.pi*physconst.rearth/physconst.gravit *
            integrate.cumtrapz(field, lev*mb_to_Pa, axis=levax, initial=0)*1E-9)
    if isinstance(field, xray.DataArray):
        result_xray = field.copy()
        result_xray.values = result
        return result_xray
    else:
        return result

#  SHOULD ALSO SET UNITS FOR EACH FIELD IN METADATA

def compute_diagnostics(run):
    '''Compute a bunch of additional diagnostics from regular CAM ouput.

    Input is an xray dataset containing a CAM simulation climatology.

    These diagnostics are all computed point by point and so should work
    on any grid.'''

    sinlat = np.sin(np.deg2rad(run.lat))
    coslat = np.cos(np.deg2rad(run.lat))
    SST = run.TS - physconst.tmelt
    #TS_global = global_mean(run.TS,run.lat)
    #SST_global = run.TS_global - const.tempCtoK
    # TOA radiation fluxes
    SWdown_toa = run.SOLIN
    OLR = run.FLNT
    ASR = run.FSNT
    OLRclr = run.FLNTC
    ASRclr = run.FSNTC
    Rtoa = ASR - OLR  # net downwelling radiation
    Rtoaclr = ASRclr - OLRclr
    ASRcld = ASR - ASRclr
    OLRcld = OLR - OLRclr
    Rtoacld = Rtoa - Rtoaclr
    SWup_toa = SWdown_toa - ASR
    ALBtoa = SWup_toa / SWdown_toa

    #  surface energy budget terms, all defined as POSITIVE UP
    #    (from ocean to atmosphere)
    LHF = run.LHFLX
    SHF = run.SHFLX
    LWsfc = run.FLNS
    LWsfc_clr = run.FLNSC
    SWsfc = -run.FSNS
    SWsfc_clr = -run.FSNSC
    SnowFlux =  ((run.PRECSC + run.PRECSL) *
                      physconst.rhoh2o * physconst.latice)
    #  we'll let the SW down term be positive down
    SWdown_sfc = run.FSDS
    SWdown_sfc_clr = run.FSDSC
    SWdown_sfc_cld = SWdown_sfc - SWdown_sfc_clr
    SWup_sfc = SWsfc + SWdown_sfc
    ALBsfc = SWup_sfc / SWdown_sfc
    # all the net surface radiation terms are defined positive up
    LWsfc_cld = LWsfc - LWsfc_clr
    SWsfc_cld = SWsfc - SWsfc_clr
    # net upward radiation from surface
    SurfaceRadiation = LWsfc + SWsfc
    SurfaceRadiation_clr = LWsfc_clr + SWsfc_clr
    SurfaceRadiation_cld = SurfaceRadiation - SurfaceRadiation_clr
    # net upward surface heat flux
    SurfaceHeatFlux = SurfaceRadiation + LHF + SHF + SnowFlux
    # net heat flux into atmosphere
    Fatmin = Rtoa + SurfaceHeatFlux

    #  hydrological cycle, all terms in  kg/m2/s or mm/s
    Evap = run.QFLX
    Precip = (run.PRECC + run.PRECL) * physconst.rhoh2o
    EminusP = Evap - Precip

    newfields = {
        'sinlat': sinlat, 'coslat': coslat,
        'SST': SST,
        'SWdown_toa': SWdown_toa, 'SWup_toa': SWup_toa, 'ALBtoa': ALBtoa,
        'OLR': OLR, 'OLRclr': OLRclr, 'OLRcld': OLRcld,
        'ASR': ASR, 'ASRclr': ASRclr, 'ASRcld': ASRcld,
        'Rtoa': Rtoa, 'Rtoaclr': Rtoaclr, 'Rtoacld': Rtoacld,
        'LHF': LHF, 'SHF': SHF,
        'LWsfc': LWsfc, 'LWsfc_clr': LWsfc_clr, 'LWsfc_cld': LWsfc_cld,
        'SWsfc': SWsfc, 'SWsfc_clr': SWsfc_clr, 'SWsfc_cld': SWsfc_cld,
        'SnowFlux': SnowFlux,
        'SWdown_sfc': SWdown_sfc, 'SWdown_sfc_clr': SWdown_sfc_clr,
            'SWdown_sfc_cld': SWdown_sfc_cld,
        'SWup_sfc': SWup_sfc, 'ALBsfc': ALBsfc,
        'SurfaceRadiation': SurfaceRadiation,
        'SurfaceRadiation_clr': SurfaceRadiation_clr,
        'SurfaceRadiation_cld': SurfaceRadiation_cld,
        'SurfaceHeatFlux': SurfaceHeatFlux, 'Fatmin': Fatmin,
        'Evap': Evap, 'Precip': Precip, 'EminusP': EminusP,
    }
    #  use the xray.Dataset.assign() method to create a new dataset
    #  with all the new fields added, and return it.
    return run.assign(**newfields)


def compute_heat_transport_xray(run):
    '''Compute heat transport diagnostics from regular CAM ouput.

    Input is an xray dataset containing a CAM simulation climatology.

    These diagnostics involve integration so probably only work on
    regular lat-lon grids (for now).'''

    HT = compute_heat_transport(run.Rtoa, run.SurfaceHeatFlux, run.EminusP)
    newfields = {
        'HT_total': HT['total'], 'HT_atm': HT['atm'], 'HT_ocean': HT['ocean'],
        'HT_latent': HT['latent'], 'HT_dse': HT['dse'],
    }
    #  use the xray.Dataset.assign() method to create a new dataset
    #  with all the new fields added, and return it.
    return run.assign(**newfields)


def compute_heat_transport(TOAfluxDown, SurfaceFluxUp, EminusP,
                           lat=None, latax=None):
    # net heat flux into atmosphere
    Fatmin = TOAfluxDown + SurfaceFluxUp

    # heat transport terms
    HT_total = inferred_heat_transport(TOAfluxDown, lat, latax)
    HT_atm = inferred_heat_transport(Fatmin, lat, latax)
    HT_ocean = inferred_heat_transport(-SurfaceFluxUp, lat, latax)
    # atm. latent heat transport from moisture imbal.
    HT_latent = inferred_heat_transport(EminusP * physconst.latvap, lat, latax)
    # dry static energy transport as residual
    HT_dse = HT_atm - HT_latent
    HT = {
        'total': HT_total, 'atm': HT_atm, 'ocean': HT_ocean,
        'latent': HT_latent, 'dse': HT_dse,
        }
    return HT

    #  still need to work on the rest!

#    run.dz = -physconst.rair * run.T / physconst.gravit * run.dp / run.p * 1.E-3  #  in km
#    run.dTdp_moistadiabat = thermo.pseudoadiabat(run.T,run.p)
#    run.dTdp,ignored = np.gradient(run.T) / run.dp
#    run.dTdp_moistanom = run.dTdp - run.dTdp_moistadiabat
#    run.dTdz = run.dTdp * run.dp / run.dz  # in K / km
#    run.dTdz_moistanom = run.dTdp_moistanom * run.dp / run.dz
#    # convert OMEGA (in Pa/s) to w (in m/s)
#    run.w = -run.omega*const.Rd/const.g*run.T/(run.p*const.mb_to_Pa)
#    # overturning mass streamfunction (in 10^9 kg/s or "mass Sverdrup")
#    run.Psi = overturning(run.V,run.lat,run.lev)
#    #  correct for mass imbalance....
#    run.V_imbal = np.trapz(run.V/run.PS, run.lev*100, axis=0)
#    run.V_bal = run.V - run.V_imbal
#    run.Psi_bal = overturning(run.V_bal,run.lat,run.lev)
#    ind700 = np.nonzero(np.abs(run.lev-700)==np.min(np.abs(run.lev-700)))  # closest vertical level to 700 mb
#    T700 = np.squeeze(run.T[ind700,:])
#    run.EIS = thermo.EIS(run.Ta,T700)
#    run.DSE = const.cp * run.T + const.g * run.Z
#    run.MSE = run.DSE + const.Lhvap * run.Q #  J / kg


#  NOTE would be better to use the xray .rename() method to change variable names
def convert_am2(run):
    '''Translate AM2 model output to the CAM naming conventions, so we can
    use the same diagnostic code.'''
    lev = run.pfull
    TS = run.t_surf
    T = run.temp
    Ta = run.t_ref # near-surface air temperature
    # TOA radiation
    SOLIN = run.swdn_toa
    FLNT = run.olr
    FLNTC = run.olr_clr
    FSNT = run.swdn_toa - run.swup_toa
    FSNTC = run.swdn_toa_clr - run.swup_toa_clr
    #  surface energy budget terms matching CAM sign and unit conventions
    LHFLX = run.evap * physconst.latvap
    SHFLX = run.shflx
    FLNS = -run.lwflx
    FLNSC = run.lwup_sfc_clr - run.lwdn_sfc_clr
    FSDS = run.swdn_sfc
    FSDSC = run.swdn_sfc_clr
    FSNS = run.swdn_sfc - run.swup_sfc
    FSNSC = run.swdn_sfc_clr - run.swup_sfc_clr
    #  snow flux
    PRECSC = run.snow_conv / physconst.rhoh2o
    PRECSL = run.snow_ls / physconst.rhoh2o
    # hydrological cycle
    QFLX = run.evap  # kg/m2/s or mm/s
    PRECC = run.prec_conv / physconst.rhoh2o  # m/s
    PRECL = run.prec_ls / physconst.rhoh2o
    # precipitable water in kg/m2
    TMQ = run.WVP
    # near-surface wind speed
    U10 = run.wind
    # Geopotential height
    Z = run.z_full
    #  moisture and clouds
    RELHUM = run.rh
    Q = run.sphum
    CLOUD = run.cld_amt_dyn
    #  surface pressure
    PS = run.ps
    #  velocity components
    U = run.ucomp
    V = run.vcomp

    newfields = {
    'lev': lev, 'TS': TS, 'T': T, 'Ta': Ta,
    'SOLIN': SOLIN, 'FLNT': FLNT, 'FLNTC': FLNTC, 'FSNT': FSNT, 'FSNTC': FSNTC,
    'LHFLX': LHFLX, 'SHFLX': SHFLX, 'FLNS': FLNS, 'FLNSC': FLNSC,
    'FSDS': FSDS, 'FSDSC': FSDSC, 'FSNS': FSNS, 'FSNSC': FSNSC,
    'PRECSC': PRECSC, 'PRECSL': PRECSL, 'PRECC': PRECC, 'PRECL': PRECL,
    'QFLX': QFLX, 'TMQ': TMQ, 'U10': U10, 'Z': Z,
    'RELHUM': RELHUM, 'Q': Q, 'CLOUD': CLOUD, 'PS': PS,
    'U': U, 'V': V,
    }
    #  use the xray.Dataset.assign() method to create a new dataset
    #  with all the new fields added, and return it.
    return run.assign(**newfields)
