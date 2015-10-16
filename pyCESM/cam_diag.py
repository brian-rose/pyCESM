# CESM analysis routines
import numpy as np
import xray
from scipy import integrate
#from climlab import thermo

import physconst

#  A re-write of the old CESM.py custom class for use with the xray package

def open_dataset(filename_or_ob, **kwargs):
    '''Convenience method to open and return an xray dataset handle,
    with precomputed CAM-specific diagnostic quantities.
    '''
    dataset = xray.open_dataset(filename_or_ob, **kwargs)
    fulldataset = compute_diagnostics(dataset)
    return fulldataset


#   NEED TO MAKE ALL THIS MORE ROBUST TO VARYING NUMBER OF DIMENSIONS
def inferred_heat_transport( energy_in, lat_deg ):
    '''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
    lat_rad = np.deg2rad( lat_deg )
    return (1E-15 * 2 * np.math.pi * physconst.rearth**2 * 
           integrate.cumtrapz(np.cos(lat_rad)*energy_in,x=lat_rad, initial=0.))

def inferred_heat_transport_xray(energy_in):
    lat_deg = energy_in.lat
    #  result as plain numpy array
    result = inferred_heat_transport(energy_in, lat_deg)
    result_xray = energy_in.copy()
    result_xray.values = result
    return result_xray

def overturning(V, lat_deg, p_mb):
    '''compute overturning mass streamfunction (SV) from meridional velocity V, latitude (degrees) and pressure levels (mb)'''
    coslat = np.cos(np.deg2rad( lat_deg ))
    return (2*np.pi*physconst.rearth/physconst.gravit * coslat *
            integrate.cumtrapz(V, p_mb*100, axis=0, initial=0)*1E-9)


#  SHOULD ALSO SET UNITS FOR EACH FIELD IN METADATA

def compute_diagnostics(run):
    '''Compute a bunch of additional diagnostics from regular CAM ouput.
    
    Input is an xray dataset containing a CAM simulation climatology'''
    
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
    
    # heat transport terms
    HT_total = inferred_heat_transport_xray(Rtoa)
    HT_atm = inferred_heat_transport_xray(Fatmin)
    HT_ocean = inferred_heat_transport_xray(-SurfaceHeatFlux)
    # atm. latent heat transport from moisture imbal.
    HT_latent = inferred_heat_transport_xray(EminusP * physconst.latvap)
    # dry static energy transport as residual
    HT_dse = HT_atm - HT_latent  

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
        'SWdown_sfc': SWdown_sfc, 'SWdown_sfc_clr': SWdown_sfc_clr,
            'SWdown_sfc_cld': SWdown_sfc_cld,
        'SWup_sfc': SWup_sfc, 'ALBsfc': ALBsfc,
        'SurfaceRadiation': SurfaceRadiation,
        'SurfaceRadiation_clr': SurfaceRadiation_clr,
        'SurfaceRadiation_cld': SurfaceRadiation_cld,
        'SurfaceHeatFlux': SurfaceHeatFlux, 'Fatmin': Fatmin,
        'Evap': Evap, 'Precip': Precip, 'EminusP': EminusP,
        'HT_total': HT_total, 'HT_atm': HT_atm, 'HT_ocean': HT_ocean,
        'HT_latent': HT_latent, 'HT_dse': HT_dse,
    }
    #  use the xray.Dataset.assign() method to create a new dataset
    #  with all the new fields added, and return it.
    return run.assign(**newfields)

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


def convert_gram(run):
    '''Translate GRaM model output to the CAM naming conventions, so we can 
    use the same diagnostic code.'''
    pass

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
    FLNSC = run.lwup_sfc_clr - run.lwdown_sfc_clr
    FSDS = run.swdown_sfc
    FSDSC = run.swdown_sfc_clr
    FSNS = run.swdown_sfc - run.swup_sfc
    FSNSC = run.swdown_sfc_clr - run.swup_sfc_clr
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
    RELHUM = run.relhum
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
