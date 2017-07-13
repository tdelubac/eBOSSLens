import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
from utils import *


def DR12Q_extractor(path='Superset_DR12Q.fits'):
    hdulist = pf.open(path)
    z_vi_DR12Q = hdulist[1].data.field('Z_VI')
    z_PCA_DR12Q = hdulist[1].data.field('Z_PIPE')
    plate_DR12Q = hdulist[1].data.field('PLATE')
    mjd_DR12Q = hdulist[1].data.field('MJD')
    fiber_DR12Q = hdulist[1].data.field('FIBERID')
    return np.transpose(np.vstack((plate_DR12Q,mjd_DR12Q,fiber_DR12Q,z_PCA_DR12Q,z_vi_DR12Q)))


def mask_QSO(ivar,z, l_width,c0,c1,Nmax):
    # Masking Lya, NV, SiIV, CIV, etc...
    l_LyA = 1215.668
    l_NV = 1240
    l_SiIV= 1400.0
    l_CIV = 1549.0
    l_HeII = 1640.0
    l_CIII = 1909.0
    l_CII = 2326.0
    l_FeII_a = 2382.765
    l_FeII_b = 2600.173
    l_MgII = 2798.0
    l_NeV = 3426.0
    l_OII = 3727
    l_NeIII = 3869
    l_Hd = 4101
    l_Hg = 4340
    l_Hb = 4861
    l_OIII_a = 4959
    l_OIII_b = 5007
    l_Ha = 6562.81

    #Mask the above emission lines of QSO's
    ivar[wave2bin((1+z)*(l_LyA -2.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_LyA +2.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_NV -0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_NV +0.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_SiIV -1.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_SiIV +1.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_CIV -2*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_CIV +2*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_HeII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_HeII +1.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_CIII -l_width),c0,c1,Nmax):wave2bin((1+z)*(l_CIII +l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_CII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_CII +0.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_FeII_a -l_width),c0,c1,Nmax):wave2bin((1+z)*(l_FeII_a +l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_FeII_b -l_width),c0,c1,Nmax):wave2bin((1+z)*(l_FeII_b +l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_MgII -2.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_MgII +2.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_NeV -0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_NeV +0.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_OII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_OII +1.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_NeIII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_NeIII +0.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_Hd -0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_Hd +0.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_Hg - 1.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_Hg + 1.5*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_Hb - l_width),c0,c1,Nmax):wave2bin((1+z)*(l_Hb + l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_OIII_a -2*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_OIII_a +2*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_OIII_b -2*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_OIII_b +2*l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(l_Ha -2*l_width),c0,c1,Nmax):wave2bin((1+z)*(l_Ha +3*l_width),c0,c1,Nmax)] = 0
    #additional Lines added after 1st results
    # Ne IV
    ivar[wave2bin((1+z)*(2427 -l_width),c0,c1,Nmax):wave2bin((1+z)*(2427 +l_width),c0,c1,Nmax)] = 0
    # Ne V
    ivar[wave2bin((1+z)*(3350 - 0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(3350 + 0.5*l_width),c0,c1,Nmax)] = 0
    # NI
    ivar[wave2bin((1+z)*(5200 -l_width),c0,c1,Nmax):wave2bin((1+z)*(5200 +l_width),c0,c1,Nmax)] = 0
    # [Fe VII]
    ivar[wave2bin((1+z)*(5721 - 0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(5721 + 0.5*l_width),c0,c1,Nmax)] = 0
    #[Fe VII]
    ivar[wave2bin((1+z)*(6087 - 0.5*l_width),c0,c1,Nmax):wave2bin((1+z)*(6087 + 0.5*l_width),c0,c1,Nmax)] = 0
    # night sky

    #ivar[wave2bin((1+z)*(6376 - l_width),c0,c1,Nmax):wave2bin((1+z)*(6376 + l_width),c0,c1,Nmax)] = 0
    # night sky
    #ivar[wave2bin((1+z)*(6307 -l_width),c0,c1,Nmax):wave2bin((1+z)*(6307 +l_width),c0,c1,Nmax)] = 0

    # SII
    ivar[wave2bin((1+z)*(6734 - l_width),c0,c1,Nmax):wave2bin((1+z)*(6734 + l_width),c0,c1,Nmax)] = 0
    # SII
    ivar[wave2bin((1+z)*(6716 - l_width),c0,c1,Nmax):wave2bin((1+z)*(6716 + l_width),c0,c1,Nmax)] = 0

    # Fe ?
    ivar[wave2bin((1+z)*(5317 -l_width),c0,c1,Nmax):wave2bin((1+z)*(5317 +l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(5691 -l_width),c0,c1,Nmax):wave2bin((1+z)*(5691 +l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(6504 - l_width),c0,c1,Nmax):wave2bin((1+z)*(6504 + l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(4490 - l_width),c0,c1,Nmax):wave2bin((1+z)*(4490 + l_width),c0,c1,Nmax)] = 0
    ivar[wave2bin((1+z)*(5080 -l_width),c0,c1,Nmax):wave2bin((1+z)*(5080 +l_width),c0,c1,Nmax)] = 0

    return ivar


def QSO_compute_FWHM(sobj, l_width):
    # Constants
    H0 = 72e3
    c = 299792458
    parsec = 3.0857e16
    # Lines
    l_o3 = [4959, 5007]
    l_Hb = 4861
    if sobj.z < 1:
        # Mask OIII
        o3Wave = np.array([[each - 2.0 * l_width, each + 2.0 * l_width]
                           for each in l_o3])
        sobj.mask(o3Wave)
        HB_flux = flux[wave2bin(4612*(1+z),c0,c1,Nmax):wave2bin(5112*(1+z),c0,c1,Nmax)]
        HB_wave = wave[wave2bin(4612*(1+z),c0,c1,Nmax):wave2bin(5112*(1+z),c0,c1,Nmax)]
        HB_weight = np.sqrt(ivar[wave2bin(4612*(1+z),c0,c1,Nmax):wave2bin(5112*(1+z),c0,c1,Nmax)])
        ### fit a line to continuum on unmasked points
        line_coeff = np.polyfit(x = HB_wave, y = HB_flux, deg = 1, w=HB_weight)
        HB_flux_r = flux[wave2bin(4812*(1+z),c0,c1,Nmax):wave2bin(4912*(1+z),c0,c1,Nmax)]
        HB_wave_r = wave[wave2bin(4812*(1+z),c0,c1,Nmax):wave2bin(4912*(1+z),c0,c1,Nmax)]
        HB_weight_r = np.sqrt(ivar[wave2bin(4812*(1+z),c0,c1,Nmax):wave2bin(4912*(1+z),c0,c1,Nmax)])
        res =  minimize(chi2Lorenz,[4862*(1+z),10,30],args=(HB_wave_r, HB_flux_r-line_coeff[0]*HB_wave_r -line_coeff[1],HB_weight_r), method='SLSQP', bounds = [(4862*(1+z)-5,4862*(1+z)+5),(1,100),(1,10000)])

        params = res.x
        FWHM = (c/1000)*2*params[1]/((1+z)*l_Hb) # km s-1
        average_flux = np.mean(flux[wave2bin(5100-40,c0,c1,Nmax):wave2bin(5100+40,c0,c1,Nmax)])
        l_times_luminosity = 5100*(1e-17)*average_flux*4*np.pi*(100*parsec*1e6*(c/H0)*quad(x12,0.0,z)[0]*(1+z))**2
    elif 6.2>z>1.5:
        HB_wave = 0.0

        #CIV
        l_CIV = 1549.0
        CIV_flux = flux[wave2bin(1300*(1+z),c0,c1,Nmax):wave2bin(1800*(1+z),c0,c1,Nmax)]
        CIV_wave = wave[wave2bin(1300*(1+z),c0,c1,Nmax):wave2bin(1800*(1+z),c0,c1,Nmax)]
        CIV_weight = np.sqrt(ivar[wave2bin(1300*(1+z),c0,c1,Nmax):wave2bin(1800*(1+z),c0,c1,Nmax)])
        ### fit a line to continuum on unmasked points
        line_coeff = np.polyfit(x = CIV_wave, y = CIV_flux, deg = 1, w=CIV_weight)
        CIV_flux_r = flux[wave2bin(1500*(1+z),c0,c1,Nmax):wave2bin(1600*(1+z),c0,c1,Nmax)]
        CIV_wave_r = wave[wave2bin(1500*(1+z),c0,c1,Nmax):wave2bin(1600*(1+z),c0,c1,Nmax)]
        CIV_weight_r = np.sqrt(ivar[wave2bin(1500*(1+z),c0,c1,Nmax):wave2bin(1600*(1+z),c0,c1,Nmax)])
        res =  minimize(chi2Lorenz,[1549*(1+z),10,10],args=(CIV_wave_r, CIV_flux_r-line_coeff[0]*CIV_wave_r -line_coeff[1],CIV_weight_r), method='SLSQP', bounds = [(1549*(1+z)-5,1549*(1+z)+5),(1,100),(1,10000)])

        params = res.x
        average_flux = 1350*np.mean(flux[wave2bin(1350-40,c0,c1,Nmax):wave2bin(1350+40,c0,c1,Nmax)])
        FWHM =  (c/1000)*2*params[1]/((1+z)*l_CIV) #km s-1
        l_times_luminosity = 1350*(1e-17)*average_flux*4*np.pi*(100*parsec*1e6*(c/H0)*quad(x12,0.0,z)[0]*(1+z))**2
    else:
        HB_wave = 0.0
        FWHM = 0.0
        l_times_luminosity = 0.0
        params = [0,0,0,0,0,0]
        line_coeff = None

    return FWHM, l_times_luminosity, HB_wave, params, line_coeff
