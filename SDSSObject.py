import os
import numpy as np
import pyfits as pf
from scipy.integrate import quad
from scipy.optimize import minimize
from utils import chi2g, chi2D


class SDSSObject():
    '''
    SDSSObject
    ==========
    A class for SDSS BOSS or eBOSS objects. Note that each fiberid will also
    be separated for better serialization
    '''
    def __init__(self, plate, mjd, fid, dataVersion, baseDir="/SCRATCH"):
        # Data version check
        if dataVersion == "v5_7_0" or dataVersion == "v5_7_2":
            flag = "BOSS"
        elif dataVersion == "v5_10_0":
            flag = "eBOSS"
        else:
            raise ValueError("Undefined data version")
        # Construct dirs for loading data
        dir1 = os.path.join(baseDir, flag, "data", dataVersion, str(plate))
        suffix = str(plate) + "-" + str(mjd) + ".fits"
        spFile = os.path.join(dir1, "spPlate-" + suffix)
        zbFile = os.path.join(dir1, dataVersion, "spZbest-" + suffix)
        zlFile = os.path.join(dir1, dataVersion, "spZline-" + suffix)
        # Match fiberid with the index in the plate
        hdulist = pf.open(zbFile)
        fiberidList = hdulist[1].data.field('FIBERID')
        hdulist.close()
        findex = np.argwhere(fiberidList == fid)[0][0]
        self.fiberid = fid
        self.mjd = mjd
        self.plate = plate
        # Load data from spPlate
        hdulist = pf.open(spFile)
        self.c0 = hdulist[0].header['coeff0']
        self.c1 = hdulist[0].header['coeff1']
        self.npix = hdulist[0].header['naxis1']
        self.wave = 10. ** (self.c0 + self.c1 * np.arange(self.npix))
        self.flux = hdulist[0].data[findex]
        self.ivar = hdulist[1].data[findex]
        hdulist.close()
        # Load data from spZbest
        hdulist = pf.open(zbFile)
        self.vdisp = hdulist[1].data.field('VDISP')[findex]
        self.synflux = hdulist[2].data[findex]
        self.RA = hdulist[1].data.field('PLUG_RA')[findex]
        self.DEC = hdulist[1].data.field('PLUG_DEC')[findex]
        self.obj_id = hdulist[1].data.field('OBJID')[findex]
        self.obj_class = hdulist[1].data.field('CLASS')[findex]
        self.obj_type = hdulist[1].data.field('OBJTYPE')[findex]
        self.z = hdulist[1].data.field('Z')[findex]
        self.zwarning = hdulist[1].data.field('ZWARNING')[findex]
        self.z_err = hdulist[1].data.field('Z_ERR')[findex]
        self.spectroflux = hdulist[1].data.field('SPECTROFLUX')[findex]
        self.rchi2 = hdulist[1].data.field('RCHI2')[findex]
        hdulist.close()
        # Load data from spZline
        hdulist = pf.open(zlFile)
        zlineList = hdulist[1].data
        hdulist.close()
        self.zline = zlineList[np.where(zlineList['fiberid'] == self.fiberid)]
        # Additional processing
        self.reduced_flux = self.flux - self.synflux
        self.nMax = len(self.flux)

    def wave2bin(self, waveLength):
        '''
        SDSSObject.wave2bin(waveLength)
        ===============================
        Return the accordinf bin of given waveLength

        Paramters:
            waveLength: given wavelength to convert (float)
        Returns:
            b: the bin of the given wavelength (same type as waveLength)
        '''
        b = int((np.log10(waveLength) / self.c1) - self.c0 / self.c1)
        if b < 0:
            b = 0
        elif b > self.nMax:
            b = self.nMax
        return b

    def mask(self, lineList):
        '''
        SDSSObject.mask(lineList)
        =========================
        Mask part of the spectra by setting ivar = 0

        Parameters:
            lineList: 2-d list with [:, 0] as start and [:, 1] as end
        Returns:
            None
        '''
        for each in lineList:
            self.ivar[self.wave2bin(each[0]): self.wave2bin(each[1])] = 0

    def nearLine(self, x0, width=10.0):
        '''
        nearLine(mjd, plate, i, x0, obj, width=10.0)
        ============================================
        Whether given wavelength is near a masked or not

        Parameters:
            x0: given wavelength
            width: default 10.0, a parameter for nearness
        Returns:
            a boolean of whether near or not
        '''
        crit = []
        crit.append(abs(self.zline['linewave'] * (1 + self.z) - x0) < width)
        crit.append(self.zline['lineew'] / self.zline['lineew_err'] > 6.0)
        crit.append(self.zline['lineew_err'] > 0)
        match = np.logical_and.reduce(crit)
        if np.any(match):
            return True
        else:
            return False

    def singletFit(self, bounds, initParam, paramLim):
        '''
        SDSSObject.singleFit(bounds, initParam, paramLim)
        ====================================================
        Fit singlet

        Parameters:
            bounds: the range that will be fitted
            initParam: initial parameter
            paramLim: limit of parameter to search for
        Returns:
            res.x: parameter value after fitting
            res.fun: chisquare of the fitting
        '''
        res = minimize(chi2g, initParam, args=(self.wave[bounds],
                                               self.reduced_flux[bounds],
                                               self.ivar[bounds]),
                       method='SLSQP', bounds=paramLim)
        return res.x, res.fun

    def doubletFit(self, bounds, initParam, paramLim):
        '''
        SDSSObject.singleFit(bounds, initParam, paramLim)
        ====================================================
        Fit doublet

        Parameters:
            bounds: the range that will be fitted
            initParam: initial parameter
            paramLim: limit of parameter to search for
        Returns:
            res.x: parameter value after fitting
            res.fun: chisquare of the fitting
        '''
        res = minimize(chi2D, initParam, args=(self.wave[bounds],
                                               self.reduced_flux[bounds],
                                               self.ivar[bounds]),
                       method='SLSQP', bounds=paramLim)
        return res.x, res.fun

    def radEinstein(self, zs):
        '''
        SDSSObject.radEinstein(z1, z2)
        ==============================
        Estimated Einstein Radius from Single Isothermal Sphere (SIS) model

        Parameters:
            zs: source redshift
        Returns:
            None
        '''
        # Necessary parameter
        H0 = 72e3
        c = 299792458.0
        OmegaM = 0.258
        OmegaL = 0.742
        vdisp = self.vdisp * 1000.0

        # Function needed for numerical computation of angular diameter distance
        def x12(z):
            return 1.0 / np.sqrt((1.0 - OmegaM - OmegaL) * (1.0 + z) *
                                 (1.0 + z) + OmegaM * (1.0 + z) ** 3.0 + OmegaL)
        # compute ang. diam. distances
        Ds = ((c / H0) * quad(x12, 0.0, zs)[0]) / (1.0 + zs)
        Dls = ((c / H0) * quad(x12, self.z, zs)[0]) / (1.0 + zs)
        # return value in arcsec
        coeff = 3600.0 * (180.0 / np.pi)
        return coeff * 4.0 * np.pi * vdisp * vdisp * Dls / (c * c * Ds)
