import os
import numpy as np
import pyfits as pf
from scipy.optimize import minimize
from utils import chi2g, chi2D


class SDSSObject():
    '''
    SDSSObject
    ==========
    A class for SDSS BOSS or eBOSS objects
    '''
    def __init__(self, mjd, plate, dataVersion, baseDir="/SCRATCH"):
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
        # Load data from spPlate
        hdulist = pf.open(spFile)
        self.c0 = hdulist[0].header['coeff0']
        self.c1 = hdulist[0].header['coeff1']
        self.npix = hdulist[0].header['naxis1']
        self.wave = 10. ** (self.c0 + self.c1 * np.arange(self.npix))
        self.flux = hdulist[0].data
        self.ivar = hdulist[1].data
        hdulist.close()
        # Load data from spZbest
        hdulist = pf.open(zbFile)
        self.vdisp = hdulist[1].data.field('VDISP')
        self.synflux = hdulist[2].data
        self.fiberid = hdulist[1].data.field('FIBERID')
        self.RA = hdulist[1].data.field('PLUG_RA')
        self.DEC = hdulist[1].data.field('PLUG_DEC')
        self.obj_id = hdulist[1].data.field('OBJID')
        self.obj_class = hdulist[1].data.field('CLASS')
        self.obj_type = hdulist[1].data.field('OBJTYPE')
        self.z = hdulist[1].data.field('Z')
        self.zwarning = hdulist[1].data.field('ZWARNING')
        self.z_err = hdulist[1].data.field('Z_ERR')
        self.spectroflux = hdulist[1].data.field('SPECTROFLUX')
        self.rchi2 = hdulist[1].data.field('RCHI2')
        hdulist.close()
        # Load data from spZline
        hdulist = pf.open(zlFile)
        self.zline = hdulist[1].data
        hdulist.close()
        # Additional processing
        self.reduced_flux = self.flux - self.synflux
        self.nMax = len(self.flux[:, 0])

    def wave2bin(self, waveLength):
        '''
        SDSSObject.wave2bin(waveLength)
        ===============================
        Return the accordinf bin of given waveLength

        Paramters:
            waveLength: given wavelength to convert (float or numpy.array)
        Returns:
            b: the bin of the given wavelength (same type as waveLength)
        '''
        b = int((np.log10(waveLength) / self.c1) - self.c0 / self.c1)
        try:
            b[b < 0] = 0
            b[b > self.nMax] = self.nMax
        except Exception:
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
            maskRange = self.wave2bin(each)
            self.ivar[:, maskRange[0]: maskRange[1]] = 0

    def singletFit(self, i, bounds, initParam, paramLim):
        res = minimize(chi2g, initParam, args=(self.wave[bounds],
                                               self.reduced_flux[i, bounds],
                                               self.ivar[i, bounds]),
                       method='SLSQP', bounds=paramLim)
        return res.x, res.fun

    def doubletFit(self, i, bounds, initParam, paramLim):
        res = minimize(chi2D, initParam, args=(self.wave[bounds],
                                               self.reduced_flux[i, bounds],
                                               self.ivar[i, bounds]),
                       method='SLSQP', bounds=paramLim)
        return res.x, res.fun
