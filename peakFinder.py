import itertools as it
import numpy as np
from utils import gauss, kernel


class peakCandidateGalGal():
    '''
    peakCandidateGalGal
    ===================
    A class for all possible candidates of peaks found in the galaxy-galaxy lensing case
    '''
    def __init__(self, x0, sn):
        self.sn = sn                            #  Original SN from Bolton 2004
        self.wavelength = x0                    # Wavelength, keep this for comp
        self.chi = 1000.0                       # Overall chisquare
        self.amp = np.array([0.0, 0.0])         # Overall amplitude
        self.var = 0.0                          # Overall var
        self.wav = x0                           # Overall wavelength
        self.isDoublet = False                  # Doublet or singlet
        # Singlet fit
        self.chiSinglet = 1000.0                # Original [1]
        self.ampSinglet = 0.0                   # Original [2]
        self.varSinglet = 0.0                   # Original [3]
        self.wavSinglet = x0                    # Original [0]
        # Doublet fit
        self.chiDoublet = 1000.0                # Original [5]
        self.ampDoublet = np.array([0.0, 0.0])  # Original [6], [7]
        self.varDoublet = 0.0                   # Original [8]
        self.wavDoublet = np.array([x0, x0])    # Original [9], [10]

    def setDoublet(self, flag):
        self.isDoublet = flag
        if flag:
            self.chi = self.chiDoublet
            self.amp = self.ampDoublet
            self.var = self.varDoublet
            self.wav = self.wavDoublet
        else:
            self.chi = self.chiSinglet
            self.amp = self.ampSinglet
            self.var = self.varSinglet
            self.wav = self.wavSinglet

    def update(self, cm):
        if self.chiDoublet >= cm * self.chiSinglet:
            self.setDoublet(False)
        elif cm * self.chiSinglet > self.chiDoublet:
            self.setDoublet(True)


class peakCandidateGalLAE(peakCandidateGalGal):
    '''
    peakCandidateGalLAE
    ===================
    A class for all possible candidates of peaks found in the galaxy-LAE lensing case
    '''
    def __init__(self, x0, sn):
        # Initialize as for classical ELG detection
        peakCandidateGalGal.__init__(x0,sn):
        # Holders for a skew-normal fit
        self.isSkew = False

        self.chiSkew = 0.0                              # Chi square for Skew fit                     
        self.ampSkew = 0.0                              # Amplitude A
        self.scaleSkew = 0.0                            # w
        self.locSkew = 0.0                              # epsilon
        self.skewSkew = 0.0                             # skewness alpha

        self.chiSkew_singlet = 0.0                      # Chi square for singlet Skew fit                     
        self.ampSkew_singlet = 0.0                      # Amplitude A
        self.scaleSkew_singlet = 0.0                    # w
        self.locSkew_singlet = 0.0                      # epsilon
        self.skewSkew_singlet = 0.0                     # skewness alpha

        self.chiSkew_doublet  = 0.0                     # Chi square for doublet Skew fit
        self.ampSkew_doublet = np.array([0.0,0.0])      # Amplitude A
        self.scaleSkew_doublet = np.array([0.0,0.0])    # w
        self.locSkew_doublet = np.array([0.0,0.0])      # epsilon
        self.skewSkew_doublet = np.array([0.0,0.0])     # skewness alpha

        self.eq_Width = 0.0								# LyA equivalent width
        self.flux = 0.0									# LyA flux
        peak.l_blue_10 = 0.0							# Wavelength at 10% peak flux on blue side
    	peak.l_red_10 = 0.0								# Wavelength at 10% peak flux on red side
    	peak.aLambda = 0.0								# A_lambda skewness indicator
    	peak.skewness = 0.0 							# 3rd order moment around peak

    def update(self, cm):
        if self.chiSkew_singlet*cm <= self.chiSkew_doublet:
            self.chiSkew = self.chiSkew_singlet                                           
            self.ampSkew = self.ampSkew_singlet             
            self.scaleSkew = self.scaleSkew_singlet              
            self.locSkew = self.locSkew_singlet
            self.skewSkew = self.skewSkew_singlet
        else:
            self.chiSkew = self.chiSkew_doublet                                     
            self.ampSkew = self.ampSkew_doublet          
            self.scaleSkew = self.scaleSkew_doublet    
            self.locSkew = self.locSkew_doublet
            self.skewSkew = self.skewSkew_doublet

        if self.chiSkew >= cm * self.chiSinglet:
            self.setSkew(False)
        elif cm * self.chiSinglet > self.chiSkew:
            self.setSkew(True)

    def setSkew(self, flag):
        self.isSkew= flag
        if flag:
            self.chi = self.chiSkew
        #else:
        #    self.setDoublet(self.chiDoublet>self.chiSinglet)


class peakCandidateQSOLAE():
    '''
    peakCandidateQSOLAE
    ===================
    A class for all possible candidates of peaks found in the QSO-LAE lensing case
    '''

    # TODO : Fill 
    def __init__(self, x0, sn):

        self.sn = sn                            # Original SN from Bolton 2004
        self.wavelength = x0                    # Peak wavelength
        self.reduced_sn = 0.0                   # Recomputed SN with QSO continuum 3order subtraction
        self.redshift = x0/1215.667 -1          # Redshift of the background emission

        self.eq_Width = 0.0						# LyA equivalent width
        self.flux = 0.0							# LyA flux
        peak.l_blue_10 = 0.0					# Wavelength at 10% peak flux on blue side
    	peak.l_red_10 = 0.0						# Wavelength at 10% peak flux on red side
    	peak.aLambda = 0.0						# A_lambda skewness indicator
    	peak.skewness = 0.0 					# 3rd order moment around peak

class peakCandidateQSOGal():
    '''
    peakCandidateQSOGal
    ===================
    A class for all possible candidates of peaks found in the QSO-ELG lensing case
    '''
    # TODO : Fill 
    def __init__(self, x0, sn):
        self.sn = sn                            # Original SN from Bolton 2004
        self.wavelength = x0                    # Peak wavelength
        self.reduced_sn = 0.0                   # Recomputed SN with QSO continuum 3order fit subtraction
        self.redshift = 0.0                     # Redshift of the background emission
        self.total_sn  = 0.0                    # Total SN of all redshifted ELG emissions (OII, OIII, Hb, Ha usually)

class peakCandidateJackpot():
    '''
    peakCandidateJackpot
    ===================
    A class for all possible candidates of peaks found in the Jackpot case (galaxy lensing 2 ELG at different z)
    '''
    def __init__(self, x0, sn):
    	self.wavelength = 0.0

        self.sn_1 = sn                     # SN 1
        self.wavelength_1 = x0             # Peak wavelength 1
        self.z_1 = 0.0                     # Redshift of the background emission 1
        self.total_sn_1  = 0.0             # Total SN of all redshifted ELG emissions (OII, OIII, Hb, Ha usually) 1

        self.sn_2 = 0.0                    # SN 2
        self.wavelength_2 = 0.0            # Peak wavelength 2
        self.z_2 = 0.0                     # Redshift of the background emission 2 
        self.total_sn_2  = 0.0             # Total SN of all redshifted ELG emissions (OII, OIII, Hb, Ha usually) 2



def combNear(pcList):
    '''
    combNear(pcList)
    ================
    Combine peakCandidate in the given list if too near

    Parameters:
        pcList: a 1-d numpy array of peakCandidate
    Returns:
        pcList: the list after combination
    '''
    k = 0
    while (k < len(pcList) - 1):
        if (abs(pcList[k].wavelength - pcList[k+1].wavelength) < 10.0):
            if pcList[k].sn < pcList[k+1].sn:
                pcList = np.delete(pcList, k)
                k = k-1
            else:
                pcList = np.delete(pcList, k+1)
                k = k-1
        k = k+1
    return pcList


def checkFore(pcList, emLines, z):
    '''
    checkFore(pcList, emLines, z)
    =============================
    Check hits are not from foreground galaxy or badly fitted QSO

    Parameters:
        pcList: a 1-d numpy array of peakCandidate
        emLines: a list with emission lines to check
        z: redshift of the object
    Returns:
        foreground_z: boolean for whether if false hits or not
    '''
    foreground_z = False
    compare = it.combinations(emLines, len(pcList))
    for group in compare:
        for k in range(len(pcList)):
            for j in range(k + 1, len(pcList)):
                crit1 = pcList[k].wavelength / group[k]
                crit2 = crit1 - pcList[j].wavelength / group[j]
                if (np.abs(crit2) < 0.01) and (0.0 < crit1 < (z - 0.05)):
                    foreground_z = True
                    break
            # Break the loop once found to save time
            if foreground_z:
                break
        if foreground_z:
            break
    return foreground_z


def bolEstm(obj, sig, width):
    '''
    peakFinder.bolEstm(obj, sig, width)
    =============================
    Performs a SN estimation at each wavelength following the max-likelihood
    approach of Bolton et al. (2012).  

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        sig: Width of the gaussian kernel
        width: Width of the convolutional window
    Returns:
        SN: The SN at each wavelength as an array. Beginning and end are filled
        with null values due to the convolution.
    '''
    NormGauss = gauss(np.linspace(-width * 0.5, width * 0.5, width), 0.0, 1.0,
                      sig ** 2.0)
    NormGauss = NormGauss / np.sum(NormGauss)
    Cj1 = np.array([np.sum(kernel(j + 0.5 * width, width, NormGauss,
                                  len(obj.wave)) * obj.reduced_flux * obj.ivar)
                    for j in range(int(len(obj.wave) - width))])
    Cj2 = np.array([np.sum(obj.ivar * kernel(j + 0.5 * width, width, NormGauss,
                                             len(obj.wave)) ** 2.0)
                    for j in range(int(len(obj.wave) - width))])
    SN = np.zeros(len(obj.wave))
    SN[int(width * 0.5): int(width * 0.5 + len(Cj1))] = Cj1 / np.sqrt(Cj2)
    return SN


def qsoContfit(obj, peak, searchLyA, window_width = 40):
    '''
    peakFinder.qsoContfit(obj, peak, searchLyA, window_width)
    =============================
    A function to fit the QSO continuum near a peak with a 3rd order polynomial 
    and subtract it. Add the new SN of the peak in its class.

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak: The inquired peak
        searchLyA: True if we search for background LAE, False for backgroung ELGs
        window_width: Half-width (Angstroms) on which the polynomial is fitted
    Returns:
        accept: True if the SN is still high after subtraction. False otherwise.
    '''
    window = np.linspace(obj.wave2bin(x0)-40,obj.wave2bin(x0)+40,81,dtype = n.int16)
    median_local = np.median(obj.reduced_flux[window])

    fit_QSO = np.poly1d(np.polyfit(x=obj.wave[window],y=obj.reduced_flux[window],deg=3, \
        w=(np.abs(obj.reduced_flux[window]-median_local)<5)*np.sqrt(obj.ivar[window])) )
    
    new_flux = obj.reduced_flux[window] - fit_QSO(obj.wave[window])

    cj1_new = np.sum(new_flux*kernel(int(len(window)/2),width,NormGauss,len(new_flux))*obj.ivar[window])
    cj2_new = np.sum(obj.ivar[window]*kernel(int(len(window)/2),width,NormGauss,len(window))**2)
    SN_fitted = cj1_new/np.sqrt(cj2_new)

    if searchLyA and SN_fitted < 6:
        return False #Reject
    elif searchLyA and SN_fitted > 6:
        peak.reduced_sn = SN_fitted
        obj.reduced_flux_QSO[window]=new_flux
    elif searchLyA == False and SN_fitted < 6:
        return False  # Reject
    elif searchLyA == False and SN_fitted > 6:
        peak.reduced_sn = SN_fitted
        obj.reduced_flux_QSO[window]=new_flux
    return True  # Accept


def backgroundELG(obj, peak, em_lines):
    '''
    peakFinder.backgroundELG(obj, peak, em_lines)
    =============================
    Determines if a SN peak in a spectra can be part of a group of ELG 
    emission lines at a certain redshift (greater than the foreground object redshift)

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak: The inquired peak
        em_lines: The set of emission tracing the background ELG 
    Returns:
        - Nothing. Updates peak attributes by recording the best emission match.
    '''

    for l in em_lines:
        test_z = peak.wavelength/l - 1.0
        if test_z > obj.z:
            quad_SN = 0.0
            for w in em_lines:
                center_bin = obj.wave2bin(w*(1+test_z))
                SN_line = np.array(SN[center_bin-2:center_bin+2])
                quad_SN += max(SN_line*(SN_line>0))**2
            quad_SN = np.sqrt(quad_SN)
            if quad_SN > peak.total_sn:
                 peak.total_sn = quad_SN
                 peak.z = test_z


def jackpotLens(obj, peak, peak_candidates em_lines, mask_width_Jackpot = 20, total_SN_threshold = 6):
    '''
    peakFinder.jackpotpLens(obj, peak, em_lines, mask_width_Jackpot = 50):
    =============================
    Determines if a peak can be attributed to a background 

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak: The inquired peak
        peak_candidates: All other high SN peak detected in the spectra
        em_lines: The set of emission tracing the background ELG 
        mask_width: Width to mask first lens features after detection
        total_SN_threshold: Thresholding to determine if additional lines provide 
        	meaningful extra SN to be considered as background ELG
    Returns:
        - Nothing. Updates peak attributes by recording the best emission match.
    '''

    first_lens = False

    for l in em_lines:
        test_z = peak.wavelength_1/l -1.0
        if test_z > obj.z+0.05:
            quad_SN_1 = 0.0
            for w in em_lines:
                if w*(1+test_z) <= 9500:
                    center_bin = obj.wave2bin(w*(1+test_z))
                    SN_line = np.array(obj.SN[center_bin-2:center_bin+2])* (not nearline(w*(1+test_z), width= mask_width_Jackpot))
                    quad_SN_1 += max(SN_line*(SN_line>0))**2
            quad_SN_1 = np.sqrt(quad_SN_1)
            if quad_SN_1 > peak.sn_1 + total_SN_threshold:
                peak.total_sn_1 = quad_SN_1
                peak.z_1 = test_z
                first_lens = True
    if first_lens:
        for peak2 in peak_candidates:
            if np.abs(peak2.wavelength_1- peak.wavelength_1)> 30:
                for l in em_lines:
                    test_z_2 = peak2.wavelength_1/l -1.0
                    if test_z_2 > obj.z+ 0.05 and abs(peak.z_1-test_z_2)> 0.05:
                        quad_SN_2 = 0.0
                        for w in em_lines:
                            if w*(1+test_z_2) < 9500:
                                center_bin = obj.wave2bin(w*(1+test_z_2))
                                SN_line = np.array(obj.SN[center_bin-2:center_bin+2])*(not nearline(w*(1+test_z_2), width = mask_width_Jackpot))
                                quad_SN_2 += max(SN_line*(SN_line>0))**2
                        quad_SN_2 = np.sqrt(quad_SN_2)
                        if quad_SN_2 > peak2.sn + total_SN_threshold:
                        	peak.wavelength_2 = peak2.wavelength_2
                        	peak.sn_2 = peak2.sn
                            peak.total_sn_2 = quad_SN_2
                            peak.z_2 = test_z_2

                            peak.wavelength = min(peak.wavelength_1, peak.wavelength_2)


def doubletO2(obj, peak, bounds, max_chi2):
    '''
    peakFinder.doubletO2(obj, peak, bounds, QSOlens, max_chi2)
    =============================
    A function to fit the investigated SN peak with a doublet Gaussian.

    Parameters:
        obj: The SDSS object/spectra at hand
        peak: The inquired peak
        bounds: Bounds to define the support of the fitted polynomial
        max_chi2: Max chi square for the fitting
    Returns:
        - Nothing. Updates the appropriate parameters in the peak class    
    '''
    x0 = peak.wavelength
    params2, chisq2 = obj.doubletFit(bounds,
                                     [1.0, 5.0, 1.0, x0 - 1.5, x0 + 1.5],
                                     [(0.1, 5.0), (1.0, 8.0), (0.1, 5.0),
                                      (x0 - 7.0, x0), (x0, x0 + 7.0)])
    # Delta OII restframe:  1.3 A  (3725.94 3727.24)
    if (chisq2 < max_chi2) and (0.5 * x0 < abs(params2[3] - params2[4]) * 3726.5 < 4.0 * x0) \
        and (not obj.obj_class.startswith('QSO')) :
        # Here above we check that the foreground is not a QSO where this procedure is useless
        peak.chiDoublet = chisq2
        peak.ampDoublet = np.array([params2[0], params2[2]])
        peak.varDoublet = params2[1]
        peak.wavDoublet = np.array([params2[3], params2[4]])

def skewFit(obj, peak, bounds, max_chi2):
    '''
    peakFinder.doubletO2(obj, peak, bounds, QSOlens, max_chi2)
    =============================
    A function to fit the investigated SN peak with a skew-normal distribution.

    Because the LyA is made from two wings, we test three model: symmetric wings,
    sharp red emission with blue tail and sharp blue emission with red tail. 

    Parameters:
        obj: The SDSS object/spectra
        peak: The inquired peak
        bounds: Bounds to define the support of the fitted polynomial
        max_chi2: Max chi square for the fitting
    Returns:
        - Nothing. Updates the appropriate parameters in the peak class    
    '''
    
    # For convenience, define
    x0 = peak.wavelength
    if obj.obj_class.startswith('QSO'):
        temp_flux = obj.reduced_flux_QSO[bounds]
    else:
        temp_flux = obj.reduced_flux[bounds]

    # Sharp blue, red tail
    init_skew = [peak.ampSinglet,0.5,2,x0]
    res_skew = minimize(chi2skew,init_skew,args=(obj.wave[bounds], temp_flux,obj.ivar[bounds]), \
     method='SLSQP', bounds = [(2,50),(0.0,10),(1,10),(x0-4,x0+4)])
    params_skew_a = res_skew.x
    chisq_skew_a = res_skew.fun

    # Sharp red, blue tail
    init_skew = [peak.ampSinglet,0.5,-2,x0]
    res_skew = minimize(chi2skew,init_skew,args=(obj.wave[bounds], temp_flux,obj.ivar[bounds]), \
        method='SLSQP', bounds = [(2,50),(0.0,10),(-10,-1),(x0-4,x0+4)])
    params_skew_b = res_skew.x
    chisq_skew_b = res_skew.fun

    # Double skew symmetric
    init_skew = [peak.ampSinglet,0.5,-2,x0, peak.ampSinglet/2.0,0.5,2,x0+8]
    res_skew = minimize(chi2skew2,init_skew,args=(obj.wave[bounds], temp_flux,obj.ivar[bounds]), \
        method='SLSQP', bounds = [(2,50),(0.0,10),(-10,-1),(x0-6,x0+6), (2,50),(0.0,10),(1,10),(x0-15,x0+15)])
    params_skew_c = res_skew.x
    chisq_skew_c = res_skew.fun

    if chisq_skew_b < chisq_skew_a:
        params_skew = params_skew_b
        chisq_skew = chisq_skew_b
    elif chisq_skew_a < chisq_skew_b:
        params_skew = params_skew_a
        chisq_skew = chisq_skew_a
    #if chisq_skew < chisq_skew_c:
    peak.chiSkew = chisq_skew
    peak.ampSkew = params_skew[0] #A
    peak.scaleSkew = params_skew[1] #w
    peak.skewSkew = params_skew[2] #a
    peak.locSkew = params_skew[3] #eps
        
        #peak[7] = params_skew[0] #A
        #peak[8] = params_skew[1] #w
        #peak[9] = params_skew[2] #a
        #peak[10] = params_skew[3] #eps
        #if chisq_skew < chi2_width:
        #    chi2_width = chisq_skew
    #else:
    peak.chiSkew_doublet = chisq_skew_c
    peak.ampSkew  = np.array([params_skew_c[0], params_skew_c[4]]) #A1,A2
    peak.scaleSkew = np.array([params_skew_c[1], params_skew_c[5]]) #w1, w2
    peak.skewSkew  = np.array([params_skew_c[2], params_skew_c[6]]) #a1, a2
    peak.locSkew = np.array([params_skew_c[3], params_skew_c[7]]) #eps1, eps2
        #peak[11] = params_skew_c[4] #A2
        #peak[12] = params_skew_c[5] #w2
        #peak[13] = params_skew_c[6] #a2
        #peak[14] = params_skew_c[7] #eps2
        #if chisq_skew_c < chi2_width:
        #    chi2_width = chisq_skew_c

    ##put back reduced flux by adding again 3rd order fit (plotting purpose)
    #if QSOlens:
    #    reduced_flux[i,window]= new_flux + fit_QSO(wave[window])




