import itertools as it
import numpy as np
from utils import gauss, kernel


class peakCandidate():
    '''
    peakCandidate
    =============
    A class for all possible candidates of the peak
    '''
    def __init__(self, x0, sn):
        self.sn = sn                            # Original [4]
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

    def update(self):
        if self.chi > self.chiDoublet > self.chiSinglet:
            self.setDoublet(False)
        elif self.chi > self.chiSinglet > self.chiDoublet:
            self.setDoublet(True)


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


def qsoContfit(obj, peak, searchLyA):
    # TODO: complete the function/parameters
    '''
    window = n.linspace(wave2bin(x0,c0,c1,Nmax)-40,wave2bin(x0,c0,c1,Nmax)+40,81,dtype = n.int16)
    median_local = n.median(reduced_flux[i,window])
    fit_QSO = n.poly1d(n.polyfit(x=wave[window],y=reduced_flux[i,window],deg=3,w=(n.abs(reduced_flux[i,window]-median_local)<5)*n.sqrt(ivar[i,window])) )
    new_flux = reduced_flux[i,window] - fit_QSO(wave[window])
    cj1_new = n.sum(new_flux*kernel(int(len(window)/2),width,NormGauss,len(new_flux))*ivar[i,window])
    cj2_new = n.sum(ivar[i,window]*kernel(int(len(window)/2),width,NormGauss,len(window))**2)
    SN_fitted = cj1_new/n.sqrt(cj2_new)
    if searchLyA and SN_fitted < 6:
        continue
    elif searchLyA and SN_fitted > 6:
        peak[19] = SN_fitted
        reduced_flux[i,window]=new_flux
    elif searchLyA == False and SN_fitted < 6:
        return False  # Reject
    elif searchLyA == False and SN_fitted > 6:
        peak[3] = SN_fitted
        reduced_flux[i,window]=new_flux
    '''
    return True  # Accept


def qsoBggal(obj, peak, em_lines):
    # TODO: complete the function/parameters
    '''
    for l in em_lines:
        test_z = peak[0]/l - 1.0
        if test_z > z[i]:
            quad_SN = 0.0
            for w in em_lines:
                center_bin = wave2bin(w*(1+test_z),c0,c1,Nmax)
                SN_line = n.array(SN[center_bin-2:center_bin+2])
                quad_SN += max(SN_line*(SN_line>0))**2
            quad_SN = n.sqrt(quad_SN)
            if quad_SN > peak[5]:
                 peak[5] = quad_SN
                 peak[4] = test_z
     '''


def jpLens(obj, peak, em_lines, mask_width_Jackpot = 50):
    # TODO: complete the function/parameters
    '''
    first_lens = False
    peak[5] = peak[6]
    peak[4] = peak[6]
    for l in em_lines:
        test_z = peak[0]/l -1.0
        if test_z > z[i]+0.05:
            quad_SN_1 = peak[6]
            for w in em_lines:
                if w*(1+test_z) < 9500:
                    center_bin = wave2bin(w*(1+test_z),c0,c1,Nmax)
                    SN_line = n.array(SN[center_bin-2:center_bin+2])*(not(nearline(w*(1+test_z), zline, fiberid[i], z[i], int(mjd), int(plate), mask_width_Jackpot)))
                    quad_SN_1 += max(SN_line*(SN_line>0))**2
            quad_SN_1 = n.sqrt(quad_SN_1)
            if quad_SN_1 > peak[6] + 6:
                peak[5] = quad_SN_1
                peak[2] = test_z
                first_lens = True
    if first_lens:
        for peak2 in peak_candidates:
            if n.abs(peak2[0]- peak[0])> 30:
                for l in em_lines:
                    test_z_2 = peak2[0]/l -1.0
                    if test_z_2 > z[i]+ 0.05 and abs(peak[2]-test_z_2)> 0.05:
                        quad_SN_2 = 0.0
                        for w in em_lines:
                            if w*(1+test_z_2) < 9500:
                                center_bin = wave2bin(w*(1+test_z_2),c0,c1,Nmax)
                                SN_line = n.array(SN[center_bin-2:center_bin+2])*(not(nearline(w*(1+test_z_2), zline, fiberid[i], z[i], int(mjd), int(plate),mask_width_Jackpot)))
                                quad_SN_2 += max(SN_line*(SN_line>0))**2
                        quad_SN_2 = n.sqrt(quad_SN_2)
                        if quad_SN_2 > peak[5] + 6:
                            peak[4] = quad_SN_2
                            peak[3] = test_z_2
    '''


def doubletO2(obj, peak, bounds, searchLyA, max_chi2):
    x0 = peak.wavelength
    params2, chisq2 = obj.doubletFit(bounds,
                                     [1.0, 5.0, 1.0, x0 - 1.5, x0 + 1.5],
                                     [(0.1, 5.0), (1.0, 8.0), (0.1, 5.0),
                                      (x0 - 7.0, x0), (x0, x0 + 7.0)])
    if ((not (searchLyA or chisq2 > max_chi2)) and
        0.5 * x0 < abs(params2[3] - params2[4]) * 3726.5 < 4.0 * x0):
        peak.chiDoublet = chisq2
        peak.ampDoublet = np.array([params2[0], params2[2]])
        peak.varDoublet = params2[1]
        peak.wavDoublet = np.array([params2[3], params2[4]])
    # TODO: complete the 2nd part of the function
    '''
    elif searchLyA and chisq2<chisq:
        peak[1] = params2[3] #x1
        peak[2] = params2[4] #x2
        peak[3] = params2[0] #amp1
        peak[4] = params2[2] #amp2
        peak[5] = params2[1] #var
        chi2_width = chisq2
        peak[16] = chisq2
    elif searchLyA:
        peak[16] = chisq2
        # Delta OII restframe:  1.3 A  (3725.94 3727.24)
    '''


def skewFit(obj, peak):
    # TODO: complete the function/parameters
    '''
                # Sharp blue, red tail
                init_skew = [params[1],0.5,2,x0]
                res_skew = minimize(chi2skew,init_skew,args=(wave[bounds], reduced_flux[i,bounds],ivar[i,bounds]), method='SLSQP', bounds = [(2,50),(0.0,10),(1,10),(x0-4,x0+4)])
                params_skew_a = res_skew.x
                chisq_skew_a = res_skew.fun
                # Double skew symmetric
                init_skew = [params[1],0.5,-2,x0, params[1]/2,0.5,2,x0+8]
                res_skew = minimize(chi2skew2,init_skew,args=(wave[bounds], reduced_flux[i,bounds],ivar[i,bounds]), method='SLSQP', bounds = [(2,50),(0.0,10),(-10,-1),(x0-6,x0+6), (2,50),(0.0,10),(1,10),(x0-15,x0+15)])
                params_skew_c = res_skew.x
                chisq_skew_c = res_skew.fun

                # Sharp red, blue tail
                init_skew = [params[1],0.5,-2,x0]
                res_skew = minimize(chi2skew,init_skew,args=(wave[bounds], reduced_flux[i,bounds],ivar[i,bounds]), method='SLSQP', bounds = [(2,50),(0.0,10),(-10,-1),(x0-4,x0+4)])
                params_skew_b = res_skew.x
                chisq_skew_b = res_skew.fun

                if chisq_skew_b < chisq_skew_a:
                    params_skew = params_skew_b
                    chisq_skew = chisq_skew_b
                elif chisq_skew_a < chisq_skew_b:
                    params_skew = params_skew_a
                    chisq_skew = chisq_skew_a
                if chisq_skew < chisq_skew_c:
                    peak[7] = params_skew[0] #A
                    peak[8] = params_skew[1] #w
                    peak[9] = params_skew[2] #a
                    peak[10] = params_skew[3] #eps
                    if chisq_skew < chi2_width:

                        chi2_width = chisq_skew
                else:
                    peak[7] = params_skew_c[0] #A1
                    peak[8] = params_skew_c[1] #w1
                    peak[9] = params_skew_c[2] #a1
                    peak[10] = params_skew_c[3] #eps1
                    peak[11] = params_skew_c[4] #A2
                    peak[12] = params_skew_c[5] #w2
                    peak[13] = params_skew_c[6] #a2
                    peak[14] = params_skew_c[7] #eps2
                    if chisq_skew_c < chi2_width:

                        chi2_width = chisq_skew_c

                peak[17] = chisq_skew
                peak[18] = chisq_skew_c

                #put back reduced flux by adding again 3rd order fit (plotting purpose)
                if QSOlens:
                    reduced_flux[i,window]= new_flux + fit_QSO(wave[window])
    '''



