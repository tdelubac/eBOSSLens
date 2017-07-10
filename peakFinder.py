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


def nearLine(mjd, plate, i, x0, obj, width=10.0):
    '''
    nearLine(mjd, plate, i, x0, obj, width=10.0)
    ============================================
    Whether given wavelength is near a masked or not

    Parameters:
        mjd: mjd
        plate: plate
        i: the fiber that is being examined
        x0: given wavelength
        obj: SDSSObject for the source
        width: default 10.0, a parameter for nearness
    Returns:
        a boolean of whether near or not
    '''
    crit = ()
    crit.append(abs(obj.zline['linewave'] * (1 + obj.z[i]) - x0) < width)
    crit.append(obj.zline['lineew'] / obj.zline['lineew_err'] > 6.0)
    crit.append(obj.zline['fiberid'] == obj.fiberid[i])
    crit.append(obj.zline['mjd'] == int(mjd))
    crit.append(obj.zline['plate'] == int(plate))
    crit.append(obj.zline['lineew_err'] > 0)
    match = np.logical_and.reduce(crit)
    if np.any(match):
        return True
    else:
        return False


def bolEstm(i, obj, sig, width):
    NormGauss = gauss(np.linspace(-width * 0.5, width * 0.5, width), 0.0, 1.0,
                      sig ** 2.0)
    NormGauss = NormGauss / np.sum(NormGauss)
    allKer = np.array([kernel(j + 0.5 * width, width, NormGauss, len(obj.wave))
                       for j in range(int(len(obj.wave)) - width)])
    Cj1 = np.array([np.sum(each * obj.reduced_flux[i, :] * obj.ivar[i, :])
                    for each in allKer])
    Cj2 = np.array([np.sum(obj.ivar[i, :] * each ** 2.0) for each in allKer])
    SN = np.zeros(len(obj.wave))
    SN[int(width * 0.5): int(width * 0.5 + len(Cj1))] = Cj1 / np.sqrt(Cj2)
    return SN
