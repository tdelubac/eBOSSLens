# Last modification: Xinlun Cheng, 2017
# Original author: Romain A. Meyer, 2016, LASTRO, EPFL
# Spectroscopic lensing systems detection tool
# Detects OII doublet for galaxy-galaxy lensing systems
# Searches for Lyman alpha emission for galaxy-LAE or QSO-LAE systems
# 'Paper' mode removes intermediate steps plots useful for analysis but
# irrelevant for publication


# Imports
import numpy as n
from SDSSObject import SDSSObject
from objFilter import qsoFilter, genFilter
from peakFinder import bolEstm, peakCandidate, combNear, checkFore, \
    qsoContfit, qsoBggal, jpLens, doubletO2, skewFit
from galSave import galSave
from qsoSave import qsoSave
from jptSave import jptSave
from lyaSave import lyaSave


# Currently unused, but might be useful parameters
'''
# Typical mask width
l_width = 15
# TODO: clean up
l_LyA = 1215.668
'''
# Set of emission lines used for lensed galaxy detection:
# OII, Hb, OIII, Ha
em_lines = n.array([3726.5, 4861.325, 4958.911, 5006.843, 6562.801])
# Waves for mask
wMask = n.array([[5570.0, 5590.0], [5880.0, 5905.0], [6285.0, 6315.0],
                 [6348.0, 6378.0], [6071.0, 6091.0], [6555.0, 6575.0],
                 [8335.0, 8355.0], [8819.0, 8839.0], [8988.0, 9008.0],
                 [9217.0, 9237.0], [9300.0, 9330.0], [9365.0, 9385.0],
                 [9425.0, 9445.0], [9777.0, 9797.0]])
# Strong emission lines should also be masked
eMask = n.array([4103.0, 4342.0, 4863.0, 6563.0, 5008.0, 4960.0, 4364.0])


def eBOSSLens(plate, mjd, fiberid, datav, searchLyA, QSOlens, Jackpot, savedir,
              datadir, max_chi2=2.5, wMask=wMask, eMask=eMask, em_lines=em_lines,
              bwidth=60.0, bsig=1.2, cMulti=1.04, doPlot=False,
              prodCrit=0.6, emWidth=10.0, snCrit=50.0):
    obj = SDSSObject(plate, mjd, fiberid, datav, datadir)
    # Mask BOSS spectra glitches + Sky
    obj.mask(wMask)
    # Filter out unwanted spectras
    accept = False
    if QSOlens:
        # TODO: complete the function
        accept = qsoFilter(obj)
    else:
        accept = genFilter(obj)
    if not accept:
        raise Exception("Rejected by filter")
    # Mask strong emission lines and sn filter
    if not (QSOlens or searchLyA or Jackpot):
        if obj.sn > snCrit:
            raise Exception("Rejected by unusual high SN ratio")
        for each in eMask:
            obj.mask(n.array([[each - emWidth, each + emWidth]]) *
                     (obj.z + 1.0))
    # Find peaks
    doublet = None
    # Bolton 2004: S/N of maximum likelihood estimator of gaussian peaks
    if searchLyA:
        width = 30.0
        sig = 2.17
    else:
        width = bwidth
        sig = bsig
    SN = bolEstm(obj, sig, width)
    # Filter out invalid sources
    if not (searchLyA or QSOlens):
        peak_candidates = n.array([peakCandidate(x0, test)
                                   for x0, test in zip(obj.wave, SN)
                                   if test > 6.0])
    else:
        return
    # TODO: peakCandidate for qso or/and lya
    '''
    elif searchLyA and QSOlens:
        peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(obj.wave,SN) if (test > 8.0 and  (l_LyA*(1+obj.z[i])+300)<x0<9500)])
    elif searchLyA == False and QSOlens:
        peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test) for x0,test in zip(obj.wave,SN) if (test>6.0 and  (l_LyA*(1+obj.z[i])+300)<x0<9500)])
    elif searchLyA == True and QSOlens == False:
        peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(obj.wave,SN) if (test>8.0 and  3600 < x0 < 4800)])
    elif Jackpot == True:
        # x0 z1 z2 Quad_SN2 SN0->Quad_SN1  free free
        peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test) for x0,test in zip(obj.wave,SN) if test>8.0])
    '''
    if len(peak_candidates) == 0:
        raise("Bolton estimator does not found any SN > 6 peaks")
    # Keep the center
    peak_candidates = combNear(peak_candidates)
    # Check hits are not from foreground galaxy or badly fitted QSO
    if QSOlens and searchLyA:
        if checkFore(peak_candidates, em_lines, obj.z):
            return
    # --------------------------------------------------------------------------
    # Search for suitable peak candidates
    for peak in peak_candidates:
        x0 = peak.wavelength
        if obj.nearLine(x0):
            continue
        x0Bin = obj.wave2bin(x0)
        bounds = n.linspace(x0Bin - 15, x0Bin + 15, 61, dtype=n.int16)
        # Fit QSO continuum and check if signal is reduced or not
        # i.e. check if line detection is produced by large features
        accept = False
        if QSOlens:
            # TODO: complete the function
            accept = qsoContfit(obj, peak, searchLyA)
            if not accept:
                continue
        # Special case: QSOlens with background galaxies
        if (not (searchLyA or Jackpot)) and QSOlens:
            # TODO: complete the function
            qsoBggal(obj, peak, em_lines)
            continue
        # Special case: Jackpot lenses
        if Jackpot:
            # TODO: complete the function
            jpLens(obj, peak, em_lines)
            continue
        # Singlet
        if searchLyA:
            init = [x0, 4.0, 6.0]
            paramLim = [(x0 - 2.0, x0 + 2.0), (1.0, 100.0), (1.0, 15.0)]
        elif not (searchLyA or QSOlens):
            init = [x0, 1.0, 2.0]
            paramLim = [(x0 - 2.0, x0 + 2.0), (0.1, 5.0), (1.0, 8.0)]
        params, chisq = obj.singletFit(bounds, init, paramLim)
        # Check for not too high chi square and save
        if not (chisq > max_chi2 or searchLyA or QSOlens):
            peak.wavSinglet = params[0]
            peak.ampSinglet = params[1]
            peak.varSinglet = params[2]
            peak.chiSinglet = chisq
        # TODO: clean up
        '''
        elif searchLyA:
            peak.wavelength = params[0]
            peak[2] = params[1]
            peak[3] = params[2]
            chi2_width = chisq
            peak[15] = chisq
        '''
        # Doublet OII
        if x0 > 3727.0 * (1.0 + obj.z) or searchLyA and (not QSOlens):
            # TODO: finish the 2nd part of the function
            doubletO2(obj, peak, bounds, searchLyA, max_chi2)
        # If looking at LAE, test a skew-normal profile as well
        if searchLyA:
            # TODO: complete the function
            skewFit(obj, peak)
    # Compare the fitting results
    if not (searchLyA or QSOlens or Jackpot):
        # Compare singlet and doublet fit within each peakCandidate
        for k in range(len(peak_candidates)):
            pk = peak_candidates[k]
            pk.update(cMulti)
        # Removing candidates that were not fitted
        peak_candidates = n.array([peak for peak in peak_candidates if
                                   peak.chi != 1000.0])
        if len(peak_candidates) == 0:
            raise Exception("Rejected since no valid peak")
        # Hard cut everything above 9200A
        peak_candidates = n.array([peak for peak in peak_candidates if
                                   peak.wavelength < 9200.0])
        if len(peak_candidates) == 0:
            raise Exception("Rejected since no peak below 9200")
        # Sorting candidates by chi square
        peak_candidates = sorted(peak_candidates, key=lambda peak: peak.chi)
        # Keeping only 5 most likely candidates
        if len(peak_candidates) > 5:
            peak_candidates = peak_candidates[0:5]
        # Find whether doublet is in the 5 most likely
        doublet_index = 0
        doublet = False
        for k in range(len(peak_candidates)):
            if peak_candidates[k].isDoublet:
                doublet_index = k
                doublet = True
                break
    # TODO: clean up
    '''
    elif (not (searchLyA or Jackpot)) and QSOlens:
        peak_candidates = n.array([peak for peak in peak_candidates if peak[5]>(1.5+peak[6])])
        if len(peak_candidates) == 0:
            continue
        peak_candidates = sorted(peak_candidates, key=lambda peak: peak[5])
        if len(peak_candidates) > 3:
            peak_candidates = peak_candidates[0:3]
    elif Jackpot:
        peak_candidates = n.array([peak for peak in peak_candidates if peak[3]>0.0])
        if len(peak_candidates) == 0:
                continue
        peak_candidates = sorted(peak_candidates, key=lambda peak: peak[5])
        if len(peak_candidates) > 3:
            peak_candidates = peak_candidates[0:3]
    '''
    # Try to infer background redshift
    if not (searchLyA or QSOlens or Jackpot):
        galSave(doublet, obj, peak_candidates, doublet_index, savedir, em_lines,
                doPlot, prodCrit)
    elif (not (searchLyA or Jackpot)) and QSOlens:
        # TODO: complete the function
        qsoSave(peak_candidates, savedir, doPlot)
    elif Jackpot:
        # TODO: complete the function
        jptSave(doPlot)
    elif searchLyA and QSOlens and (not Jackpot):
        # TODO: complete the function
        lyaSave(doPlot)
