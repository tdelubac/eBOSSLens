# Last modification: Xinlun Cheng, 2017
# Original author: Romain A. Meyer, 2016, LASTRO, EPFL
# Spectroscopic lensing systems detection tool
# Detects OII doublet for galaxy-galaxy lensing systems
# Searches for Lyman alpha emission for galaxy-LAE or QSO-LAE systems
# 'Paper' mode removes intermediate steps plots useful for analysis but
# irrelevant for publication


# Imports
import numpy as n
import os
from SDSSObject import SDSSObject
from objFilter import qsoFilter, genFilter
from peakFinder import bolEstm, peakCandidateGalGal, peakCandidateQSOGal, \
    peakCandidateQSOLAE, peakCandidateQSOGal, peakCandidateJackpot, \
    combNear, checkFore, qsoContfit, backgroundELG, jackpotLens,  \
    doubletO2, skewFit
from galSave import galSave
from qsoSave import qsoSave
from jptSave import jptSave
from lyaSave import lyaSave

# Set of emission lines used for lensed galaxy detection:
# OII, Hb, OIII, Ha
em_lines = n.array([3726.5, 4861.325, 4958.911, 5006.843, 6562.801])
# Sky pollution Wavebands for masking
wMask = n.array([[5570.0, 5590.0], [5880.0, 5905.0], [6285.0, 6315.0],
                 [6348.0, 6378.0]])

def eBOSSLens(plate, mjd, fiberid, datav, searchLyA, QSOlens, Jackpot, savedir,
              datadir, max_chi2=4.0, wMask=wMask, em_lines=em_lines,
              bwidth=30.0, bsig=1.2, cMulti=1.04, doPlot=False,
              prodCrit=1000.0, QSO_line_width = 15, minSN = 6.0, threshold_SN = 6.0,
              jackpotwidth = 20):
    '''
    eBOSSLens
    =============
    The Spectroscopic Lens Finder algorithm.
  
    Parameters:
        obj: A SDSSObject instance to inspect
        datav: 
        searchLyA: True if looking for background LAE. False for background ELGs
        QSOlens: True if looking for QSO as foreground lens. False for foreground ELGs
        Jackpot: True if looking for ELGs lensing **2** background ELGs. Overrides the 2 previous Booleans
        savedir: Directory where to save detections. (Note for QSO Lenses, DR12Q must be stored here)
        datadir: Directory where the .fits files (spPlate, spZBest, spzLine) of the SDSS spectra are stored
        max_chi2: Max chi^2 for the doublet fitting in the case Galaxy-ELG lensing
        wMask: Sky lines / Known defects to mask in the spectra
        em_lines: The emission lines to look at for background ELG search
        bwidth: Width of window to peform (Bolton 2004) SN convolution
        bsig: Width of the Gaussian kernel for (Bolton 2004) SN search
        cMulti: Multiple to set threshold on non-gaussian fits (chi^2_{ng} > cMulti*chi^2_{g}). Must be >= 1
        doPlot: True to plot directly the detections. False to save the plots as images
        prodCrit:
        QSO_line_width: Typical half-width of QSO strong emission to mask
        minSN: Minimal SN to detect promising peaks
        threshold_SN: Threshold to select potential ELG emissions
        jackpotwidth: Width of masking on first lensed features before search for second background source
    Returns:
        returns: Nothing. Raises appropriate exceptions when candidates are discarded.
            Detections are saved (plot+attributes in txt file) in a child directory of 'savedir'
    '''

    if cMulti < 1:
        raise Exception('cMulti must be greater than 1 for proper chi2 comparisons.')

    # Define LyA wavelength (Angstroms)
    l_LyA = 1215.668

    # TODO : Externalize the SDSSObject creation
    obj = SDSSObject(plate, mjd, fiberid, datav, datadir)

    accept = False
    if QSOlens:
        # Filter out unwanted spectras
        DR12Q = DR12Q_extractor(path=os.path.join(savedir, 'Superset_DR12Q.fits'))
        accept = qsoFilter(obj, DR12Q)
    else:
        accept = genFilter(obj)
        obj.mask(wMask)
    #if not accept:
    #    raise Exception("Rejected by filter")

    # Mask BOSS spectra glitches + Sky
    obj.mask(wMask)
    if QSOlens:
        obj.mask(mask_QSO(QSO_line_width))

    # Find peaks
    doublet = None
    # Bolton 2004: S/N of maximum likelihood estimator of gaussian peaks
    if searchLyA:
        width = bwidth
        sig = bsig 
    else:
        width = bwidth
        sig = bsig
    obj.SN = bolEstm(obj, sig, width)

    # Select the max likelihood peak depending on the case
    if Jackpot:
        peak_candidates = np.array([peakCandidateJackpot(x0, test) \
            for x0, test in zip(obj.wave, obj.SN) if test > 6.0])
    elif not (searchLyA or QSOlens):
        peak_candidates = np.array([peakCandidateGalGal(x0, test) \
            for x0, test in zip(obj.wave, obj.SN) if test > 6.0])
    elif (not(QSOlens) and searchLyA):
        peak_candidates = np.array([peakCandidateGalLAE(x0, test) \
            for x0, test in zip(obj.wave, obj.SN) if test > 8.0])
    elif (not(searchLyA) and QSOlens):
        peak_candidates = np.array([peakCandidateQSOGal(x0, test) \
            for x0, test in zip(obj.wave, obj.SN) if (test > 8.0 and  (l_LyA*(1+obj.z[i])+300)<x0<9500)])
    elif (searchLyA and QSOlens): 
        peak_candidates = np.array([peakCandidateQSOLAE(x0, test) \
            for x0, test in zip(obj.wave, obj.SN) if (test > 8.0 and  (l_LyA*(1+obj.z[i])+300)<x0<9500)])
    else:    
        raise Exception('Error: Foreground/background objects boolean combinations not found.')

    # TODO: delete if above is working
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
            accept = qsoContfit(obj, peak, searchLyA)
            if not accept:
                raise Exception('Rejected: Peak detection due to incorrect \
                    QSO continuum removal')

        # Special case: QSOlens with background galaxies
        if (not (searchLyA or Jackpot)) and QSOlens:
            backgroundELG(obj, peak, em_lines)
            continue
        # Special case: Jackpot lenses
        if Jackpot:
            jackpotLens(obj, peak, peak_candidates, em_lines, jackpotwidth, threshold_SN)
            continue

        # Singlet
        if searchLyA:
            init = [x0, 4.0, 6.0]
            paramLim = [(x0 - 2.0, x0 + 2.0), (1.0, 100.0), (1.0, 15.0)]
        elif not (searchLyA or QSOlens):
            init = [x0, 1.0, 2.0]
            paramLim = [(x0 - 2.0, x0 + 2.0), (0.1, 5.0), (1.0, 8.0)]
        #Perform the single Gaussian fit    
        params, chisq = obj.singletFit(bounds, init, paramLim)
        # Check for not too high chi square and save
        if not (chisq > max_chi2 or searchLyA or QSOlens):
            peak.wavSinglet = params[0]
            peak.ampSinglet = params[1]
            peak.varSinglet = params[2]
            peak.chiSinglet = chisq
        elif searchLyA:
            peak.wav_g = params[0]
            peak.amp_g = params[1]
            peak.var_g = params[2]
            peak.chi_g = chisq

        # Doublet OII
        if x0 > 3727.0 * (1.0 + obj.z) or searchLyA and (not QSOlens):
            doubletO2(obj, peak, bounds, max_chi2)
        # If looking at LAE, test a skew-normal profile as well
        if searchLyA and (not QSOlens):
            skewFit(obj, peak)

    # Compare the fitting results
    if not (QSOlens or Jackpot):
        # Compare singlet and doublet fit within each peakCandidate
        # If searching for Gal-LAE, compares also to skew fit
        for k in range(len(peak_candidates)):
            pk = peak_candidates[k]
            pk.update(cMulti)
        # Removing candidates that were not fitted
        peak_candidates = n.array([peak for peak in peak_candidates if
                                   (peak.chi != 1000.0)])
        if len(peak_candidates) == 0:
            raise Exception("Rejected since no peak candidates")
        # Sorting candidates by chi square
        peak_candidates = sorted(peak_candidates, key=lambda peak: peak.chi)
        # Keeping only 5 most likely candidates
        if len(peak_candidates) > 5:
            peak_candidates = peak_candidates[0:5]
        # Check that doublet is in the 5 most likely remaining
        if (not searchLyA):
            doublet_index = 0
            doublet = False
            for k in range(len(peak_candidates)):
                if peak_candidates[k].isDoublet:
                    doublet_index = k
                    doublet = True
                    break
        # Conversely, check that skew fit is in the 5 most likely remaining
        elif:
            skew_index = 0
            skew = False
            for k in range(len(peak_candidates)):
                if peak_candidates[k].isSkew:
                    skew_index = k
                    skew = True
                    break

    '''
    QSO-ELG case: We retain only the significant boost in SN when taking 
    multiple emissions into account compared to a single peak
    '''
    elif (not (searchLyA or Jackpot)) and QSOlens:
        peak_candidates = np.array([peak for peak in peak_candidates  \
            if peak.total_sn >(1.5+peak.sn)])
        if len(peak_candidates) == 0:
            raise Exception('Candidate rejected. SN peak is not part of a \
            series of background ELG emissions.')
        # Keep only 3 best candidates since useless duplicates often happen.     
        peak_candidates = sorted(peak_candidates, key=lambda peak: peak.total_sn)
        if len(peak_candidates) > 3:
            peak_candidates = peak_candidates[0:3]
    '''
    Jackpot case: We retain only the significant boost in SN when taking 
    multiple emissions into account compared to a single peak
    '''       
    elif Jackpot:
        peak_candidates = n.array([peak for peak in peak_candidates if peak.z_2 > 0.0])
        if len(peak_candidates) == 0:
            raise Exception("Rejected since no peak retained")
        peak_candidates = sorted(peak_candidates, key=lambda peak: peak.total_sn_1 + peak.total_sn_2)
        # Keep only 3 best candidates since useless duplicates often happen.   
        if len(peak_candidates) > 3:
            peak_candidates = peak_candidates[0:3]

    if len(peak_candidates) == 0:
        raise Exception("Rejected since no peak retained")
    # Check that at least 1 candidate is below 9500 Angstrom cut
    below_9500 = False
    for peak in peak_candidates:
        if peak.wavelength < 9500.0:
            below_9500 = True
            break
    if not below_9500:
        raise Exception("Rejected since no peak below 9500 A")

    # Try to infer background redshift and if successful, save the detection

    if not (searchLyA or QSOlens or Jackpot):
        galSave(doublet, obj, peak_candidates, doublet_index, savedir, em_lines,
                doPlot, prodCrit)
    elif (not (searchLyA or Jackpot)) and QSOlens:
        # TODO: complete the function QSOGAL case
        qsoSave(obj, peak_candidates, savedir,em_lines) 
    elif searchLyA and (not Jackpot):
        lyaSave(obj, peak_candidates, savedir,em_lines, threshold_SN, QSOlens, paper_mode)
    elif Jackpot:
        # TODO: complete the function
        jptSave(doPlot)
