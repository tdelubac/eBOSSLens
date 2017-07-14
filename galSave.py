import os
import itertools as it
import numpy as np
from utils import make_sure_path_exists, gauss
from SDSSObject import SDSSObject
# Matplotlib trick
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


def galSave(doublet, obj, peak_candidates, doublet_index, savedir, em_lines,
            prodCrit=10.0):
    detection = False
    preProd = 1.0
    nxtProd = 1.0
    if doublet:
        if len(peak_candidates):
            preProd, nxtProd = compSpec(obj, peak_candidates[doublet_index])
            if preProd > prodCrit or nxtProd > prodCrit:
                raise Exception("Rejected by comparing to other fibers")
        z_s = peak_candidates[doublet_index].wavelength / 3727.24 - 1.0
        detection = _doubletSave(obj, z_s, peak_candidates, doublet_index,
                                 savedir)
        if len(peak_candidates):
            detection = _dblmultSave(obj, z_s, peak_candidates, savedir,
                                     detection, em_lines)
    elif len(peak_candidates) > 1:
        detection = _multletSave(obj, peak_candidates, savedir, em_lines)
    peaks = []
    for k in range(len(peak_candidates)):
        peak = peak_candidates[k]
        if k == doublet_index and doublet:
            peaks.append([peak.wavDoublet[0], peak.ampDoublet[0],
                          peak.varDoublet])
            peaks.append([peak.wavDoublet[1], peak.ampDoublet[1],
                          peak.varDoublet])
        else:
            peaks.append([peak.wavSinglet, peak.ampSinglet, peak.varSinglet])
    peak_number = len(peak_candidates)
    if (peak_number > 1 or doublet) and detection:
        fit = 0.0
        for k in np.arange(len(peaks)):
            fit = fit + gauss(obj.wave, x_0=peaks[k][0], A=peaks[k][1],
                              var=peaks[k][2])
        plotGalaxyLens(doublet, obj, savedir, peak_candidates, preProd,
                       nxtProd, doublet_index, fit)


def _doubletSave(obj, z_s, peak_candidates, doublet_index, savedir):
    score = 0.0
    detection = False
    fileD = open(os.path.join(savedir, 'candidates_doublet.txt'), 'a')
    if z_s > obj.z + 0.05:
        detection = True
        score += peak_candidates[doublet_index].chi
        fileD.write(str(obj.radEinstein(z_s)) + " " + str(score) +
                    " " + str(z_s) + " " + str(obj.RA) + " " +
                    str(obj.DEC) + " " + str(obj.plate) + " " +
                    str(obj.mjd) + " " + str(obj.fiberid) + " " +
                    str(peak_candidates[doublet_index].wavDoublet[0]) +
                    "\n")
    fileD.close()
    return detection


def _dblmultSave(obj, z_s, peak_candidates, savedir, detection, em_lines):
    confirmed_lines = []
    score = 0.0
    det = detection
    fileM = open(os.path.join(savedir, 'candidates_DM.txt'), 'a')
    # Generating all combinations of lines from above list to compare with
    # candidates
    temp = [peak for peak in peak_candidates if peak.chi != peak.chiDoublet]
    compare = em_lines[1: 5]
    if z_s > obj.z + 0.05:
        for peak in temp:
            for line in compare:
                if abs(peak.wavelength/line - 1.0 - z_s) < 0.01:
                    det = True
                    confirmed_lines.append(line)
                    score += peak.chiDoublet
    if confirmed_lines != []:
        fileM.write(str(obj.radEinstein(z_s)) + " " + str(score) + " " +
                    str(z_s) + " " + str(obj.RA) + " " + str(obj.DEC) +
                    " " + str(obj.plate) + " " + str(obj.mjd) + " " +
                    str(obj.fiberid) + " " + str(confirmed_lines) + "\n")
    fileM.close()
    return det


def _multletSave(obj, peak_candidates, savedir, em_lines):
    confirmed_lines = []
    score = 0.0
    detection = False
    compare = it.combinations(em_lines, len(peak_candidates))
    fileM = open(os.path.join(savedir, 'candidates_multi.txt'), 'a')
    for group in compare:
        for k in range(len(peak_candidates)):
            for j in range(k + 1, len(peak_candidates)):
                crit1 = peak_candidates[k].wavelength / group[k]
                crit2 = crit1 - peak_candidates[j].wavelength / group[j]
                if abs(crit2) < 0.01 and crit1 - 1.05 > obj.z:
                    detection = True
                    z_s = peak_candidates[k].wavelength / group[k] - 1.0
                    confirmed_lines.append([group,
                                            peak_candidates[k].wavelength /
                                            group[k] - 1.0])
                    score += peak_candidates[j].chi ** 2.0 + \
                        peak_candidates[k].chi ** 2.0
    if confirmed_lines != []:
        fileM.write(str(obj.radEinstein(z_s)) + " " + str(score) + " " +
                    str(z_s) + " " + str(obj.RA) + " " + str(obj.DEC) +
                    " " + str(obj.plate) + " " + str(obj.mjd) + " " +
                    str(obj.fiberid) + " " + str(confirmed_lines) + "\n")
    fileM.close()
    return detection


def plotGalaxyLens(doublet, obj, savedir, peak_candidates, preProd, nxtProd,
                   doublet_index, fit):
    if not doublet:
        ax = plt.subplot(1, 1, 1)
        plt.title('RA=' + str(obj.RA) + ', Dec=' + str(obj.DEC) + ', Plate=' +
                  str(obj.plate) + ', Fiber=' + str(obj.fiberid) +
                  ', MJD=' + str(obj.mjd) + '\n$z=' + str(obj.z) + ' \pm' +
                  str(obj.z_err) + '$, Class=' + str(obj.obj_class))
        ax.plot(obj.wave, obj.reduced_flux, 'k')
        plt.xlabel('$Wavelength\, (Angstroms)$')
        plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$')
        ax.plot(obj.wave, fit, 'r')
        make_sure_path_exists(savedir + '/plots/')
        plt.savefig(savedir + '/plots/' + str(obj.plate) + '-' + str(obj.mjd) +
                    '-' + str(obj.fiberid) + '.png')
        plt.close()
    # If doublet, plot in two different windows
    else:
        # Plot currently inspecting spectra
        x_doublet = np.mean(peak_candidates[doublet_index].wavDoublet)
        bounds = np.linspace(obj.wave2bin(x_doublet) - 10,
                             obj.wave2bin(x_doublet) + 10, 21, dtype=np.int16)
        f = open(os.path.join(savedir, 'doublet_ML.txt'), 'a')
        f.write(str(obj.plate) + ' ' + str(obj.mjd) + ' ' + str(obj.fiberid) +
                ' ' + str(preProd) + ' ' + str(nxtProd) + ' ' +
                str(obj.reduced_flux[bounds]) + "\n")
        f.close()
        plt.figure(figsize=(14, 6))
        ax1 = plt.subplot2grid((1, 3), (0, 0), colspan=2)
        plt.suptitle('RA=' + str(obj.RA) + ', Dec=' + str(obj.DEC) +
                     ', Plate=' + str(obj.plate) + ', Fiber='+str(obj.fiberid) +
                     ', MJD=' + str(obj.mjd) + '\n$z=' + str(obj.z) + ' \pm' +
                     str(obj.z_err) + '$, Class=' + str(obj.obj_class))
        ax2 = plt.subplot2grid((1, 3), (0, 2))
        ax1.plot(obj.wave[10:-10], obj.reduced_flux[10:-10], 'k')
        ax1.plot(obj.wave, fit, 'r')
        ax1.set_xlabel('$\lambda \, [\AA]$ ')
        ax1.set_ylabel(
            '$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$')
        ax2.set_xlabel('$\lambda \, [\AA]$ ')
        ax2.locator_params(tight=True)
        ax2.set_xlim([peak_candidates[doublet_index].wavelength - 30.0,
                      peak_candidates[doublet_index].wavelength + 30.0])
        ax2.plot(obj.wave, obj.reduced_flux, 'k')
        ax2.plot(obj.wave, fit, 'r')
        ax2.set_ylim([-5, 10])
        ax2.vlines(x=obj.zline['linewave'] * (1.0 + obj.z), ymin=-10, ymax=10,
                   colors='g', linestyles='dashed')
        ax1.set_xlim([np.min(obj.wave), np.max(obj.wave)])
        # Plot previous one
        if obj.fiberid != 1:
            objPre = SDSSObject(obj.plate, obj.mjd, obj.fiberid - 1,
                                obj.dataVersion, obj.baseDir)
            ax2.plot(objPre.wave, objPre.reduced_flux, 'b')
        # Plot next one
        if obj.fiberid != 1000:
            objNxt = SDSSObject(obj.plate, obj.mjd, obj.fiberid + 1,
                                obj.dataVersion, obj.baseDir)
            ax2.plot(objNxt.wave, objNxt.reduced_flux, 'g')
        # Save to file
        make_sure_path_exists(os.path.join(savedir, 'plots'))
        plt.savefig(os.path.join(savedir, 'plots', str(obj.plate) + '-' +
                                 str(obj.mjd) + '-' + str(obj.fiberid) +
                                 '.png'))
        plt.close()


def compSpec(obj, peak, width=2.0):
    bounds = np.arange(obj.wave2bin(peak.wavDoublet.min() - width *
                                    np.sqrt(peak.varDoublet)),
                       obj.wave2bin(peak.wavDoublet.max() + width *
                                    np.sqrt(peak.varDoublet)), 1.0, dtype=int)
    nowFlux = obj.reduced_flux[bounds]
    nowLeng = np.sqrt(np.dot(nowFlux, nowFlux))
    if obj.fiberid != 1:
        objPre = SDSSObject(obj.plate, obj.mjd, obj.fiberid - 1,
                            obj.dataVersion, obj.baseDir)
        preFlux = objPre.reduced_flux[bounds]
        preLeng = np.sqrt(np.dot(preFlux, preFlux))
        preProd = np.dot(preFlux, nowFlux) / (preLeng * nowLeng)
    else:
        preProd = 0.0
    if obj.fiberid != 1000:
        objNxt = SDSSObject(obj.plate, obj.mjd, obj.fiberid + 1,
                            obj.dataVersion, obj.baseDir)
        nxtFlux = objNxt.reduced_flux[bounds]
        nxtLeng = np.sqrt(np.dot(nxtFlux, nxtFlux))
        nxtProd = np.dot(nxtFlux, nowFlux) / (nxtLeng * nowLeng)
    else:
        nxtProd = 0.0
    return preProd, nxtProd


def fitcSpec(obj, peak, width=2.0):
    bounds = np.arange(obj.wave2bin(peak.wavDoublet.min() - width *
                                    np.sqrt(peak.varDoublet)),
                       obj.wave2bin(peak.wavDoublet.max() + width *
                                    np.sqrt(peak.varDoublet)), 1.0, dtype=int)
    initP = [peak.ampDoublet[0], peak.varDoublet, peak.ampDoublet[1],
             peak.wavDoublet[3], peak.wavDoublet[4]]
    limP = [(0.1, 5.0), (1.0, 8.0), (0.1, 5.0),
            (peak.ampDoublet[0] - width * np.sqrt(peak.varDoublet),
             peak.ampDoublet[0] + width * np.sqrt(peak.varDoublet)),
            (peak.ampDoublet[1] - width * np.sqrt(peak.varDoublet),
             peak.ampDoublet[1] - width * np.sqrt(peak.varDoublet))]
    if obj.fiberid != 1:
        objPre = SDSSObject(obj.plate, obj.mjd, obj.fiberid - 1,
                            obj.dataVersion, obj.baseDir)
        resp, preChi2 = objPre.doubletFit(bounds, initP, limP)
    else:
        preChi2 = 1000.0
    if obj.fiberid != 1000:
        objNxt = SDSSObject(obj.plate, obj.mjd, obj.fiberid + 1,
                            obj.dataVersion, obj.baseDir)
        resp, nxtChi2 = objNxt.doubletFit(bounds, initP, limP)
    else:
        nxtChi2 = 1000.0
    return [preChi2, nxtChi2]
