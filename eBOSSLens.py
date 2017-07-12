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
from resultSave import galSave
from utils import gauss
from utils_Gal import plot_GalaxyLens


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
                 [6348.0, 6378.0]])
# Save dir
savedir = "../Meyer2016"


def eBOSSLens(plate, mjd, fiberid, searchLyA, QSOlens, Jackpot, max_chi2=4.0,
              wMask=wMask, em_lines=em_lines, bwidth=30.0, bsig=1.2):
    obj = SDSSObject(plate, mjd, fiberid, "v5_7_0", "../SCRATCH")
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
        bounds = n.linspace(x0Bin - 15, x0Bin + 15, 31, dtype=n.int16)
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
            #eq_Width = quad(gauss,x0-200,x0+200,args=(params[0],params[1],params[2]))
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
            pk.update()
            if pk.wavelength > 9200.0:
                pk.setDoublet(False)
        # Removing candidates that were not fitted
        peak_candidates = n.array([peak for peak in peak_candidates if
                                   (peak.chi != 1000.0)])
        if len(peak_candidates) == 0:
            raise Exception("Rejected since no candidates")
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
    elif searchLyA == False and QSOlens == True and Jackpot == False:
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
    if len(peak_candidates) == 0:
        raise Exception("Rejected since all is lost")
    # Check that at least 1 candidate is below 9200 Angstrom cut
    below_9500 = False
    for peak in peak_candidates:
        if peak.wavelength < 9500.0:
            below_9500 = True
            break
    if not below_9500:
        raise Exception("Rejected since no below 9200")
    # Try to infer background redshift
    detection = False
    if not (searchLyA or QSOlens or Jackpot):
        detection = galSave(doublet, obj, peak_candidates, doublet_index,
                            savedir, em_lines)
    # TODO: clean up
    '''
    elif searchLyA == False and QSOlens == True and Jackpot==False:
        for k in range(len(peak_candidates)):
            z_backgal = peak_candidates[k][4]
            peak = peak_candidates[k]
            # compute OII,OIII flux (of lensed galaxy)
            temp_fluxes_OII = n.zeros(5)
            temp_fluxes_OIII = n.zeros(5)
            dwave = n.array([wave[gen_i+1]-wave[gen_i] if gen_i<len(wave)-1 else 0 for gen_i in range(len(wave))])
            if z_backgal < 1:
                for j in range(4,9):
                    temp_bounds = n.linspace(wave2bin((1+z_backgal)*3727,c0,c1,Nmax)-j,wave2bin((1+z_backgal)*3727,c0,c1,Nmax)+j,2*j+1,dtype = n.int16)
                    temp_fluxes_OII[j-4] = n.sum((flux[i,temp_bounds]-synflux[i,temp_bounds])*dwave[temp_bounds])
                    temp_bounds = n.linspace(wave2bin((1+z_backgal)*5007,c0,c1,Nmax)-j,wave2bin((1+z_backgal)*5007,c0,c1,Nmax)+j,2*j+1,dtype = n.int16)
                    temp_fluxes_OIII[j-4] = n.sum((flux[i,temp_bounds]-synflux[i,temp_bounds])*dwave[temp_bounds])
            OII_flux = n.median(temp_fluxes_OII)
            OIII_flux = n.median(temp_fluxes_OIII)
            fileQSO = open(savedir +  '/candidates_QSO.txt','a')
            fileQSO.write('\n' + str([RA[i], DEC[i], int(plate), int(mjd), fiberid[i], z[i], peak[0],peak[4],peak[5],peak[6],spectroflux[i,1], spectroflux[i,3], OII_flux, OIII_flux ]))
            fileQSO.close()
            plot_QSOGal(k=k,RA = RA[i],DEC= DEC[i],plate = int(plate), mjd = int(mjd), fiberid = fiberid[i],z=z[i], z_backgal= z_backgal,flux=flux[i,:],wave=wave,synflux=synflux[i,:],ivar= ivar[i,:], \
                reduced_flux = reduced_flux[i,:], c0=c0,c1=c1,Nmax=Nmax,show = plot_show,topdir=topdir, savedir=savedir, HB_wave = HB_wave , params_beta=params_beta, line_coeff =line_coeff)
            # plot with +-1 AA for paper
            if paper:
                plot_QSOGal(k=k+len(peak_candidates),RA = RA[i],DEC= DEC[i],plate = int(plate), mjd = int(mjd), fiberid = fiberid[i],z=z[i], z_backgal= z_backgal+0.0005,flux=flux[i,:],wave=wave,synflux=synflux[i,:],ivar= ivar[i,:], \
                reduced_flux = reduced_flux[i,:], c0=c0,c1=c1,Nmax=Nmax,show = plot_show,topdir=topdir, savedir=savedir, HB_wave = HB_wave , params_beta=params_beta, line_coeff =line_coeff)
                plot_QSOGal(k=k+len(peak_candidates)+1,RA = RA[i],DEC= DEC[i],plate = int(plate), mjd = int(mjd), fiberid = fiberid[i],z=z[i], z_backgal= z_backgal-0.0005,flux=flux[i,:],wave=wave,synflux=synflux[i,:],ivar= ivar[i,:], \
                reduced_flux = reduced_flux[i,:], c0=c0,c1=c1,Nmax=Nmax,show = plot_show,topdir=topdir, savedir=savedir, HB_wave = HB_wave , params_beta=params_beta, line_coeff =line_coeff)
    elif Jackpot == True:
        k = 0
        for peak in peak_candidates:
            fileJ = open(topdir + savedir +  '/candidates_Jackpot.txt','a')
            fileJ.write('\n' + str([RA[i], DEC[i], int(plate), int(mjd), fiberid[i], z[i], peak[2],peak[3],peak[4],peak[5],peak[6],spectroflux[i,1], spectroflux[i,3]]))
            fileJ.close()
            k += 1
            plot_Jackpot(RA= RA[i],DEC=DEC[i],plate =int(plate), mjd=int(mjd), fiberid=fiberid[i], z=z[i],wave=wave, flux =flux[i,:], synflux = synflux[i,:],topdir=topdir,savedir=savedir ,peak = peak, show = plot_show, counter = k)
    elif searchLyA==True and QSOlens==True and Jackpot == False:
        if QSOlens:
            fileLyA = open(topdir + savedir +  '/candidates_QSO_LyA.txt','a')
        else:
            fileLyA = open(topdir + savedir +  '/candidates_LAE.txt','a')
        n_peak = 1
        for peak in peak_candidates:
            x0 = peak[0]
            z_O2 = x0/3727.24 - 1.0
            # compute SN of Ha, Hb, OIII to conclude if LyA candidate is or not OII at low-redshift
            SNlines = 0
            em_lines = n.array([4861.325,5006.843,6562.801])
            for l in em_lines:
                center_bin = wave2bin(l*(1+z_O2),c0,c1,Nmax)
                SNlines += max(SN[center_bin-2:center_bin+2])**2
            SNlines = n.sqrt(SNlines)
            if SNlines < 100:
                if n.abs(peak[9]) > n.abs(peak[13]):
                    skewness = peak[9]
                else:
                    skewness = peak[13]
                params = peak[1:6]
                params_skew = peak[7:15]
                if not(n.any(params) and n.any(params_skew)):
                    continue
                bounds = n.linspace(wave2bin(x0,c0,c1,Nmax)-15,wave2bin(x0,c0,c1,Nmax)+15,31,dtype = n.int16)
                #compute equivalent width before manipulating the spectra
                dwave = n.array([wave[gen_i+1]-wave[gen_i] if gen_i<len(wave)-1 else 0 for gen_i in range(len(wave))])
                eq_Width = n.sum((flux[i,bounds]/synflux[i,bounds]-1)*dwave[bounds])
                # compute LyA flux
                temp_fluxes = n.zeros(5)
                for j in range(4,9):
                    temp_bounds = n.linspace(wave2bin(x0,c0,c1,Nmax)-j,wave2bin(x0,c0,c1,Nmax)+j,2*j+1,dtype = n.int16)
                    temp_fluxes[j-4] = n.sum((flux[i,temp_bounds]-synflux[i,temp_bounds])*dwave[temp_bounds])
                lyA_flux = n.median(temp_fluxes)
                #compute skewness indicator (on reduced flux but without 3rd order fit (which is meant for bad fit QSO and not present in final candidates)
                I = n.sum(reduced_flux[i,bounds])
                xmean = n.sum(reduced_flux[i,bounds]*wave[bounds])/I
                sigma2 = n.sum(reduced_flux[i,bounds]*(wave[bounds] - xmean)**2)/I
                S = n.sum(reduced_flux[i,bounds]*(wave[bounds] - xmean)**3)/(I*n.sign(sigma2)*n.abs(sigma2)**1.5)
                local_wave = wave[bounds]
                local_flux = reduced_flux[i,bounds]
                peak_index = local_flux.argmax()
                F = 0.1*local_flux[peak_index]
                #Find blue and red 10% peak flux
                f = 10*F
                k = peak_index
                while f > F and 0<k:
                    k = k-1
                    f = local_flux[k]
                a = (local_flux[k+1]-local_flux[k])/(local_wave[k+1]-local_wave[k])
                b = local_flux[k] - a*local_wave[k]
                l_blue_10 = (F-b)/a
                k = peak_index
                f = 10*F
                while f > F and k<len(local_flux)-1:
                    k = k+1
                    f = local_flux[k]
                a = (local_flux[k]-local_flux[k-1])/(local_wave[k]-local_wave[k-1])
                b = local_flux[k-1] - a*local_wave[k-1]
                l_red_10 = (F-b)/a
                Sw = S*(l_red_10-l_blue_10)
                a_lambda = (l_red_10-local_wave[peak_index])/(local_wave[peak_index]-l_blue_10)
                # save the parameters
                fileLyA.write('\n' + str([peak[0],peak[6],z[i], SNlines ,RA[i], DEC[i], int(plate), int(mjd), fiberid[i],params[0],params[1],params[2],params[3],params[4],
                                params_skew[0],params_skew[1],params_skew[2],params_skew[3],params_skew[4],params_skew[5],params_skew[6],params_skew[7],peak[15],peak[16],peak[17],peak[18], skewness, S, Sw, a_lambda,rchi2[i],peak[19],spectroflux[i,1], spectroflux[i,3],lyA_flux]))
                # Make the graph
                plot_QSOLAE(RA= RA[i],DEC = DEC[i],z=z[i],flux=flux[i,:],wave=wave,synflux=synflux[i,:],x0= x0, ivar = ivar[i,:], reduced_flux = reduced_flux[i,:],window=window,peak =peak,
                    params = params,params_skew=params_skew, topdir = topdir, savedir = savedir, n_peak = n_peak, plate = int(plate), mjd = int(mjd), fiberid = fiberid[i],c0=c0,c1=c1,Nmax=Nmax, show = plot_show, paper = paper, QSOlens = QSOlens)
            n_peak = n_peak +1
        fileLyA.close()
        '''
    # Save surviving candidates (Galaxy-Galaxy case)
    peaks = []
    if not (searchLyA or QSOlens or Jackpot):
        for k in range(len(peak_candidates)):
            peak = peak_candidates[k]
            if k == doublet_index and doublet:
                # Save doublet gaussians
                peaks.append([peak.wavDoublet[0], peak.ampDoublet[0],
                              peak.varDoublet])
                peaks.append([peak.wavDoublet[1], peak.ampDoublet[1],
                              peak.varDoublet])
            else:
                # Save emission line
                peaks.append([peak.wavSinglet, peak.ampSinglet,
                              peak.varSinglet])
        peak_number = len(peak_candidates)
        # Graphs OII doublet
        if (peak_number > 1 or doublet) and below_9500 and detection:
            # Computing total fit of all peaks
            fit = 0.0
            for k in n.arange(len(peaks)):
                fit = fit + gauss(obj.wave, x_0=peaks[k][0], A=peaks[k][1],
                                  var=peaks[k][2])
            plot_GalaxyLens(doublet, obj, savedir, peak_candidates,
                            doublet_index, fit)
