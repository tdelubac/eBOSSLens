# Last modification: Xinlun Cheng, 2017
# Original author: Romain A. Meyer, 2016, LASTRO, EPFL
# Spectroscopic lensing systems detection tool
# Detects OII doublet for galaxy-galaxy lensing systems
# Searches for Lyman alpha emission for galaxy-LAE or QSO-LAE systems
# 'Paper' mode removes intermediate steps plots useful for analysis but
# irrelevant for publication


# show graphs or not? (do not use show on PC clusters!!!)
# Imports
import argparse
import itertools as it
import numpy as n
from scipy.optimize import minimize
from SDSSObject import SDSSObject
from peakFinder import bolEstm, peakCandidate, combNear, checkFore, nearLine
from utils import gauss, kernel


# Parse argument from command line
parser = argparse.ArgumentParser()
parser.add_argument("chi2", type=float)
parser.add_argument("width", type=float)
parser.add_argument("sig", type=float)
parser.add_argument("sdir", type=str)
parser.add_argument("dataVersion", type=str)
args = parser.parse_args()
# Never show plot
plot_show = False
# Operation mode
searchLyA = False
QSOlens = False
paper = True
Jackpot = False
# Data version
dataVersion = args.dataVersion
# Max chi2 for gaussian/ doublet fitting OII
setWidth = args.width
setSig = args.sig
max_chi2 = args.chi2
# Set topdir,savedir:
topdir = '..'
savedir = args.sdir
baseDir = "/SCRATCH"
# Give file in [plate mjd] format of plates you want to inspect
plates_list = 'list_QSOGal.txt'
# Waves for mask
wMask = n.array([[5570.0, 5590.0], [5880.0, 5905.0], [6285.0, 6315.0],
                 [6348.0, 6378.0]])
# Set of emission lines used for lensed galaxy detection: OII, Hb, OIII,
# OIII, Ha
em_lines = n.array([3726.5, 4861.325, 4958.911, 5006.843, 6562.801])
# Typical mask width
l_width = 15
# ------------------------------------------------------------------------------
# ------------------------- INITIALIZATION -------------------------------------
plate_mjd = [line.strip().split() for line in
             open(topdir + savedir + plates_list)]
l_LyA = 1215.668

# Needed for QSO lenses (QSO in DR12 with proper classificaiton)
# DR12Q = DR12Q_extractor(path='./Superset_DR12Q.fits')

# Counters used to check the numbers of candidates at key step
goodObjs = 0
counter2 = 0
counter3 = 0
counter4 = 0

# Loop over plates
for j in n.arange(len(plate_mjd)):
    # Loading the data
    flux = 0
    mjd = plate_mjd[j][1]
    plate = plate_mjd[j][0]
    obj = SDSSObject(mjd, plate, dataVersion, baseDir=baseDir)
    # Mask BOSS spectra glitches + Sky
    obj.mask(wMask)
    # Loop over fibers
    for i in n.arange(0, len(flux[:, 0])):
        # TODO: clean up
        if QSOlens:
            '''
            redshift_warning = False
            # Take only QSOs
            if obj.obj_type[i].startswith('sky') or \
                    obj.obj_type[i].startswith('SPECTROPHOTO_STD') or \
                    (not obj.obj_type[i] == 'QSO   '):
                continue
            # Reject badly fitted QSOs
            if obj.zwarning[i] != 0 or obj.rchi2[i] > 10 or obj.z_err[i] < 0:
                continue
            # ?
            index = [x for x in DR12Q if int(x[0]) == int(plate) and
                     int(x[1]) == int(mjd) and int(x[2]) == int(fiberid[i])]
            if len(index) > 0:
                if n.abs(index[0][3] - obj.z[i]) < 0.005 and \
                        n.abs(index[0][4] - obj.z[i]) > 0.1:
                    continue
            else:
                redshift_warning = True

            # Before masking, compute the FWHM of CIV or HBeta depending on redshift:
            FWHM, l_times_luminosity, HB_wave, params_beta, line_coeff = \
                QSO_compute_FWHM(ivar = ivar[i,:],flux = flux[i,:], wave = wave,c0=c0,c1=c1,Nmax=Nmax,z =z[i],l_width = l_width)

            M_BH = 10**(6.91 + n.log10(n.sqrt(5100*l_times_luminosity/1e44)*(FWHM/1000)**2)) # Masses solaires
            sigma_host = 200*10**((n.log10(M_BH) - 7.92)/3.93) ### km s-1


            ivar[i,:] = mask_QSO(ivar=ivar[i,:],z=z[i], l_width = l_width, c0 = c0 , c1=c1,Nmax=Nmax)
            '''
        # End TODO

        else:
            # ?
            if obj.obj_class[i] == 'STAR  ' or obj.obj_class[i] == 'QSO   ' or \
                    obj.obj_type[i] == 'SKY             ' or \
                    obj.obj_type[i] == 'SPECTROPHOTO_STD':
                continue
        # Good objects
        goodObjs = goodObjs + 1
        # Find peaks
        peaks = []
        doublet = None
        # Bolton 2004: S/N of maximum likelihood estimator of gaussian peaks
        if searchLyA:
            width = 30.0
            sig = 2.17
        else:
            width = setWidth
            sig = setSig
        NormGauss = gauss(n.linspace(-width * 0.5, width * 0.5, width), 0.0,
                          1.0, sig ** 2.0)
        NormGauss = NormGauss / n.sum(NormGauss)
        SN = bolEstm(i, obj, width, sig)
        '''
        if searchLyA and QSOlens:
            peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(obj.wave,SN) if (test > 8.0 and  (l_LyA*(1+obj.z[i])+300)<x0<9500)])
        elif searchLyA == False and QSOlens:
            peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test) for x0,test in zip(obj.wave,SN) if (test>6.0 and  (l_LyA*(1+obj.z[i])+300)<x0<9500)])
        elif searchLyA == True and QSOlens == False:
            peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(obj.wave,SN) if (test>8.0 and  3600 < x0 < 4800)])
        elif Jackpot == True:
            # x0 z1 z2 Quad_SN2 SN0->Quad_SN1  free free
            peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test) for x0,test in zip(obj.wave,SN) if test>8.0])
        elif not (searchLyA or QSOlens):
            peak_candidates = n.array([(x0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(obj.wave,SN) if test>6.0])
        else:
            continue
        '''
        # Filter out invalid sources
        if not (searchLyA or QSOlens):
            peak_candidates = n.array([peakCandidate(x0, test)
                for x0, test in zip(obj.wave, SN) if test > 6.0])
        # TODO: clean up
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
        else:
            continue
        # Keep the center
        peak_candidates = combNear(peak_candidates)
        # Check hits are not from foreground galaxy or badly fitted QSO
        if QSOlens and searchLyA:
            if checkFore(peak_candidates, em_lines, obj.z[i]):
                continue
        # Search for suitable peak candidates
        for peak in peak_candidates:
            x0 = peak.wavSinglet
            if nearLine(int(mjd), int(plate), i, x0, obj):
                continue
            x0Bin = obj.wave2bin(x0)
            bounds = n.linspace(x0Bin - 15, x0Bin + 15, 31, dtype=n.int16)
            # Fit QSO continuum and check if signal is reduced or not
            # i.e. check if line detection is produced by large features
            # TODO: clean up
            '''
            if QSOlens:
                window = n.linspace(x0Bin - 40, x0Bin + 40, 81, dtype=n.int16)
                median_local = n.median(obj.reduced_flux[i, window])
                weight = n.abs(obj.reduced_flux[i, window] - median_local) < 5
                fit_QSO = n.poly1d(n.polyfit(x=obj.wave[window],
                                             y=obj.reduced_flux[i, window],
                                             deg=3,
                                             w=weight *
                                             n.sqrt(obj.ivar[i, window])))
                new_flux = obj.reduced_flux[i, window] - \
                    fit_QSO(obj.wave[window])
                cj1_new = n.sum(new_flux * kernel(int(len(window) / 2), width,
                                                  NormGauss, len(new_flux)) *
                                obj.ivar[i, window])
                cj2_new = n.sum(obj.ivar[i, window] *
                                kernel(int(len(window) / 2), width, NormGauss,
                                       len(window)) ** 2.0)
                SN_fitted = cj1_new / n.sqrt(cj2_new)
                if SN_fitted < 6:
                    continue
                else:
                    peak.snFitted = SN_fitted
                    obj.reduced_flux[i, window] = new_flux
            # Special case: QSOlens with background galaxies
            if (not (searchLyA or Jackpot)) and QSOlens:
                for l in em_lines:
                    test_z = peak.wavelength / l - 1.0
                    if test_z > obj.z[i]:
                        quad_SN = 0.0
                        for w in em_lines:
                            center_bin = obj.wave2bin(w * (1.0 + test_z))
                            SN_line = n.array(SN[center_bin - 2:center_bin + 2])
                            quad_SN += max(SN_line * (SN_line > 0)) ** 2.0
                        quad_SN = n.sqrt(quad_SN)
                        if quad_SN > peak.sn:
                            peak.sn = quad_SN
                            peak.z = test_z
                continue
            # Special case: Jackpot
            mask_width_Jackpot = 50
            if Jackpot:
                first_lens = False
                peak[5] = peak.sn
                peak[4] = peak.sn
                for l in em_lines:
                    test_z = x0 / l - 1.0
                    if test_z > obj.z[i] + 0.05:
                        quad_SN_1 = peak.sn
                        for w in em_lines:
                            if w * (1.0 + test_z) < 9500:
                                center_bin = obj.wave2bin(w * (1.0 + test_z))
                                SN_line = n.array(
                                    SN[center_bin - 2:center_bin + 2]) * \
                                    (not(nearline(w*(1+test_z), zline, fiberid[i], z[i], int(mjd), int(plate), mask_width_Jackpot)))
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
                continue
            '''
            # Singlet
            if searchLyA:
                init = [x0, 4.0, 6.0]
                paramLim = [(x0 - 2.0, x0 + 2.0), (1.0, 100.0), (1.0, 15.0)]
            elif not (searchLyA or QSOlens):
                init = [x0, 1.0, 2.0]
                paramLim = [(x0 - 2.0, x0 + 2.0), (0.1, 5.0), (1.0, 8.0)]
            params, chisq = obj.singletFit(i, bounds, init, paramLim)
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
            #Doublet OII
            if (x0 > 3727.0 * (1.0 + obj.z[i]) or searchLyA and not QSOlens:
                params2, chisq2 = obj.doubletFit(i, bounds,
                                                 [1.0, 5.0, 1.0, x0 - 1.5, x0 + 1.5],
                                                 [(0.1, 5.0), (1.0, 8.0), (0.1, 5.0),
                                                  (x0 - 7.0, x0), (x0, x0 + 7.0)])
                if  (not searchLyA and
                     0.5 * x0 < abs(params2[3] - params2[4]) * 3726.5 < 4.0 * x0 and
                     not (chisq2 > max_chi2)):
                    peak.chiDoublet = chisq2
                    peak.ampDoublet = np.array([params2[0], params2[2]])
                    peak.varDoublet = params2[1]
                    peak.wavDoublet = np.array([params2[3], params2[4]])
                # TODO: clean up
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
            # TODO: clean up
            '''
            # If looking at LAE, test a skew-normal profile as well
            if searchLyA:
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
        counter2 = counter2 + 1;
        if not (searchLyA or QSOlens or Jackpot):
            # Compare singlet and doublet fit within each peakCandidate
            for k in range(len(peak_candidates)):
                peak = peak_candidates[k]
                peak.update()
                if peak.wavelength > 9200.0:
                    peak.setDoublet(False)
            # Removing candidates that were not fitted : params still 0
            peak_candidates = n.array([peak for peak in peak_candidates if (peak.chiDoublet == peak.chiSinglet != 1000.0)])
            if len(peak_candidates) == 0:
                continue
            counter3 = counter3 + 1;
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
            continue
        # Check that at least 1 candidate is below 9200 Angstrom cut, if not, go to next fiber
        below_9500 = False
        for peak in peak_candidates:
            if peak.wavelength < 9500.0:
                below_9500 = True
                break
        if not below_9500:
            continue
        counter4 = counter4 + 1;
        #Try to infer background redshift
        detection = False
        score = 0.0
        if not (searchLyA or QSOlens or Jackpot):
                if doublet:
                    fileD = open(topdir + savedir +  '/candidates_doublet.txt','a')
                    z_s = peak_candidates[doublet_index].wavelength / 3727.24 - 1.0
                    if (z_s > obj.z[i]+0.05):
                        detection = True
                        score = peak_candidates[doublet_index].chi
                        fileD.write('\n' +str([radEinstein(z[i],z_s,vdisp[i]*1000), score,z_s, RA[i], DEC[i], int(plate), int(mjd), fiberid[i],peak_candidates[doublet_index][9]]))
                    fileD.close()
                    if len(peak_candidates):
                        fileDM = open(topdir + savedir + '/candidates_DM.txt','a')
                        confirmed_lines = []
                        # Generating all combinations of lines from above list to compare with candidates
                        temp = [peak for peak in peak_candidates if peak[1]!=peak[5]]
                        compare = em_lines[1:5]
                        if z_s > obj.z[i] + 0.05 :
                            for peak in temp:
                                for line in compare:
                                    if (abs(peak[0]/line -1 - z_s) < 0.01):
                                        detection = True
                                        confirmed_lines.append(line)
                                        score+= peak[5]
                        if (confirmed_lines != []):
                            fileDM.write('\n'+str([radEinstein(z[i],z_s,vdisp[i]*1000), score,z_s, RA[i], DEC[i], int(plate), int(mjd), fiberid[i],confirmed_lines]))
                        fileDM.close()
                elif len(peak_candidates) > 1:
                    compare = it.combinations(em_lines,len(peak_candidates))
                    confirmed_lines = []
                    fileM = open(topdir + savedir +'/candidates_multi.txt','a')
                    for group in compare:
                        for k in range(len(peak_candidates)):
                            for j in range(k+1,len(peak_candidates)):
                                if ( abs(peak_candidates[k][0]/group[k] - peak_candidates[j][0]/group[j]) < 0.01 and peak_candidates[k][0]/group[k]-1.0 > (z[i] + 0.05) ):
                                    detection = True
                                    z_s = peak_candidates[k][0]/group[k]-1.0
                                    confirmed_lines.append([group, peak_candidates[k][0]/group[k]-1.0])
                                    score+= peak_candidates[j][4]**2+peak_candidates[k][4]**2
                    if (confirmed_lines != []):
                        fileM.write('\n'+str([radEinstein(z[i],z_s,vdisp[i]*1000), score, z_s, RA[i], DEC[i], int(plate), int(mjd), fiberid[i],confirmed_lines]))
                    fileM.close()
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


                fileQSO = open(topdir + savedir +  '/candidates_QSO.txt','a')
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
        if searchLyA == False and QSOlens == False and Jackpot == False:
            for k in range(len(peak_candidates)):
                peak = peak_candidates[k]
                if (k == doublet_index and doublet ==True):
                    # Save doublet gaussians
                    peaks.append([peak[9], peak[6], peak[8]])
                    peaks.append([peak[10], peak[7], peak[8]])
                else:
                    # Save emission line
                    peaks.append([peak[0], peak[2], peak[3]])
        peak_number = len(peak_candidates)
        #Graphs OII doublet
        if ((peak_number>1 or doublet==True) and below_9500 and detection and searchLyA==False and QSOlens==False and Jackpot == False):

            #Computing total fit of all peaks
            fit=0
            for k in n.arange(len(peaks)):
                fit = fit + gauss(wave, x_0 = peaks[k][0] , A=peaks[k][1], var=peaks[k][2])

            plot_GalaxyLens(doublet = doublet,RA = RA[i],DEC = DEC[i],plate = plate,fiberid=fiberid[i],mjd=mjd,z=z[i],z_err =z_err[i], \
                obj_class = obj_class[i],wave = wave, reduced_flux=reduced_flux[i,:], zline = zline, fit = fit,topdir = topdir, savedir = savedir, peak_candidates =peak_candidates, \
                doublet_index = doublet_index, c0 = c0,c1 = c1,Nmax = Nmax, show = plot_show)
