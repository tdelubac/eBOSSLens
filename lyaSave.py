import os
import numpy as np
from SDSSObject import SDSSObject
# Matplotlib trick
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
from utils import SDSSname, make_sure_path_exists

def foregroundELG(obj, peak, em_lines): 
    '''
    lyaSave.foregroundELG(obj, peak, em_lines)
    ==========================================
    Test if LAE emission cannot be attributed to foreground OII

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak: The inquired peak
        em_lines: The rest-frame wavelength of ELG emission lines
    Returns:
        SNlines: the squared sum SN of the foreground em_lines
    '''
    x0 = peak.wavelength
    z_O2 = x0/3726.5 - 1.0
    # compute SN of Ha, Hb, OIII to test if LyA is not OII at low-z
    SNlines = 0
    for l in em_lines:
        center_bin = obj.wave2bin(l*(1+z_O2))
        SNlines += max(obj.SN[center_bin-2:center_bin+2])**2
    return np.sqrt(SNlines)

def lyaFlux_eqWidth(obj, peak):
    '''
    lyaSave.lyaFlux_eqWidth(obj, peak)
    ====================================
    Compute apparent LyA flux and eq_width

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak: The inquired peak
    Returns:
        - Nothing. Updates the peak attributes.
    '''
    x0 = peak.wavelength
    bounds = np.linspace(obj.wave2bin(x0)-15,obj.wave2bin(x0)+15,31,dtype = np.int16)
    #compute equivalent width 
    dwave = np.array([obj.wave[gen_i+1]-obj.wave[gen_i] 
        if gen_i<len(obj.wave)-1 else 0 for gen_i in range(len(obj.wave))])
    eq_Width = np.sum((obj.flux[bounds]/obj.synflux[bounds]-1)*dwave[bounds])
    # compute LyA flux
    temp_fluxes = np.zeros(5)
    for j in range(4,9):
        temp_bounds = np.linspace(obj.wave2bin(x0)-j,obj.wave2bin(x0)+j,2*j+1,dtype = np.int16)
        temp_fluxes[j-4] = np.sum((obj.flux[temp_bounds]-obj.synflux[temp_bounds])*dwave[temp_bounds])
    flux = np.median(temp_fluxes)

    peak.eq_Width = eq_Width
    peak.flux =  flux

def lyaSkewness(obj,peak):
    '''
    lyaSave.lyaSkewness(obj, peak)
    ====================================
    Compute skewness indicators (S and a_lambda)

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak: The inquired peak
    Returns:
        - Nothing. Updates the peak attributes.
    '''
    x0 = peak.wavelength
    bounds = np.linspace(obj.wave2bin(x0)-15,obj.wave2bin(x0)+15,31,dtype = np.int16)

    # Skewness indicator (3rd order moment)
    I = np.sum(obj.reduced_flux[bounds])
    xmean = np.sum(obj.reduced_flux[bounds]*obj.wave[bounds])/I
    sigma2 = np.sum(obj.reduced_flux[bounds]*(obj.wave[bounds] - xmean)**2)/I
    S = np.sum(obj.reduced_flux[bounds]*(obj.wave[bounds] - xmean)**3)/(I*np.sign(sigma2)*(np.abs(sigma2)**1.5))

    # A_lambda (Rhoads et al. 2003)
    local_wave = obj.wave[bounds]
    local_flux = obj.reduced_flux[bounds]
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
    
    a_lambda = (l_red_10-local_wave[peak_index])/(local_wave[peak_index]-l_blue_10)

    peak.l_blue_10 = l_blue_10
    peak.l_red_10 = l_red_10
    peak.aLambda = a_lambda    

def lyaSave(obj, peak_candidates,savedir,em_lines, threshold_SN, QSOlens, paper_mode):
    '''
    lyaSave.lyaSave(obj, peak)
    =============================
    Performs specific LAE contamination test and saves the detections

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak_candidates: The peak candidates on this object
        QSOlens: Boolean if foreground lens is QSO or not
        savedir: Directory to save the plots and peak informations
        em_lines: The rest-frame wavelength of ELG emission lines
        threshold_SN: To reject LyA in favor of foreground OII 
        QSOlens: If foreground object is QSO
        paper_mode: Simple plots or not. See plotQSOLAE()
    Returns:
        Nothing.
    '''
    if QSOlens:
        fileLyA = open(os.path.join(savedir, 'candidates_QSOLAE.txt'), 'a')
    else:
        fileLyA = open(os.path.join(savedir, 'candidates_GALLAE.txt'), 'a')

    n = 0   
    for peak in peak_candidates:
        # Test if Foreground is likely or not
        if (foregroundELG(obj,peak,em_lines) > peak.sn + threshold_SN):
            raise Exception('Rejected as attributed to foreground ELG (OII)')
        else:
            # Compute apparent LyA flux and eq_width
            lyaFlux_eqWidth(obj,peak)
            # Compute skewness indicators
            lyaSkewness(obj,peak)
            # Save and plot
            plot_QSOLAE(obj,peak,n, QSOlens, savedir, paper_mode)
            n+=1
            # RA DEC plate mjd fiber 
            fileLyA.write('\n' + str(obj.RA) + " " + str(obj.DEC) +
                    " " + str(obj.plate) + " " + str(obj.mjd) + " " +
                    str(obj.fiberid) + str(obj.spectroflux[1]) + " " +
                    str(obj.spectroflux[3]) + " " + str(obj.z)+ " " +str(obj.rchi2) + 
                    " " + str(peak.wavelength) + " " + str(peak.sn) + 
                    " " + str(peak.reduced_sn) + " " + str(peak.redshift) 
                    + " " + str(peak.eq_Width) + " " + str(peak.flux) +
                    " " + str(peak.l_blue_10) + " " + str(peak.l_red_10) + 
                    " " + str(peak.aLambda) + " " + str(peak.skewness)) 
    fileLyA.close()

def plot_QSOLAE(obj,peak, n_peak, QSOlens, savedir, paper_mode= True):
    '''
    lyaSave.plot_QSOLAE(obj,peak, n_peak, QSOlens, savedir, paper_mode= True)
    =============================
    Plot the detection for QSOLAE case

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak: Detected LyA peak
        savedir: Directory to save plot
        n_peak: Number of peak for this spectrum
        QSOlens: If foreground lens is QSO or not
        paper_mode: True:Produces simple plots. False:more detailed
    Returns:
        Nothing. 
    '''
    
    make_sure_path_exists(savedir +'/plots/')

    # Create and save graph
    fontP = FontProperties()
    fontP.set_size('medium')
    plt.figure(figsize=(12,3))
    plt.suptitle(SDSSname(obj.RA,obj.DEC)+'\n'+'RA='+str(obj.RA)+
        ', Dec='+str(obj.DEC) +', $z_{QSO}='+'{:03.3}'.format(obj.z)+ '$')
    plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}$')

    x0 = peak.wavelength

    if paper_mode:
        gs = gridspec.GridSpec(1,3)

        smoothed_flux = np.array([np.mean(obj.flux[ii-2:ii+3]) for ii 
            in range(len(obj.flux)) if (ii>4 and ii<len(obj.flux)-4)])

        p1 = plt.subplot(gs[0,:2])
        #p1.plot(wave,  flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
        p1.plot(obj.wave[5:-4], smoothed_flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
        p1.plot(obj.wave, obj.synflux, 'r', label = 'PCA fit', drawstyle='steps-mid')
        p1.set_ylim(np.min(obj.synflux)-3, np.max(obj.synflux)+3)

        p1.vlines(x = peak.wavelength,ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
        box = p1.get_position()
        p1.set_position([box.x0,box.y0+0.06,box.width,box.height*0.85])
        plt.ylabel('Flux [$10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}$]')
        plt.xlabel('Observed wavelength [$\AA$]')
        p2 = plt.subplot(gs[0,2:3])

        p2.plot(obj.wave, obj.flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
        p2.plot(obj.wave, obj.synflux, 'r', label = 'PCA fit', drawstyle='steps-mid')

        window = np.linspace(obj.wave2bin(x0)-20,obj.wave2bin(x0)+20,41,dtype = np.int16)
        p2.set_ylim(np.min(obj.flux[window]), np.max(obj.flux[window])+0.5)
        p2.legend(loc='upper right', bbox_to_anchor = (1.3,1.1), ncol = 1, prop=fontP)
        box = p2.get_position()
        p2.set_position([box.x0,box.y0+0.06,box.width*0.9,box.height*0.85])
        x1 = int(x0/10.)*10
        plt.xticks([x1-10,x1,x1+10,x1+20])
        p2.set_xlim(x0-15,x0+25)
        plt.xlabel('Observed wavelength [$\AA$]')

    else:

        gs = gridspec.GridSpec(2,2)
        p1 = plt.subplot(gs[0,:2])
        p1.plot(obj.wave,  obj.flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
        p1.plot(obj.wave, obj.synflux, 'r', label = 'PCA fit', drawstyle='steps-mid')

        ''' Structure to make masks appear
        p1.fill_between(obj.wave,np.min(obj.synflux)-10,np.max(obj.synflux)+10, \
        where=np.logical_and(6348<wave,wave<6378),facecolor='c', alpha=0.2)
        '''    
        p1.legend(loc='upper right', bbox_to_anchor = (1,1), ncol = 1, prop=fontP)
        box = p1.get_position()
        p1.set_position([box.x0,box.y0,box.width,box.height])

        p1.set_ylim(np.min(obj.synflux)-3, np.max(obj.synflux)+3)
        plt.ylabel('$f_{\lambda}\, [10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')
        if QSOlens == False:
            p1.set_xlim(3600,6000)

        window = np.linspace(obj.wave2bin(x0)-40,obj.wave2bin(x0)+40,81,dtype = np.int16)
        if QSOlens:
            p3 = plt.subplot(gs[1,:1])
            p3.plot(obj.wave,  obj.flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
            p3.plot(obj.wave, obj.synflux, 'r', label = 'PCA fit', drawstyle='steps-mid')
            p3.set_xlim(np.min(obj.wave[window]),np.max(obj.wave[window]))
            p3.set_ylim(np.min(obj.synflux[window])-1, np.max(obj.flux[window])+1)
            box = p3.get_position()

            p3.set_position([box.x0,box.y0,box.width,box.height])
            plt.ylabel('$f_{\lambda}\, [10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')
            p3.legend(loc='upper right', bbox_to_anchor = (1,1), ncol = 1,prop=fontP)

            p3.locator_params(axis='x',nbins=6)

            median_local = np.median(obj.reduced_flux[window])
            fit_QSO = np.poly1d(np.polyfit(x=obj.wave[window],y=obj.reduced_flux[window],deg=3, \
                w=(np.abs(obj.reduced_flux[window]-median_local)<5)*np.sqrt(obj.ivar[window])) )
            p4 = plt.subplot(gs[1,1:2])
            p4.plot(obj.wave[window], fit_QSO(obj.wave[window]), '-m',label = 'Order 3 fit')
            box = p4.get_position()
            p4.set_xlim(np.min(obj.wave[window]),np.max(obj.wave[window]))

            p4.set_position([box.x0,box.y0,box.width,box.height])
            p4.plot(obj.wave[window], obj.reduced_flux_QSO[window],'k', label = 'Reduced flux', drawstyle='steps-mid')
            p4.legend(loc='upper right', bbox_to_anchor = (1,1), ncol = 1,prop=fontP)

            p4.locator_params(axis='x',nbins=6)
        else:
            p3 = plt.subplot(gs[1,:2])
            p3.plot(obj.wave,  obj.flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
            p3.plot(obj.wave, obj.synflux, 'r', label = 'PCA fit', drawstyle='steps-mid')
            p3.legend(prop=fontP)
            p3.set_xlim(x0-50,x0+60)
            p3.set_ylim(np.min(obj.synflux[bounds])-2, np.max(obj.flux[bounds])+3)
            plt.ylabel('$f_{\lambda}\, [10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$', fontsize=18)

    plt.savefig(savedir +'/plots/'+SDSSname(obj.RA,obj.DEC)+ '-' + str(obj.plate) + '-' 
        + str(obj.mjd) + '-' + str(obj.fiberid) + '-' + str(n_peak)+ '.eps', format = 'eps', dpi = 2000)

    plt.close()
