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

    # TODO: complete the function/parameters/returns
    '''
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


#def plot_QSOLAE(RA,DEC,z,flux,wave,synflux,x0,ivar, reduced_flux,window,peak,params,params_skew, topdir, savedir, n_peak, plate, mjd, fiberid,c0,c1,Nmax, show = False, paper=True, QSOlens = True):
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
    
    #if show ==False:
    #    mpl.use('Agg')
    
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

        ##### Old code to plot skew fit
        #p2 = plt.subplot(gs[2,:2])
        #if QSOlens:
            #p2.plot(wave[window], reduced_flux[window]-fit_QSO(wave[window]),'k', label = 'Reduced flux', drawstyle='steps-mid')
        #else:
            #p2.plot(wave[window], reduced_flux[window],'k', label = 'Reduced flux')
        #if 0.0<peak[16]<peak[15]:
            #p2.plot(wave,gauss2(x=wave,x1=params[0],x2=params[1],A1=params[2],A2=params[3],var=params[4]),'g', label = r'$\chi_D^2 = $' + '{:.4}'.format(peak[16]))
        #else:
            #p2.plot(wave,gauss(x=wave, x_0=params[0], A=params[1], var=params[2]),'r', label = r'$\chi_G^2 = $' + '{:.4}'.format(peak[15]) )
        #if 0.0<peak[17]<peak[18]:
            #p2.plot(wave,skew(x=wave,A = params_skew[0], w=params_skew[1], a=params_skew[2], eps=params_skew[3]), 'b', label =r'$\chi_S^2 = $' + '{:.4}'.format(peak[17]))
        #else:
            #p2.plot(wave,skew2(x=wave,A1 = params_skew[0], w1=params_skew[1], a1=params_skew[2], eps1 = params_skew[3], A2 = params_skew[4], w2=params_skew[5], a2=params_skew[6], eps2=params_skew[7]), 'c',label= r'$\chi_{S2}^2 = $' + '{:.4}'.format(peak[18]))
        #box = p2.get_position()
        #p2.set_position([box.x0,box.y0,box.width*0.9,box.height])
        #p2.legend(loc='upper right', bbox_to_anchor = (1.2,1), ncol = 1,prop=fontP)
        #plt.xlabel('$Wavelength\, [\AA]$',fontsize = 18)
        #p2.set_xlim(np.min(wave[window]),np.max(wave[window]))
        #if QSOlens:
            #p2.set_ylim(-1, np.max(reduced_flux[window]-fit_QSO(wave[window])+1))

    plt.savefig(savedir +'/plots/'+SDSSname(obj.RA,obj.DEC)+ '-' + str(obj.plate) + '-' 
        + str(obj.mjd) + '-' + str(obj.fiberid) + '-' + str(n_peak)+ '.eps', format = 'eps', dpi = 2000)

    plt.close()
