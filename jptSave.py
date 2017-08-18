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

def jptSave(obj, peak_candidates, savedir, em_lines):
    '''
    jptSave.jptSave(obj, peak_candidates)
    ==========================================
    Saves Jackpot lens candidates (data+plots)

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak_candidates: The inquired peaks
        savedir: Directory to save the plots/data
        em_lines: Rest frame ELG emission lines
    Returns:
        - Nothing. Checks candidates and save data/plots if OK.
    '''
    k = 0

    fileJpt = open(os.path.join(savedir, 'candidates_Jackpot.txt'), 'a')
    for peak in peak_candidates:
        fileJpt.write('\n' + str(obj.RA) + " " + str(obj.DEC) +
            " " + str(obj.plate) + " " + str(obj.mjd) + " " +
            str(obj.fiberid)+ " " + str(obj.spectroflux[1]) + " " +
            str(obj.spectroflux[3]) + " " + str(obj.z)+ " " + str(obj.rchi2) + 
            " " + str(peak.wavelength_1) + " " + str(peak.sn_1) + 
            " " + str(peak.total_sn_1) + " " + str(peak.z_1) +
            " " + str(peak.wavelength_2) + " " + str(peak.sn_2) + 
            " " + str(peak.total_sn_2) + " " + str(peak.z_2)) 

        plot_Jackpot(obj,peak,em_lines,savedir,k)
        k += 1
    fileJpt.close()
    #plot_Jackpot(RA= RA[i],DEC=DEC[i],plate =int(plate), mjd=int(mjd), fiberid=fiberid[i], z=z[i],wave=wave, flux =flux[i,:], synflux = synflux[i,:],topdir=topdir,savedir=savedir ,peak = peak, show = plot_show, counter = k)


def plot_Jackpot(obj, peak,em_lines, savedir, counter):
    '''
    jptSave.plot_Jackpot(obj, peak,em_lines, savedir, counter)
    =========================================================
    Plots Jackpot lens candidates

    Parameters:
        obj: The SDSS object/spectra on which applied the subtraction
        peak_candidates: The inquired peaks
        savedir: Directory to save the plots/data
        em_lines: Rest frame ELG emission lines
        counter: To keep track of each candidates per spectra
    Returns:
        - Nothing. Create and save the plot.
    '''
    
    fontP = FontProperties()
    fontP.set_size('medium')
    plt.suptitle(SDSSname(obj.RA,obj.DEC)+'\n'+'RA='+str(obj.RA)+
        ', Dec='+str(obj.DEC) +', $z_{QSO}='+'{:03.3}'.format(obj.z)+ '$')

    gs = gridspec.GridSpec(1,4)
    p1 = plt.subplot(gs[0,:4])

    smoothed_flux = np.array([np.mean(obj.flux[ii-2:ii+3]) 
        for ii in range(len(obj.flux)) if (ii>4 and ii<len(obj.flux)-4)])

    p1.plot(obj.wave[5:-4], smoothed_flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
    #p1.plot(wave,  flux, 'k', label = 'BOSS Flux')
    p1.plot(obj.wave, obj.synflux, 'r', label = 'PCA fit')
    box = p1.get_position()
    p1.set_position([box.x0,box.y0+0.02,box.width*0.9,box.height])
    p1.set_ylim(np.min(obj.synflux)-3, np.max(obj.synflux)+3)
    p1.vlines(x = em_lines*(1+peak.z_1),ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
    p1.vlines(x = em_lines*(1+peak.z_2),ymin= -100,ymax= 100,colors= 'b',linestyles='dashed')

    p1.legend(loc='upper right', bbox_to_anchor = (1.2,1), ncol = 1, prop=fontP)
    p1.set_xlim(3500,10500)
    plt.ylabel('Flux [$10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')

    make_sure_path_exists(savedir +'/plots/')

    plt.savefig(savedir +'/plots/'+SDSSname(obj.RA,obj.DEC)+ '-' + str(obj.plate) 
        + '-' + str(obj.mjd) + '-' + str(obj.fiberid) + '-'+str(counter) +'.png')
    plt.close()

