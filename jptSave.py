def jptSave():
    # TODO: complete the function/parameters/returns
    '''
    k = 0
    for peak in peak_candidates:
        fileJ = open(topdir + savedir +  '/candidates_Jackpot.txt','a')
        fileJ.write('\n' + str([RA[i], DEC[i], int(plate), int(mjd), fiberid[i], z[i], peak[2],peak[3],peak[4],peak[5],peak[6],spectroflux[i,1], spectroflux[i,3]]))
        fileJ.close()
        k += 1
    plot_Jackpot(RA= RA[i],DEC=DEC[i],plate =int(plate), mjd=int(mjd), fiberid=fiberid[i], z=z[i],wave=wave, flux =flux[i,:], synflux = synflux[i,:],topdir=topdir,savedir=savedir ,peak = peak, show = plot_show, counter = k)
    '''


def plot_Jackpot(obj, peak, savedir, counter):
    # TODO: complete the function
    '''
    em_lines = np.array([3726.5,4861.325,4958.911,5006.843,6562.801])
    if show==False:
        mpl.use('Agg')
    fontP = FontProperties()
    fontP.set_size('medium')
    plt.suptitle(SDSSname(RA,DEC)+'\n'+'RA='+str(RA)+', Dec='+str(DEC) +', $z_{QSO}='+'{:03.3}'.format(z)+ '$')

    gs = gridspec.GridSpec(1,4)
    p1 = plt.subplot(gs[0,:4])

    smoothed_flux = np.array([np.mean(flux[ii-2:ii+3]) for ii in range(len(flux)) if (ii>4 and ii<len(flux)-4)])

    p1.plot(wave[5:-4], smoothed_flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
    #p1.plot(wave,  flux, 'k', label = 'BOSS Flux')
    p1.plot(wave, synflux, 'r', label = 'PCA fit')
    box = p1.get_position()
    p1.set_position([box.x0,box.y0+0.02,box.width*0.9,box.height])
    p1.set_ylim(np.min(synflux)-3, np.max(synflux)+3)
    p1.vlines(x = em_lines*(1+peak[2]),ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
    p1.vlines(x = em_lines*(1+peak[3]),ymin= -100,ymax= 100,colors= 'b',linestyles='dashed')

    p1.legend(loc='upper right', bbox_to_anchor = (1.2,1), ncol = 1, prop=fontP)
    p1.set_xlim(3500,10500)
    plt.ylabel('Flux [$10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')


    make_sure_path_exists(topdir + savedir +'/plots/')
    #plt.show()
    plt.savefig(topdir + savedir +'/plots/'+SDSSname(RA,DEC)+ '-' + str(plate) + '-' + str(mjd) + '-' + str(fiberid) + '-'+str(counter) +'.png')
    plt.close()
    '''
