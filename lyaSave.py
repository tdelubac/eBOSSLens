def lyaSave():
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


def plot_QSOLAE(RA,DEC,z,flux,wave,synflux,x0,ivar, reduced_flux,window,peak,params,params_skew, topdir, savedir, n_peak, plate, mjd, fiberid,c0,c1,Nmax, show = False, paper=True, QSOlens = True):
    # TODO: clean up
    '''
    if show ==False:
        mpl.use('Agg')
    # Create and save graph
    fontP = FontProperties()
    fontP.set_size('medium')
    plt.figure(figsize=(12,3))
    plt.suptitle(SDSSname(RA,DEC)+'\n'+'RA='+str(RA)+', Dec='+str(DEC) +', $z_{QSO}='+'{:03.3}'.format(z)+ '$')
    plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}$')

    if paper:
        gs = gridspec.GridSpec(1,3)

        smoothed_flux = np.array([np.mean(flux[ii-2:ii+3]) for ii in range(len(flux)) if (ii>4 and ii<len(flux)-4)])

        p1 = plt.subplot(gs[0,:2])
        #p1.plot(wave,  flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
        p1.plot(wave[5:-4], smoothed_flux, 'k', label = 'eBOSS Flux', drawstyle='steps-mid')
        p1.plot(wave, synflux, 'r', label = 'PCA fit', drawstyle='steps-mid')
        p1.set_ylim(np.min(synflux)-3, np.max(synflux)+3)

        p1.vlines(x = x0,ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
        box = p1.get_position()
        p1.set_position([box.x0,box.y0+0.06,box.width,box.height*0.85])
        plt.ylabel('Flux [$10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}$]')
        plt.xlabel('Observed wavelength [$\AA$]')
        p2 = plt.subplot(gs[0,2:3])

        p2.plot(wave,  flux, 'k', label = 'eBOSS Flux', drawstyle='steps-mid')
        p2.plot(wave, synflux, 'r', label = 'PCA fit', drawstyle='steps-mid')

        p2.set_ylim(np.min(flux[window]), np.max(flux[window])+0.5)
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
        p1.plot(wave,  flux[:], 'k', label = 'BOSS Flux', drawstyle='steps-mid')
        p1.plot(wave, synflux[:], 'r', label = 'PCA fit', drawstyle='steps-mid')
        #p1.plot(wave,ivar[:]-15, 'g', label = 'var$^{-1}-15$')

        p1.fill_between(wave,np.min(synflux[:])-10,np.max(synflux[:])+10,where=(ivar[:]==0),facecolor='k', alpha=0.2)
        p1.fill_between(wave,np.min(synflux[:])-10,np.max(synflux[:])+10,where=np.logical_and(5560<wave, wave<5600),facecolor='c', alpha=0.2)
        p1.fill_between(wave,np.min(synflux[:])-10,np.max(synflux[:])+10,where=np.logical_and(5880<wave, wave<5905),facecolor='c', alpha=0.2)
        p1.fill_between(wave,np.min(synflux[:])-10,np.max(synflux[:])+10,where=np.logical_and(6285<wave, wave<6315),facecolor='c', alpha=0.2)
        p1.fill_between(wave,np.min(synflux[:])-10,np.max(synflux[:])+10,where=np.logical_and(6348<wave,wave<6378),facecolor='c', alpha=0.2)

        p1.legend(loc='upper right', bbox_to_anchor = (1,1), ncol = 1, prop=fontP)
        box = p1.get_position()
        p1.set_position([box.x0,box.y0,box.width,box.height])

        p1.set_ylim(np.min(synflux[:])-3, np.max(synflux[:])+3)
        plt.ylabel('$f_{\lambda}\, [10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')
        if QSOlens == False:
            p1.set_xlim(3600,6000)


        window = np.linspace(wave2bin(x0,c0,c1,Nmax)-40,wave2bin(x0,c0,c1,Nmax)+40,81,dtype = np.int16)
        if QSOlens:
            p3 = plt.subplot(gs[1,:1])
            p3.plot(wave,  flux[:], 'k', label = 'BOSS Flux', drawstyle='steps-mid')
            p3.plot(wave, synflux[:], 'r', label = 'PCA fit', drawstyle='steps-mid')
            p3.set_xlim(np.min(wave[window]),np.max(wave[window]))
            p3.set_ylim(np.min(synflux[window])-1, np.max(flux[window])+1)
            box = p3.get_position()

            p3.set_position([box.x0,box.y0,box.width,box.height])
            plt.ylabel('$f_{\lambda}\, [10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')
            p3.legend(loc='upper right', bbox_to_anchor = (1,1), ncol = 1,prop=fontP)

            p3.locator_params(axis='x',nbins=6)

            median_local = np.median(reduced_flux[window])
            fit_QSO = np.poly1d(np.polyfit(x=wave[window],y=reduced_flux[window],deg=3,w=(np.abs(reduced_flux[window]-median_local)<5)*np.sqrt(ivar[window])) )
            p4 = plt.subplot(gs[1,1:2])
            p4.plot(wave[window], fit_QSO(wave[window]), '-m',label = 'Order 3 fit')
            box = p4.get_position()
            p4.set_xlim(np.min(wave[window]),np.max(wave[window]))

            p4.set_position([box.x0,box.y0,box.width,box.height])
            p4.plot(wave[window], reduced_flux[window],'k', label = 'Reduced flux', drawstyle='steps-mid')
            p4.legend(loc='upper right', bbox_to_anchor = (1,1), ncol = 1,prop=fontP)

            p4.locator_params(axis='x',nbins=6)
        else:
            p3 = plt.subplot(gs[1,:2])
            p3.plot(wave,  flux[:], 'k', label = 'BOSS Flux', drawstyle='steps-mid')
            p3.plot(wave, synflux[:], 'r', label = 'PCA fit', drawstyle='steps-mid')
            p3.legend(prop=fontP)
            p3.set_xlim(peak[0]-50,peak[0]+60)
            p3.set_ylim(np.min(synflux[bounds])-2, np.max(flux[bounds])+3)
            plt.ylabel('$f_{\lambda}\, [10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$', fontsize=18)



        ##### Old code to
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

    make_sure_path_exists(topdir + savedir +'/plots/')
    plt.savefig(topdir + savedir +'/plots/'+SDSSname(RA,DEC)+ '-' + str(plate) + '-' + str(mjd) + '-' + str(fiberid) + '-' + str(n_peak)+ '.eps', format = 'eps', dpi = 2000)

    plt.close()
    '''
