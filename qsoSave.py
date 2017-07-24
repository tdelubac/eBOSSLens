def qsoSave():
    # TODO: complete the function/parameters/returns
    '''
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
                reduced_flux = reduced_flux[i,:], c0=c0,c1=c1,Nmax=Nmax,show = plot_show, savedir=savedir, HB_wave = HB_wave , params_beta=params_beta, line_coeff =line_coeff)
            # plot with +-1 AA for paper
            plot_QSOGal(k=k+len(peak_candidates),RA = RA[i],DEC= DEC[i],plate = int(plate), mjd = int(mjd), fiberid = fiberid[i],z=z[i], z_backgal= z_backgal+0.0005,flux=flux[i,:],wave=wave,synflux=synflux[i,:],ivar= ivar[i,:], \
                reduced_flux = reduced_flux[i,:], c0=c0,c1=c1,Nmax=Nmax,show = plot_show, savedir=savedir, HB_wave = HB_wave , params_beta=params_beta, line_coeff =line_coeff)
            plot_QSOGal(k=k+len(peak_candidates)+1,RA = RA[i],DEC= DEC[i],plate = int(plate), mjd = int(mjd), fiberid = fiberid[i],z=z[i], z_backgal= z_backgal-0.0005,flux=flux[i,:],wave=wave,synflux=synflux[i,:],ivar= ivar[i,:], \
                reduced_flux = reduced_flux[i,:], c0=c0,c1=c1,Nmax=Nmax,show = plot_show, savedir=savedir, HB_wave = HB_wave , params_beta=params_beta, line_coeff =line_coeff)
    '''


def plotQSOGal(obj, peak, savedir, hbwave, paramsBeta, lineCoeff):
    # TODO: complete the function/parameters/returns
    '''
    em_lines = np.array([3726.5,4861.325,4958.911,5006.843,6562.801])

    if show ==False:
        mpl.use('Agg')
    fontP = FontProperties()
    fontP.set_size('medium')

    plt.suptitle(SDSSname(RA,DEC)+'\n'+'RA='+str(RA)+', Dec='+str(DEC) +', $z_{QSO}='+'{:03.3}'.format(z)+ '$')


    gs = gridspec.GridSpec(2,4)
    p1 = plt.subplot(gs[0,:4])


    smoothed_flux = np.array([np.mean(flux[ii-2:ii+3]) for ii in range(len(flux[:])) if (ii>4 and ii<len(flux[:])-4)])


    p1.plot(wave[5:-4], smoothed_flux, 'k', label = 'BOSS Flux', drawstyle='steps-mid')
    #p1.plot(wave,  flux[:], 'k', label = 'BOSS Flux')
    p1.plot(wave, synflux[:], 'r', label = 'PCA fit')

    if z<1 and show == True:

        p1.plot(HB_wave, lorentz(HB_wave, params_beta[0],params_beta[1],params_beta[2]) + HB_wave*line_coeff[0] + line_coeff[1], '--g')
    box = p1.get_position()

    p1.set_position([box.x0,box.y0+0.02,box.width*0.9,box.height])

    p1.set_ylim(np.min(synflux[:])-3, np.max(synflux[:])+3)

    p1.vlines(x = em_lines*(1+z_backgal),ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
    p1.legend(loc='upper right', bbox_to_anchor = (1.2,1), ncol = 1, prop=fontP)
    p1.set_xlim(3500,10500)
    plt.ylabel('Flux [$10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')

    p2 = plt.subplot(gs[1,:1])
    p2.vlines(x = em_lines*(1+z_backgal),ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
    loc_flux = flux[wave2bin((1+z_backgal)*(3727-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(3727+10),c0,c1,Nmax)]
    p2.plot(wave[wave2bin((1+z_backgal)*(3727-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(3727+10),c0,c1,Nmax)],loc_flux,'k', label = 'OII', drawstyle='steps-mid')
    p2.plot(wave[wave2bin((1+z_backgal)*(3727-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(3727+10),c0,c1,Nmax)],synflux[wave2bin((1+z_backgal)*(3727-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(3727+10),c0,c1,Nmax)],'r', label = 'OII', drawstyle='steps-mid')
    if loc_flux != []:

        p2.set_ylim(np.min(loc_flux)-1,np.max(loc_flux)+1)

    plt.title('[OII] 3727')
    p2.set_xlim((1+z_backgal)*(3727-10),(1+z_backgal)*(3727+10))
    x1 = int((1+z_backgal)*3727)
    plt.xticks([x1-15,x1,x1+15])
    plt.ylabel('Flux [$10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}]$')

    #If Ha is below 9500 A, show it
    if z>0.44:
        p3 = plt.subplot(gs[1,1:4])
    else:
        p3 = plt.subplot(gs[1,1:3])
    p3.vlines(x = em_lines*(1+z_backgal),ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
    loc_flux = flux[wave2bin((1+z_backgal)*(4861-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(5007+10),c0,c1,Nmax)]
    p3.plot(wave[wave2bin((1+z_backgal)*(4861-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(5007+10),c0,c1,Nmax)],loc_flux,'k', label = 'OIII, Hb', drawstyle='steps-mid')
    p3.plot(wave[wave2bin((1+z_backgal)*(4861-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(5007+10),c0,c1,Nmax)],synflux[wave2bin((1+z_backgal)*(4861-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(5007+10),c0,c1,Nmax)],'r', label = 'OIII, Hb', drawstyle='steps-mid')
    if loc_flux != []:

        p3.set_ylim(np.min(loc_flux)-1,np.max(loc_flux)+1)

    plt.title(r'H$\beta$,[OIII] 4959, [OIII] 5007')
    plt.xlabel(r'Observed wavelength [$\AA$]')
    p3.set_xlim((1+z_backgal)*(4861-10),(1+z_backgal)*(5007+10))
    x1 = int((1+z_backgal)*4862/10.)*10

    if x1<7600:
        plt.xticks([x1,x1+50 , x1+100 ,  x1 +150 ,x1+200])
    else:
        plt.xticks([x1,x1+50 , x1+100 ,  x1 +150 ,x1+200, x1+ 250])


    box = p3.get_position()
    p3.set_position([box.x0+0.02,box.y0,box.width*0.9,box.height])

    if z<0.44:
        p4 = plt.subplot(gs[1,3:4])
        p4.vlines(x = em_lines*(1+z_backgal),ymin= -100,ymax= 100,colors= 'g',linestyles='dashed')
        loc_flux = flux[wave2bin((1+z_backgal)*(6562-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(6562+10),c0,c1,Nmax)]
        p4.plot(wave[wave2bin((1+z_backgal)*(6562-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(6562+10),c0,c1,Nmax)],loc_flux,'k', label = 'Ha', drawstyle='steps-mid')
        p4.plot(wave[wave2bin((1+z_backgal)*(6562-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(6562+10),c0,c1,Nmax)],synflux[wave2bin((1+z_backgal)*(6562-10),c0,c1,Nmax) :wave2bin((1+z_backgal)*(6562+10),c0,c1,Nmax)],'r', label = 'Ha', drawstyle='steps-mid')
        if loc_flux != []:

            p4.set_ylim(np.min(loc_flux)-1,np.max(loc_flux)+1)
        plt.title(r'H$\alpha$')
        p4.set_xlim((1+z_backgal)*(6562-10),(1+z_backgal)*(6562+10))
        x1 = int((1+z_backgal)*6562)
        if x1 < 9900:
            plt.xticks([x1-10,x1,x1+10], [str(x1-10),str(x1),str(x1+10)])
        else:
            plt.xticks([x1-10,x1,x1+10], [str(x1-10),'',str(x1+10)])
    make_sure_path_exists(savedir +'/plots/')
    plt.savefig(savedir +'/plots/'+SDSSname(RA,DEC)+ '-' + str(plate) + '-' + str(mjd) + '-' + str(fiberid) + '-'+str(k+1) +'.png')

    if show:
        plt.show()
    else:
        plt.close()
    '''



