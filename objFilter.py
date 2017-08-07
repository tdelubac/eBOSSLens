def qsoFilter(obj, DR12Q, rchi2_max):
    '''
    objFilter.qsoFilter(obj)
    ========================
    Filter out sky or calibration fiber and non-QSO targets
    Filters QSOs with redshifts not matching DR12Q (Pâris et al. 2016)
    
    Parameters:
        obj: A SDSSObject instance to select
        DR12Q: The DR12Q catalog, already loaded
        rchi2_max: Maximum acceptable chi2 on redshift
    Returns:
        True for accept, False for reject
    '''

    # Take only QSOs
    if obj.obj_type[i].startswith('sky') or \
            obj.obj_type[i].startswith('SPECTROPHOTO_STD') or \
            (not obj.obj_class.startswith('QSO') ):
        raise Exception("Rejected: Not a proper QSO target")
    # Reject badly fitted QSOs in BOSS pipeline
    if obj.zwarning[i] != 0 or obj.rchi2[i] > rchi2_max or obj.z_err[i] < 0:
        raise Exception("Rejected: Bad redsfhit fitting")
    # Search for matching QSOs in DR12Q
    matching_QSOs = [x for x in DR12Q if int(x[0]) == int(self.plate) and
             int(x[1]) == int(self.mjd) and int(x[2]) == int(self.fid)]
    # Check that the redshifts match         
    if len(index) > 0:
        if n.abs(index[0][3] - obj.z[i]) < 0.005 and \
                n.abs(index[0][4] - obj.z[i]) > 0.1:
            continue
    else:
        raise Exception("Rejected: .fits redshift does not match DR12Q tabulated redshift")

    # Before masking, compute the FWHM of CIV or HBeta depending on redshift:
    #FWHM, l_times_luminosity, HB_wave, params_beta, line_coeff = \
    #    QSO_compute_FWHM(ivar = ivar[i,:],flux = flux[i,:], wave = wave,c0=c0,c1=c1,Nmax=Nmax,z =z[i],l_width = l_width)

    #M_BH = 10**(6.91 + n.log10(n.sqrt(5100*l_times_luminosity/1e44)*(FWHM/1000)**2)) # Masses solaires
    #sigma_host = 200*10**((n.log10(M_BH) - 7.92)/3.93) ### km s-1


    #ivar[i,:] = mask_QSO(ivar=ivar[i,:],z=z[i], l_width = l_width, c0 = c0 , c1=c1,Nmax=Nmax)

    return True


def genFilter(obj):
    '''
    objFilter.genFilter(obj)
    ========================
    Filter out sky, star, qso, or calibration fiber
    
    Parameters:
        obj: A SDSSObject instance to select
    Returns:
        True for accept, False for reject
    '''
    if obj.obj_class.startswith('STAR') or obj.obj_class.startswith('QSO') or \
                    obj.obj_type.startswith('SKY') or \
                    obj.obj_type.startswith('SPECTROPHOTO_STD'):
        raise Exception("Rejected: Not a proper galaxy target")
    else:
        return True
