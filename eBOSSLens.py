# Removes the fitted template from the spectra
# Fits the peaks with a Gaussian profile
# Determines OII doublet and lensed object redshift
# Estimates the Eintein radius of the lensing system

# Imports
import numpy as n
import pyfits as pf
import matplotlib as mpl
#mpl.use('Agg')
from matplotlib import pyplot as plt
import math
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.integrate import quad
from scipy import special as sp
from scipy import interpolate
import datetime
import copy
from matplotlib import rcParams
import os
import errno
import itertools as it
#-----------------------------------------------------------------------------------------------------
# Operation mode
searchLyA = True
QSOlens = False
testMode  = False

#-----------------------------------------------------------------------------------------------------
# Constants: flat Universe is assumed, OmegaLambda = 0.7, OmegaM = 0.3
H0 = 72e3 #m s-1 Mpc-1
c = 299792.458 #m s-1
OmegaM = 0.3
OmegaL = 0.7
l_LyA = 1215.668 #Angstroms
#-----------------------------------------------------------------------------------------------------
# Function definitions

# Generate a Gaussian around x_0 with amplitude A and variance var
def gauss(x,x_0,A,var):
	y = A*n.exp( (-(x-x_0)**2) / (2*var) )
	return y
	
def skew(x,A,w,a,eps):
	phi = 0.5*(1+sp.erf(a*(x-eps)/(w*n.sqrt(2))))
	return A*2*gauss(x,eps,1/n.sqrt(2*n.pi),w**2)*phi/w

#Reduced Chi square for one gaussian
def chi2g(params, xdata, ydata, ivar, x0):
	return n.sum(ivar*(ydata - gauss(x=xdata, x_0=params[0], A=params[1], var=params[2]))**2)/(len(xdata)-len(params)-1)
#Reduced Chi square for Doublet
def chi2D(params, xdata, ydata, ivar):
	return n.sum(ivar*(ydata - gauss(x=xdata, x_0=params[3], A=params[0], var=params[1])-gauss(x=xdata, x_0=params[4], A=params[2], var=params[1]))**2)/(len(xdata)-len(params) -1)
#Reduced Chi square for skew profile
def chi2skew(params, xdata, ydata, ivar):
	return n.sum(ivar*(ydata - skew(x=xdata,A = params[0], w=params[1], a=params[2], eps = params[3]))**2)/(len(xdata)-len(params)-1)
#Reduced Chi square for  double skew profile
def chi2skew2(params, xdata, ydata, ivar):
	return n.sum(ivar*(ydata - skew(x=xdata,A = params[0], w=params[1], a=params[2], eps = params[3]) - skew(x=xdata, A = params[4], w = params[5], a=params[6], eps=params[7]))**2)/(len(xdata)-len(params)-1)


# Check if x0 is near any emission line redshifted by z
def nearline(x0, zline, fiberid, z, mjd, plate):
	match1 = n.logical_and(abs(zline['linewave']*(1+z) -x0) < 10, zline['lineew']/zline['lineew_err'] > 6)
	match2 = n.logical_and(zline['fiberid']==fiberid,zline['mjd']==int(mjd))
	match3 = n.logical_and(zline['plate']==int(plate), zline['lineew_err']>0)
	match4 = n.logical_and(match1,n.logical_and(match2,match3))
	if (n.sum(match4)>0):
		return True
	else:
		return False

#Gaussian kernel used in first feature search (Bolton et al.,2004 method)	
def kernel(j,width,NormGauss,length):
	ker = n.zeros(length)
	ker[j-int(width*0.5):j+int(width*0.5)] = NormGauss
	return ker
	
# Estimated Einstein Radius from Single Isothermal Sphere (SIS) model
def radEinstein(z1,z2,vdisp):
	#compute ang. diam. distances
	Ds = ((c/H0)*quad(x12,0.0,z2)[0])/(1+z2)
	Dls = ((c/H0)*quad(x12,z1,z2)[0])/(1+z2)
	### return value in arcsec
	coeff = 3600*(180/n.pi)
	return coeff*4.0*n.pi*vdisp*vdisp*Dls/(c*c*Ds)

#Function needed for numerical computation of angular diameter distance
def x12(z):
	return 1.0/n.sqrt((1-OmegaM-OmegaL)*(1+z)*(1+z) + OmegaM*(1+z)**3 + OmegaL)
	
	
#Convert wavelength to bin number
def wave2bin(x,c0,c1,Nmax):
	b = int((n.log10(x)-c0)/c1)
	if b <= 0:
		return 0
	elif b >= Nmax:
		return Nmax
	else:
		return b
#Give BOSS approximated resolution as a function of wavelength
def resolution(x):
	if 4000<x<5800:
		a = (2000-1400)/(5800-4000)
		b = 1400-a*4000
		return a*x+b
	elif 5800<x<6200:
		a = (1900-2000)/(6200-5800)
		b = 2000-a*5800
		return a*x+b
	elif 6200<x<9400:
		a = (2600-1900)/(9400-6200)
		b = 2600-a*9400
		return a*x+b
	else:
		return 2500

#Prepare the flux in the BOSS bins starting from MC template/any datapoints array
def template_stretch(template_x, template_y,x0,A,B,eps):
	template_y = template_y*A
	m = n.mean(template_x) 
	template_x = B*(template_x -m) + m + eps
	
	sigma = x0/resolution(x0)
	gaussian_kernel = gauss(template_x,x_0=x0,A=1/n.sqrt(sigma*2*n.pi),var=sigma**2)
	template_y = n.convolve(template_y, gaussian_kernel,mode='same')
	return template_x, template_y
	
def chi2template(params,xdata,ydata, template_x, template_y, x0, ivar):
	template_x,template_y = template_stretch(template_x, template_y, x0, params[0],params[1],params[2])
	nearest_y = [template_y[(n.abs(template_x - w)).argmin()] for w in xdata]
	return n.sum(ivar*(ydata - nearest_y)**2)/(len(xdata)-len(params)-1)


# Check if a path exists, if not make it
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
#-----------------------------------------------------------------------------------------------------
#Set topdir,savedir:
topdir = '..'
savedir = '/LRG-LAE'
## import list of plates to analyze
plate_mjd = [line.strip().split() for line in open(topdir + '/fits_files/plates_Shu2016.txt')]

plate = 0
fiberid = [0]
if searchLyA == True and QSOlens == True:
	f = open(topdir + savedir + '/candidates_QSO_LyA.txt','a')
	f.close
elif searchLyA == True and QSOlens == False:
	f = open(topdir + savedir + '/candidates_LyA.txt','a')
	f.close
elif searchLyA == False:
	f = open(topdir + savedir + '/candidates_doublet.txt','a')
	f.close() 
	f = open(topdir + savedir + '/candidates_DM.txt','a')
	f.close()
	f = open(topdir + savedir + '/candidates_multi.txt','a')
	f.close()

#Set of emission lines used for lensed galaxy detection OII, Hb, OIII, OIII, Ha
em_lines = n.array([3726.5,4861.325,4958.911,5006.843,6562.801])
countDoublet = 0
countDM = 0
countMulti = 0

counter1 = 0
counter2 = 0
counter3 = 0
counter4 = 0

#Loop over plates
for j in n.arange(len(plate_mjd)):
	#print "initialization plate " , plate_mjd[j][0], "\n",
	
	# Initialization
	flux = 0
	plate = plate_mjd[j][0]
	mjd = plate_mjd[j][1]
	
	# Pick your plate/mjd and read the data:
	plate = plate_mjd[j][0]
	spfile = topdir + '/fits_files/spPlate-' + str(plate) + '-' + str(mjd) + '.fits'
	zbfile = topdir + '/fits_files/spZbest-' + str(plate) + '-' + str(mjd) + '.fits'
	zlfile = topdir + '/fits_files/spZline-' + str(plate) + '-' + str(mjd) + '.fits'
	hdulist = pf.open(spfile)
	c0 = hdulist[0].header['coeff0']
	c1 = hdulist[0].header['coeff1']
	npix = hdulist[0].header['naxis1']
	wave = 10.**(c0 + c1 * n.arange(npix))
	# Following commented-out bit was needed for some of the early redux:
	#bzero = hdulist[0].header['bzero']
	bunit = hdulist[0].header['bunit']
	flux = hdulist[0].data	
	ivar = hdulist[1].data
	ivar_copy = copy.deepcopy(ivar)
	hdulist.close()
	hdulist = 0
	hdulist = pf.open(zbfile)
	vdisp = hdulist[1].data.field('VDISP')
	vdisp_err = hdulist[1].data.field('VDISP_ERR')
	synflux = hdulist[2].data
	fiberid = hdulist[1].data.field('FIBERID')
	RA = hdulist[1].data.field('PLUG_RA')
	DEC = hdulist[1].data.field('PLUG_DEC')
	obj_class = hdulist[1].data.field('CLASS')
	obj_type = hdulist[1].data.field('OBJTYPE')
	z = hdulist[1].data.field('Z')
	zwarning = hdulist[1].data.field('ZWARNING')
	z_err = hdulist[1].data.field('Z_ERR')
	rchi2 = hdulist[1].data.field('RCHI2')
	hdulist.close()
	hdulist = 0
	hdulist = pf.open(zlfile)
	zline = hdulist[1].data
	hdulist.close()
	hdulist = 0
	##### PlotSpec ends here 
	#-----------------------------------------------------------------------------------------------------
	
	reduced_flux = n.array(flux - synflux)
	sqrtivar=copy.deepcopy(ivar)
	
	Nmax = len(flux[0,:])
	
	
	# Masks atmosphere
	#ivar[:,542:551] = 0 # Hg line
	#ivar[:,868:873] = 0 # Hg line
	#ivar[:,1847:1852] = 0 # Hg line
	#ivar[:,1938:1944] = 0 # OI line
	#ivar[:,2022:2031] = 0 # lines
	#ivar[:,2175:2185] = 0 # lines
	#ivar[:,423:428] = 0 # Ca absorption
	#ivar[:,460:467] = 0 # Ca absorption
	
	#Mask 5580 A galaxy spectra
	ivar[:,wave2bin(5570,c0,c1,Nmax): wave2bin(5590,c0,c1,Nmax)] = 0
	ivar[:,wave2bin(5880,c0,c1,Nmax): wave2bin(5905,c0,c1,Nmax)] = 0
	

	
	startTime = datetime.datetime.now()

	# Loop over objects
	for i in  n.arange(len(flux[:,0])):
	#for i in n.array([788,820,702,721,148,346,638,651,754,340,377,609,696,884,377,614,940,59,253,540]):
		if QSOlens:
			# Take only QSOs
			if (not(obj_class[i] == 'QSO   ') or obj_type[i] =='SKY             ' or 'SPECTROPHOTO_STD'==obj_type[i]): 
				continue
			# Reject badly fitted QSOs
			if (zwarning[i] != 0 or rchi2[i]>3 or z_err[i] < 0):
				continue
			# Masking Lya, NV, SiIV, CIV, etc...
			l_NV = 1240
			l_SiIV= 1400.0
			l_CIV = 1549.0
			l_HeII = 1640.0
			l_CIII = 1909.0
			l_CII = 2326.0
			l_FeII_a = 2382.765
			l_FeII_b = 2600.173
			l_MgII = 2798.0
			l_NeV = 3426.0
			l_OII = 3727
			l_NeIII = 3869
			l_Hd = 4101
			l_Hg = 4340
			l_Hb = 4861
			l_OIII_a = 4959
			l_OIII_b = 5007
			l_width = 20
			#Mask the above emission lines of QSO's
			ivar[i,wave2bin((1+z[i])*(l_LyA -2.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_LyA +2.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_NV -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_NV +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_SiIV -1.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_SiIV +1.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_CIV -2*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_CIV +2*l_width),c0,c1,Nmax)] = 0 
			ivar[i,wave2bin((1+z[i])*(l_HeII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_HeII +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_CIII -l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_CIII +l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_CII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_CII +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_FeII_a -l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_FeII_a +l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_FeII_b -l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_FeII_b +l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_MgII -1.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_MgII +1.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_NeV -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_NeV +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_OII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_OII +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_NeIII -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_NeIII +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_Hd -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_Hd +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1.2+z[i])*(l_Hg - l_width),c0,c1,Nmax):wave2bin((1.2+z[i])*(l_Hg + l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_Hb - l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_Hb + l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_OIII_a -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_OIII_a +0.5*l_width),c0,c1,Nmax)] = 0
			ivar[i,wave2bin((1+z[i])*(l_OIII_b -0.5*l_width),c0,c1,Nmax):wave2bin((1+z[i])*(l_OIII_b +0.5*l_width),c0,c1,Nmax)] = 0
		else:
			if (obj_class[i] == 'STAR  ' or obj_class[i] == 'QSO   ' or obj_type[i] =='SKY             ' or 'SPECTROPHOTO_STD'==obj_type[i]): 
				continue
				
		
		peaks = []
		peak_number = len(peaks)
		doublet = None
		
		counter1 = counter1 +1;
		
		sqrtivar[i,:] = n.sqrt(ivar[i,:])
		
		### Bolton 2004: S/N of maximum likelihood estimator of gaussian peaks
		if searchLyA == True:
			width = 30.0
			sig = 2
		elif searchLyA == False:
			width = 30.0
			sig = 2.4
		## Prepare normalized gaussian
		NormGauss = gauss(n.linspace(-width*0.5,width*0.5,width),0.0,1.0,sig)
		NormGauss = NormGauss/n.sum(NormGauss)
		Cj1 = n.array([n.sum(reduced_flux[i,:]*kernel(j+0.5*width,width,NormGauss,len(wave))*ivar[i,:]) for j in range(int(len(wave)-width))])
		Cj2 = n.array([n.sum(ivar[i,:]*kernel(j+0.5*width,width,NormGauss,len(wave))**2) for j in range(int(len(wave)-width))])
		SN = n.zeros(len(wave))
		SN[width*0.5:len(wave)-width*0.5] = Cj1/n.sqrt(Cj2)
		
		if searchLyA == True and QSOlens:
			peak_candidates = n.array([(x0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(wave,SN) if (test>12.0 and  5600<x0)])
		elif searchLyA == True and QSOlens == False:
			peak_candidates = n.array([(x0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(wave,SN) if (test>8.0 and  3600<x0<4800)])
		elif searchLyA == False and QSOlens == False:
			peak_candidates = n.array([(x0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(wave,SN) if test>6.0])
		else:
			continue
		# (Legend key:) wavelength(x0) chisq_doublet amp_gauss var_gauss S/N chisq_doublet amp1_doublet amp2_doublet var_doublet x1 x2
		
		#Keep only center of candidate peaks
		k = 0
		while (k < (len(peak_candidates)-1)):
			if (abs(peak_candidates[k][0] - peak_candidates[k+1][0]) < 10):
				if peak_candidates[k][4] < peak_candidates[k+1][4]:
					peak_candidates = n.delete(peak_candidates,k, axis=0)
					k = k-1
				else:
					peak_candidates = n.delete(peak_candidates,k+1,axis = 0)
					k = k-1					
			k = k+1
		chisq = 10000
		chisq2 = 10000
		chisq_skew = 10000
		#Search for suitable peak candidates
		for peak in peak_candidates:
			x0 = peak[0]
			if nearline(x0, zline, fiberid[i], z[i], int(mjd), int(plate)):
				continue
 			
 			bounds = n.linspace(wave2bin(x0,c0,c1,Nmax)-20,wave2bin(x0,c0,c1,Nmax)+20,41,dtype = n.int16)

			#Single Line: Gaussian fit around x_0
			if searchLyA == True:
				init = [x0,4,6]
				res =  minimize(chi2g,init,args=(wave[bounds], reduced_flux[i,bounds],ivar[i,bounds],x0), method='SLSQP', bounds = [(x0-2,x0+2),(1,100),(1,15)])
			elif searchLyA == False:
				init = [x0,1,2]
				res =  minimize(chi2g,init,args=(wave[bounds], reduced_flux[i,bounds],ivar[i,bounds],x0), method='SLSQP', bounds = [(x0-2,x0+2)(0.1,5),(1,8)])
			params = res.x
			chisq = res.fun
			#residue_squared=((reduced_flux[i,:] - gauss(x=wave, x_0=params[0], A=params[1], var=params[2]))**2)*ivar[i,:]
			#chisq =  n.sum(residue_squared[bounds])/(len(bounds)-len(params) -1)
			chigauss1 = chisq
			#Check for not too high chi square and save
			if (not(chisq > 2.0) or searchLyA):
				peak[1] = chisq
				peak[2] = params[1]
				peak[3] = params[2]
				
			#Doublet OII: Gaussian fit around x_0
			if (x0 > 3727.0*(1+z[i]) or searchLyA): 

				res2 = minimize(chi2D,[1.0,5,1.0,x0-1.5,x0+1.5],args=(wave, reduced_flux[i,:],ivar[i,:]), method='SLSQP', bounds = [(0.1,5),(1,8),(0.1,5),(x0-7,x0),(x0,x0+7)])
				params2 = res2.x
				residue_squared=((reduced_flux[i,:] - gauss(x=wave, x_0=params2[3], A=params2[0], var=params2[1]) - gauss(x=wave, x_0=params2[4], A=params2[2], var=params2[1]))**2)*ivar[i,:]
				chisq2 = n.sum(residue_squared[bounds])/(len(bounds)-len(params2) -1)
				
				if  0.5*x0/3726.5>abs(params2[3]-params2[4])>2.1*x0/3726.5 and not(chisq2 > 2):	
					peak[5] = chisq2
					peak[6] = params2[0] #amp1
					peak[7] = params2[2] #amp2
					peak[8] = params2[1] #var
					peak[9] = params2[3] #x1
					peak[10] = params2[4] #x2
					print 'Valid doublet'
					
					
					
				# Delta OII restframe: 	1.3 A  (3725.94 3727.24)
					
				
			# If looking at LAE, test a skew-normal profile as well 
			if searchLyA:
				# Sharp blue, red tail
				init_skew = [params[1],0.5,2,x0]
				res_skew = minimize(chi2skew,init_skew,args=(wave, reduced_flux[i,:],ivar[i,:]), method='SLSQP', bounds = [(2,100),(0.0,10),(1,10),(x0-4,x0+4)])
				params_skew_a = res_skew.x
				residue_squared=((reduced_flux[i,:] - skew(x=wave,A = params_skew_a[0], w=params_skew_a[1], a=params_skew_a[2], eps=params_skew_a[3]) )**2)*ivar[i,:]

				chisq_skew_a = n.sum(residue_squared[bounds])/(len(bounds)-len(params_skew_a) -1)
				# Double Sharp blue, red tail
				#init_skew = [params[1],0.5,2,x0, params[1]/2,0.5,2,x0+8]
				#res_skew = minimize(chi2skew2,init_skew,args=(wave, reduced_flux[i,:],ivar[i,:]), method='SLSQP', bounds = [(2,100),(0.0,10),(1,10),(x0-4,x0+4), (1,10),(0.0,10),(1,10),(x0+4,x0+10)])
				#params_skew_c = res_skew.x
				#residue_squared=((reduced_flux[i,:] - skew(x=wave,A = params_skew_c[0], w=params_skew_c[1], a=params_skew_c[2], eps=params_skew_c[3]) -  skew(x=wave,A = params_skew_c[4], w=params_skew_c[5], a=params_skew_c[6], eps=params_skew_c[7]))**2)*ivar[i,:]
				#bounds = n.linspace(wave2bin(x0,c0,c1,Nmax)-10,wave2bin(x0,c0,c1,Nmax)+10,21,dtype = n.int16)
				#chisq_skew_c = n.sum(residue_squared[bounds])/(len(bounds)-len(params_skew_c) -1)
				# Sharp red, blue tail
				init_skew = [params[1],0.5,-2,x0]
				res_skew = minimize(chi2skew,init_skew,args=(wave, reduced_flux[i,:],ivar[i,:]), method='SLSQP', bounds = [(2,100),(0.0,10),(1,10),(x0-4,x0+4)])
				params_skew_b = res_skew.x
				residue_squared =((reduced_flux[i,:] - skew(x=wave,A = params_skew_b[0], w=params_skew_b[1], a=params_skew_b[2], eps=params_skew_b[3]) )**2)*ivar[i,:]
				chisq_skew_b = n.sum(residue_squared[bounds])/(len(bounds)-len(params_skew_b) -1)
				
				LyAdoublet = False
				if chisq_skew_b < chisq_skew_a: #and chisq_skew_b < chisq_skew_c:
					params_skew = params_skew_b
					chisq_skew = chisq_skew_b
				elif chisq_skew_a < chisq_skew_b: #and chisq_skew_a < chisq_skew_c:
					params_skew = params_skew_a
					chisq_skew = chisq_skew_a
				#elif chisq_skew_c < chisq_skew_a and chisq_skew_c < chisq_skew_b:
				#	params_skew = params_skew_c
				#	chisq_skew = chisq_skew_c
				#	LyAdoublet = True
				
				# Testing A.Verhamme MCMC LAE profiles
				
				template_keys = n.array(['./LAE_TEMPLATE/spec_V0_2N17_B40_D0_E150_F500.dat','./LAE_TEMPLATE/spec_V0_2N18_B40_D0_E150_F500.dat'])
										#'./LAE_TEMPLATE/spec_V0_2N19_B40_D0_E150_F500.dat','./LAE_TEMPLATE/spec_V0_2N17_B40_D0_E100_F150.dat', 
										#'./LAE_TEMPLATE/spec_V0_2N18_B40_D0_E100_F150.dat','./LAE_TEMPLATE/spec_V0_2N19_B40_D0_E100_F150.dat', 
										#'./LAE_TEMPLATE/spec_V0_2N20_B40_D0_E100_F150.dat','./LAE_TEMPLATE/spec_V0_2N21_B40_D0_E100_F150.dat', 
										#'./LAE_TEMPLATE/spec_V150_2N17_B40_D0_E100_F150.dat','./LAE_TEMPLATE/spec_V150_2N18_B40_D0_E100_F150.dat',
										#'./LAE_TEMPLATE/spec_V150_2N19_B40_D0_E100_F150.dat','./LAE_TEMPLATE/spec_V150_2N20_B40_D0_E100_F150.dat', 
										#'./LAE_TEMPLATE/spec_V150_2N21_B40_D0_E100_F150.dat','./LAE_TEMPLATE/spec_V300_2N19_B40_D0_E100_F150.dat',
										#'./LAE_TEMPLATE/spec_V300_2N20_B40_D0_E100_F150.dat','./LAE_TEMPLATE/spec_V300_7N20_B40_D0_E100_F150.dat'])
				best_template_key = ''
				best_template_chisq = 100000000
				best_template_params = [1,1,0]
				init = [1,1,0]
				for key in template_keys:
					data = n.loadtxt('.'+ key)
					template_x = data[:,0]
					template_y = data[:,1]
					
					z_LAE = x0/l_LyA-1
					template_x = (1+z_LAE)*l_LyA*(template_x/c +1 )
					template_y = 5*template_y/n.max(template_y)
					
					res =  minimize(chi2template,init,args=(wave[bounds], reduced_flux[i,bounds], template_x, template_y, x0, ivar[i,bounds]), method='SLSQP', bounds = [(0.,10000),(0.,1000),(-5,+5)])
					params = res.x
					chisq_template = res.fun
					print res.x
					if chisq_template < best_template_chisq:
						best_template_chisq = chisq_template
						best_template_key = key
						best_template_params = params
						
				print chisq, chisq2, chisq_skew, chisq_template		
				
				ax = plt.subplot(1,1,1)
				plt.title('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))
				ax.plot(wave, reduced_flux[i,:],'k')
				plt.xlabel('$Wavelength\, (Angstroms)$')
				plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
				ax.plot(wave,gauss(x=wave, x_0=params2[3], A=params2[0], var=params2[1]) + gauss(x=wave, x_0=params2[4], A=params2[2], var=params2[1]),'g', label = r'$\chi_D^2 = $' + '{:.4}'.format(chisq2))
				#ax.plot(wave,gauss(x=wave, x_0=params[0], A=params[1], var=params[2]),'r', label = r'$\chi_G^2 = $' + '{:.4}'.format(chisq) )
				ax.plot(wave,skew(x=wave,A = params_skew[0], w=params_skew[1], a=params_skew[2], eps=params_skew[3]), 'b', label = r'$\chi_S^2 = $' + '{:.4}'.format(chisq_skew))
				best_template = n.loadtxt('.'+ best_template_key)
				data_x = best_template[:,0]
				data_y = best_template[:,1]
				plot_x, plot_y = template_stretch(data_x,data_y,x0, best_template_params[0],best_template_params[1],best_template_params[2])
				print plot_x
				ax.plot(plot_x,plot_y,'r')
				ax.legend()			
				ax.set_xlim(3600,6000)	
				plt.show()
		
		counter2 = counter2 + 1;					
		#Finding peak with lowest chi square for doublet and see if it is better fitted by single line or not
		doublet_index = 0
		chi2saved = 1000.0
		# Find the doublet index
		if searchLyA == False:
			for k in range(len(peak_candidates)):
				peak = peak_candidates[k]
				if (peak[1]>peak[5]>0 and peak[5]< chi2saved and peak[0]<9200) :
					peak[1] = peak[5]
					chi2saved = peak[5]
					doublet = True
					doublet_index = k
		#elif searchLyA:
		#	check chi^2 gauss, doublet, LAE with delta preset
		# 	condition on separation of doublet
		#	if abs(peak[] - peak[]) > 
		#		max(gauss,lyA)
		#	else 
		#		max(gauss,doublet,lyA)
		
		
		#Removing candidates that were not fitted : params still 0
		peak_candidates = n.array([peak for peak in peak_candidates if (peak[2]!=0 or peak[5]==peak[1]!=0)])
		if len(peak_candidates) == 0:
			continue
		
		counter3 = counter3+1;
		#Sorting candidates by chi square
		peak_candidates = sorted(peak_candidates, key=lambda peak: peak[1])
		
		# Keeping only 5 most likely candidates
		if len(peak_candidates) > 5:
			peak_candidates = peak_candidates[0:5]
		#find again doublet index
		found = False
		for k in range(len(peak_candidates)):
			if (peak_candidates[k][5] == chi2saved):
				doublet_index = k
				found = True
		if found==False:
			doublet = False
			
		# Check that at least 1 candidate is below 9200 Angstrom cut, if not, go to next fiber
		below_9200 = False
		for peak in peak_candidates:
			if peak[0] < 9200:
				below_9200 = True
		if below_9200 ==False:
			continue
		
		counter4 = counter4+1;
		#Try to infer background redshift
		detection = False
		score = 0.0
		
		#print counter1, counter2, counter3, counter4 
		if (doublet == True and searchLyA == False):
			fileD = open(topdir + savedir +  '/candidates_doublet.txt','a')
			z_s = peak_candidates[doublet_index][0]/3727.24 - 1.0
			if (z_s > z[i]+0.05):
				detection = True
				score += peak_candidates[doublet_index][5]
				if testMode == False:
					fileD.write('\n' +str([radEinstein(z[i],z_s,vdisp[i]*1000), score,z_s, RA[i], DEC[i], int(plate), int(mjd), fiberid[i],peak_candidates[doublet_index][9]]))
				fileD.close()
			if len(peak_candidates):	
				fileDM = open(topdir + savedir + '/candidates_DM.txt','a')
				confirmed_lines = []
				#Generating all combinations of lines from above list to compare with candidates
				temp = [peak for peak in peak_candidates if peak[1]!=peak[5]]
				compare = em_lines[1:5]
				if z_s > z[i]+0.05 :
					for peak in temp:
						for line in compare:
							if ( abs(peak[0]/line -1 - z_s) < 0.01):
								detection = True
								confirmed_lines.append(line)
								score+= peak[5]
					if (confirmed_lines != [] and testMode == False):
						fileDM.write('\n'+str([radEinstein(z[i],z_s,vdisp[i]*1000), score,z_s, RA[i], DEC[i], int(plate), int(mjd), fiberid[i],confirmed_lines]))
					fileDM.close()
		elif (doublet != True and len(peak_candidates) > 1 and searchLyA == False ):
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
			if (confirmed_lines != [] and  testMode == False):
				fileM.write('\n'+str([radEinstein(z[i],z_s,vdisp[i]*1000), score, z_s, RA[i], DEC[i], int(plate), int(mjd), fiberid[i],confirmed_lines]))
			fileM.close()
		elif searchLyA:
			if QSOlens:
				fileLyA = open(topdir + savedir +  '/candidates_QSO_LyA.txt','a')
			else:
				fileLyA = open(topdir + savedir +  '/candidates_LyA.txt','a')
			for peak in peak_candidates:
				z_O2 = peak[0]/3727.24 - 1.0
				# compute SN of Ha, Hb, OIII to conclude if LyA candidate is or not OII at low-redshift
				SNlines = 0
				em_lines = n.array([4861.325,5006.843,6562.801])
				for l in em_lines:
					center_bin = wave2bin(l*(1+z_O2),c0,c1,Nmax)
					SNlines += max(SN[center_bin-2:center_bin+2])**2
				SNlines = n.sqrt(SNlines)
				
				if SNlines < 5:
					fileLyA.write('\n' + str([peak[0],peak[4],z[i], SNlines ,RA[i], DEC[i], int(plate), int(mjd), fiberid[i], chisq, chisq2, chisq_skew]))
					detection = True
			fileLyA.close()

		# Save surviving candidates
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
		#print "plate ", plate, " fiber ", fiberid[i], " peak ", peak_number,  "\n",		
	
		######################################
		#TEST MODE
		#if testMode: 
			#print 'Doublet: ', doublet,'9200: ', below_9200, 'Detection: ', detection
			#detection = True
			#doublet = True
			#below_9200 = True
			#print peak_candidates
		######################################
		
		
		#Graphs
		if ((peak_number>1 or doublet==True) and below_9200 and detection and searchLyA==False):
			#Computing total fit of all peaks
			fit=0
			for k in n.arange(len(peaks)):
				fit = fit + gauss(wave, x_0 = peaks[k][0] , A=peaks[k][1], var=peaks[k][2])

			if doublet != True: 
				ax = plt.subplot(1,1,1)
				plt.title('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))
				ax.plot(wave, reduced_flux[i,:],'k')
				plt.xlabel('$Wavelength\, (Angstroms)$')
				plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
				ax.plot(wave,fit,'r')
				make_sure_path_exists(topdir + savedir +'/plots/')
				if testMode:
					plt.show()
				else:
					plt.savefig(topdir + savedir +'/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid[i]) + '.png')
				plt.close()
	
			# If doublet, plot in two different windows
			if doublet==True:
				x_doublet = 0.5*(peak_candidates[doublet_index][9]+ peak_candidates[doublet_index][10])
				bounds = n.linspace(wave2bin(x_doublet,c0,c1,Nmax)-10,wave2bin(x_doublet,c0,c1,Nmax)+10,21,dtype = n.int16)
				f = open(topdir + savedir +'/doublet_ML.txt','a')
				f.write('\n' +  str(plate) + ' ' + str(mjd) + ' ' + str(fiberid[i]) + ' ' + str(reduced_flux[i,bounds]))
				f.close()
				
				plt.figure(figsize=(14,6))
				ax1 = plt.subplot2grid((1,3), (0,0), colspan=2)
				plt.suptitle('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))
				ax2 = plt.subplot2grid((1,3), (0,2))
				ax1.plot(wave[10:-10], reduced_flux[i,10:-10],'k')
				ax1.plot(wave,fit,'r')
				ax1.set_xlabel('$\lambda \, [\AA]$ ')
				ax1.set_ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
				ax2.set_xlabel('$\lambda \, [\AA]$ ')
				ax2.locator_params(tight=True)
				
				ax2.set_xlim([peak_candidates[doublet_index][0]-30,  peak_candidates[doublet_index][0]+30])
				ax2.plot(wave, reduced_flux[i,:],'k')
				ax2.plot(wave,fit,'r')
				ax2.set_ylim([-5,10])
				
				ax2.vlines(x = zline['linewave']*(1+z[i]),ymin= -10,ymax= 10,colors= 'g',linestyles='dashed')
				#ax1.vlines(x = zline['linewave']*(1+z[i]),ymin= -10,ymax= 10,colors= 'g',linestyles='dashed')
				ax1.set_xlim([n.min(wave),n.max(wave)])
				
				make_sure_path_exists(topdir + savedir +'/plots/')
				if testMode:
					plt.show()
				else:
					plt.savefig(topdir + savedir +'/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid[i]) + '.png')
				plt.close()
		elif searchLyA == True and QSOlens == False and detection:
			ax = plt.subplot(1,1,1)
			plt.title('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))
			ax.plot(wave, reduced_flux[i,:],'k')
			plt.xlabel('$Wavelength\, (Angstroms)$')
			plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
			ax.set_xlim(3600,4800)
			maxAmp = 0
			for peak in peaks:
				if peak[1] > maxAmp:
					maxAmp = peak[1]
				ax.vlines(x = peak[0],ymin= -50,ymax= 50,colors= 'r',linestyles='dashed')
			ax.set_ylim(-5, int(maxAmp+10))
 			make_sure_path_exists(topdir + savedir +'/plots/')
			if testMode:
				plt.show()
			else:
				plt.savefig(topdir + savedir +'/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid[i]) + '.png')
			plt.close()
		elif searchLyA == True and QSOlens == True and detection:
			p1 = plt.subplot2grid((2,1), (0,0))
			p1.plot(wave,  flux[i,:], 'k', label = 'BOSS flux')
			p1.plot(wave, synflux[i,:], 'r', label = 'Fitted PCA template')
			plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2}  \AA^{-1}$)')
			p1.legend()
			p1.set_xlim(3600,6000)
			
			plt.suptitle('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))

			p2 = plt.subplot2grid((2,1), (1,0))
			p2.plot(wave, reduced_flux[i,:],'k')
			p2.plot(wave, ivar[i,:]>0,'g')
			plt.xlabel('$Wavelength\, (Angstroms)$')
			plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
			p2.set_xlim(3600,6000)
			peaks_list = n.array(peaks)
			p2.vlines(x = peaks_list[:,0],ymin= -50,ymax= 50,colors= 'r',linestyles='dashed')
			p2.set_ylim(-5, 2+max(reduced_flux[i,wave2bin(3600,c0,c1,Nmax):wave2bin(6000,c0,c1,Nmax)])) 
 			make_sure_path_exists(topdir + savedir +'/plots/')
			if testMode:
				plt.show()
			else:
				plt.savefig(topdir + savedir +'/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid[i]) + '.png')
			plt.close()
		

			
#print 'Time taken ', (datetime.datetime.now() - startTime)

