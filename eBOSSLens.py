# Removes the fitted template from the spectra
# Fits the peaks with a Gaussian profile
# Based on PlotSpec.py

# TO DO
# Sort objects by # significant peaks


# Imports
import numpy as n
import pyfits as pf
import matplotlib as mpl
from matplotlib import pyplot as p
import math
from scipy.optimize import leastsq
from scipy.optimize import minimize
import datetime
import copy
from matplotlib import rcParams
#rcParams['text.usetex'] = True
import os
import errno
import itertools as it

#-----------------------------------------------------------------------------------------------------
# Function definitions

# Generate a Gaussian around x_0 with amplitude A and variance var
def gauss(x,x_0,A,var):
	y = A*n.exp( (-(x-x_0)**2) / (2*var) )
	return y

# The function whose square is to be minimised.
# params = list of parameters tuned to minimise function.
# Further arguments:
# xdata = xaxis
# ydata = observed data
# yerr = variance of y
def func(params, xdata, ydata, sqrtivar, x0, llimit, ulimit):
	for i in n.arange(len(llimit)):
		if ((params[i] < llimit[i]) or (params[i] > ulimit[i])):
			return (ydata - gauss(x=xdata, x_0=x0, A=params[0], var=params[1]))*sqrtivar + 10 
	return (ydata - gauss(x=xdata, x_0=x0, A=params[0], var=params[1]))*sqrtivar

def func2(params, xdata, ydata, sqrtivar, x0, llimit, ulimit):
	if ( not(10 > abs(params[3]) > 1.3) or  not(ulimit[0]>params[0] > llimit[0]) or not(ulimit[1]>params[1] > llimit[1]) or not(ulimit[0] > params[2] > llimit[0]) ):
			return (ydata - gauss(x=xdata, x_0=x0, A=params[0], var=params[1])-gauss(x=xdata, x_0=x0-params[3], A=params[2], var=params[1]))*sqrtivar + 10 
	return (ydata - gauss(x=xdata, x_0=x0, A=params[0], var=params[1])-gauss(x=xdata, x_0=x0-params[3], A=params[2], var=params[1]))*sqrtivar

def chi2g(params, xdata, ydata, ivar, x0):
	return sum(ivar*(ydata - gauss(x=xdata, x_0=x0, A=params[0], var=params[1]))**2)/(len(xdata)-len(params)-1)


def chi2(params, xdata, ydata, ivar, x0):
	return sum(ivar*(ydata - gauss(x=xdata, x_0=x0, A=params[0], var=params[1])-gauss(x=xdata, x_0=x0-params[3], A=params[2], var=params[1]))**2)/(len(xdata)-len(params) -1)

# Check if x0 is near any emission line redshifted by z
def nearline(x0, zline, fiberid, z, mjd, plate):
	match1 = n.logical_and(abs(zline['linewave']*(1+zline['linez']) -x0) < 3*zline['lineew'], zline['linearea']!=0)
	match2 = n.logical_and(zline['fiberid']==fiberid,zline['mjd']==int(mjd))
	match3 = n.logical_and(zline['plate']==int(plate), zline['linearea']/zline['linearea_err'] > 2)
	match4 = n.logical_and(match1,n.logical_and(match2,match3))
	if (n.sum(match4)>0):
		return True
	else:
		return False

# Check if a path exists, if not make it
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

#-----------------------------------------------------------------------------------------------------
# Set topdir:
topdir = '..'

#print " "*60, "\ropening platelist.fits\r",
#hdulist = pf.open('/fits_files/plates-dr12.fits')
#platelist = hdulist[1].data
#hdulist.close()
#hdulist = 0

## import list of plates to analyze
print " "*20, "importing list of plates",
#plate_mjd = [line.strip() for line in open('Stripe82.platelist.txt')]
plate_mjd = [line.strip().split() for line in open(topdir + '/fits_files/test_mjd.txt')]

#for i in n.arange(len(plate_mjd)):
	#for k in n.arange(len(platelist['plate'])):
		#if (int(plate_mjd[i]) == platelist['plate'][k]):
			#plate_mjd[i] = [int(plate_mjd[i]) , platelist['mjd'][k]]
			#break

plate = 0
fiberid = [0]

f = open(topdir + '/candidates.txt','a')
f.write('Type RA DEC plate mjd fiber peak_wavelength peak_amp peak_amp_err peak_width peak_width_err  peak_number\n')
f.close()

#Set of emission lines used for lensed galaxy detection OII, Hb, OIII, OIII, Ha
em_lines = n.array([3726.5,4861.325,4958.911,5006.843,6562.801])


#Loop over plates
for j in n.arange(len(plate_mjd)):
	#print ' '*60, '\r', "initialization plate " , plate_mjd[j][0], "\r",
	print "initialization plate " , plate_mjd[j][0], "\n",
	
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
	synflux = hdulist[2].data
	#zstruc = hdulist[1].data
	fiberid = hdulist[1].data.field('FIBERID')
	RA = hdulist[1].data.field('PLUG_RA')
	DEC = hdulist[1].data.field('PLUG_DEC')
	obj_class = hdulist[1].data.field('CLASS')
	z = hdulist[1].data.field('Z')
	zwarning = hdulist[1].data.field('ZWARNING')
	z_err = hdulist[1].data.field('Z_ERR')
	hdulist.close()
	hdulist = 0
	hdulist = pf.open(zlfile)
	zline = hdulist[1].data
	hdulist.close()
	hdulist = 0
	##### PlotSpec ends here 
	#-----------------------------------------------------------------------------------------------------

	reduced_flux = n.array(flux - synflux)
	detected = []
	sqrtivar=copy.deepcopy(ivar)
	
	#Upper and lower limit on amplitude and variance of peaks
	llimit = [0.1, 2]
	ulimit = [20, 15]
	
	# Masks atmosphere
	ivar[:,542:551] = 0 # Hg line
	ivar[:,868:873] = 0 # Hg line
	ivar[:,1847:1852] = 0 # Hg line
	ivar[:,1938:1944] = 0 # OI line
	ivar[:,2022:2031] = 0 # lines
	ivar[:,2175:2185] = 0 # lines
	ivar[:,423:428] = 0 # Ca absorption
	ivar[:,460:467] = 0 # Ca absorption
	
	#Mask 5580 A galaxy spectrums
	ivar[:,int((n.log10(5570)-c0)/c1) : int((n.log10(5585)-c0)/c1)] = 0
	
	startTime = datetime.datetime.now()

	# Loop over objects
	for i in n.arange(len(flux[:,0])):
		if (obj_class[i] == 'STAR  ' or obj_class[i] == 'QSO   ' or zwarning[i]!=0):
			continue
		peaks = []
		peaks_err = []
		peak_number = len(peaks)
		searchpeaks = True
		doublet = None
		
		sqrtsaved = n.sqrt(ivar[i,:])
		sqrtivar[i,:] = n.sqrt(ivar[i,:])
		
		# Initial guess
		init = [3,5]
		# Initial guess OII doublet
		init2 = [1,4,1,2.6]
		
		#peak candidates: a priori search if nearline foreground emission line or not
		peak_candidates = n.array([(x0,0.0,0.0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(wave,reduced_flux[i,:]* sqrtivar[i,:]) if test>5])
		# Legend key x0: wavelength chisq_doublet amp_gauss var_gauss err_amp err_var S/N chisq_doublet amp1_doublet amp2_doublet var_doublet err1_doublet err2_doublet err_var_doublet
		#Keep only center of candidate peak
		k = 0
		while (k < (len(peak_candidates)-1)):
			if (abs(peak_candidates[k][0] - peak_candidates[k+1][0]) < 5):
				if peak_candidates[k][6] < peak_candidates[k+1][6]:
					peak_candidates = n.delete(peak_candidates,k, axis=0)
					k = k-1
				else:
					peak_candidates = n.delete(peak_candidates,k+1,axis = 0)
					k = k-1					
			k = k+1

		for peak in peak_candidates:
			x0 = peak[0]
			if nearline(x0, zline, fiberid[i], z[i], int(mjd), int(plate)):
				continue
			cov=0
			#Single Line: Gaussian fit around x_0
			params, cov, infodict, mesg, ier = leastsq(func, init,
				 args=(wave, reduced_flux[i,:], sqrtivar[i,:], x0, llimit, ulimit),
				 full_output=True)
			res =  minimize(chi2g,init,args=(wave, reduced_flux[i,:],ivar[i,:],x0), method='TNC', bounds = [(0.2,10),(2,10)])
			params = res.x
			residue_squared=((reduced_flux[i,:] - gauss(x=wave, x_0=x0, A=params[0], var=params[1]))**2)*ivar[i,:]
			chisq =  sum(residue_squared)/(len(wave)-len(params) -1)
			#Check for S/N > 4 fit
			if not(cov is None or abs(params[0]/math.sqrt(cov[0,0])) < 4):
				peak[1] = chisq
				peak[2] = params[0]
				peak[3] = params[1]
				peak[4] = cov[0,0]
				peak[5] = cov[1,1]	
				
			#Doublet OII: Gaussian fit around x_0
			if (x0 > 3727*(1+z[i])):
				init2= [peak[6],5,peak[6],-3]
				#params, cov, infodict, mesg, ier = leastsq(func2, init2,
				#	args=(wave, reduced_flux[i,:], sqrtivar[i,:], x0, llimit, ulimit),
				#	full_output=True)
				res = minimize(chi2,init2,args=(wave, reduced_flux[i,:],ivar[i,:],x0), method='TNC', bounds = [(0.2,10),(2,10),(0.2,10),(-3,3)])
				params = res.x
				residue_squared=((reduced_flux[i,:] - gauss(x=wave, x_0=x0, A=params[0], var=params[1]) -gauss(x=wave, x_0=x0-params[3], A=params[2], var=params[1]))**2)*ivar[i,:]
				chisq = sum(residue_squared)/(len(wave)-len(params) -1)
				
				
			
				#if (not(cov is None)): #or abs(params[0]/math.sqrt(cov[0,0])) < 4 or abs(params[2]/math.sqrt(cov[2,2])) < 4):
				if abs(params[3])>1.5:	
					peak[7] = chisq
					peak[8] = params[0] #amp1
					peak[9] = params[1] #var
					peak[10] = params[2] #amp2
					peak[11] = params[3] #delta x
					peak[12] = 1.0 #cov[0,0] #err amp1
					peak[13] = 1.0#cov[2,2] #err var
					peak[14] = 1.0#cov[1,1] #err amp2
						
		#Finding peak with highest chi square for doublet and see if it is better fitted by single line or not
		doublet_index = 0
		chi2saved = 1000.0
		# Find the doublet index
		for k in range(len(peak_candidates)):
			peak = peak_candidates[k]
			if (peak[1]>peak[7]>0 and peak[7]< chi2saved):
				peak[1] = peak[7]
				chi2saved = peak[7]
				doublet = True
				doublet_index = k		
		
		#Removing candidates that were not fitted : params still 0
		peak_candidates = n.array([peak for peak in peak_candidates if (peak[2]!=0 or peak[7]==peak[1]!=0)])
		if len(peak_candidates) == 0:
			continue
		
		#Sorting candidates by chi square
		peak_candidates = sorted(peak_candidates, key=lambda peak: peak[1])
		
		# Keeping only 5 most likely candidates
		if len(peak_candidates) > 5:
			peak_candidates = peak_candidates[0:5]
		#find again doublet index
		found = False
		for k in range(len(peak_candidates)):
			if (peak_candidates[k][7] == chi2saved):
				doublet_index = k
				found = True
		if found==False:
			doublet = False
			
		# Check that at least 1 candidate is below 9000 Angstrom cut, if not, go to next fiber
		below_9000 = False
		for peak in peak_candidates:
			if peak[0] < 9000:
				below_9000 = True
		if below_9000 ==False:
			continue
			
		#Try to infer background redshift
		#Generating all combinations of lines from above list to compare with candidates
		detection = False
		if (doublet == True and len(peak_candidates)<2):
			z_background = peak_candidates[doublet_index][0]/3727.24 - 1.0
			if (z_background > z[i]+0.05):
				detection = True
				print '(Doublet) Lensed object at z1 = ', z_background 
		elif (doublet == True and len(peak_candidates)>1):
			z_background = peak_candidates[doublet_index][0]/3727.24 - 1.0
			temp = [peak for peak in peak_candidates if peak[1]!=peak[7]]
			compare = em_lines[1:5]
			if z_background > z[i]+0.05 :
				for peak in temp:
					for line in compare:
						if ( abs(peak[0]/line -1 - z_background) < 0.01):
							detection = True
							print '(Doublet+Multi) Lensed object at z1 = ', z_background, 'em_line: ', line 				
		elif (doublet != True and len(peak_candidates) > 1 ):
			compare = it.combinations(em_lines,len(peak_candidates))
			for group in compare:
				for k in range(len(peak_candidates)):
					for j in range(k+1,len(peak_candidates)):
						if ( abs(peak_candidates[k][0]/group[k] - peak_candidates[j][0]/group[j]) < 0.01 and peak_candidates[k][0]/group[k]-1.0 > (z[i] + 0.05) ):
							detection = True
							print '(Multi Lines) Lensed object at z = ', peak_candidates[k][0]/group[k] -1.0, 'em_lines: ', group[k], group[j]
							
		# Save surviving candidates
		for k in range(len(peak_candidates)):
			peak = peak_candidates[k]
			x0_saved = peak[0] #peak wavelength
			#Redefine variables for clarity
			chisq_saved = peak[1]
			params = [peak[2],peak[3]] #amplitude, variance
			err_amp = peak[4]
			err_var = peak[5]
			if (k == doublet_index and doublet ==True):
				#Specify special variables for one of two peaks
				params = [peak[8],peak[9]]
				err_amp = peak[12]
				err_var = peak[13]
				# Save doublet gaussians
				peaks.append([peak[0], peak[8], peak[9]])
				peaks_err.append([math.sqrt(peak[12]*chisq_saved), math.sqrt(peak[13]*chisq_saved)])
				peaks.append([peak[0]-peak[11], peak[10], peak[9]])
				peaks_err.append([math.sqrt(peak[14]*chisq_saved), math.sqrt(peak[13]*chisq_saved)])
			else:
				# Save emission line 
				peaks.append([x0_saved, params[0], params[1]])
				peaks_err.append([math.sqrt(err_amp*chisq_saved), math.sqrt(err_var*chisq_saved)])
			# Set ivar = 0 for points around peak/emission line found
			delta = 2*int((math.log(1 + math.sqrt(params[1])/x0_saved)/math.log(10)) / c1)
			if (delta == 0):
				delta=1
			center = int(((math.log(x0_saved)/math.log(10))-c0)/c1)
			ivar[i][max(center-delta,0):min(center+delta,len(sqrtivar[i,:])-1)+1] = 0
		
		peak_number = len(peak_candidates)
		print "plate ", plate, " fiber ", fiberid[i], " peak ", peak_number,  "\n",		
			
		#Graphs
		if ((peak_number>1 or doublet==True) and below_9000 and detection):
			#p.subplot(2,1,1)
			#p.plot(wave, flux[i,:], 'k', hold=False)
			#p.plot(wave, synflux[i,:], 'g', hold=True)
			#p.plot(wave, ivar[i,:], 'r', hold=True)
			p.title('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+
				', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+
				'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))
			#p.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
			
			ax = p.subplot(1,1,1)
			#p.plot(wave, reduced_flux[i,:]*sqrtsaved,'k', hold=False)
			p.plot(wave, reduced_flux[i,:],'k', hold=True)
			#p.errorbar(wave, reduced_flux[i,:], yerr= 1.0/sqrtivar[i,:])
			p.xlabel('$Wavelength\, (Angstroms)$')
			p.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
			#Save candidate coordinates and peaks
			fit=0
			for k in n.arange(len(peaks)):
				detected.append([RA[i], DEC[i], int(plate), int(mjd), fiberid[i],
						peaks[k][0], peaks[k][1], math.sqrt(peaks[k][2]), peaks_err[k][0], n.sqrt(peaks_err[k][1]), peak_number])
				fit = fit + gauss(wave, x_0 = peaks[k][0] , A=peaks[k][1], var=peaks[k][2])
				#p.annotate(str(int(peaks[k][1]/peaks_err[k][0]))+' sig \n|',
				#		xy=(peaks[k][0], reduced_flux[i,int(((math.log(peaks[k][0])/math.log(10))-c0)/c1)]),
				#		xytext=(peaks[k][0], reduced_flux[i,int(((math.log(peaks[k][0])/math.log(10))-c0)/c1)] + 3),
				#		ha='center', fontsize=10)
			p.plot(wave,fit,'r',hold=True)
			ax.set_ylim(ymin = -0.5, ymax = 5)
			print ' '*60, '\r', "saving figure\r",
			make_sure_path_exists(topdir + '/plots/')
			p.savefig(topdir + '/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid[i]) + '.png')
			#p.show()
	
	print 'Time taken ', (datetime.datetime.now() - startTime)
	detected = sorted(detected, key = lambda obj : obj[10], reverse = True)
	f = open(topdir + '/candidates.txt','a')
	for item in detected:
		if (doublet == True and len(peak_candidates)<2):
			f.write('Doublet: ' +  "\n" + str(item))
		elif (doublet == True and len(peak_candidates)>1):
			f.write('D + M : ' +  "\n" + str(item))
		else:
			f.write('Multi: ' +  "\n" + str(item))
	f.close()
