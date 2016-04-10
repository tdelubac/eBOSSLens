# Removes the fitted template from the spectra
# Fits the peaks with a Gaussian profile
# Based on PlotSpec.py

# TO DO
# Sort objects by # significant peaks
# plots centered on doublets

# Imports
import numpy as n
import pyfits as pf
import matplotlib as mpl
from matplotlib import pyplot as p
import math
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import datetime
import copy
from matplotlib import rcParams
import os
import errno
import itertools as it

#-----------------------------------------------------------------------------------------------------
# Function definitions

# Generate a Gaussian around x_0 with amplitude A and variance var
def gauss(x,x_0,A,var):
	y = A*n.exp( (-(x-x_0)**2) / (2*var) )
	return y

#Chi square for one gaussian
def chi2g(params, xdata, ydata, ivar, x0):
	return sum(ivar*(ydata - gauss(x=xdata, x_0=x0, A=params[0], var=params[1]))**2)/(len(xdata)-len(params)-1)
#Chi square for Doublet
def chi2D(params, xdata, ydata, ivar):
	return sum(ivar*(ydata - gauss(x=xdata, x_0=params[3], A=params[0], var=params[1])-gauss(x=xdata, x_0=params[4], A=params[2], var=params[1]))**2)/(len(xdata)-len(params) -1)

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

#Gaussian kernel used in first feature search (Bolton et al.,2004 method)	
def kernel(j,width,NormGauss,length):
	ker = n.zeros(length)
	ker[j-int(width*0.5):j+int(width*0.5)] = NormGauss
	return ker
	
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

## import list of plates to analyze
print " "*20, "importing list of plates",
plate_mjd = [line.strip().split() for line in open(topdir + '/fits_files/test_mjd.txt')]

plate = 0
fiberid = [0]

f = open(topdir + '/candidates.txt','a')
f.write('Type RA DEC plate mjd fiber peak_wavelength peak_amp peak_width peak_number\n')
f.close()

#Set of emission lines used for lensed galaxy detection OII, Hb, OIII, OIII, Ha
em_lines = n.array([3726.5,4861.325,4958.911,5006.843,6562.801])
countDoublet = 0
countDM = 0
countMulti = 0

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
	
	# Masks atmosphere
	ivar[:,542:551] = 0 # Hg line
	ivar[:,868:873] = 0 # Hg line
	ivar[:,1847:1852] = 0 # Hg line
	ivar[:,1938:1944] = 0 # OI line
	ivar[:,2022:2031] = 0 # lines
	ivar[:,2175:2185] = 0 # lines
	ivar[:,423:428] = 0 # Ca absorption
	ivar[:,460:467] = 0 # Ca absorption
	
	#Mask 5580 A galaxy spectra
	ivar[:,int((n.log10(5570)-c0)/c1) : int((n.log10(5585)-c0)/c1)] = 0
	
	startTime = datetime.datetime.now()

	# Loop over objects
	for i in n.arange(len(flux[:,0])):
		if (obj_class[i] == 'STAR  ' or obj_class[i] == 'QSO   '): # or zwarning[i]!=0):
			continue
		peaks = []
		peak_number = len(peaks)
		searchpeaks = True
		doublet = None
		
		sqrtsaved = n.sqrt(ivar[i,:])
		sqrtivar[i,:] = n.sqrt(ivar[i,:])
		
		### Bolton 2004: S/N of maximum likelihood estimator of gaussian peaks
		width = 30.0
		sig = 2.0
		## Prepare normalized gaussian
		NormGauss = gauss(n.linspace(-width*0.5,width*0.5,width),0.0,1.0,sig)
		NormGauss = NormGauss/sum(NormGauss)
		Cj1 = n.array([sum(reduced_flux[i,:]*kernel(j+0.5*width,width,NormGauss,len(wave))*ivar[i,:]) for j in range(int(len(wave)-width))])
		Cj2 = n.array([sum(ivar[i,:]*kernel(j+0.5*width,width,NormGauss,len(wave))**2) for j in range(int(len(wave)-width))])
		SN = n.zeros(len(wave))
		SN[width*0.5:len(wave)-width*0.5] = Cj1/n.sqrt(Cj2)
		
		peak_candidates = n.array([(x0,0.0,0.0,0.0,test,0.0,0.0,0.0,0.0,0.0,0.0) for x0,test in zip(wave,SN) if test>6.0])
		# Legend key x0: wavelength(x0) chisq_doublet amp_gauss var_gauss S/N chisq_doublet amp1_doublet amp2_doublet var_doublet x1 x2

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
			
		#Search for suitable peak candidates
		for peak in peak_candidates:
			x0 = peak[0]
			if nearline(x0, zline, fiberid[i], z[i], int(mjd), int(plate)):
				continue
			#Single Line: Gaussian fit around x_0
			init = [1,2]
			res =  minimize(chi2g,init,args=(wave, reduced_flux[i,:],ivar[i,:],x0), method='SLSQP', bounds = [(0.1,5),(1,8)])
			params = res.x
			residue_squared=((reduced_flux[i,:] - gauss(x=wave, x_0=x0, A=params[0], var=params[1]))**2)*ivar[i,:]
			chisq =  sum(residue_squared)/(len(wave)-len(params) -1)
			chigauss1 = chisq
			#Check for not to high chi square and save
			if not(chisq > 2.0):
				peak[1] = chisq
				peak[2] = params[0]
				peak[3] = params[1]
				
			#Doublet OII: Gaussian fit around x_0
			if (x0 > 3727.0*(1+z[i])):

				res2 = minimize(chi2D,[1.0,5,1.0,x0-1.5,x0+1.5],args=(wave, reduced_flux[i,:],ivar[i,:]), method='SLSQP', bounds = [(0.1,5),(1,8),(0.1,5),(x0-7,x0),(x0,x0+7)])
				params2 = res2.x
				residue_squared=((reduced_flux[i,:] - gauss(x=wave, x_0=params2[3], A=params2[0], var=params2[1]) - gauss(x=wave, x_0=params2[4], A=params2[2], var=params2[1]))**2)*ivar[i,:]
				chisq2 = sum(residue_squared)/(len(wave)-len(params2) -1)
				
				if abs(params2[3]-params2[4])>1.5:	
					peak[5] = chisq2
					peak[6] = params2[0] #amp1
					peak[7] = params2[2] #amp2
					peak[8] = params2[1] #var
					peak[9] = params2[3] #x1
					peak[10] = params2[4] #x2
									
		#Finding peak with lowest chi square for doublet and see if it is better fitted by single line or not
		doublet_index = 0
		chi2saved = 1000.0
		# Find the doublet index
		for k in range(len(peak_candidates)):
			peak = peak_candidates[k]
			if (peak[1]>peak[5]>0 and peak[5]< chi2saved):
				peak[1] = peak[5]
				chi2saved = peak[5]
				doublet = True
				doublet_index = k		
		
		#Removing candidates that were not fitted : params still 0
		peak_candidates = n.array([peak for peak in peak_candidates if (peak[2]!=0 or peak[5]==peak[1]!=0)])
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
			if (peak_candidates[k][5] == chi2saved):
				doublet_index = k
				found = True
		if found==False:
			doublet = False
			
		# Check that at least 1 candidate is below 9000 Angstrom cut, if not, go to next fiber
		below_9000 = False
		for peak in peak_candidates:
			if peak[0] < 9500:
				below_9000 = True
		if below_9000 ==False:
			continue
			
		#Try to infer background redshift
		#Generating all combinations of lines from above list to compare with candidates
		detection = False
		if (doublet == True):
			z_background = peak_candidates[doublet_index][0]/3727.24 - 1.0
			if (z_background > z[i]+0.05):
				detection = True
				print '(Doublet) Lensed object at z1 = ', z_background 
			if len(peak_candidates):	
				temp = [peak for peak in peak_candidates if peak[1]!=peak[5]]
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
			if (k == doublet_index and doublet ==True):
				# Save doublet gaussians
				peaks.append([peak[9], peak[6], peak[8]])
				peaks.append([peak[10], peak[7], peak[8]])
			else:
				# Save emission line 
				peaks.append([x0_saved, params[0], params[1]])
			# Set ivar = 0 for points around peak/emission line found
			delta = 2*int((math.log(1 + math.sqrt(params[1])/x0_saved)/math.log(10)) /c1)
			if (delta == 0):
				delta=1
			center = int(((math.log(x0_saved)/math.log(10))-c0)/c1)
			ivar[i][max(center-delta,0):min(center+delta,len(sqrtivar[i,:])-1)+1] = 0
		
		peak_number = len(peak_candidates)
		print "plate ", plate, " fiber ", fiberid[i], " peak ", peak_number,  "\n",		
			
		#Graphs
		if ((peak_number>1 or doublet==True) and below_9000 and detection):
			p.title('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+
				', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+
				'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))
			ax = p.subplot(1,1,1)
			#p.plot(wave, reduced_flux[i,:]*sqrtsaved,'k', hold=False)
			p.plot(wave, reduced_flux[i,:],'k', hold=True)
			p.xlabel('$Wavelength\, (Angstroms)$')
			p.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
			#Save candidate coordinates and peaks
			fit=0
			for k in n.arange(len(peaks)):
				detected.append([RA[i], DEC[i], int(plate), int(mjd), fiberid[i],
						peaks[k][0], peaks[k][1], math.sqrt(peaks[k][2]), peak_number])
				fit = fit + gauss(wave, x_0 = peaks[k][0] , A=peaks[k][1], var=peaks[k][2])

			p.plot(wave,fit,'r',hold=True)
			ax.set_ylim(ymin = -0.5, ymax = 5)
			print ' '*60, '\r', "saving figure\r",
			make_sure_path_exists(topdir + '/plots/')
			p.savefig(topdir + '/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid[i]) + '.png')
			p.show()
			f = open(topdir + '/candidates.txt','a')
			for item in detected:
				if (doublet == True and len(peak_candidates)>1):
					f.write('\n D + M : ' +  "\n" + str(item))
					countDoublet +=1.0/len(peak_candidates)
				elif (doublet == True and len(peak_candidates)<2):
					f.write('\n Doublet: ' +  "\n" + str(item))
					countDoublet +=1 
				elif (doublet != True and len(peak_candidates) > 1) :
					f.write('\n Multi: ' +  "\n" + str(item))
					countMulti =+1.0/len(peak_candidates)
			f.close()
	
	print 'Time taken ', (datetime.datetime.now() - startTime)
	detected = sorted(detected, key = lambda obj : obj[8], reverse = True)

f = open(topdir + '/candidates.txt','a')
f.write('\n Doublet: ' + str(countDoublet) + ' D+M: ' + str(countDM) + ' Multi: ' + str(countMulti)	)
f.close()
