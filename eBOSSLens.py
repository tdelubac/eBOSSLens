# Removes the fitted template from the spectra
# Fits the peaks with a Gaussian profile
# Based on PlotSpec.py

# TO DO
# Sort objects by # significant peaks				Not here
# Save objects that have at least 1 peak below 9000 A		XXX
# Plots in 1 folder, plate-mjd-fiber.png			XXX
# Reject 2 sigma around emission lines instead of 10 A		XXX




# Imports
import numpy as n
import pyfits as pf
import matplotlib as mpl
from matplotlib import pyplot as p
import math
from scipy.optimize import leastsq
import datetime
import copy
from matplotlib import rcParams
#rcParams['text.usetex'] = True
import os
import errno



#-----------------------------------------------------------------------------------------------------
#
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


# Check if x0 is near any emission line redshifted by z
def nearline(x0, zline, fiberid, z, mjd, plate):
	t = datetime.datetime.now()
	
	match1 = n.logical_and((x0 - zline['linewave']*(1+z)) < 2*zline['lineew'], zline['linearea']!=0)
	match2 = n.logical_and(zline['fiberid']==fiberid,zline['mjd']==int(mjd))
	match3 = n.logical_and(zline['plate']==int(plate), zline['linearea']/zline['linearea_err'] > 2)
	match4 = n.logical_and(match1,n.logical_and(match2,match3))
	if (n.sum(match4)>0):
		#test = True 
		print datetime.datetime.now() - t
		return True
	else:
		#test = False
		print datetime.datetime.now() - t
		return False
	#for i in n.arange(len(zline['fiberid'])):
		#if ((x0 - zline['linewave'][i]*(1+z)) < 2*zline['lineew'][i] and zline['linearea'][i]!=0 and 
				#zline['fiberid'][i]==fiberid and zline['mjd'][i]==int(mjd) and
				#zline['plate'][i]==int(plate) and zline['linearea'][i]/zline['linearea_err'][i] > 2):
			#print test == True, '\n'
			#return True	
	#print test==False, '\n'
	#return False


# Check if a path exists, if not make it
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

#-----------------------------------------------------------------------------------------------------

#print " "*60, "\ropening platelist.fits\r",
#hdulist = pf.open('platelist.fits')
#platelist = hdulist[1].data
#hdulist.close()
#hdulist = 0

## import list of plates to analyze
#print " "*60, "\rimporting list of plates\r",
#plate_mjd = [line.strip() for line in open('Stripe82.platelist.txt')]
plate_mjd = [line.strip() for line in open('test_mjd.txt')]

#for i in n.arange(len(plate_mjd)):
	#for k in n.arange(len(platelist['plate'])):
		#if (int(plate_mjd[i]) == platelist['plate'][k]):
			#plate_mjd[i] = [int(plate_mjd[i]) , platelist['mjd'][k]]
			#break


plate = 0
fiberid = [0]
peak_number = 0
i = 0
# Set topdir:
topdir = '.'

f = open(topdir + '/candidates.txt','a')
f.write('RA DEC plate mjd fiber peak_wavelength peak_amp peak_amp_err peak_width peak_width_err\n')
f.close()


#Loop over plates
for j in n.arange(len(plate_mjd)):
	#print ' '*60, '\r', "initialization plate " , plate_mjd[j][0], "\r",
	print "initialization plate " , plate_mjd[j][0], "\n",
	
	# Initialization
	flux = 0
	
	# Pick your plate/mjd and read the data:
	#plate = plate_mjd[j][0]
	#mjd = plate_mjd[j][1]
	plate = 4191
	mjd = 55444
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
	z_err = hdulist[1].data.field('Z_ERR')
	hdulist.close()
	hdulist = 0
	hdulist = pf.open(zlfile)
	zline = hdulist[1].data
	hdulist.close()
	hdulist = 0
	##### PlotSpec ends here 
	
	reduced_flux = flux - synflux
	detected = []
	detected_err = []
	sqrtivar=copy.deepcopy(ivar)
	
	llimit = [0, 10]
	ulimit = [1000, 2500]
	
	
	# Masks
	ivar[:,542:551] = 0 # Hg line
	ivar[:,868:873] = 0 # Hg line
	ivar[:,1847:1852] = 0 # Hg line
	ivar[:,1938:1944] = 0 # OI line
	ivar[:,2022:2031] = 0 # lines
	ivar[:,2175:2185] = 0 # lines
	ivar[:,423:428] = 0 # Ca absorption
	ivar[:,460:467] = 0 # Ca absorption
	
	
	startTime = datetime.datetime.now()
	#print ' '*60, '\r', "initialization done\r",
	i = 0
	# Loop over objects
	for i in n.arange(len(flux[:,0])):
		if (obj_class[i] == 'STAR  ' or obj_class[i] == 'QSO   '):
			continue
		peaks = []
		peaks_err = []
		peak_number = len(peaks)
		below_9000 = False
		
		
		while (peak_number > -1):
			#print ' '*60, '\r', "plate ", plate, " fiber ", fiberid[i], " peak ", peak_number+1, "\r",
			print "plate ", plate, " fiber ", fiberid[i], " peak ", peak_number, "\n",
			#for k in n.arange(len(ivar[i,:])):
			#	sqrtivar[i,k]=math.sqrt(ivar[i,k])
			sqrtivar = n.sqrt(ivar)
			### Test for new solution sqrtivar
			#print '\n', sum(sqrtivar[i,:] == sqrtivar2[i,:])== len(sqrtivar[i,:])
			
			
			# Initial guess
			init = [1, 30]
			chisq_saved=1000
			x0_saved=0
			sig_points = []
			for k in n.arange(len(wave)):
				if (abs(reduced_flux[i,k]*sqrtivar[i,k]) > 5):
					sig_points.append(k)
				
			
			# Loop over all data points
			for k in n.arange(len(sig_points)):
				x0=wave[sig_points[k]]
				#print ' '*60, '\r', "plate ", plate, " fiber ", fiberid[i], " peak ", peak_number+1, "| ", len(sig_points)-k, " points left\r",
				#print "plate ", plate, " fiber ", fiberid[i], " peak ", peak_number+1, "| ", len(sig_points)-k, " points left\n",
				
				
				cov=0
				params, cov, infodict, mesg, ier = leastsq(func, init,
					 args=(wave, reduced_flux[i,:] * (ivar[i,:]>0), sqrtivar[i,:], x0, llimit, ulimit),
					 full_output=True)
	
				residue_squared=(((reduced_flux[i,:] - gauss(x=wave, x_0=x0, A=params[0], var=params[1]))**2)*ivar[i,:])/(len(wave)-3)
				chisq = sum(residue_squared)
				
				if (chisq<chisq_saved):
					x0_saved = x0
					chisq_saved = chisq
					#print chisq, x0, 'saved'
					
			# Gaussian fit around x0_saved
			params, cov, infodict, mesg, ier = leastsq(func, init, 
					args=(wave, reduced_flux[i,:], sqrtivar[i,:], x0_saved, llimit, ulimit), 
					full_output=True)
			
			if (cov is None or abs(chisq_saved*params[0]/math.sqrt(cov[0,0])) < 4):
				if (peak_number == 0):
					#print 'No peak found'
					peak_number = -2
					continue
				else:
					peak_number = -1
					continue

			
			if (nearline(x0_saved, zline, fiberid[i], z[i], int(mjd), int(plate))):
				delta = 2*int((math.log(1 + math.sqrt(params[1])/x0_saved)/math.log(10)) / c1)
				if (delta == 0):
					delta=1
				center = int(((math.log(x0_saved)/math.log(10))-c0)/c1)
				ivar[i][max(center-delta,0):min(center+delta,len(sqrtivar[i,:])-1)+1] = 0
				continue
			if (x0_saved < 9000):
				below_9000 = True
			
			#print 'Peak number', peak_number+1
			
			peaks.append([x0_saved, params[0], params[1]])
			#chisq = sum((((reduced_flux[i,:] - gauss(x=wave, x_0=x0, A=params[0], var=params[1]))**2)*ivar[i,:]))/(len(wave)-3)
			peaks_err.append([math.sqrt(cov[0,0]*chisq_saved), math.sqrt(cov[1,1]*chisq_saved)])
			
			
			#print 'x0 = ', x0_saved
			#print 'A = ', params[0], 'err ', math.sqrt(cov[0,0]*chisq)
			#print 'var = ', params[1], 'sigma = ', math.sqrt(params[1]), 'err ', math.sqrt(math.sqrt(cov[1,1]*chisq))
			#print 'chi_sq = ', chisq
			#print 'cov = ', cov * chisq
			
			# Set ivar = 0 for points around peaks
			delta = 2*int((math.log(1 + math.sqrt(params[1])/x0_saved)/math.log(10)) / c1)
			if (delta == 0):
				delta=1
			
			center = int(((math.log(x0_saved)/math.log(10))-c0)/c1)
			ivar[i][max(center-delta,0):min(center+delta,len(sqrtivar[i,:])-1)+1] = 0
			
			peak_number = peak_number+1
		
		
		
		p.subplot(2,1,1)
		p.plot(wave, flux[i,:] * (ivar_copy[i,:]>0), 'k', hold=False)
		p.plot(wave, synflux[i,:], 'g', hold=True)
		p.title('RA='+str(RA[i])+', Dec='+str(DEC[i])+', Plate='+str(plate)+
			', Fiber='+str(fiberid[i])+', MJD='+str(mjd)+
			'\n$z='+str(z[i])+' \pm'+str(z_err[i])+'$, Class='+str(obj_class[i]))
		p.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
		
		p.subplot(2,1,2)
		p.plot(wave, reduced_flux[i,:] *(ivar_copy[i,:]>0),'k', hold=False)
		for k in n.arange(len(ivar_copy[i,:])):
			sqrtivar[i,k]=math.sqrt(ivar_copy[i,k])
		#p.errorbar(wave, reduced_flux[i,:], yerr= 1.0/sqrtivar[i,:])
		p.xlabel('$Wavelength\, (Angstroms)$')
		p.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$)')
		
		fit=0
		if (len(peaks)>1 and below_9000):
			for k in n.arange(len(peaks)):
				detected.append([RA[i], DEC[i], int(plate), int(mjd), fiberid[i],
						peaks[k][0], peaks[k][1], math.sqrt(peaks[k][2])])
				detected_err.append([peaks_err[k][0], math.sqrt(peaks_err[k][1])])
				fit = fit + gauss(wave, x_0 = peaks[k][0] , A=peaks[k][1], var=peaks[k][2])
				p.annotate(str(int(peaks[k][1]/peaks_err[k][0]))+' sig \n|',
						xy=(peaks[k][0], reduced_flux[i,int(((math.log(peaks[k][0])/math.log(10))-c0)/c1)]),
						xytext=(peaks[k][0], reduced_flux[i,int(((math.log(peaks[k][0])/math.log(10))-c0)/c1)] + 3),
						ha='center', fontsize=10)
			p.plot(wave,fit,'g',hold=True)
			print ' '*60, '\r', "saving figure\r",
			make_sure_path_exists(topdir + '/plots/')
			p.savefig(topdir + '/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid[i]) + '.png')
			#p.show()
	
	print 'Time taken ', (datetime.datetime.now() - startTime)
	
	f = open(topdir + '/candidates.txt','a')
	for i in n.arange(len(detected)):
		f.write(str(detected[i][0])+ ' '+ str(detected[i][1])+ ' '+ str(detected[i][2])+ ' '+ 
			str(detected[i][3])+ ' '+ str(detected[i][4])+' '+ str(detected[i][5])+ ' '+ 
			str(detected[i][6])+ ' '+ str(detected_err[i][0])+' ' + 
			str(detected[i][7])+ ' '+ str(detected_err[i][1])+ '\n')
	f.close()