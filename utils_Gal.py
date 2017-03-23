import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
from utils import *

def plot_GalaxyLens(doublet,RA,DEC,plate,fiberid,mjd,z,z_err,obj_class,wave, reduced_flux, zline,fit,topdir,savedir, peak_candidates, doublet_index, c0,c1,Nmax, show = False):
	if show==False:
		mpl.use('Agg')
		
	if doublet != True: 
		ax = plt.subplot(1,1,1)
		plt.title('RA='+str(RA)+', Dec='+str(DEC)+', Plate='+str(plate)+', Fiber='+str(fiberid)+', MJD='+str(mjd)+'\n$z='+str(z)+' \pm'+str(z_err)+'$, Class='+str(obj_class))
		ax.plot(wave, reduced_flux[:],'k')
		plt.xlabel('$Wavelength\, (Angstroms)$')
		plt.ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$')
		ax.plot(wave,fit,'r')
		make_sure_path_exists(topdir + savedir +'/plots/')
		plt.savefig(topdir + savedir +'/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid) + '.png')
		if show:
			plt.show()
		else:
			plt.close()
	
	# If doublet, plot in two different windows
	if doublet==True:
		x_doublet = 0.5*(peak_candidates[doublet_index][9]+ peak_candidates[doublet_index][10])
		bounds = np.linspace(wave2bin(x_doublet,c0,c1,Nmax)-10,wave2bin(x_doublet,c0,c1,Nmax)+10,21,dtype = np.int16)
		f = open(topdir + savedir +'/doublet_ML.txt','a')
		f.write('\n' +  str(plate) + ' ' + str(mjd) + ' ' + str(fiberid) + ' ' + str(reduced_flux[bounds]))
		f.close()
		
		plt.figure(figsize=(14,6))
		ax1 = plt.subplot2grid((1,3), (0,0), colspan=2)
		plt.suptitle('RA='+str(RA)+', Dec='+str(DEC)+', Plate='+str(plate)+', Fiber='+str(fiberid)+', MJD='+str(mjd)+'\n$z='+str(z)+' \pm'+str(z_err)+'$, Class='+str(obj_class))
		ax2 = plt.subplot2grid((1,3), (0,2))
		ax1.plot(wave[10:-10], reduced_flux[10:-10],'k')
		ax1.plot(wave,fit,'r')
		ax1.set_xlabel('$\lambda \, [\AA]$ ')
		ax1.set_ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$')
		ax2.set_xlabel('$\lambda \, [\AA]$ ')
		ax2.locator_params(tight=True)
		
		ax2.set_xlim([peak_candidates[doublet_index][0]-30,  peak_candidates[doublet_index][0]+30])
		ax2.plot(wave, reduced_flux[:],'k')
		ax2.plot(wave,fit,'r')
		ax2.set_ylim([-5,10])
		
		ax2.vlines(x = zline['linewave']*(1+z),ymin= -10,ymax= 10,colors= 'g',linestyles='dashed')
		#ax1.vlines(x = zline['linewave']*(1+z),ymin= -10,ymax= 10,colors= 'g',linestyles='dashed')
		ax1.set_xlim([np.min(wave),np.max(wave)])
		
		make_sure_path_exists(topdir + savedir +'/plots/')
		plt.savefig(topdir + savedir +'/plots/' + str(plate) + '-' + str(mjd) + '-' + str(fiberid) + '.png')
		if show:
			plt.show()
		else:
			plt.close()
	return
	
def plot_Jackpot(RA,DEC,plate,fiberid,mjd,z,wave, flux, synflux,topdir,savedir,peak, show ,counter ): 
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
