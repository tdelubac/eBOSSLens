import os
import numpy as np
import matplotlib.pyplot as plt
from pysynphot import observation
from pysynphot import spectrum
from SDSSObject import SDSSObject
from peakFinder import bolEstm


def lineAtx(x0, ax, s, ymin, ymax, ytext):
    ax.vlines(x=x0, ymin=ymin, ymax=ymax, colors='g', linestyles='dashed')
    ax.text(x0 + 2.0, ytext, s, rotation=90)


def plotSpec(obj, ax, xlim, ylim):
    ax.plot(obj.wave[10: -10], obj.reduced_flux[10: -10])
    ax.set_xlabel('$\lambda \, [\AA]$ ')
    ax.set_ylabel('$f_{\lambda}\, (10^{-17} erg\, s^{-1} cm^{-2} Ang^{-1}$')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plotSn(obj, ax, s, xlim):
    ax.plot(obj.wave[10: -10], s[10: -10])
    ax.set_xlabel('$\lambda \, [\AA]$ ')
    ax.set_ylabel("SN Ratio")
    ax.set_xlim(xlim)
    ax.set_ylim([-5, 15])


def rebin_spec(wave, specin, wnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin, keepneg=True)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wnew, force='taper')
    return obs.binflux


def findPeak(obj, waveHint, sn, width=10.0):
    tmp = [obj.wave2bin(waveHint - 10.0), obj.wave2bin(waveHint + 10.0)]
    try:
        res = obj.wave[tmp[0]: tmp[1]][np.argmax(sn[tmp[0]: tmp[1]])]
    except Exception:
        res = 0.0
    return res


def plotter(obj, savedir, o2wave=None):
    o3wave2 = 65535
    width = 60.0
    sn = bolEstm(obj, 1.2, width)
    # Decide layout
    lo = ()
    if o2wave is not None:
        lo = (2, 7)
        plt.figure(figsize=(35, 15))
    else:
        lo = (2, 2)
        plt.figure(figsize=(10, 15))
    plt.suptitle('RA=' + str(obj.RA) + ', Dec=' + str(obj.DEC) +
                 ', Plate=' + str(obj.plate) + ', Fiber='+str(obj.fiberid) +
                 ', MJD=' + str(obj.mjd) + '\n$z=' + str(obj.z) + ' \pm' +
                 str(obj.z_err) + '$, Class=' + str(obj.obj_class))
    # Reduced flux overall
    ax1 = plt.subplot2grid(lo, (0, 0), colspan=2)
    xl = [np.min(obj.wave), np.max(obj.wave)]
    plotSpec(obj, ax1, xl, [1.1 * np.min(obj.reduced_flux),
                            1.1 * np.max(obj.reduced_flux)])
    # SN overall
    ax11 = plt.subplot2grid(lo, (1, 0), colspan=2)
    plotSn(obj, ax11, sn, xl)
    if o2wave is not None:
        # Infer other lines
        zsource = o2wave / 3727.0 - 1.0
        hbwave = 4853.0 * (1.0 + zsource)
        o3wave = 4960.0 * (1.0 + zsource)
        o3wave2 = 5008.0 * (1.0 + zsource)
        hawave = 6563.0 * (1.0 + zsource)
        # Annotate on main plot
        yt = 0.9 * np.max(obj.reduced_flux)
        lineAtx(o2wave, ax1, "OII 3727", -100, 100, yt)
        lineAtx(hbwave, ax1, "Hbeta 4853", -100, 100, yt)
        lineAtx(o3wave, ax1, "OIII 4960", -100, 100, yt)
        lineAtx(o3wave2, ax1, "OIII 5008", -100, 100, yt)
        lineAtx(hawave, ax1, "Halpha 6563", -100, 100, yt)
        lineAtx(o2wave, ax11, "OII 3727", -5, 15, 13.5)
        lineAtx(hbwave, ax11, "Hbeta 4853", -5, 15, 13.5)
        lineAtx(o3wave, ax11, "OIII 4960", -5, 15, 13.5)
        lineAtx(o3wave2, ax11, "OIII 5008", -5, 15, 13.5)
        lineAtx(hawave, ax11, "Halpha 6563", -5, 15, 13.5)
        # OII reduced flux detail
        o2xl = [o2wave - 30.0, o2wave + 30.0]
        ax2 = plt.subplot2grid(lo, (0, 2))
        plotSpec(obj, ax2, o2xl, [-5, 10])
        lineAtx(o2wave, ax2, "OII 3727", -5, 10, 9)
        # OII SN detail
        ax12 = plt.subplot2grid(lo, (1, 2))
        plotSn(obj, ax12, sn, o2xl)
        lineAtx(o2wave, ax12, "OII 3727", -5, 15, 13.5)
        # Hbeta reduced flux detail
        hbxl = [hbwave - 30.0, hbwave + 30.0]
        hbfound = findPeak(obj, hbwave, sn)
        ax3 = plt.subplot2grid(lo, (0, 3))
        plotSpec(obj, ax3, hbxl, [-5, 10])
        lineAtx(hbwave, ax3, "Hbeta 4853", -5, 10, 9)
        lineAtx(hbfound, ax3, "Found", -5, 10, 9)
        # Hbeta SN detail
        ax13 = plt.subplot2grid(lo, (1, 3))
        plotSn(obj, ax13, sn, hbxl)
        lineAtx(hbwave, ax13, "Hbeta 4853", -5, 15, 13.5)
        lineAtx(hbfound, ax13, "Found", -5, 15, 13.5)
        # O3 reduced flux detail
        o3xl = [o3wave - 30.0, o3wave + 30.0]
        o3found = findPeak(obj, o3wave, sn)
        ax4 = plt.subplot2grid(lo, (0, 4))
        plotSpec(obj, ax4, o3xl, [-5, 10])
        lineAtx(o3wave, ax4, "OIII 4960", -5, 10, 9)
        lineAtx(o3found, ax4, "Found", -5, 10, 9)
        # O3 SN detail
        ax14 = plt.subplot2grid(lo, (1, 4))
        plotSn(obj, ax14, sn, o3xl)
        lineAtx(o3wave, ax14, "OIII 4960", -5, 15, 13.5)
        lineAtx(o3found, ax14, "Found", -5, 15, 13.5)
        # O3 reduced flux detail
        o3xl = [o3wave2 - 30.0, o3wave2 + 30.0]
        o3found = findPeak(obj, o3wave2, sn)
        ax5 = plt.subplot2grid(lo, (0, 5))
        plotSpec(obj, ax5, o3xl, [-5, 10])
        lineAtx(o3wave2, ax5, "OIII 5008", -5, 10, 9)
        lineAtx(o3found, ax5, "Found", -5, 10, 9)
        # O3 SN detail
        ax15 = plt.subplot2grid(lo, (1, 5))
        plotSn(obj, ax15, sn, o3xl)
        lineAtx(o3wave2, ax15, "OIII 5008", -5, 15, 13.5)
        lineAtx(o3found, ax15, "Found", -5, 15, 13.5)
        # Halpha reduced flux detail
        haxl = [hawave - 30.0, hawave + 30.0]
        hafound = findPeak(obj, hawave, sn)
        ax6 = plt.subplot2grid(lo, (0, 6))
        plotSpec(obj, ax6, haxl, [-5, 10])
        lineAtx(hawave, ax6, "Halpha 6563", -5, 10, 9)
        lineAtx(hafound, ax6, "Found", -5, 10, 9)
        # O3 SN detail
        ax16 = plt.subplot2grid(lo, (1, 6))
        plotSn(obj, ax16, sn, haxl)
        lineAtx(hawave, ax16, "Halpha 6563", -5, 15, 13.5)
        lineAtx(hafound, ax16, "Found", -5, 15, 13.5)
    plt.savefig(os.path.join(savedir, str(obj.plate) + '-' + str(obj.mjd) +
                             '-' + str(obj.fiberid) + '.png'))
    plt.close()
    w = obj.wave / (1.0 + zsource)
    sn = sn / sn[obj.wave2bin(o2wave)] * 10.0
    print(sn[obj.wave2bin(o2wave)])
    if o3wave2 < 9200.0:
        return [w, sn]
    else:
        return []


pmf = np.loadtxt("toPlot.txt")
bn = np.arange(3000.0, 7000.0, 0.5)
nsn = []
for each in pmf:
    obj = SDSSObject(int(each[0]), int(each[1]), int(each[2]), 'v5_7_0',
                     '../SCRATCH')
    o2w = each[3]
    if each[3] == 0.0:
        o2w = None
    wsn = plotter(obj, '../PlotCheck', o2w)
    if wsn != []:
        nsn.append(rebin_spec(wsn[0], wsn[1], bn))
nsn = np.array(nsn)
y = np.mean(nsn, axis=0)
plt.plot(bn, y)
plt.xlabel("Wavelength")
plt.ylabel("SN Ratio")
plt.vlines(3727.0, np.min(y), np.max(y), linestyles='dashed')
plt.vlines(4853.0, np.min(y), np.max(y), linestyles='dashed')
plt.vlines(4960.0, np.min(y), np.max(y), linestyles='dashed')
plt.vlines(5008.0, np.min(y), np.max(y), linestyles='dashed')
plt.vlines(6563.0, np.min(y), np.max(y), linestyles='dashed')
plt.show()
