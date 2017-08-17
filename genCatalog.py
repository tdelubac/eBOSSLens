import argparse
import os
import numpy as np
from astropy.io import fits


# Commandline argument parser
parser = argparse.ArgumentParser()
parser.add_argument("dataDir", help="Directory for data", type=str)
parser.add_argument("saveFile", help="Output file", type=str)
args = parser.parse_args()
dataDir = args.dataDir
saveFile = args.saveFile
# Generate catalog
er = []
sc = []
zs = []
ra = []
dec = []
z = []
plate = []
mjd = []
fid = []
w = []
snp = []
snn = []
sns = []
o3w = []
o3s = []
for each in os.listdir(dataDir):
    tmpDir = os.path.join(dataDir, each)
    if os.path.isdir(tmpDir):
        dDir = os.path.join(tmpDir, "candidates_doublet.txt")
        try:
            temp = np.loadtxt(dDir)
            er.extend(temp[:, 0])
            sc.extend(temp[:, 1])
            zs.extend(temp[:, 2])
            ra.extend(temp[:, 3])
            dec.extend(temp[:, 4])
            z.extend(temp[:, 5])
            plate.extend(temp[:, 6])
            mjd.extend(temp[:, 7])
            fid.extend(temp[:, 8])
            w.extend(temp[:, 9])
            snp.extend(temp[:, 10])
            snn.extend(temp[:, 11])
            sns.extend(temp[:, 12])
            o3w.extend(temp[:, 13])
            o3s.extend(temp[:, 14])
        except Exception:
            pass
# Save to file
c0 = fits.Column(name="einsteinRadius", array=np.array(er), format='D')
c1 = fits.Column(name='score', array=np.array(sc), format='D')
c2 = fits.Column(name='zsource', array=np.array(zs), format='D')
c3 = fits.Column(name='ra', array=np.array(ra), format='D')
c4 = fits.Column(name='dec', array=np.array(dec), format='D')
c5 = fits.Column(name='z', array=np.array(z), format='D')
c6 = fits.Column(name='plate', array=np.array(plate), format='D')
c7 = fits.Column(name='mjd', array=np.array(mjd), format='D')
c8 = fits.Column(name='fiberid', array=np.array(fid), format='D')
c9 = fits.Column(name='wavelength', array=np.array(w), format='D')
ca = fits.Column(name='snPrevFiber', array=np.array(snp), format='D')
cb = fits.Column(name='snNextFiber', array=np.array(snn), format='D')
cc = fits.Column(name='snSpectra', array=np.array(sns), format='D')
cd = fits.Column(name='o3FoundWave', array=np.array(o3w), format='D')
ce = fits.Column(name='o3FoundSig', array=np.array(o3s), format='D')
t = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, ca, cb,
                                   cc, cd, ce])
t.writeto(saveFile, overwrite=True)
