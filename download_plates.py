import os
import urlparse
import urllib
import numpy as n


def create_directory(dire):
    '''
    create_directory(dire)
    ======================
    Input: dire: a relative directory that needs to be created
    Output: None
    '''
    direlast = os.path.dirname(dire)
    if not os.path.exists(direlast):
        os.makedirs(direlast)


def download_plates(is570):
    # Download_list.txt should be written plain as two columns:
    # "plate_number" "mjd"
    plate_mjd = [line.strip().split() for line in open('download_list.txt')]
    # Update to dr12
    # Add 5.7.0 and 5.7.2
    if is570:
        flag = "v5_7_0"
    else:
        flag = "v5_7_2"
    baseDir = "../SCRATCH/BOSS/data/" + flag
    baseUrl = "http://data.sdss3.org/sas/dr12/sdss/spectro/redux/" + flag + "/"
    for i in n.arange(len(plate_mjd)):
        # Construct file
        fileSfx = plate_mjd[i][0] + "-" + plate_mjd[i][1] + ".fits"
        sppFile = plate_mjd[i][0] + "/spPlate-" + fileSfx
        spbFile = plate_mjd[i][0] + "/" + flag + "/spZbest-" + fileSfx
        splFile = plate_mjd[i][0] + "/" + flag + "/spZline-" + fileSfx
        # Construct saving directories
        sppDir = os.path.join(baseDir, sppFile)
        spbDir = os.path.join(baseDir, spbFile)
        splDir = os.path.join(baseDir, splFile)
        create_directory(sppDir)
        create_directory(spbDir)
        create_directory(splDir)
        # Construct urls
        sppUrl = urlparse.urljoin(baseUrl, sppFile)
        spbUrl = urlparse.urljoin(baseUrl, spbFile)
        splUrl = urlparse.urljoin(baseUrl, splFile)
        # download spPlate
        urllib.URLopener()
        print("Downloading " + str(sppUrl))
        urllib.urlretrieve(sppUrl, sppDir)
        # download spZbest
        print("Downloading " + str(spbUrl))
        urllib.urlretrieve(spbUrl, spbDir)
        # download spZline
        print("Downloading " + str(splUrl))
        urllib.urlretrieve(splUrl, splDir)


if __name__ == "__main__":
    download_plates(True)
