import urllib
import pyfits as pf
import numpy as n


dir = "http://data.sdss3.org/sas/dr10/sdss/spectro/redux/v5_5_12/"

# Download_list.txt should be written plain as two columns: "plate_number" "mjd"

plate_mjd = [line.strip().split() for line in open('download_list.txt')]

			
for i in n.arange(len(plate_mjd)):
	print "Downloading spPlate ", plate_mjd[i][1], "\r",
	# download spPlate
	u = urllib.URLopener()
	urllib.urlretrieve( dir + plate_mjd[i][0] + "/spPlate-" + plate_mjd[i][0] + "-"+ plate_mjd[i][1] + ".fits", "./fits_files/spPlate-" + plate_mjd[i][0] + "-"+ plate_mjd[i][1] + ".fits")

	print "Downloading spZbest ", plate_mjd[i][0], "\r",
	# download spZbest
	urllib.urlretrieve( dir + plate_mjd[i][0] + "/v5_5_12/spZbest-" + plate_mjd[i][0] + "-" + plate_mjd[i][1] + ".fits", "./fits_files/spZbest-" + plate_mjd[i][0] + "-"+ plate_mjd[i][1] + ".fits")

	print "Downloading spZline ", plate_mjd[i][0], "\r",
	# download spZline
	urllib.urlretrieve( dir + plate_mjd[i][0] + "/v5_5_12/spZline-" + plate_mjd[i][0] + "-"+ plate_mjd[i][1] + ".fits", "./fits_files/spZline-" + plate_mjd[i][0] + "-"+ plate_mjd[i][1] + ".fits")

	
