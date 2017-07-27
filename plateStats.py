import os
import logging
import logging.config
import numpy as np


def warningReport(w, wLevel, prefix):
    bnw = 5
    bn = np.arange(3000.0, 9200.0, bnw)
    histg = np.histogram(w, bins=bn)
    for each, eachW in zip(histg[0], histg[1][:-1]):
        if each > wLevel * 1000.0:
            logging.getLogger("root")
            logging.warning(prefix + " waverange: " + str(int(eachW)) + "-" +
                            str(int(eachW + bnw)) + " " + str(each) + " hits")


def plateStats(plate, mjd, savedir, wLevel=0.01):
    logging.config.fileConfig("plateStats.conf")
    sd = os.path.join(savedir, str(plate) + "-" + str(mjd),
                      "candidates_doublet.txt")
    try:
        allData = np.loadtxt(sd)
        z = allData[:, 5]
        w = allData[:, 9]
        warningReport(w, wLevel, str(int(plate)) + " " + str(int(mjd)) +
                      " Obs frame")
        warningReport(w / (1.0 + z), wLevel, str(int(plate)) + " " +
                      str(int(mjd)) + " Rst frame")
    except Exception as reason:
        print(str(reason))
