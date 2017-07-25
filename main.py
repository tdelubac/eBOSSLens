import argparse
import os
import multiprocessing as mul
import numpy as np
from eBOSSLens import eBOSSLens
from utils import make_sure_path_exists


def para_return(func, params, num_thread=4):
    res_list = []
    p = mul.Pool(processes=num_thread)
    for each in params:
        def _log_result(res):
            res_list.append(res)
        p.apply_async(func, each, callback=_log_result)
    p.close()
    p.join()
    return res_list


def lensFinder(plate, mjd, fiberid, datav, datadir, savedir, lya, qso, jpt,
               bwidth, bsig, maxchi2, doplot):
    sd = os.path.join(savedir, str(plate) + "-" + str(mjd))
    make_sure_path_exists(sd)
    try:
        eBOSSLens(plate, mjd, fiberid, datav, lya, qso, jpt, sd, datadir,
                  max_chi2=maxchi2, bwidth=bwidth, bsig=bsig, doPlot=doplot)
    except Exception as reason:
        text = str(plate) + " " + str(mjd) + " " + str(fiberid) + " " + \
            str(reason)
        print(text)


def setArg(argument, defValue):
    if argument is not None:
        return argument
    else:
        return defValue


if __name__ == "__main__":
    '''
    # Command line argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("pmfile", help="The file for plate and mjd input",
                        type=str)
    parser.add_argument("-c", help="Max accepted chisquare", type=float)
    parser.add_argument("-w", help="Width for peak detection", type=float)
    parser.add_argument("-s", help="Sig for peak detection", type=float)
    parser.add_argument("--dataversion", help="Data version string", type=str)
    parser.add_argument("--datadir", help="Data directory string", type=str)
    parser.add_argument("--savedir", help="Save directory string", type=str)
    parser.add_argument("--lya", help="LyA search", action="store_true")
    parser.add_argument("--qso", help="QSO search", action="store_true")
    parser.add_argument("--jpt", help="Jackpot search", action="store_true")
    args = parser.parse_args()
    # Storing command line argument
    pmfile = args.pmfile
    maxchi2 = setArg(args.c, 2.5)
    bwidth = setArg(args.w, 60.0)
    bsig = setArg(args.s, 1.2)
    dataversion = setArg(args.dataversion, 'v5_7_0')
    datadir = setArg(args.datadir, '/SCRATCH')
    savedir = setArg(args.savedir, '../FullSearch')
    lya = setArg(args.lya, False)
    qso = setArg(args.qso, False)
    jpt = setArg(args.jpt, False)
    # Read and begin
    platemjd = np.loadtxt(pmfile, dtype=int)
    for each in platemjd:
        fiberid = np.arange(1, 1001, 1, dtype=int)
        args = []
        for fid in fiberid:
            args.append((each[0], each[1], fid, dataversion, datadir, savedir,
                         lya, qso, jpt, bwidth, bsig, maxchi2,))
        res = para_return(lensFinder, args, 12)
    '''
    # Uncomment below and comment above to debug
    # lensFinder(4388, 55536, 69, datadir='/SCRATCH')
    # lensFinder(4391, 55866, 179, datadir='/SCRATCH')
    lensFinder(4201, 55443, 26, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    '''
    lensFinder(4398, 55946, 515, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4398, 55946, 737, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4398, 55946, 883, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7238, 56660, 171, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6717, 56397, 501, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4321, 55504, 613, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6741, 56394, 920, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6807, 56429, 59, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6807, 56429, 553, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4397, 55921, 266, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4397, 55921, 339, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4397, 55921, 412, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4394, 55924, 49, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4394, 55924, 243, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4394, 55924, 984, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4391, 55866, 42, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4391, 55866, 235, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6739, 56393, 272, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4388, 55536, 333, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6933, 56398, 201, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6933, 56398, 259, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6933, 56398, 548, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6933, 56398, 562, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6933, 56398, 655, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6933, 56398, 708, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6933, 56398, 740, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 65, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 84, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 128, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 170, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 203, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 412, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 414, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7184, 56629, 452, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6804, 56447, 45, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4393, 55944, 364, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7032, 56471, 981, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4387, 55534, 252, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4387, 55534, 474, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4199, 55481, 876, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4199, 55481, 999, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7235, 56603, 849, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6932, 56397, 81, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6932, 56397, 322, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6932, 56397, 349, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6932, 56397, 707, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6932, 56397, 861, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(6932, 56397, 897, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(7237, 56662, 245, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    lensFinder(4349, 55803, 775, 'v5_7_0', '../SCRATCH', '../PlotCheck', False, False, False, 60.0, 1.2, 2.5, True)
    '''
