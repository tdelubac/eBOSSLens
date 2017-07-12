import itertools as it


def galSave(doublet, obj, peak_candidates, doublet_index, savedir, em_lines):
    confirmed_lines = []
    score = 0.0
    detection = False
    if doublet:
        fileD = open(savedir + '/candidates_doublet.txt', 'a')
        z_s = peak_candidates[doublet_index].wavelength / 3727.24 - 1.0
        if z_s > obj.z + 0.05:
            detection = True
            score += peak_candidates[doublet_index].chi
            fileD.write(str(obj.radEinstein(z_s)) + " " + str(score) +
                        " " + str(z_s) + " " + str(obj.RA) + " " +
                        str(obj.DEC) + " " + str(obj.plate) + " " +
                        str(obj.mjd) + " " + str(obj.fiberid) + " " +
                        str(peak_candidates[doublet_index].wavDoublet[0]) +
                        "\n")
        fileD.close()
        if len(peak_candidates):
            fileM = open(savedir + '/candidates_DM.txt', 'a')
            # Generating all combinations of lines from above list to
            # compare with candidates
            temp = [peak for peak in peak_candidates
                    if peak.chi != peak.chiDoublet]
            compare = em_lines[1: 5]
            if z_s > obj.z + 0.05:
                for peak in temp:
                    for line in compare:
                        if abs(peak.wavelength/line - 1.0 - z_s) < 0.01:
                            detection = True
                            confirmed_lines.append(line)
                            score += peak.chiDoublet
    elif len(peak_candidates) > 1:
        compare = it.combinations(em_lines, len(peak_candidates))
        fileM = open(savedir + '/candidates_multi.txt', 'a')
        for group in compare:
            for k in range(len(peak_candidates)):
                for j in range(k + 1, len(peak_candidates)):
                    crit1 = peak_candidates[k].wavelength / group[k]
                    crit2 = crit1 - peak_candidates[j].wavelength / group[j]
                    if abs(crit2) < 0.01 and crit1 - 1.05 > obj.z:
                        detection = True
                        z_s = peak_candidates[k].wavelength / group[k] - 1.0
                        confirmed_lines.append([group,
                                                peak_candidates[k].wavelength /
                                                group[k] - 1.0])
                        score += peak_candidates[j].chi ** 2.0 + \
                            peak_candidates[k].chi ** 2.0
    if confirmed_lines != []:
        fileM.write(str(obj.radEinstein(z_s)) + " " + str(score) + " " +
                    str(z_s) + " " + str(obj.RA) + " " + str(obj.DEC) +
                    " " + str(obj.plate) + " " + str(obj.mjd) + " " +
                    str(obj.fiberid) + " " + str(confirmed_lines) + "\n")
        fileM.close()
    return detection
