import sys
import numpy as np
import cooler as cool
import math
from os.path import exists as pexists
from os.path import join as pjoin
from os import mkdir
import pyBigWig

import parameters as pars


def make_train_data(hic_cooler, boundaries, windows_bed=None, found_boundaries_bed=None):
    """

    :param hic_cooler: Cooler file with HIC interaction matrix
    :param boundaries: dictionary with lists of TAD-boundaries for every chromosome
    :param windows_bed: file in which the windows will be written
    :param found_boundaries_bed: windows in which found boundaries will be written

    """

    binsize = hic_cooler.binsize

    numberbins = []  # list of number of bins for each chromosome

    for chromsize in hic_cooler.chromsizes:
        numberbins.append(math.ceil(chromsize / binsize))

    windows = windows_bed
    found_boundaries = found_boundaries_bed

    for chr in hic_cooler.chromnames:
        if (chr not in pars.TRAINCHORMS):
            continue
        print('Chromname: ' + chr)
        if chr in boundaries.keys():
            chrbounds = boundaries[chr]
        else:
            continue

        offset = hic_cooler.offset(chr)

        count = 0
        for bound in chrbounds:
            if (bound < pars.window_size / 2):
                continue
            wstart = int(bound - int(pars.window_size / 2)) + offset
            wend = int(bound + math.ceil(pars.window_size / 2)) + offset
            if windows_bed:
                windows.write(f"{chr}\t{(wstart-offset)*binsize}\t{(wend-offset-1)*binsize}\t{chr}_POS{count}\t1000\n")
            if found_boundaries_bed:
                found_boundaries.write(f"{chr}\t{int((bound-1)*binsize)}\t{int((bound+1)*binsize)}\n")

            count += 1

        # define first window start and end

        wstart = offset
        wend = wstart + pars.window_size
        boundarystart = 0

        count = 0

        # go trough hicmatrix and build submatrix with corresponding boundaries flag
        while wend <= hic_cooler.extent(chr)[1]:
            # while next boundary is before the start of the next window, got to the next boundary
            while chrbounds[boundarystart] < wstart - offset:
                if chrbounds[boundarystart + 1] >= len(chrbounds):
                    break
                else:
                    boundarystart += 1

            label = 0
            # checks if boundary in center area detected
            for bound in boundaries[chr][boundarystart:]:
                # break if next boundary after end of the window
                if bound > wend - offset:
                    break

                #TODO: export in additional function? Make tests?
                centersize = 2 * pars.detection_range + 1  # define center in which boundary is detected
                if wstart - offset + int((pars.window_size - centersize) / 2) <= bound < wend - offset - round(
                        (pars.window_size - centersize) / 2):
                    label = 1
                if label == 1:
                    break  # end for-loop if already found a boundary inside center

            if label != 1:
                #submat_neg.append(hic_matrix[wstart:wend, wstart:wend])
                if windows_bed:
                    windows.write(f"{chr}\t{(wstart - offset) * binsize}\t{(wend - offset - 1) * binsize}\t"
                                  f"{chr}_NEG{count}\t0\n")

            wstart, wend = wstart + pars.window_size - pars.overlap_size, wend + pars.window_size - pars.overlap_size
            count += 1




def make_val_and_test_data(hic_cooler, boundaries, windows_bed = None, found_boundaries_bed = None):
    """

    :param hic_cooler: Cooler file with HIC interaction matrix
    :param boundaries: dictionary with lists of TAD-boundaries for every chromosome
    :param windows_bed: file in which the windows will be written
    :param found_boundaries_bed: windows in which found boundaries will be written
    """

    numberbins = []  # list of number of bins for each chromosome
    binsize = hic_cooler.binsize


    for chromsize in hic_cooler.chromsizes:
        numberbins.append(math.ceil(chromsize / binsize))

    windows = windows_bed
    found_boundaries = found_boundaries_bed

    for chr in hic_cooler.chromnames:
        if chr in pars.TRAINCHORMS:
            continue

        print('Chromname: ' + chr)
        if chr in boundaries.keys():  # try to get the boundaries of the chromosome if chromosome has boundaries
            chrbounds = boundaries[chr]
        else:  # if no boundaries, skip chromosome
            continue

        offset = hic_cooler.offset(chr)

        wstart = offset
        wend = wstart + pars.window_size
        boundarystart = 0

        poscount = 0
        negcount = 0

        # go trough hicmatrix and build submatrix with corresponding boundaries flag
        while wend <= hic_cooler.extent(chr)[1]:
            # while next boundary is before the start of the next window, got to the next boundary
            while chrbounds[boundarystart] < wstart - offset:
                if chrbounds[boundarystart + 1] >= len(chrbounds):
                    break
                else:
                    boundarystart += 1
            label = 0

            # checks if boundary in center area detected
            for bound in boundaries[chr][boundarystart:]:
                if bound > wend - offset:
                    break
                centersize = 2 * pars.detection_range + 1  # computed center in which boundary is recognized
                if wstart - offset + int((pars.window_size - centersize) / 2) <= bound < wend - offset - round(
                        (pars.window_size - centersize) / 2):
                    label = 1
                if label == 1:
                    #submat_pos.append(hic_matrix[wstart:wend, wstart:wend])
                    if windows_bed:
                        windows.write(f"{chr}\t{(wstart-offset) * pars.binsize}\t{(wend - offset - 1) * pars.binsize}\t"
                                      f"{chr}_POS{poscount}\t1000\n")
                    if found_boundaries_bed:
                        found_boundaries.write(f"{chr}\t{int((bound - 1) * pars.binsize)}\t"
                                               f"{int((bound + 1) * pars.binsize)}\n")
                    poscount += 1
            if (label != 1):
                #submat_neg.append(hic_matrix[wstart:wend, wstart:wend])
                if windows_bed:
                    windows.write(f"{chr}\t{(wstart - offset) * pars.binsize}\t{(wend - offset - 1) * pars.binsize}\t"
                                  f"{chr}_NEG{negcount}\t0\n")
                negcount += 1
            wstart, wend = wstart + pars.window_size - pars.overlap_size, wend + pars.window_size - pars.overlap_size


def save_as_npz(submats_pos, submats_neg, hic_cooler):
    submats_pos_train = []
    submats_pos_val = []
    submats_pos_test = []

    submats_neg_train = []
    submats_neg_val = []
    submats_neg_test = []

    for chr in hic_cooler.chromnames:
        if chr in pars.TRAINCHORMS:
            submats_pos_train += submats_pos[chr]
            submats_neg_train += submats_neg[chr]
        elif chr in pars.VALCHROMS:
            submats_pos_val += submats_pos[chr]
            submats_neg_val += submats_neg[chr]
        elif chr in pars.TESTCHROMS:
            submats_pos_test += submats_pos[chr]
            submats_neg_test += submats_neg[chr]

    submats_pos_train = np.array(submats_pos_train)
    submats_pos_val = np.array(submats_pos_val)
    submats_pos_test = np.array(submats_pos_test)
    submats_neg_train = np.array(submats_neg_train)
    submats_neg_val = np.array(submats_neg_val)
    submats_neg_test = np.array(submats_neg_test)

    np.savez(pars.interaction_matrices_pos_npz,
             train=submats_pos_train, val=submats_pos_val, test=submats_pos_test)
    np.savez(pars.interaction_matrices_neg_npz,
             train=submats_neg_train, val=submats_neg_val, test=submats_neg_test)

def use_A_compartment(TAD_boundaries, PCAs):
    """

    :param TAD_boundaries: a dict with key:chromosomes and value:list of chromosomes
    :param PCAs: a list with tuples of PCA.bigwig files and the chromosomes for which it has the compartments
    :return: a new list of boundaries, with boundaries inside B-compartments removed
    """

    new_bounds = {}

    for PCA in PCAs:
        bw = pyBigWig.open(PCA[0])
        for chrom in PCA[1]:
            new_bounds[chrom] = []
            del_bounds = []
            A_comp = pars.A_comp[chrom]
            sign = 1
            if A_comp == "neg":
                sign = -1
            for bound in TAD_boundaries[chrom]:
                try: boundval = bw.stats(chrom, int(bound*pars.binsize), min(int((bound+1)*pars.binsize), bw.chroms(chrom)))
                except:
                    continue
                if sign * boundval[0] >= 0:
                    new_bounds[chrom].append(bound)
                else:
                    del_bounds.append(bound)
    return new_bounds






def prepare_data(filepath_hicmatrix, filepath_TAD_domains, write_windows=False, usehicpca = False):
    """Prepares the data (submatrices, TAD-boundaries, windows) for later usage

    :param filepath_hicmatrix: String with path to the HIC-File in cooler format
    :param filepath_TAD_domains: String with path to the TAD domains in BED format
    :param write_windows: Flag that indicates if the new windows and boundaries should be saved
    :return: None
    """
    filepath_hicmatrix = filepath_hicmatrix
    filepath_bound = filepath_TAD_domains
    write_windows = write_windows
    window_size = pars.window_size
    overlap = pars.overlap_size
    # binsize = pars.binsize
    windows = None

    if (overlap >= window_size):
        sys.exit("overlap must be smaller than window size")

    # make directories to save preparations
    if not pexists("preparations"):
        mkdir("preparations")
    if not pexists(pjoin("preparations", pars.save_preparation_id)):
        mkdir(pjoin("preparations", pars.save_preparation_id))
    if not pexists(pjoin(pjoin("preparations", pars.save_preparation_id), "windows.bed")) or not pexists(
            pjoin(pjoin("preparations", pars.save_preparation_id), "foundboundaries.bed")):
        write_windows = True  # sets write_windows to true because there are no windows yet

    # begin bedfile with window boundaries
    if (write_windows):
        windows = open(pjoin(pjoin("preparations", pars.save_preparation_id), "windows.bed"), "w")
        found_boundaries = open(pjoin(pjoin("preparations", pars.save_preparation_id), "foundboundaries.bed"), "w")

    # load hicmatrix
    c = cool.Cooler(filepath_hicmatrix)
    m = c.matrix(balance = False)
    binsize = c.binsize

    boundaries = {}

    # extract boundaries from TAD_domains.bed, save in lists all boundaries per chrom
    with open(filepath_bound, 'r') as f:
        for line in f:  # get bins of boundaries
            cont = line.strip().split()
            if cont[0] not in boundaries:  # make new entry for chrom in dict if chrom not in dict
                boundaries[cont[0]] = [int(cont[1]) / binsize]  # Add first boundary
            boundaries[cont[0]].append((int(cont[2]) / binsize))  # append every end boundary of TAD

    boundaries = use_A_compartment(boundaries, [["data/PC1.bigwig", pars.PC1], ["data/PC2.bigwig", pars.PC2]])


    windows.write("chrom\tchromStart\tchromEnd\tname\tscore\n")
    found_boundaries.write("chrom\tchromStart\tchromEnd\tname\tscore\n")

    # gets all submatrices for validation and training datasets
    make_val_and_test_data(c, boundaries, windows, found_boundaries)
    make_train_data(c, boundaries, windows, found_boundaries)

    if (write_windows):
        windows.close()
        found_boundaries.close()


if __name__ == "__main__":
    prepare_data(pars.hic_matrix, pars.TAD_domains, True)
