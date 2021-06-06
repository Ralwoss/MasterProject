import numpy as np
import os.path as op
import cooler as cool
import math
import pickle as pick
from os.path import exists as pexists
from os.path import join as pjoin
from os import mkdir
import parameters as pars


def make_train_data(hic_cooler, boundaries, write_windows, windows_bed, found_boundaries_bed):
    offset = 0  # how far start point is moved because of chromosome concatenating in hicmatrix
    numberbins = []  # list of number of bins for each chromosome
    hic_matrix = hic_cooler.matrix(balance=False)

    for chromsize in hic_cooler.chromsizes:
        numberbins.append(math.ceil(chromsize / pars.binsize))

    windows = None
    if (write_windows):
        windows = open(pjoin(pjoin("preparations", pars.save_preparation_id), "windows.bed"), "w")
        windows.write("chrom\tchromStart\tchromEnd\n")
        found_boundaries = open(pjoin(pjoin("preparations", pars.save_preparation_id), "foundboundaries.bed"),
                                "w")
        found_boundaries.write("chrom\tchromStart\tchromEnd\n")

    submats_pos = {}
    submats_neg = {}

    chrcount = 0

    for number in numberbins:
        chrom = hic_cooler.chromnames[chrcount]
        if(chrom not in pars.TRAINCHORMS):
            chrcount += 1
            if chrcount >= len(hic_cooler.chromnames):
                return submats_pos, submats_neg
            continue

        print('Chromname: ' + chrom)
        try:  # try to get the boundaries of the chromosome if chromosome has boundaries
            chrbounds = boundaries[chrom]
        except:  # if no boundaries, skip chromosome
            chrcount += 1
            continue

        # new start and end position depending on last positions, window_size and overlap sizes
        start = offset
        end = start + pars.window_size
        startboundary = 0

        # for saving submats of negative and positive windows
        submat_neg = []
        submat_pos = []

        # constructs window for every boundary with boundary in the center (positive windows)
        for bound in chrbounds:
            if (bound < pars.window_size / 2):
                continue
            start = int(bound - int(pars.window_size / 2)) + offset
            end = int(bound + round(pars.window_size / 2)) + offset
            submat_pos.append(hic_matrix[start:end, start:end])
            if (write_windows):
                windows.write(chrom + "\t" + str((start) * pars.binsize) + "\t" + str((end - 1) * pars.binsize) + "\n")
                found_boundaries.write(
                    chrom + "\t" + str(int(bound - 1) * pars.binsize) + "\t" + str(int(bound + 1) * pars.binsize) + "\n")

        start = offset
        end = start + pars.window_size

        # go trough hicmatrix and build submatrix with corresponding boundaries flag
        while end <= number + offset:

            # while next boundary is before the start of the next window, got to the next boundary
            while chrbounds[startboundary] < start - offset:
                # tries if there is another boundary, if not break
                try:
                    chrbounds[startboundary + 1]
                    startboundary += 1

                except:
                    break
            label = 0

            # checks if boundary in center area detected
            for bound in boundaries[chrom][startboundary:]:
                if bound > end - offset:
                    break
                centersize = 2 * pars.detection_range + 1
                if bound >= start - offset + int((pars.window_size - centersize) / 2) and bound < end - offset- round(
                        (pars.window_size - centersize) / 2):
                    label = 1
                if label == 1:
                    break
            if (label != 1):
                submat_neg.append(hic_matrix[start:end, start:end])
                if (write_windows):
                    windows.write(
                        chrom + "\t" + str((start - offset) * pars.binsize) + "\t" + str((end - offset) * pars.binsize) + "\n")

            start, end = start + pars.window_size - pars.overlap_size, end + pars.window_size - pars.overlap_size
        offset = offset + number  # computes new offset for chromosome
        submats_pos[chrom] = submat_pos
        submats_neg[chrom] = submat_neg
        chrcount += 1
        print("Number of positive windows in " + chrom + ": " + str(len(submat_pos)))
        print("Number of negative windows in " + chrom + ": " + str(len(submat_neg)))
        print()
    return submats_pos, submats_neg

def make_val_and_test_data(hic_cooler, boundaries, write_windows, windows_bed, found_boundaries_bed):
    offset = 0  # how far start point is moved because of chromosome concatenating in hicmatrix
    numberbins = []  # list of number of bins for each chromosome
    hic_matrix = hic_cooler.matrix(balance=False)

    for chromsize in hic_cooler.chromsizes:
        numberbins.append(math.ceil(chromsize / pars.binsize))

    windows = None
    if (write_windows):
        windows = open(pjoin(pjoin("preparations", pars.save_preparation_id), "windows.bed"), "w")
        windows.write("chrom\tchromStart\tchromEnd\n")
        found_boundaries = open(pjoin(pjoin("preparations", pars.save_preparation_id), "foundboundaries.bed"),
                                "w")
        found_boundaries.write("chrom\tchromStart\tchromEnd\n")

    submats_pos = {}
    submats_neg = {}

    chrcount = 0

    for number in numberbins:
        chrom = hic_cooler.chromnames[chrcount]

        if chrom in pars.TRAINCHORMS:
            chrcount += 1
            if chrcount >= len(hic_cooler.chromnames):
                return submats_pos, submats_neg
            continue

        print('Chromname: ' + chrom)
        try:  # try to get the boundaries of the chromosome if chromosome has boundaries
            chrbounds = boundaries[chrom]
        except:  # if no boundaries, skip chromosome
            chrcount += 1
            continue

        # new start and end position depending on last positions, window_size and overlap sizes
        start = offset
        end = start + pars.window_size
        startboundary = 0

        # for saving submats of negative and positive windows
        submat_neg = []
        submat_pos = []

        """"# constructs window for every boundary with boundary in the center (positive windows)
        for bound in chrbounds:
            if (bound < pars.window_size / 2):
                continue
            start = int(bound - int(pars.window_size / 2)) + offset
            end = int(bound + round(pars.window_size / 2)) + offset
            submat_pos.append(hic_matrix[start:end, start:end])
            if (write_windows):
                windows.write(chrom + "\t" + str((start) * pars.binsize) + "\t" + str((end - 1) * pars.binsize) + "\n")
                found_boundaries.write(
                    chrom + "\t" + str(int(bound - 1) * pars.binsize) + "\t" + str(int(bound + 1) * pars.binsize) + "\n")"""

        start = offset
        end = start + pars.window_size

        # go trough hicmatrix and build submatrix with corresponding boundaries flag
        while end <= number + offset:

            # while next boundary is before the start of the next window, got to the next boundary
            while chrbounds[startboundary] < start - offset:
                # tries if there is another boundary, if not break
                try:
                    chrbounds[startboundary + 1]
                    startboundary += 1

                except:
                    break
            label = 0

            # checks if boundary in center area detected
            for bound in boundaries[chrom][startboundary:]:
                if bound > end - offset:
                    break
                centersize = 2 * pars.detection_range + 1
                if bound >= start - offset + int((pars.window_size - centersize) / 2) and bound < end - offset - round(
                        (pars.window_size - centersize) / 2):
                    label = 1
                if label == 1:
                    submat_pos.append(hic_matrix[start:end, start:end])
                    if (write_windows):
                        windows.write(
                            chrom + "\t" + str((start) * pars.binsize) + "\t" + str((end - 1) * pars.binsize) + "\n")
                        found_boundaries.write(
                            chrom + "\t" + str(int(bound - 1) * pars.binsize) + "\t" + str(
                                int(bound + 1) * pars.binsize) + "\n")

            if (label != 1):
                submat_neg.append(hic_matrix[start:end, start:end])
                if (write_windows):
                    windows.write(
                        chrom + "\t" + str((start - offset) * pars.binsize) + "\t" + str((end - offset) * pars.binsize) + "\n")

            start, end = start + pars.window_size - pars.overlap_size, end + pars.window_size - pars.overlap_size
        offset = offset + number  # computes new offset for chromosome
        submats_pos[chrom] = submat_pos
        submats_neg[chrom] = submat_neg
        chrcount += 1
        print("Number of positive windows in " + chrom + ": " + str(len(submat_pos)))
        print("Number of negative windows in " + chrom + ": " + str(len(submat_neg)))
        print()
    return submats_pos, submats_neg



def prepare_data(filepath_hicmatrix, filepath_TAD_domains, write_windows = False):
    filepath_cool = filepath_hicmatrix
    filepath_bound = filepath_TAD_domains
    write_windows = write_windows
    window_size = pars.window_size
    overlap = pars.overlap_size
    binsize = pars.binsize
    windows = None

    if(overlap >= window_size):
        print("overlap must be smaller than window size")
        return

    offset = 0 #how far start point is moved because of chromosome concatenating in hicmatrix


    #make directories to save preparations
    if not pexists("preparations"):
        mkdir("preparations")
    if not pexists(pjoin("preparations", pars.save_preparation_id)):
        mkdir(pjoin("preparations", pars.save_preparation_id))
    if not pexists(pjoin(pjoin("preparations", pars.save_preparation_id), "windows.bed")) or not pexists(pjoin(pjoin("preparations", pars.save_preparation_id), "foundboundaries.bed")):
        write_windows = True


    #begin bedfile with window boundaries
    if(write_windows):
        windows = open(pjoin(pjoin("preparations", pars.save_preparation_id), "windows.bed"), "w")
        windows.write("chrom\tchromStart\tchromEnd\n")
        found_boundaries = open(pjoin(pjoin("preparations", pars.save_preparation_id), "foundboundaries.bed"), "w")
        found_boundaries.write("chrom\tchromStart\tchromEnd\n")


    #load hicmatrix
    c = cool.Cooler(filepath_cool)
    """numberbins = [] # list of number of bins for each chromosome
    hic_matrix = c.matrix(balance=False)



    #compute nr of bins per chromosome
    for chromsize in c.chromsizes:
        numberbins.append(math.ceil(chromsize / binsize))



    #dictionaries for submats/boundaries (key:chromosome, value:list of submats/boundaries)
    submats_pos = {}
    submats_neg = {}"""
    boundaries = {}


    #extract boundaries from TAD_domains.bed
    with open(filepath_bound, 'r') as f:
        for line in f: #get boundaries bins
            cont = line.strip().split()
            if cont[0] not in boundaries: #make new entry in dict if not existing
                boundaries[cont[0]] = [int(cont[1]) / binsize]
            boundaries[cont[0]].append((int(cont[2]) / binsize)) #append every end boundary of TAD

    """chrcount = 0"""

    submats_pos, submats_neg = make_val_and_test_data(c, boundaries, write_windows, windows, found_boundaries)
    a, b = make_train_data(c, boundaries, write_windows, windows, found_boundaries)
    submats_pos.update(a), submats_neg.update(b)

    #go through every chromosome and compute submatrices
    """for number in numberbins:
        chrom = c.chromnames[chrcount]

        print('Chromname: ' + chrom)
        try:#try to get the boundaries of the chromosome if chromosome has boundaries
            chrbounds = boundaries[chrom]
        except:#if no boundaries, skip chromosome
            chrcount+=1
            continue

        
        #new start and end position depending on last positions, window_size and overlap sizes
        start = offset
        end = start + window_size
        startboundary = 0

        #for saving submats of negative and positive windows
        submat_neg = []
        submat_pos = []

        # constructs window for every boundary with boundary in the center (positive windows)
        for bound in chrbounds:
            if (bound < window_size / 2):
                continue
            start = int(bound - int(window_size / 2)) + offset
            end = int(bound + round(window_size / 2)) + offset
            submat_pos.append(hic_matrix[start:end, start:end])
            if (write_windows):
                windows.write(chrom + "\t" + str((start)*binsize) + "\t" + str((end-1)*binsize) + "\n")
                found_boundaries.write(chrom + "\t" + str(int(bound-1)*binsize) + "\t" + str(int(bound+1)*binsize) + "\n")

        start = offset
        end = start + window_size


        #go trough hicmatrix and build submatrix with corresponding boundaries flag
        while end <= number + offset:

            #while next boundary is before the start of the next window, got to the next boundary
            while chrbounds[startboundary] < start - offset:
                #tries if there is another boundary, if not break
                try:
                    chrbounds[startboundary+1]
                    startboundary += 1

                except:
                    break
            label = 0



            #checks if boundary in center area detected
            for bound in boundaries[chrom][startboundary:]:
                if bound > end - offset:
                    break
                centersize = 2 * parameters.detection_range + 1
                if bound >=  start + int((window_size - centersize)/2) and bound < end - round((window_size-centersize)/2):
                    label = 1
                if label == 1:
                    break
            if(label != 1):
                submat_neg.append(hic_matrix[start:end, start:end])
                if (write_windows):
                    windows.write(chrom + "\t" + str((start - offset)*binsize) + "\t" + str((end - offset)*binsize) + "\n")

            start, end = start + window_size - overlap, end + window_size - overlap
        offset = offset + number #computes new offset for chromosome
        submats_pos[chrom] = submat_pos
        submats_neg[chrom] = submat_neg
        chrcount += 1
        print("Number of positive windows in " + chrom + ": " + str(len(submat_pos)))
        print("Number of negative windows in " + chrom + ": " + str(len(submat_neg)))
        print()"""


    #save interaction matrices and corresponding labels
    pick.dump(submats_pos, open('preparations/'+ pars.save_preparation_id +'/InteractionMatricesPos' , 'wb'))
    pick.dump(submats_neg, open('preparations/' + pars.save_preparation_id + '/InteractionMatricesNeg', 'wb'))
    #pick.dump(labels, open('preparations/' + parameters.save_preparation_id + '/labels', 'wb'))
    if(write_windows):
        windows.close()
        found_boundaries.close()

if __name__ == "__main__":
    prepare_data(pars.hic_matrix, pars.TAD_domains, True)