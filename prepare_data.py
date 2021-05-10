import numpy as np
import os.path as op
import cooler as cool
import math
import pickle as pick
from os.path import exists as pexists
from os.path import join as pjoin
from os import mkdir
import parameters


def prepare_data(filepath_hicmatrix, filepath_TAD_domains, write_windows = False):
    filepath_cool = filepath_hicmatrix
    filepath_bound = filepath_TAD_domains
    write_windows = write_windows
    window_size = parameters.window_size
    overlap = parameters.overlap_size
    binsize = parameters.binsize
    windows = None

    if(overlap >= window_size):
        print("overlap must be smaller than window size")

    offset = 0 #how far start point is moved because of chromosome concatenating in hicmatrix
    #make directories to save preparations
    if not pexists("preparations"):
        mkdir("preparations")
    if not pexists(pjoin("preparations", parameters.save_preparation_id)):
        mkdir(pjoin("preparations", parameters.save_preparation_id))
    if not pexists(pjoin(pjoin("preparations", parameters.save_preparation_id), "windows.bed")) or not pexists(pjoin(pjoin("preparations", parameters.save_preparation_id), "foundboundaries.bed")):
        write_windows = True


    #begin bedfile with window boundaries
    if(write_windows):
        windows = open(pjoin(pjoin("preparations", parameters.save_preparation_id), "windows.bed"), "w")
        windows.write("chrom\tchromStart\tchromEnd\n")
        found_boundaries = open(pjoin(pjoin("preparations", parameters.save_preparation_id), "foundboundaries.bed"), "w")
        found_boundaries.write("chrom\tchromStart\tchromEnd\n")


    #load hicmatrix
    c = cool.Cooler(filepath_cool)
    numberbins = [] # list of number of bins for each chromosome
    arr = c.matrix(balance=False)



    #compute nr of bins per chromosome
    for chromsize in c.chromsizes:
        numberbins.append(math.ceil(chromsize / binsize))



    #dictionaries for submats/boundaries (key:chromosome, value:list of submats/boundaries)
    submats = {}
    boundaries = {}
    #extract boundaries from TAD_domains.bed
    with open(filepath_bound, 'r') as f:
        for line in f: #get boundaries bins
            cont = line.strip().split()
            if cont[0] not in boundaries:
                boundaries[cont[0]] = [int(cont[1]) / binsize]
            boundaries[cont[0]].append((int(cont[2]) / binsize))

    chrcount = 0
    labels = {}
    #go through every chromosome and compute submatrices
    for number in numberbins:
        chrom = c.chromnames[chrcount]

        print('Chromname: ' + chrom)
        try:#try to get the boundaries of the chromosome if chromosome has boundaries
            chrbounds = boundaries[chrom]
        except:#if no boundaries, skip chromosome
            chrcount+=1
            continue

        labels[chrom] = []
        #new start and end position depending on last positions, window_size and overlap sizes
        start = offset
        end = start + window_size
        submat = []
        startboundary = 0
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
            searchboundary = startboundary
            while(chrbounds[searchboundary] >= start - offset and chrbounds[searchboundary] <= end -offset):
            #checks if the next boundary is in the current window center and set labels
                if start - offset + int((end-start)/2) == chrbounds[startboundary]:
                    label = 1
                    if(write_windows):
                        found_boundaries.write(chrom+"\t"+str((start - offset + int((end-start)/2)-1)*binsize) + "\t" + str((start - offset + int((end-start)/2)+1)*binsize) + "\n")
                    break
                try:
                    chrbounds[searchboundary+1]
                    searchboundary +=1
                except:
                    break
            submat.append(arr[start:end, start:end])
            labels[chrom].append(label)
            if (write_windows):
                windows.write(chrom + "\t" + str((start - offset)*binsize) + "\t" + str((end - offset)*binsize) + "\n")

            start, end = start + window_size - overlap, end + window_size - overlap
        offset = offset + number #computes new offset for chromosome
        submats[chrom] = submat
        chrcount += 1
        print()


    #save interaction matrices and corresponding labels
    pick.dump(submats, open('preparations/'+ parameters.save_preparation_id +'/InteractionMatrices' , 'wb'))
    pick.dump(labels, open('preparations/' + parameters.save_preparation_id + '/labels', 'wb'))
    if(write_windows):
        windows.close()
        found_boundaries.close()

if __name__ == "__main__":
    prepare_data(parameters.hic_matrix, parameters.TAD_domains, True)