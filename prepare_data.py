import numpy as np
import os.path as op
import cooler as cool
import math
import pickle as pick
from os.path import exists as pexists
from os.path import join as pjoin
from os import mkdir


def prepare_data(filepath_hicmatrix, filepath_TAD_domains,  binsize, window_size, overlap_size, heatmap=True, write_windows = False):
    filepath_cool = filepath_hicmatrix
    filepath_bound = filepath_TAD_domains
    write_windows = write_windows
    step = window_size
    overlap = overlap_size
    binsize = binsize
    windows = None
    save_preparation_id = str(window_size)+'_'+ str(overlap_size)+'_'+ str(binsize)
    offset = 0 #how far start point is moved because of chromosome concatenating in hicmatrix
    #make directories to save preparations
    if not pexists("preparations"):
        mkdir("preparations")
    if not pexists(pjoin("preparations", save_preparation_id)):
        mkdir(pjoin("preparations", save_preparation_id))
    if not pexists(pjoin(pjoin("preparations", save_preparation_id), "heatmaps")):
        mkdir(pjoin(pjoin("preparations", save_preparation_id), "heatmaps"))
    if not pexists(pjoin(pjoin("preparations", save_preparation_id), "windows.bed")):
        write_windows = True

    #begin bedfile with window boundaries
    if(write_windows):
        windows = open(pjoin(pjoin("preparations", save_preparation_id), "windows.bed"), "w")
        windows.write("chrom\tchromStart\tchromEnd\n")


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
                boundaries[cont[0]] = []
            boundaries[cont[0]].append((int(cont[1]) / binsize, int(cont[2]) / binsize))

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
        #new start and end position depending on last positions, step and overlap sizes
        start = offset
        end = start + step
        submat = []
        bcount = 0
        #go trough hicmatrix and build submatrix with corresponding boundaries flag
        while end <= number + offset:

            #while next boundary is before the start of the next window, got to the next boundary
            while chrbounds[bcount][1] < start - offset:
                #tries if there is another boundary, if not break
                try:
                    chrbounds[bcount+1][0]
                    bcount += 1

                except:
                    break

            #checks if the next boundary is in the current window and set labels
            if start - offset <= chrbounds[bcount][0] < end - offset \
                    or start - offset <= chrbounds[bcount][1] < end - offset:
                labels[chrom].append(1)
            else:
                labels[chrom].append(0)
            submat.append(arr[start:end, start:end])
            if (write_windows):
                windows.write(chrom + "\t" + str((start - offset)*binsize) + "\t" + str((end - offset)*binsize) + "\n")
            start, end = start + step - overlap, end + step - overlap
        offset = offset + number #computes new offset for chromosome
        submats[chrom] = submat
        chrcount += 1
        print()


    #save interaction matrices and corresponding labels
    pick.dump(submats, open('preparations/'+ save_preparation_id +'/InteractionMatrices' , 'wb'))
    pick.dump(labels, open('preparations/' + save_preparation_id + '/labels', 'wb'))
    if(write_windows):
        windows.close()

if __name__ == "__main__":
    prepare_data("data/GM12878-MboI-allreps-filtered.10kb.cool", "data/GM12878-MboI-allreps-filtered-TAD-domains.bed",
                 10000, 25, 5)