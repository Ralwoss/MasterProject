import numpy as np
import os.path as op
import cooler as cool
import math
import h5py
import pickle as pick

filepath_cool = 'data/GM12878-MboI-allreps-filtered.10kb.cool'
filepath_bound = 'data/GM12878-MboI-allreps-filtered-TAD-domains.bed'
step = 25
overlap = 5
binsize = 10000
offset = 0

c = cool.Cooler(filepath_cool)
numberbins = []
arr = c.matrix(balance=False)

for chromsize in c.chromsizes:
    numberbins.append(math.ceil(chromsize / binsize))

submats = {}
boundaries = {}
with open(filepath_bound, 'r') as f:
    for line in f:
        cont = line.strip().split()
        if cont[0] not in boundaries:
            boundaries[cont[0]] = []
        boundaries[cont[0]].append((int(cont[1]) / binsize, int(cont[2]) / binsize))


chrcount = 0
labels = {}
for number in numberbins:
    chrom = c.chromnames[chrcount]
    print('Chromname: ' + chrom)
    #print(c.chromsizes)
    try:
        chrbounds = boundaries[chrom]
    except:
        continue

    #print("Chrbounds: " + str(chrbounds))
    labels[chrom] = []
    start = offset
    end = start + step
    submat = []
    bcount = 0

    while end <= number + offset:
        while chrbounds[bcount][1] < start - offset:
            try:
                chrbounds[bcount+1][0]
                bcount += 1

            except:
                break
        if start - offset <= chrbounds[bcount][0] < end - offset \
                or start - offset <= chrbounds[bcount][1] < end - offset:
            labels[chrom].append(1)
        else:
            labels[chrom].append(0)
        submat.append(arr[start:end, start:end])
        start, end = start + step - overlap, end + step - overlap
    offset = offset + number
    #print(submat[0:5])
    submats[chrom] = submat
    chrcount += 1
    #print(chrbounds)
    print()
    #print(labels)

pick.dump(submats, open('InteractionMatrices', 'wb'))
pick.dump(labels, open('labels', 'wb'))


# arr = c.matrix(balance=False,sparse=True)
# shape = arr.shape
# print(shape)
# submatrices = arr[0:step, 0:step]
#
# start = step - overlap
# end = start + step
# while end <= arr.shape[0]:
#     submatrices = np.concatenate((submatrices, arr[start:end, start:end]),axis=0)
#     start, end = start + step-overlap, end + step - overlap
# print(submatrices.shape)

# print(c.matrix(balance=False)[0:1000, 0:1000])


# def main():
#    return


# if __name__ == '__main__' :
#    main()
