import tensorflow as tf
import numpy as np
import load_data as load
import parameters
import cooler as cool

TRAINCHORMS = parameters.TRAINCHORMS
VALCHROMS = parameters.VALCHROMS
TESTCHROMS = parameters.TESTCHROMS
#TODO: split function for each kind of set


def build_training_set (cooler,windows_file = parameters.windows_bed, chroms = parameters.TRAINCHORMS):
    x = []
    y = []
    matrix = cooler.matrix(balance=False)
    with open(windows_file, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count % 1000 == 0:
                print(f"{count} Lines of parsed")
            cont = line.strip().split()
            if cont[0] not in chroms:
                continue
            offset = cooler.offset(cont[0])
            binsize = cooler.binsize
            x.append(matrix[offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize),
                     offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize)])
            y.append(int(cont[4])//1000)

    return x,y







# split data in training, validation and testset, still hardcoded split of chromosomes
def build_tvt(verbose=False):
    # load data
    interaction_matrices_pos, interaction_matrices_neg = load.load_data(parameters.interaction_matrices_pos, parameters.interaction_matrices_neg)
    # construct training set
    xtrain = []
    ytrain = []
    for chr in TRAINCHORMS:
        xtrain = xtrain + interaction_matrices_pos[chr]+interaction_matrices_neg[chr]
        ytrain = ytrain + len(interaction_matrices_pos[chr]) * [1] + len(interaction_matrices_neg[chr]) * [0]
    xtrain = np.array(xtrain)
    ytrain = np.array(ytrain)
    if verbose:
        zeros = ytrain[ytrain == 0]
        print("Percent zeros train: " + str(len(zeros) / len(ytrain)))

    # construct validation set
    xval = []
    yval = []
    for chr in VALCHROMS:
        xval = xval + interaction_matrices_pos[chr] + interaction_matrices_neg[chr]
        yval = yval + len(interaction_matrices_pos[chr]) * [1] + len(interaction_matrices_neg[chr]) * [0]
    xval = np.array(xval)
    yval = np.array(yval)

    if verbose:
        zeros = yval[yval == 0]
        print("Percent zeros validation: " + str(len(zeros) / len(yval)))

    # construct test set
    xtest = []
    ytest = []
    for chr in TESTCHROMS:
        xtest = xtest + interaction_matrices_pos[chr] + interaction_matrices_neg[chr]
        ytest = ytest + len(interaction_matrices_pos[chr]) * [1] + len(interaction_matrices_neg[chr]) * [0]
    xtest = np.array(xtest)
    ytest = np.array(ytest)

    if(verbose):
        print(xtrain[0][9])

    if verbose:
        zeros = ytest[ytest == 0]
        print("Percent zeros validation: " + str(len(zeros) / len(ytest)))
    return (xtrain, ytrain), (xval, yval), (xtest, ytest)