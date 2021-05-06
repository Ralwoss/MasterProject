import tensorflow as tf
import numpy as np
import load_data as load
import parameters

TRAINCHORMS = ['chr5', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
VALCHROMS = ['chr1', 'chr3', 'chr7']
TESTCHROMS = ['chr2', 'chr4', 'chr6']

#split data in training, validation and testset, still hardcoded split of chromosomes
def build_tvt(verbose=False):
    # load data
    data, labels = load.load_data(parameters.interaction_matrices, parameters.labels)
    # construct training set
    xtrain = []
    ytrain = []
    for chr in TRAINCHORMS:
        xtrain = xtrain + data[chr]
        ytrain = ytrain + labels[chr]
    xtrain = np.array(xtrain)
    ytrain = np.array(ytrain)
    if verbose:
        zeros = ytrain[ytrain == 0]
        print("Percent zeros train: " + str(len(zeros) / len(ytrain)))

    # construct validation set
    xval = []
    yval = []
    for chr in VALCHROMS:
        xval = xval + data[chr]
        yval = yval + labels[chr]
    xval = np.array(xval)
    yval = np.array(yval)

    if verbose:
        zeros = yval[yval == 0]
        print("Percent zeros validation: " + str(len(zeros) / len(yval)))

    # construct test set
    xtest = []
    ytest = []
    for chr in TESTCHROMS:
        xtest = xtest + data[chr]
        ytest = ytest + labels[chr]
    xtest = np.array(xtest)
    ytest = np.array(ytest)

    if verbose:
        zeros = ytest[ytest == 0]
        print("Percent zeros validation: " + str(len(zeros) / len(ytest)))
    return (xtrain, ytrain), (xval, yval), (xtest, ytest)