import tensorflow as tf
import numpy as np
import load_data as load
import parameters

TRAINCHORMS = parameters.TRAINCHORMS
VALCHROMS = parameters.VALCHROMS
TESTCHROMS = parameters.TESTCHROMS

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

    if(verbose):
        print(xtrain[0][9])

    if verbose:
        zeros = ytest[ytest == 0]
        print("Percent zeros validation: " + str(len(zeros) / len(ytest)))
    return (xtrain, ytrain), (xval, yval), (xtest, ytest)