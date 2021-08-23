import numpy as np
import parameters

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
                print(f"{count} Lines parsed")
            cont = line.strip().split()
            if cont[0] not in chroms:
                continue
            offset = cooler.offset(cont[0])
            binsize = cooler.binsize
            x.append(matrix[offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize)+1,
                     offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize)+1])
            y.append(int(cont[4])//1000)
    x, y = np.array(x), np.array(y)

    return x,y

def build_validation_set (cooler,windows_file = parameters.windows_bed, chroms = parameters.VALCHROMS):
    x = []
    y = []
    matrix = cooler.matrix(balance=False)
    with open(windows_file, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count % 1000 == 0:
                print(f"{count} Lines parsed")
            cont = line.strip().split()
            if cont[0] not in chroms:
                continue
            offset = cooler.offset(cont[0])
            binsize = cooler.binsize
            x.append(matrix[offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize)+1,
                     offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize)+1])
            y.append(int(cont[4])//1000)
    x, y = np.array(x), np.array(y)
    return x,y

def build_test_set (cooler,windows_file = parameters.windows_bed, chroms = parameters.TESTCHROMS):
    x = []
    y = []
    matrix = cooler.matrix(balance=False)
    with open(windows_file, 'r') as f:
        count = 0
        for line in f:
            count += 1
            if count % 1000 == 0:
                print(f"{count} Lines parsed")
            cont = line.strip().split()
            if cont[0] not in chroms:
                continue
            offset = cooler.offset(cont[0])
            binsize = cooler.binsize
            x.append(matrix[offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize)+1,
                     offset+(int(cont[1])//binsize):offset+(int(cont[2])//binsize)+1])
            y.append(int(cont[4])//1000)
    x, y = np.array(x), np.array(y)

    return x,y







# split data in training, validation and testset, still hardcoded split of chromosomes
def build_tvt(cooler,windows_file = parameters.windows_bed, trainchroms = parameters.TRAINCHORMS,
              valchroms = parameters.VALCHROMS, testchroms = parameters.TESTCHROMS, verbose=False):
    # construct training set
    xtrain, ytrain = build_training_set(cooler,windows_file, trainchroms)
    if verbose:
        zeros = ytrain[ytrain == 0]
        print("Percent zeros train: " + str(len(zeros) / len(ytrain)))

    # construct validation set
    xval, yval = build_validation_set(cooler,windows_file,valchroms)

    if verbose:
        zeros = yval[yval == 0]
        print("Percent zeros validation: " + str(len(zeros) / len(yval)))

    # construct test set
    xtest, ytest = build_test_set(cooler,windows_file,testchroms)

    if verbose:
        zeros = ytest[ytest == 0]
        print("Percent zeros validation: " + str(len(zeros) / len(ytest)))
    return (xtrain, ytrain), (xval, yval), (xtest, ytest)