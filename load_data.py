import numpy as np
import pickle as pick

def load_data(pos_matrices, neg_matrices, verbose = False):
    pos_matrices = pick.load(open(pos_matrices, 'rb'))
    neg_matrices = pick.load(open(neg_matrices, 'rb'))


    return pos_matrices,neg_matrices
