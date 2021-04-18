import numpy as np
import pickle as pick

def load_data(interaction_matrices, labels, verbose = False):
    a = pick.load(open(interaction_matrices, 'rb'))
    b = pick.load(open(labels, 'rb'))
    return a,b
