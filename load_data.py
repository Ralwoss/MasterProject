import numpy as np
import parameters as pars
import pickle as pick

def load_data(pos_matrices_npz, neg_matrices_npz, verbose = False):
    pos_matrices_npz = np.load(pos_matrices_npz,allow_pickle=True)
    neg_matrices_npz = np.load(neg_matrices_npz, allow_pickle=True)


    return pos_matrices_npz,neg_matrices_npz


if __name__ == "__main__":
    load_data('preparations/' + pars.save_preparation_id + '/InteractionMatricesPosNp.npz',
              'preparations/' + pars.save_preparation_id + '/InteractionMatricesNegNp.npz'
              )