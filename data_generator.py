import tensorflow as tf
from tensorflow.keras.utils import Sequence
import numpy as np
import parameters as pars


class dataGenerator(Sequence):
    def __init__(self, pos_matrices_npz, neg_matrices_npz, batch_size=32, dim = (pars.window_size, pars.window_size),
                 n_classes = 2, shuffle = True, balanced_data = False):

        self.dim = dim
        self.batch_size = batch_size
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.balanced_data = balanced_data


        pos_matrices_npz = np.load(pos_matrices_npz)
        neg_matrices_npz = np.load(neg_matrices_npz)

        self.pos_train_matrices = pos_matrices_npz["train"]
        self.neg_train_matrices = neg_matrices_npz["train"]

        self.on_epoch_end()


    def on_epoch_end(self):
        self.indexes = np.arange(len(self.pos_train_matrices)+len(self.neg_train_matrices))
        #print(self.indexes)
        if self.shuffle:
            np.random.shuffle(self.indexes)

    def __getitem__(self, index):
        #generate indexes for this batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        X,y = self.__data_generation(indexes)

        return X,y

    def __len__(self):
        return int(np.floor((len(self.pos_train_matrices)+len(self.neg_train_matrices)) / self.batch_size))

    def __data_generation(self, indexes):
        'Generates data containing batch_size samples'  # X : (n_samples, *dim, n_channels)
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty((self.batch_size), dtype=int)

        for i in np.arange(len(indexes)):
            if indexes[i] < len(self.pos_train_matrices):
                X[i, ] = self.pos_train_matrices[indexes[i]]
                y[i] = 1
            else :
                X[i, ] = self.neg_train_matrices[indexes[i]-len(self.pos_train_matrices)]
                y[i] = 0
        return X, y
