import tensorflow as tf
from tensorflow.keras.utils import Sequence
import numpy as np
import parameters as pars


class dataGenerator(Sequence):
    def __init__(self, pos_matrices_npz, neg_matrices_npz, batch_size=32, dim=(pars.window_size, pars.window_size),
                 n_classes=2, shuffle=True, balance_method=None):

        self.dim = dim
        self.batch_size = batch_size
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.balance_method = balance_method

        pos_matrices_npz = np.load(pos_matrices_npz)
        neg_matrices_npz = np.load(neg_matrices_npz)

        self.pos_train_matrices = pos_matrices_npz["train"]
        self.neg_train_matrices = neg_matrices_npz["train"]

        self.on_epoch_end()

    def on_epoch_end(self):
        if str(self.balance_method).lower() == "oversampling":
            # for smaller sample set multiply size to get more indexes for sampling
            maxsize = max(len(self.pos_train_matrices), len(self.neg_train_matrices))
            if len(self.pos_train_matrices) < maxsize:  # if less positive samples
                proportion = int(np.ceil(maxsize / len(self.pos_train_matrices)))
                self.indexes_pos = np.arange(len(self.pos_train_matrices) * proportion)
                self.indexes_neg = np.arange(len(self.neg_train_matrices))
            elif len(self.neg_train_matrices) < maxsize:  # if less negative samples
                proportion = int(np.ceil(maxsize / len(self.neg_train_matrices)))
                self.indexes_pos = np.arange(len(self.pos_train_matrices))
                self.indexes_neg = np.arange(len(self.neg_train_matrices) * proportion)
            else:  # if equal positive and negative samples
                self.indexes_pos = np.arange(len(self.pos_train_matrices))
                self.indexes_neg = np.arange(len(self.neg_train_matrices))

        elif str(self.balance_method).lower() == "undersampling":
            self.indexes_pos = np.arange(len(self.pos_train_matrices))
            self.indexes_neg = np.arange(len(self.neg_train_matrices))
        else:
            self.indexes = np.arange(len(self.pos_train_matrices) + len(self.neg_train_matrices))

        # print(self.indexes)
        if self.shuffle:
            if str(self.balance_method).lower() == "oversampling" \
                    or str(self.balance_method).lower() == "undersampling":
                np.random.shuffle(self.indexes_pos)
                np.random.shuffle(self.indexes_neg)
            else:
                np.random.shuffle(self.indexes)

    def __getitem__(self, index):
        # generate indexes for this batch
        # division by two to get half the number of indexes for pos and neg index lists
        X, y = None, None
        if str(self.balance_method).lower() == "oversampling":
            indexes_pos = self.indexes_pos[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            indexes_neg = self.indexes_neg[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            X, y = self.__data_generation_posneg(indexes_pos, indexes_neg)
        elif str(self.balance_method).lower() == "undersampling":
            indexes_pos = self.indexes_pos[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            indexes_neg = self.indexes_neg[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            X, y = self.__data_generation_posneg(indexes_pos, indexes_neg)
        else:  # no data balancing
            indexes = self.indexes[(index * self.batch_size): (index + 1) * self.batch_size]
            X, y = self.__data_generation(indexes)
        return X, y

    def __len__(self):
        if str(self.balance_method).lower() == "oversampling":
            return int(np.floor(2 * max(len(self.pos_train_matrices), len(self.neg_train_matrices)) / self.batch_size))
        elif str(self.balance_method).lower() == "undersampling":
            return int(np.floor(2 * min(len(self.pos_train_matrices), len(self.neg_train_matrices)) / self.batch_size))
        else:
            return int(np.floor((len(self.pos_train_matrices) + len(self.neg_train_matrices)) / self.batch_size))

    def __data_generation_posneg(self, indexes_pos, indexes_neg):
        'Generates data containing batch_size samples'  # X : (n_samples, *dim, n_channels)
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty(self.batch_size, dtype=int)

        if str(self.balance_method).lower() == "oversampling":
            for i in np.arange(len(indexes_pos)):
                X[2 * i,] = self.pos_train_matrices[indexes_pos[i] % len(self.pos_train_matrices)]
                y[2 * i] = 1
                X[2 * i + 1,] = self.neg_train_matrices[indexes_neg[i] % len(self.neg_train_matrices)]
                y[2 * i + 1] = 0
        elif str(self.balance_method).lower() == "undersampling":
            for i in np.arange(len(indexes_pos)):
                X[2 * i,] = self.pos_train_matrices[indexes_pos[i]]
                y[2 * i] = 1
                X[2 * i + 1,] = self.neg_train_matrices[indexes_neg[i]]
                y[2 * i + 1] = 0
        return X, y

    def __data_generation(self, indexes):
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty(self.batch_size, dtype=int)

        for i in np.arange(len(indexes)):
            if indexes[i] < len(self.pos_train_matrices):
                X[i, ] = self.pos_train_matrices[indexes[i]]
                y[i] = 1
            else:
                X[i, ] = self.neg_train_matrices[indexes[i] - len(self.pos_train_matrices)]
                y[i] = 0
        return X, y
        

"""
# Try to load submatrices when needed in dataGenerator. Takes too long
class dataGenerator(Sequence):
    def __init__(self, cooler, windows, chromosomes, batch_size=32, dim=(pars.window_size, pars.window_size),
                 n_classes=2, shuffle=True, balance_method=None):

        self.cooler = cooler
        self.matrix = cooler.matrix(balance = False)
        self.windows = windows
        self.chromosomes = chromosomes
        self.dim = dim
        self.batch_size = batch_size
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.balance_method = balance_method

        self.X_pos, self.X_neg = self.prepare_Xs(cooler, windows, chromosomes)

        self.on_epoch_end()

    def on_epoch_end(self):
        if str(self.balance_method).lower() == "oversampling":
            # for smaller sample set multiply size to get more indexes for sampling
            maxsize = max(len(self.X_pos), len(self.X_neg))
            if len(self.X_pos) < maxsize:  # if less positive samples
                proportion = int(np.ceil(maxsize / len(self.X_pos)))
                self.indexes_pos = np.arange(len(self.X_pos) * proportion)
                self.indexes_neg = np.arange(len(self.X_neg))
            elif len(self.X_neg) < maxsize:  # if less negative samples
                proportion = int(np.ceil(maxsize / len(self.X_neg)))
                self.indexes_pos = np.arange(len(self.X_pos))
                self.indexes_neg = np.arange(len(self.X_neg) * proportion)
            else:  # if equal positive and negative samples
                self.indexes_pos = np.arange(len(self.X_pos))
                self.indexes_neg = np.arange(len(self.X_neg))

        elif str(self.balance_method).lower() == "undersampling":
            self.indexes_pos = np.arange(len(self.X_pos))
            self.indexes_neg = np.arange(len(self.X_neg))
        else:
            self.indexes = np.arange(len(self.X_pos) + len(self.X_neg))

        # print(self.indexes)
        if self.shuffle:
            if str(self.balance_method).lower() == "oversampling" \
                    or str(self.balance_method).lower() == "undersampling":
                np.random.shuffle(self.indexes_pos)
                np.random.shuffle(self.indexes_neg)
            else:
                np.random.shuffle(self.indexes)



    def __getitem__(self, index):
        # generate indexes for this batch
        # division by two to get half the number of indexes for pos and neg index lists
        X, y = None, None
        if str(self.balance_method).lower() == "oversampling":
            indexes_pos = self.indexes_pos[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            indexes_neg = self.indexes_neg[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            X, y = self.__data_generation_posneg(indexes_pos, indexes_neg)
        elif str(self.balance_method).lower() == "undersampling":
            indexes_pos = self.indexes_pos[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            indexes_neg = self.indexes_neg[int((index * self.batch_size) / 2):int(((index + 1) * self.batch_size) / 2)]
            X, y = self.__data_generation_posneg(indexes_pos, indexes_neg)
        else:  # no data balancing
            indexes = self.indexes[(index * self.batch_size): (index + 1) * self.batch_size]
            X, y = self.__data_generation(indexes)
        return X, y

    def __len__(self):
        if str(self.balance_method).lower() == "oversampling":
            return int(np.floor(2 * max(len(self.X_pos), len(self.X_neg)) / self.batch_size))
        elif str(self.balance_method).lower() == "undersampling":
            return int(np.floor(2 * min(len(self.X_pos), len(self.X_neg)) / self.batch_size))
        else:
            return int(np.floor((len(self.X_pos) + len(self.X_neg)) / self.batch_size))

    def __data_generation_posneg(self, indexes_pos, indexes_neg):
        'Generates data containing batch_size samples'  # X : (n_samples, *dim, n_channels)
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty(self.batch_size, dtype=int)

        if str(self.balance_method).lower() == "oversampling":
            for i in np.arange(len(indexes_pos)):
                bounds = self.X_pos[indexes_pos[i] % len(self.X_pos)]
                X[2 * i,] = self.matrix[bounds[0]:bounds[1], bounds[0]:bounds[1]]
                y[2 * i] = 1
                bounds = self.X_neg[indexes_neg[i] % len(self.X_neg)]
                X[2 * i + 1,] = self.matrix[bounds[0]:bounds[1], bounds[0]:bounds[1]]
                y[2 * i + 1] = 0
        elif str(self.balance_method).lower() == "undersampling":
            for i in np.arange(len(indexes_pos)):
                X[2 * i,] = self.pos_train_matrices[indexes_pos[i]]
                y[2 * i] = 1
                X[2 * i + 1,] = self.neg_train_matrices[indexes_neg[i]]
                y[2 * i + 1] = 0
        return X, y

    def __data_generation(self, indexes):
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty(self.batch_size, dtype=int)

        for i in np.arange(len(indexes)):
            if indexes[i] < len(self.pos_train_matrices):
                X[i, ] = self.pos_train_matrices[indexes[i]]
                y[i] = 1
            else:
                X[i, ] = self.neg_train_matrices[indexes[i] - len(self.pos_train_matrices)]
                y[i] = 0
        return X, y

    def prepare_Xs(self, cooler, windows, chromosomes):
        matrix = cooler.matrix(balance=False)
        X_pos = []
        X_neg = []
        with open(windows, 'r') as f:
            count = 0
            for line in f:
                count += 1
                if count % 1000 == 0:
                    print(f"{count} Lines parsed")
                cont = line.strip().split()
                if cont[0] not in chromosomes:
                    continue
                offset = cooler.offset(cont[0])
                binsize = cooler.binsize
                if int(cont[4]):
                    X_pos.append((offset + (int(cont[1]) // binsize), offset + (int(cont[2]) // binsize) + 1))
                else:
                    X_neg.append((offset + (int(cont[1]) // binsize), offset + (int(cont[2]) // binsize) + 1))
        return X_pos, X_neg
"""