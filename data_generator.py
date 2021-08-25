from tensorflow.keras.utils import Sequence
import numpy as np
import parameters as pars
import random

# Try to load submatrices when needed in dataGenerator. Takes too long
class dataGenerator(Sequence):
    #TODO: write/update functionality for no balance method
    def __init__(self, cooler, windows, chromosomes, batch_size=32, dim=(pars.window_size*pars.window_size,),
                 n_classes=2, shuffle=True, balance_method=None):

        self.cooler = cooler
        self.matrix = cooler.matrix(balance = False) #do not try balance = True, does not work
        self.windows = windows
        self.chromosomes = chromosomes
        self.dim = dim
        self.batch_size = batch_size
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.balance_method = balance_method

        self.X_pos_prep, self.X_neg_prep = self.prepare_Xs(cooler,windows,chromosomes)

        self.on_epoch_end()
        self.X_pos, self.X_neg = self.make_Xs()

    def on_epoch_end(self):
        if str(self.balance_method).lower() == "oversampling":
            # for smaller sample set multiply size to get more indexes for sampling
            maxsize = max(len(self.X_pos_prep), len(self.X_neg_prep))
            if len(self.X_pos_prep) < maxsize:  # if less positive samples
                proportion = int(np.ceil(maxsize / len(self.X_pos_prep)))
                self.indexes_pos = np.arange(len(self.X_pos_prep) * proportion)
                self.indexes_neg = np.arange(len(self.X_neg_prep))
            elif len(self.X_neg_prep) < maxsize:  # if less negative samples
                proportion = int(np.ceil(maxsize / len(self.X_neg_prep)))
                self.indexes_pos = np.arange(len(self.X_pos_prep))
                self.indexes_neg = np.arange(len(self.X_neg_prep) * proportion)
            else:  # if equal positive and negative samples
                self.indexes_pos = np.arange(len(self.X_pos_prep))
                self.indexes_neg = np.arange(len(self.X_neg_prep))

        elif str(self.balance_method).lower() == "undersampling":
            self.indexes_pos = np.arange(len(self.X_pos_prep))
            self.indexes_neg = np.arange(len(self.X_neg_prep))
        else:  #balance method = None
            self.indexes = np.arange(len(self.X_pos_prep) + len(self.X_neg_prep))

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
        #print(X.shape)
        return X, y

    def __len__(self):
        if str(self.balance_method).lower() == "oversampling":
            return int(np.floor(2 * max(len(self.X_pos_prep), len(self.X_neg_prep)) / self.batch_size))
        elif str(self.balance_method).lower() == "undersampling":
            return int(np.floor(2 * min(len(self.X_pos_prep), len(self.X_neg_prep)) / self.batch_size))
        else:
            return int(np.floor((len(self.X_pos) + len(self.X_neg)) / self.batch_size))

    def __data_generation_posneg(self, indexes_pos, indexes_neg):
        'Generates data containing batch_size samples'  # X : (n_samples, *dim, n_channels)
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty(self.batch_size, dtype=int)

        if str(self.balance_method).lower() == "oversampling":
            for i in np.arange(len(indexes_pos)):
                X[2 * i,] = self.X_pos[indexes_pos[i] % len(self.X_pos)]
                y[2 * i] = 1
                X[2 * i + 1] = self.X_neg[indexes_neg[i] % len(self.X_neg)]
                y[2 * i + 1] = 0
        elif str(self.balance_method).lower() == "undersampling":
            for i in np.arange(len(indexes_pos)):
                X[2 * i,] = self.X_pos[indexes_pos[i]]
                y[2 * i] = 1
                X[2 * i + 1,] = self.X_neg[indexes_neg[i]]
                y[2 * i + 1] = 0
        return X, y

    def __data_generation(self, indexes):
        #TODO write function new
        X = np.empty((self.batch_size, *self.dim))
        y = np.empty(self.batch_size, dtype=int)

        for i in np.arange(len(indexes)):
            if indexes[i] < len(self.X_pos):
                X[i, ] = self.X_pos[indexes[i]]
                y[i] = 1
            else:
                X[i, ] = self.X_neg[indexes[i] - len(self.X_pos)]
                y[i] = 0
        return X, y

    def prepare_Xs(self, cooler, windows, chromosomes):

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
        print("finished parsing all Lines")
        return X_pos, X_neg

    def make_Xs(self):
        #TODO: move random addition to prepare Xs?

        X_pos = []
        X_neg = []
        X_pos_bounds, X_neg_bounds = self.X_pos_prep, self.X_neg_prep
        if (self.balance_method != None):
            max_pos = len(self.indexes_pos)
            max_neg = len(self.indexes_neg)
            count = 0

            while count < max_pos:
                random_addition = random.randint(-pars.detection_range,pars.detection_range)
                bound = X_pos_bounds[count % len(X_pos_bounds)]
                if not (bound[0] + random_addition < 0 or bound[1] + random_addition >= len(self.matrix)):
                    bound = (bound[0]+random_addition, bound[1]+random_addition)
                X_pos.append((self.matrix[bound[0]:bound[1], bound[0]:bound[1]]).flatten())
                count += 1
                if count % 1000 == 0:
                    print(f"{count} positive submatrices build")
            print("finished building all positive submatrices")
            count = 0
            while count < max_neg:
                random_addition = random.randint(-pars.detection_range, pars.detection_range)
                bound = X_neg_bounds[count % len(X_neg_bounds)]
                if not (bound[0] + random_addition < 0 or bound[1] + random_addition >= len(self.matrix)):
                    bound = (bound[0]+random_addition, bound[1]+random_addition)
                X_neg.append((self.matrix[bound[0]:bound[1], bound[0]:bound[1]]).flatten())
                count += 1
                if count % 1000 == 0:
                    print(f"{count} negative submatrices build")
            print("finished building all negative submatrices")
        else:
            count = 0
            for bounds in X_pos_bounds:
                random_addition = random.randint(-pars.detection_range, pars.detection_range)
                if not (bounds[0] + random_addition < 0 or bounds[1] + random_addition >= len(self.matrix)):
                    bounds = (bounds[0] + random_addition, bounds[1] + random_addition)
                X_pos.append((self.matrix[bounds[0]:bounds[1], bounds[0]:bounds[1]]).flatten())
                count += 1
                if count % 1000 == 0:
                    print(f"{count} positive submatrices build")

            print("finished building all positive submatrices")
            count = 0
            for bounds in X_neg_bounds:
                random_addition = random.randint(-pars.detection_range, pars.detection_range)
                if not (bounds[0] + random_addition < 0 or bounds[1] + random_addition >= len(self.matrix)):
                    bounds = (bounds[0] + random_addition, bounds[1] + random_addition)
                X_neg.append((self.matrix[bounds[0]:bounds[1], bounds[0]:bounds[1]]).flatten())
                count += 1
                if count % 1000 == 0:
                    print(f"{count} negative submatrices build")
            print("finished building all negative submatrices")

        return np.array(X_pos), np.array(X_neg)



