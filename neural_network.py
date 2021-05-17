# -*- coding: utf-8 -*-
import load_data as load
import build_tvt_sets
import tensorflow as tf
import numpy as np
import parameters



verbose = True



def oversample(more_matrices, less_matrices):
    ids = np.arange(len(less_matrices))
    oversampled_matrices = less_matrices[np.random.choice(ids, len(more_matrices))]
    return oversampled_matrices

def undersample(more_matrices, less_matrices):
    ids = np.arange(len(less_matrices))
    undersampled_matrices = more_matrices[np.random.choice(ids, len(less_matrices))]
    return undersampled_matrices

def build_network():
    balanceData = 3 #method to balance data: 0 - none; 1 - class weights; 2 - oversampling; 3 - undersampling
    METRICS = [

          tf.keras.metrics.TruePositives(name='tp'),

        tf.keras.metrics.FalsePositives(name='fp'),

        tf.keras.metrics.TrueNegatives(name='tn'),

        tf.keras.metrics.FalseNegatives(name='fn'),

        tf.keras.metrics.BinaryAccuracy(name='accuracy'),

        tf.keras.metrics.Precision(name='precision'),

        tf.keras.metrics.Recall(name='recall'),

        tf.keras.metrics.AUC(name='auc')
    ]
    #build training, validation and test datasets
    tvt = build_tvt_sets.build_tvt(verbose=True)

    xtrain, ytrain = tvt[0]
    xval, yval = tvt[1]
    xtest, ytest = tvt[2]



    #compute class weights
    try:
        zerocount, onecount = np.bincount(ytrain)
    except:
        print("Not enough classes with/without boundaries to learn")
        return

    weights = {0: 1, 1: 1}
    initbias = 0

    if(balanceData == 1):
        weightzero = len(ytrain) / (2 * zerocount)
        weightone = len(ytrain) / (2*onecount)

        weights = {0:weightzero, 1:weightone}
        if verbose: print("Weights: " + str(weights))
    elif (balanceData == 2):
        print()

        if(zerocount > onecount):

            oversampled_positive_submatrices = oversample(xtrain[ytrain==0], xtrain[ytrain==1])

            xtrain = np.append(xtrain[ytrain==0],oversampled_positive_submatrices, axis = 0)
            ytrain = np.array(len(oversampled_positive_submatrices) * [0] + len(oversampled_positive_submatrices) * [1])



        elif(onecount > zerocount):
            oversampled_negative_submatrices = oversample((xtrain[ytrain==1], xtrain[ytrain==0]))
            xtrain = np.append(xtrain[ytrain == 1], oversampled_negative_submatrices, axis=0)
            ytrain = np.array(len(oversampled_negative_submatrices) * [1] + len(oversampled_negative_submatrices) * [0])

    elif (balanceData == 3):
        print()

        if (zerocount > onecount):
            undersampled_negative_submatrices = undersample(xtrain[ytrain == 0], xtrain[ytrain == 1])

            xtrain = np.append(xtrain[ytrain == 1], undersampled_negative_submatrices, axis=0)
            ytrain = np.array(len(undersampled_negative_submatrices) * [1] + len(undersampled_negative_submatrices) * [0])
            print(len(ytrain))

        elif (onecount > zerocount):
            undersampled_positive_submatrices = undersample((xtrain[ytrain == 1], xtrain[ytrain == 0]))
            xtrain = np.append(xtrain[ytrain == 0], undersampled_positive_submatrices, axis=0)
            ytrain = np.array(len(undersampled_positive_submatrices) * [0] + len(undersampled_positive_submatrices) * [1])

    ids = np.arange(len(xtrain))
    np.random.shuffle(ids)
    xtrain = xtrain[ids]
    ytrain = ytrain[ids]

    initbias=tf.keras.initializers.Constant(np.log([onecount/zerocount]))

    #build a model
    model = tf.keras.Sequential([
                                 tf.keras.layers.Flatten(input_shape=(parameters.window_size,parameters.window_size)),
                                 tf.keras.layers.BatchNormalization(),
                                 tf.keras.layers.Dense(256, activation='elu'),
                                 tf.keras.layers.Dropout(0.2),
                                 tf.keras.layers.Dense(256, activation='elu'),
                                 tf.keras.layers.Dropout(0.2),
                                 tf.keras.layers.Dense(1, activation='sigmoid', bias_initializer=initbias)
    ])

    #compile the model
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=METRICS)
    #fit the model
    model.fit(x=xtrain, y=ytrain, epochs = 100, class_weight=weights)

    #evaluate the validation set
    model.evaluate(xval, yval, verbose=2)

    #save the model
    model.save(parameters.model)

if (__name__ == "__main__"):
    build_network()