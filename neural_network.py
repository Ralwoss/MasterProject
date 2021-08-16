# -*- coding: utf-8 -*-
import sys

import load_data as load
import build_tvt_sets
import tensorflow as tf
import numpy as np
import cooler as cool

import parameters
import parameters as pars
import data_generator

import load_model

verbose = True

"""
def oversample(more_matrices, less_matrices):
    ids = np.arange(len(less_matrices))
    oversampled_matrices = less_matrices[np.random.choice(ids, len(more_matrices))]
    return oversampled_matrices

def undersample(more_matrices, less_matrices):
    ids = np.arange(len(less_matrices))
    undersampled_matrices = more_matrices[np.random.choice(ids, len(less_matrices))]
    return undersampled_matrices
"""

def build_network(balanceData):
    balanceData = balanceData  #method to balance data: 0 - none; 1 - class weights; 2 - oversampling; 3 - undersampling
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
    """
    xtrain, ytrain = build_tvt_sets.build_training_set()



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
    """
    # build a model
    model = tf.keras.Sequential([
                                 tf.keras.layers.Flatten(input_shape=(pars.window_size,pars.window_size)),
                                 tf.keras.layers.BatchNormalization(),
                                 tf.keras.layers.Dense(256, activation='elu'),
                                 tf.keras.layers.Dropout(0.2),
                                 tf.keras.layers.Dense(256, activation='elu'),
                                 tf.keras.layers.Dropout(0.2),
                                 tf.keras.layers.Dense(1, activation='sigmoid')
    ])

    """model = tf.keras.models.Sequential([
        tf.keras.layers.Flatten(input_shape=(pars.window_size, pars.window_size)),
        tf.keras.layers.Dense(128, activation='relu'),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(1, activation='sigmoid', bias_initializer=initbias)
    ])"""

    #compile the model
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=METRICS)
    #fit the model
    #model.fit(x=xtrain, y=ytrain, epochs = 10, class_weight=weights, batch_size=32)

    dg = data_generator.dataGenerator(cool.Cooler(pars.hic_matrix), pars.windows_bed, pars.TRAINCHORMS, balance_method="oversampling")

    #a = dg.__getitem__(0)

    model.fit(dg, epochs=100)

    #save the model
    model.save(pars.model)

def evaluate_network(model, detailed = False, results_dir = pars.results_dir):
    # evaluate the validation set
    hic_cooler = cool.Cooler(parameters.hic_matrix)
    if not detailed:
        xval, yval = build_tvt_sets.build_validation_set(hic_cooler, parameters.windows_bed)
        model.evaluate(xval, yval, verbose=2)
    else:
        fp_path = f"{results_dir}false_positives.txt"
        fn_path = f"{results_dir}false_negatives.txt"
        windows = parameters.windows_bed
        chroms = parameters.VALCHROMS
        binsize = hic_cooler.binsize
        matrix = hic_cooler.matrix(balance=False)
        with open(windows, 'r') as f, open(fp_path, 'w') as fp_file, open(fn_path,'w') as fn_file:
            fn_file.write(f"chrom\tchromStart\tchromEmd\taccurate_result\n")
            fp_file.write(f"chrom\tchromStart\tchromEmd\taccurate_result\n")
            count = 0
            for line in f:
                count += 1
                if count % 1000 == 0:
                    print(f"{count} Lines parsed")
                cont = line.strip().split()
                if cont[0] not in chroms:
                    continue
                offset = hic_cooler.offset(cont[0])
                pos = (offset + (int(cont[1]) // binsize), offset + (int(cont[2]) // binsize) + 1)
                submatrix = np.expand_dims(matrix[pos[0]:pos[1], pos[0]:pos[1]], axis=0)
                result_acc = model(submatrix)
                result = 0 if result_acc[0,0] < 0.5 else 1
                if (result == 0 and int(cont[4]) == 1000):
                    fn_file.write(f"{cont[0]}\t{cont[1]}\t{cont[2]}\t{result_acc[0,0]}\n")
                elif (result == 1 and int(cont[4]) == 0):
                    fp_file.write(f"{cont[0]}\t{cont[1]}\t{cont[2]}\t{result_acc[0,0]}\n")




if (__name__ == "__main__"):
    build_network(pars.balanceData)

    #evaluate_network(load_model.load_model(pars.model), detailed=True)
    evaluate_network(load_model.load_model(pars.model), detailed=False)

    #resultx, resulty = build_tvt_sets.build_training_set(cool.Cooler(parameters.hic_matrix))
    #print(resulty)

