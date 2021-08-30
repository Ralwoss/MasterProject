# -*- coding: utf-8 -*-

import build_tvt_sets
import tensorflow as tf
import numpy as np
import cooler as cool
import matplotlib.pyplot as plt

import parameters
import parameters as pars
import data_generator

import load_model

from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay

verbose = True


def build_network(balance_method = None):
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


    # build a model
    model = tf.keras.Sequential([
                                 tf.keras.Input(shape=(pars.window_size ** 2, )),
                                 tf.keras.layers.BatchNormalization(),
                                 tf.keras.layers.Dense(256, activation='relu'),
                                 tf.keras.layers.Dropout(0.2),
                                 tf.keras.layers.Dense(256, activation='relu'),
                                 tf.keras.layers.Dropout(0.2),
                                 tf.keras.layers.Dense(1, activation='sigmoid')
    ], name=pars.string_model(balance_method)[6:-3])

    #compile the model
    model.compile(optimizer='adam',
                  loss='binary_crossentropy',
                  metrics=METRICS)

    #fit the model
    dg = data_generator.dataGenerator(cool.Cooler(pars.hic_matrix), pars.windows_bed, pars.TRAINCHORMS,
                                      balance_method=balance_method)

    model.fit(dg, epochs=100)

    #save the model
    model.save(pars.string_model(balance_method))

def evaluate_network(model, detailed = False, results_dir = pars.results_dir, build_confusion_matrix_b = False):
    # evaluate the validation set
    hic_cooler = cool.Cooler(parameters.hic_matrix)
    windows = parameters.windows_bed
    chroms = parameters.VALCHROMS
    binsize = hic_cooler.binsize
    matrix = hic_cooler.matrix(balance=False)
    if not detailed:
        xval = []
        yval = []
        with open(parameters.windows_bed, 'r') as f:
            count = 0
            for line in f:
                count +=1
                if count % 1000 == 0:
                    print(f"{count} Lines parsed")
                cont = line.strip().split()
                if cont[0] not in chroms:
                    continue
                offset = hic_cooler.offset(cont[0])
                pos = (offset + (int(cont[1]) // binsize), offset + (int(cont[2]) // binsize) + 1)
                submatrix = matrix[pos[0]:pos[1], pos[0]:pos[1]].flatten()
                xval.append(submatrix)
                yval.append(int(cont[4])/1000)
        xval = np.array(xval)
        yval = np.array(yval)
        model.evaluate(xval, yval, verbose=2)
        if build_confusion_matrix_b:
            build_confusion_matrix(model, xval, yval)
    else:
        fp_path = f"{results_dir}false_positives.txt"
        fn_path = f"{results_dir}false_negatives.txt"
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
                submatrix = np.expand_dims(matrix[pos[0]:pos[1], pos[0]:pos[1]].flatten(), axis=0)
                #print(submatrix.ndim)
                result_acc = model(submatrix)
                result = 0 if result_acc[0,0] < 0.5 else 1
                if (result == 0 and int(cont[4]) == 1000):
                    fn_file.write(f"{cont[0]}\t{cont[1]}\t{cont[2]}\t{result_acc[0,0]}\n")
                elif (result == 1 and int(cont[4]) == 0):
                    fp_file.write(f"{cont[0]}\t{cont[1]}\t{cont[2]}\t{result_acc[0,0]}\n")

def build_confusion_matrix(model, xval, yval):

    predict = model.predict(xval)

    cm = confusion_matrix(yval, np.rint(predict), normalize=None)

    disp = ConfusionMatrixDisplay(cm)

    disp.plot(cmap=plt.cm.Blues)

    plt.savefig(f"{model.name}.png")


if (__name__ == "__main__"):
    balance_method = "oversampling"
    build_network(balance_method)

    #evaluate_network(load_model.load_model(pars.string_model(balance_method)), detailed=True)
    evaluate_network(load_model.load_model(pars.string_model(balance_method)), detailed=False, build_confusion_matrix_b = True)

    #resultx, resulty = build_tvt_sets.build_training_set(cool.Cooler(parameters.hic_matrix))
    #print(resulty)

