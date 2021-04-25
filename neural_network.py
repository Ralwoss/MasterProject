# -*- coding: utf-8 -*-
import load_data as load
import build_tvt_sets
import tensorflow as tf
import numpy as np


# set split of chromosomes
TRAINCHORMS = ['chr5', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
VALCHROMS = ['chr1', 'chr3', 'chr7']
TESTCHROMS = ['chr2', 'chr4', 'chr6']
verbose = False

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
tvt = build_tvt_sets.build_tvt()

xtrain, ytrain = tvt[0]
xval, yval = tvt[1]
xtest, ytest = tvt[2]



#compute class weights
zerocount, onecount = np.bincount(ytrain)
weightzero = len(ytrain) / (2 * zerocount)
weightone = len(ytrain) / (2*onecount)

weights = {0:weightzero, 1:weightone}
if verbose: print(weights)

initbias=tf.keras.initializers.Constant(np.log([onecount/zerocount]))

#build a model
model = tf.keras.Sequential([
                             tf.keras.layers.Flatten(input_shape=(5,5)),
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
model.save("model/5_0_10kb_batch_normalized.h5")