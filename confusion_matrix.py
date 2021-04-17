import load_model
import build_tvt_sets
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay

def build_confusion_matrix(model):
    valset = build_tvt_sets.build_tvt()[1]

    model = load_model.load_model("5_0_10kb_batch_normalized.h5")
    labels = valset[1]

    predict = model.predict(valset[0])

    cm = confusion_matrix(valset[1], np.rint(predict), normalize=None)

    disp = ConfusionMatrixDisplay(cm)

    disp.plot(cmap=plt.cm.Blues)

    plt.show()