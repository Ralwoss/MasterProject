import numpy
import pickle as pick
import matplotlib
import matplotlib.pyplot as plt
import os
from os.path import join as pjoin

import numpy as np


def plotHeatmap(dir, chromosome="all"):
    matrices = {}
    if chromosome == "all":
        for file in os.listdir(dir):
            matrix = np.array(pick.load(open(pjoin(dir, file),'rb')))
            print(matrix)
            plt.imshow(matrix,cmap="RdYlBu_r")
            plt.show()
    else:
        try:
            matrix = np.array(pick.load(open(pjoin(dir, chromosome), 'rb')))
            print(matrix)
            plt.imshow(matrix, cmap="RdYlBu_r")
            plt.show()
        except:
            print("Chromosome not found")




if __name__ == "__main__":
    plotHeatmap(pjoin(pjoin("preparations", "10_0_10000"), "heatmaps"), chromosome="chr20")