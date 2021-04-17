import numpy as np
import pickle as pick

def load_data(interaction_matrices, labels, verbose = False):
    a = pick.load(open(interaction_matrices, 'rb'))
    b = pick.load(open(labels, 'rb'))
    return a,b
#
# data, labels = load_data()
# x = 0
# for key in data.keys():
#     x += len(data[key])
#     print(key + ": " + str(len(data[key])))
#
# print("sum: " + str(x))
# print("train: " + str(x*0.6)) # rest
# print("validation: " + str(x*0.2)) # 1, 3, 7
# print("test: " + str(x*0.2)) # 2, 4, 6