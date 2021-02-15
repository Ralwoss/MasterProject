import numpy as np
import pickle as pick

a = pick.load(open('InteractionMatrices', 'rb'))
print(a)
b = pick.load(open('labels', 'rb'))
print(b)

