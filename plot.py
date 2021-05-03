import numpy
import pickle as pick
import matplotlib
import matplotlib.pyplot as plt
import os
from pygenometracks import plotTracks

def main (args = None):
    print(args)
    plotTracks.main(args)

if __name__ == "__main__":
    main()
