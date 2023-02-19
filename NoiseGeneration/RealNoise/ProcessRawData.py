import os
from tkinter import *
from tkinter import filedialog

import numpy as np
import pandas as pd


def main():
    """ Initialize Files"""
    root = Tk()
    dirname = filedialog.askdirectory(parent=root, initialdir="/", title='Pick a directory')
    root.destroy()
    dt = 0.0001  # 10 kHz sampling rate
    NoiseSequence = []
    for filename in sorted(os.listdir(dirname)):
        print(filename)
        df = pd.read_csv(dirname + "/" + filename, sep=",", header=None, usecols=[1])
        NoiseSequence = NoiseSequence + df.stack().tolist()
    NoiseSequence = np.array(NoiseSequence)
    np.savetxt("../Data/RealNoiseSequence.csv", NoiseSequence, delimiter=",")
    print("File Created")
    return NoiseSequence
