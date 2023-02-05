import levy
import numpy as np
import pandas as pd
from tkinter import *
from tkinter import filedialog
import matplotlib.pyplot as plt

import os
root = Tk()
dirname = filedialog.askdirectory(parent=root,initialdir="/",title='Pick a directory')
root.destroy()
dt = 0.0001  # 10 kHz sampling rate
NoiseSequence = []
for filename in os.listdir(dirname):
    print(filename)
    df = pd.read_csv(dirname+"/"+filename, sep=",", header=None,usecols=[1])
    test = df.to_numpy()
    AlphaStablePara = levy.fit_levy(test)
    print(AlphaStablePara)
    print("pause")
NoiseSequence = np.array(NoiseSequence)
NoiseSequence.tofile("NoiseSequence.csv", sep=",")

AlphaStablePara = levy.fit_levy(NoiseSequence)
print(AlphaStablePara)
