import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal

from NoiseGeneration.RealNoise import ProcessRawData


def resampeld(TargetSampleRate):
    print("Reading Start")
    try:
        NoiseArray = pd.read_csv("../Data/RealNoiseSequence.csv", delimiter=",", engine="c").to_numpy().reshape(-1)[
                     :12000001]  # 20 minutes
    except FileNotFoundError:
        ProcessRawData.main()
        NoiseArray = pd.read_csv("../Data/RealNoiseSequence.csv", delimiter=",", engine="c").to_numpy().reshape(-1)[
                     :12000001]  # 20 minutes
    print("Reading Complete")
    """Real Sample Rate"""
    dt = 0.0001  # 10 kHz sampling rate
    fs = 1 / dt
    maxt = len(NoiseArray) / fs
    TargetSampleSize = TargetSampleRate * maxt
    print(int(np.floor(TargetSampleSize) + 1))
    resampledNoise = signal.resample(NoiseArray, int(np.ceil(TargetSampleSize)))
    RealTime = np.arange(0, maxt, dt)
    TargetTime = np.arange(0, maxt, 1 / TargetSampleRate)
    plt.plot(RealTime, NoiseArray, linewidth=3, alpha=0.5)
    plt.plot(TargetTime, resampledNoise)
    plt.show()
    return resampledNoise
