import matplotlib.pyplot as plt
import numpy
import numpy as np
import pandas as pd
import scipy.stats as ss
from scipy.stats import levy_stable
from tkinter import *
from tkinter import filedialog
import os

# """ Initialize Files"""
# root = Tk()
# dirname = filedialog.askdirectory(parent=root,initialdir="/",title='Pick a directory')
# root.destroy()
# dt = 0.0001  # 10 kHz sampling rate
# NoiseSequence = []
# for filename in sorted(os.listdir(dirname)):
#     print(filename)
#     df = pd.read_csv(dirname+"/"+filename, sep=",", header=None,usecols=[1])
#     NoiseSequence = NoiseSequence+df.stack().tolist()
# NoiseSequence = np.array(NoiseSequence)
# np.savetxt("NoiseSequence.csv", NoiseSequence, delimiter=",")
# print("File Created")


"""ReAnalysis Data"""
dt = 0.0001  # 10 kHz sampling rate
fs = 1 / dt
maxt = 20 * 60
freq_range = np.arange(0, 1 / dt + 1 / maxt, 1 / maxt)
print("Reading Start")
NoiseArray = pd.read_csv("NoiseSequence.csv", delimiter=",", engine="c").to_numpy().reshape(-1)[:12000001]
print("Reading Complete")
timeSequence = np.arange(0, len(NoiseArray) * dt, dt)[:12000001]
# print(ss.levy_stable.fit(NoiseArray, f1=0))
fig, ax = plt.subplots(3, 1)
fig.suptitle('RealNoise', fontsize=16)
ax[0].plot(timeSequence, NoiseArray / max(abs(NoiseArray)))
ax[1].hist(NoiseArray, density=True, bins=100, histtype='bar', alpha=0.5)
""" Manually fitted levy stable PDF """
alpha = 1.85
beta = 0.1
x = np.linspace(levy_stable.ppf(0.01, alpha, beta, loc=0, scale=0.007),
                levy_stable.ppf(0.99, alpha, beta, loc=0, scale=0.007), 100)
ax[1].plot(x, levy_stable.pdf(x, alpha, beta, loc=0, scale=0.007),
           'r-', lw=5, alpha=0.6, label='levy_stable pdf')
ax[1].set_xlim([x[0], x[-1]])
Seg_Length_mult = len(freq_range) // 1024
FFTLength = Seg_Length_mult * 1024
WindowLength = len(freq_range) / 2
OverlapLength = len(freq_range) / 5
ax[2].psd(NoiseArray, NFFT=FFTLength, Fs=fs, window=numpy.blackman(FFTLength), noverlap=OverlapLength)
ax[2].set_xscale("log")
ax[2].set_xlim([1e-2, 1 / dt])

"""Apply frequency modification """
gradient = -30
testx = freq_range
testx = np.where(testx == 0, 1e-8, testx)
testx = np.where(testx > 105, 105, testx)
filterIndex = np.where(freq_range > 0.1)[0][0]
y = gradient * np.log10(testx)
y = np.sqrt(10 ** (y / 10))
y = np.where(freq_range < 0.1, y[filterIndex], y)
sym_y = np.zeros(len(y), dtype=float)
halfindex = int(np.ceil(len(y) / 2))
sym_y[:halfindex] = y[:halfindex]
sym_y[halfindex:] = y[:halfindex - 1][::-1]
GenNoise = ss.levy_stable.rvs(alpha, beta, loc=0, scale=0.007, size=len(timeSequence))
fftGenNoise = np.fft.fft(GenNoise, len(GenNoise))
# sym_y = sym_y/np.sqrt(np.mean(sym_y**2))
fftmodified = fftGenNoise * sym_y
adjustedGenNoise = np.real(np.fft.ifft(fftmodified))
ax[0].plot(timeSequence, adjustedGenNoise / max(abs(adjustedGenNoise)), alpha=0.5)
ax[1].hist(adjustedGenNoise, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
ax[2].psd(adjustedGenNoise, NFFT=FFTLength, Fs=fs, window=numpy.blackman(FFTLength), noverlap=OverlapLength, alpha=0.5)
plt.show()
