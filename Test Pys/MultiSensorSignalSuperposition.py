from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

print("Start")
"""Prompt user to select two csv files"""
file1 = filedialog.askopenfilename()
file2 = filedialog.askopenfilename()
"""Read the csv files into numpy arrays"""
data1 = np.genfromtxt(file1, delimiter=',')
print("Data1 loaded")
data2 = np.genfromtxt(file2, delimiter=',')
T1, S2, S3 = data1[:, 0], data1[:, 1], data1[:, 2]
T2, S1, S4 = data2[:, 0], data2[:, 1], data2[:, 2]
InputArray = np.array([S1, S2, S3, S4])
InputLabel = np.array([0, 1, 2, 3])
DistancetoSensor0 = np.array([0, 500, 1100, 1400])
dt = T1[1] - T1[0]
wavespeed = 1000
print("Data loaded")

"""signal filtering"""


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y


cutoff = 25
fs = 1000
order = 6
filteredS1 = butter_lowpass_filter(S1, cutoff, fs, order)
filteredS2 = butter_lowpass_filter(S2, cutoff, fs, order)
filteredS3 = butter_lowpass_filter(S3, cutoff, fs, order)
filteredS4 = butter_lowpass_filter(S4, cutoff, fs, order)
print("Signal filtering done")

"""Compute and plot filteredS1 and filteredS2 cross-correlation"""
plt.figure(1)
corr = np.correlate(filteredS2, filteredS3, mode='full')
time = np.linspace(-T1[-1], T1[-1], len(corr))
plt.plot(time, corr)
plt.title("Cross-correlation between filteredS1 and filteredS2")
plt.xlabel("Time (s)")
plt.ylabel("Cross-correlation")

"""postive roll == delayed signal"""
"""negative roll == advanced signal"""

rolledS1 = np.roll(filteredS1, int(1000))
rolledS4 = np.roll(filteredS4, -int(600))
adjustedS2 = filteredS2 + rolledS1
adjustedS3 = filteredS3 + rolledS4

# """Compute and plot adjustedS2 and adjustedS3 cross-correlation"""
# plt.figure(2)
# corr = np.correlate(adjustedS2, adjustedS3, mode='full')
# time = np.linspace(-T1[-1], T1[-1], len(corr))
# plt.plot(time, corr)
# plt.title("Cross-correlation between adjustedS2 and adjustedS3")
# plt.xlabel("Time (s)")
# plt.ylabel("Cross-correlation")
# plt.show()


"""Compute and plot adjustedS2 auto correlation"""
plt.figure(2)
corr = np.correlate(adjustedS2, adjustedS2, mode='full')
time = np.linspace(-T1[-1], T1[-1], len(corr))
plt.plot(time, corr)
plt.title("Auto-correlation between adjustedS2")
plt.xlabel("Time (s)")
plt.ylabel("Auto-correlation")

plt.show()
