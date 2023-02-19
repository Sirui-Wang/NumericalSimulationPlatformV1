import math

import PassiveImagingFramework.RealNoiseFiles.GenerateNoise as GenerateNoise
import matplotlib.pyplot as plt
import numpy
import numpy as np
import pyexcel
from tqdm import tqdm

np.random.seed(1107)
D = 0.3
A = (np.pi * D ** 2) / 4
ro = 1000
a = 1000
Z = ro * a / A
j = 1j
df = 0.1
maxf = 1000
FreqRange = np.arange(0, maxf + df, df)
omega = 2 * np.pi * FreqRange
k = omega / a
PipeLength = 3100
h1 = 300
h2 = 300
fault_length = 150
L = 400
dy = 1
alpha = 4.49666e-5
# blockage_alpha = 0.0018473
Y_SUM = np.zeros(len(omega))
R_Freq = np.array(np.genfromtxt("ReflectionCoefficient.csv", delimiter=",", dtype=complex))
T_Freq = np.array(np.genfromtxt("TransmissionCoefficient.csv", delimiter=",", dtype=complex))
R = np.real(np.fft.ifft(R_Freq, (len(R_Freq))))
T = np.real(np.fft.ifft(T_Freq, (len(T_Freq))))
time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)
time_rotated = time - (1 / df) / 2

Pert_Sum = np.zeros(len(time))
Pert_Sum_withNoise = np.zeros(len(time))
"""Y1, -1550 ~ -300"""
for y in tqdm(range(int(-PipeLength / 2), -h1, dy)):
    # for y in tqdm(range(-590, -h1, dy)):
    """Sensor 1"""
    # SgnFunction = (1 / (j * omega))
    Sensor1Direct = ((A * a * np.exp(((alpha * a + j * omega) * (h1 + y)) / a)) / 2)
    Sensor1Scattered = ((A * R_Freq * a * np.exp(-(alpha * a + j * omega) * (h1 + 2 * L - y) / a)) / 2)
    Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
    Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
    """Sensor 2 - Complex Conjugate"""
    Sensor2Direct = ((-A * a * np.exp((h1 - y) * (omega * j - alpha * a)/a)) / 2)
    Sensor2Scattered = -(A * R_Freq[::-1] * a * np.exp(-(h1 + y - 2 * L) * (omega * j - a * alpha) / a)) / 2
    Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
    Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
    """Noise"""
    Noise = GenerateNoise.Gen(len(Sensor1SignalTime), "Levy", [1.9, 0, 0, 1])
    CCNoise = np.correlate(Noise, Noise, mode="same")
    """Frequency domain multiplication - time domain cross correlation"""
    Y1CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
    Y1CC_Time = np.real(np.fft.ifft(Y1CC_Freq, (len(Y1CC_Freq)))) * -1
    # print(math.ceil(len(Y1CC_Time) / 2) - 1)
    Y1CC_Time = numpy.roll(Y1CC_Time, math.ceil(len(Y1CC_Time) / 2) - 1)
    Y1CC_Time_Noise = np.convolve(Y1CC_Time, CCNoise, mode="same")
    # """Debug Verification"""
    # plt.figure(1)
    # plt.plot(time, Sensor1SignalTime,label="sensor1")
    # plt.plot(time, Sensor2SignalTime, label="sensor2")
    # plt.legend()
    # plt.figure(3)
    # plt.plot(time, R)
    # plt.plot(time, T)
    # plt.figure(2)
    # plt.plot(time_rotated, Y1CC_Time)
    # plt.show()
    Pert_Sum += Y1CC_Time
    Pert_Sum_withNoise += Y1CC_Time_Noise
"""Y2, -300 ~ -50"""
for y in tqdm(range(-h1, h2, dy)):
    """Sensor 1"""
    # SgnFunction = (1 / (j * omega))
    Sensor1Direct = ((A * a * np.exp((-(h1 + y) * (omega * j + alpha * a)) / a)) / 2)
    Sensor1Scattered = ((A * R_Freq * a * np.exp((-(h1 + 2 * L - y) * (omega * j + alpha * a)) / a)) / 2)
    Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
    Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
    """Sensor 2 - Complex Conjugate"""
    Sensor2Direct = (-(A * a * np.exp(((h2 - y) * (omega * j - alpha * a)) / a)) / 2)
    Sensor2Scattered = (-(A * R_Freq[::-1] * a * np.exp(-(h2 + y - 2 * L) * (omega * j - alpha * a) / a)) / 2)
    Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
    Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
    """Noise"""
    Noise = GenerateNoise.Gen(len(Sensor1SignalTime), "Levy", [1.9, 0, 0, 1])
    CCNoise = np.correlate(Noise, Noise, mode="same")
    """Frequency domain multiplication - time domain cross correlation"""
    Y2CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
    Y2CC_Time = np.real(np.fft.ifft(Y2CC_Freq, (len(Y2CC_Freq)))) * -1
    # print(math.ceil(len(Y2CC_Time) / 2) - 1)
    Y2CC_Time = numpy.roll(Y2CC_Time, math.ceil(len(Y2CC_Time) / 2) - 1)
    Y2CC_Time_Noise = np.convolve(Y2CC_Time, CCNoise, mode="same")
    # """Debug Verification"""
    # plt.figure(1)
    # plt.plot(time, Sensor1SignalTime)
    # plt.plot(time, Sensor2SignalTime)
    # plt.figure(3)
    # plt.plot(time, R)
    # plt.plot(time, T)
    # plt.figure(2)
    # plt.plot(time_rotated, Y2CC_Time)
    # plt.show()
    Pert_Sum += Y2CC_Time
    Pert_Sum_withNoise += Y2CC_Time_Noise
"""Y3, -50 ~ 300"""
for y in tqdm(range(h2, L, dy)):
    """Sensor 1"""
    Sensor1Direct = (A * a * np.exp(-((h1 + y) * (omega * j + alpha * a)) / (a))) / (2)
    Sensor1Scattered = A * R_Freq * a * np.exp(-(h1 + 2 * L - y) * (omega * j + alpha * a) / a) / 2
    Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
    Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
    """Sensor 2 - Complex Conjugate"""
    Sensor2Direct = -(A * a * np.exp(-(h2 - y) * (omega * j - alpha * a) / a)) / (2)
    Sensor2Scattered = -(A * R_Freq[::-1] * a * np.exp(-(h2 + y - 2 * L) * (omega * j - alpha * a) / a)) / (2)
    Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
    Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
    """Noise"""
    Noise = GenerateNoise.Gen(len(Sensor1SignalTime), "Levy", [1.9, 0, 0, 1])
    CCNoise = np.correlate(Noise, Noise, mode="same")
    """Frequency domain multiplication - time domain cross correlation"""
    Y3CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
    Y3CC_Time = np.real(np.fft.ifft(Y3CC_Freq, (len(Y3CC_Freq)))) * -1
    # print(math.ceil(len(Y3CC_Time) / 2) - 1)
    Y3CC_Time = numpy.roll(Y3CC_Time, math.ceil(len(Y3CC_Time) / 2) - 1)
    Y3CC_Time_Noise = np.convolve(Y3CC_Time, CCNoise, mode="same")
    # """Debug Verification"""
    # plt.figure(1)
    # plt.plot(time, Sensor1SignalTime)
    # plt.plot(time, Sensor2SignalTime)
    # plt.figure(3)
    # plt.plot(time, R)
    # plt.plot(time, T)
    # plt.figure(2)
    # plt.plot(time_rotated, Y3CC_Time)
    # plt.show()
    Pert_Sum += Y3CC_Time
    Pert_Sum_withNoise += Y3CC_Time_Noise
"""Y4, 300 ~ 1550"""
for y in tqdm(range(L, int(PipeLength / 2) - fault_length, dy)):
    """Sensor 1"""
    Sensor1Direct = (A * T_Freq * a * np.exp(-((h1 + y) * (omega * j + alpha * a)) / (a))) / (2)
    Sensor1Scattered = np.zeros(len(omega), dtype=complex)
    Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
    Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
    """Sensor 2 - Complex Conjugate"""
    Sensor2Direct = -(A * T_Freq[::-1]* a * np.exp(-((h2 - y) * (omega * j - alpha * a)) / a)) / 2
    Sensor2Scattered = np.zeros(len(omega), dtype=complex)
    Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
    Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
    """Noise"""
    Noise = GenerateNoise.Gen(len(Sensor1SignalTime), "Levy", [1.9, 0, 0, 1])
    CCNoise = np.correlate(Noise, Noise, mode="same")
    """Frequency domain multiplication - time domain cross correlation"""
    Y4CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
    Y4CC_Time = np.real(np.fft.ifft(Y4CC_Freq, (len(Y4CC_Freq)))) * -1
    # print(math.ceil(len(Y4CC_Time) / 2) - 1)
    Y4CC_Time = numpy.roll(Y4CC_Time, math.ceil(len(Y4CC_Time) / 2) - 1)
    Y4CC_Time_Noise = np.convolve(Y4CC_Time, CCNoise, mode="same")
    # """Debug Verification"""
    # plt.figure(1)
    # plt.plot(time, Sensor1SignalTime,label="sensor1")
    # plt.plot(time, Sensor2SignalTime, label="sensor2")
    # plt.legend()
    # plt.figure(3)
    # plt.plot(time, R)
    # plt.plot(time, T)
    # plt.figure(2)
    # plt.plot(time_rotated, Y4CC_Time)
    # plt.show()
    Pert_Sum += Y4CC_Time
    Pert_Sum_withNoise += Y4CC_Time_Noise

SaveDict = {}
SaveTime = np.column_stack(
    (time, Pert_Sum, Pert_Sum_withNoise))
SaveDict["Result"] = SaveTime.tolist()
pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name="PassiveImagingResult_Excel.xlsx")
pyexcel.free_resources()
plt.figure("No Noise")
plt.plot(time_rotated, Pert_Sum)
plt.figure("With Noise")
plt.plot(time_rotated, Pert_Sum_withNoise)
plt.show()
