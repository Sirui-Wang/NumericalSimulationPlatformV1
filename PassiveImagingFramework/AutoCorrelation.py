import math

import matplotlib.pyplot as plt
import numpy
import numpy as np
import pyexcel
from tqdm import tqdm

D = 0.3
A = (np.pi * D ** 2) / 4
ro = 1000
a = 1000
Z = ro * a / A
j = 1j
df = 0.01
maxf = 1000
FreqRange = np.arange(df, maxf + df, df)
omega = 2 * np.pi * FreqRange
k = omega / a
PipeLength = 3100
h1 = 300
h2 = 150
fault_length = 150
L1 = -50
dy = 1
alpha = 4.49666e-5
blockage_alpha = 0.0018473
Y_SUM = np.zeros(len(omega))
R_Freq = np.array(np.genfromtxt("ReflectionCoefficient.csv", delimiter=",", dtype=complex))
T_Freq = np.array(np.genfromtxt("TransmissionCoefficient.csv", delimiter=",", dtype=complex))
R = np.real(np.fft.ifft(R_Freq, (len(R_Freq))))
T = np.real(np.fft.ifft(T_Freq, (len(T_Freq))))
time = np.arange(1 / maxf, (1 / df) + 1 / maxf, 1 / maxf)
time_rotated = time - (1 / df) / 2

Pert_Sum = np.zeros(len(time))
"""Y1, -1550 ~ -300"""
for y in tqdm(range(int(-PipeLength / 2), -h1, dy)):
    # for y in tqdm(range(-590, -h1, dy)):
    """Sensor 1"""
    SgnFunction = (1 / (j * omega))
    Sensor1Direct = ((A * a * np.exp(((alpha * a + j * omega) * (h1 + y)) / a)) / 2)
    Sensor1Scattered = ((A * R_Freq * a * np.exp(-(alpha * a + j * omega) * (h1 + 2 * L1 - y) / a)) / 2)
    Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
    Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
    """Sensor 1 - Complex Conjugate"""
    Sensor1DirectConjugate = (-(A * a * np.exp((-(j * omega-alpha * a ) * (h1 + y)) / a)) / 2)
    Sensor1ScatteredConjugate = (-(A * R_Freq[::-1] * a * np.exp((j * omega-alpha * a ) * (h1 + 2 * L1 - y) / a)) / 2)
    Sensor1ConjugateSignalFreq = Sensor1DirectConjugate + Sensor1ScatteredConjugate
    Sensor1ConjugateSignalTime = np.real(np.fft.ifft(Sensor1ConjugateSignalFreq, (len(Sensor1ConjugateSignalFreq))))
    # Sensor2SignalTime_NonConjugate = np.real(np.fft.ifft(Sensor2Direct_NonCongugate, (len(Sensor2Direct_NonCongugate))))
    """Frequency domain multiplication - time domain cross correlation"""
    Y1CC_Freq = Sensor1SignalFreq * Sensor1ConjugateSignalFreq
    Y1CC_Time = np.real(np.fft.ifft(Y1CC_Freq, (len(Y1CC_Freq)))) * -1
    # print(math.ceil(len(Y1CC_Time) / 2) - 1)
    Y1CC_Time = numpy.roll(Y1CC_Time, math.ceil(len(Y1CC_Time) / 2) - 1)
    """Debug Verification"""
    # plt.figure(1)
    # plt.plot(time, Sensor1SignalTime,label="sensor1")
    # plt.plot(time, Sensor1ConjugateSignalTime, label="sensor2")
    # plt.legend()
    # plt.figure(3)
    # plt.plot(time, R)
    # plt.plot(time, T)
    # plt.figure(2)
    # plt.plot(time_rotated, Y1CC_Time)
    # plt.show()
    Pert_Sum += Y1CC_Time
"""Y2, -300 ~ -50"""
for y in tqdm(range(-h1, L1, dy)):
    """Sensor 1"""
    SgnFunction = (1 / (j * omega))
    Sensor1Direct = ((A * a * np.exp((-(h1 + y) * (omega * j + alpha * a)) / a)) / 2)
    Sensor1Scattered = ((A * R_Freq * a * np.exp((-(h1 + 2 * L1 - y) * (omega * j + alpha * a)) / a)) / 2)
    Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
    Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
    """Sensor 1 - Complex Conjugate"""
    Sensor1DirectConjugate = (-(A * a * np.exp(((h1 + y) * (omega * j - alpha * a)) / a)) / 2)
    Sensor1ScatteredConjugate = (-(A * R_Freq[::-1] * a * np.exp(((h1 + 2 * L1 - y) * (omega * j - alpha * a)) / a)) / 2)
    Sensor1ConjugateSignalFreq = Sensor1DirectConjugate + Sensor1ScatteredConjugate
    Sensor1ConjugateSignalTime = np.real(np.fft.ifft(Sensor1ConjugateSignalFreq, (len(Sensor1ConjugateSignalFreq))))
    """Frequency domain multiplication - time domain cross correlation"""
    Y2CC_Freq = Sensor1SignalFreq * Sensor1ConjugateSignalFreq
    Y2CC_Time = np.real(np.fft.ifft(Y2CC_Freq, (len(Y2CC_Freq)))) * -1
    # print(math.ceil(len(Y2CC_Time) / 2) - 1)
    Y2CC_Time = numpy.roll(Y2CC_Time, math.ceil(len(Y2CC_Time) / 2) - 1)
    """Debug Verification"""
    # plt.figure(1)
    # plt.plot(time, Sensor1SignalTime)
    # plt.plot(time, Sensor1ConjugateSignalTime)
    # plt.figure(3)
    # plt.plot(time, R)
    # plt.plot(time, T)
    # plt.figure(2)
    # plt.plot(time_rotated, Y2CC_Time)
    # plt.show()
    Pert_Sum += Y2CC_Time

"""Y3, -50 ~ end"""
for y in tqdm(range(L1, int(PipeLength/2)-fault_length, dy)):
    """Sensor 1"""
    Sensor1Direct = (A * T_Freq * a * np.exp(-((h1 + y) * (omega * j + alpha * a)) / (a))) / (2)
    Sensor1Scattered = np.zeros(len(omega), dtype=complex)
    Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
    Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
    """Sensor 2 - Complex Conjugate"""
    Sensor1DirectConjugate = -(A * T_Freq[::-1] * a * np.exp(((h1 + y) * (omega * j - alpha * a)) / (a))) / (2)
    Sensor1ScatteredConjugate = np.zeros(len(omega), dtype=complex)
    Sensor1ConjugateSignalFreq = Sensor1DirectConjugate + Sensor1ScatteredConjugate
    Sensor1ConjugateSignalTime = np.real(np.fft.ifft(Sensor1ConjugateSignalFreq, (len(Sensor1ConjugateSignalFreq))))
    """Frequency domain multiplication - time domain cross correlation"""
    Y3CC_Freq = Sensor1SignalFreq * Sensor1ConjugateSignalFreq
    Y3CC_Time = np.real(np.fft.ifft(Y3CC_Freq, (len(Y3CC_Freq)))) * -1
    # print(math.ceil(len(Y3CC_Time) / 2) - 1)
    Y3CC_Time = numpy.roll(Y3CC_Time, math.ceil(len(Y3CC_Time) / 2) - 1)
    """Debug Verification"""
    # plt.figure(1)
    # plt.plot(time, Sensor1SignalTime)
    # plt.plot(time, Sensor1ConjugateSignalTime)
    # plt.figure(3)
    # plt.plot(time, R)
    # plt.plot(time, T)
    # plt.figure(2)
    # plt.plot(time_rotated, Y3CC_Time)
    # plt.show()
    Pert_Sum += Y3CC_Time

SaveDict = {}
SaveTime = np.column_stack(
    (time, Pert_Sum))
SaveDict["Result"] = SaveTime.tolist()
pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name="AutoCorrelationResult_Excel.xlsx")
pyexcel.free_resources()
plt.plot(time_rotated, Pert_Sum)
plt.show()
