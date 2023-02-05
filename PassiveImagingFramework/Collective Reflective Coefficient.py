import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

ro = 1000
""" Pipe 1 """
D1 = 0.300
Area1 = np.pi * D1 ** 2
a1 = 1000
Z1 = ro * a1 / Area1
""" Pipe 2 """
D2 = 0.1
Area2 = np.pi * D2 ** 2
a2 = 1000
Z2 = ro * a2 / Area2

df = 0.1
maxf = 1000

print("Z1 =", Z1)
print("Z2 =", Z2)

frequencyRange = np.arange(df, maxf + df, df)
SaveMatrix = np.zeros((7, len(frequencyRange)), dtype=complex)

t = 0
x1 = 0
L = 150
x2 = x1 + L
print("x1 =", x1, ";", "x2 =", x2)

for f in tqdm(range(0, len(frequencyRange))):
    omega = 2 * np.pi * frequencyRange[f]
    k1 = -omega / a1
    k2 = -omega / a2
    j = 1j
    Pin = 1
    """
    P[0] = Prf1
    P[1] = Ptr1
    P[2] = Prf2
    P[3] = Ptr2
    """
    A = [[np.exp(j * (-k1 * x1 - omega * t)), -np.exp(j * (k2 * x1 - omega * t)), -np.exp(j * (-k2 * x1 - omega * t)),
          0],
         [np.exp(j * (-k1 * x1 - omega * t)) / Z1, np.exp(j * (k2 * x1 - omega * t)) / Z2,
          -np.exp(j * (-k2 * x1 - omega * t)) / Z2, 0],
         [0, np.exp(j * (k2 * x2 - omega * t)), np.exp(j * (-k2 * x2 - omega * t)),
          -np.exp(j * (k1 * (x2 - L) - omega * t))],
         [0, np.exp(j * (k2 * x2 - omega * t)) / Z2, -np.exp(j * (-k2 * x2 - omega * t)) / Z2,
          -np.exp(j * (k1 * (x2 - L) - omega * t)) / Z1]]

    B = [[-Pin * np.exp(j * (k1 * x1 - omega * t))],
         [(Pin * np.exp(j * ((k1 * x1 - omega * t)))) / Z1],
         [0],
         [0]]

    P = np.linalg.solve(A, B)
    Prf1 = P[0][0]
    Ptr2 = P[3][0]
    Prf1Alt = -((Z2 - Z1) * (Z2 + Z1) * (
                (np.exp(2 * j * x1 * (k1 + k2))) - (np.exp(2 * j * k2 * x2 + 2 * j * k1 * x1)))) / \
              ((((Z2 ** 2) - 2 * Z1 * Z2 + (Z1 ** 2)) * np.exp(2 * j * k2 * x2)) - (
                          ((Z2 ** 2) + 2 * Z1 * Z2 + (Z1 ** 2)) * np.exp(2 * j * k2 * x1)))
    Ptr2Alt = -(4 * Z1 * Z2 * np.exp(j * k2 * (x1 + x2))) / \
              (((Z2 ** 2) - 2 * Z1 * Z2 + (Z1 ** 2)) * np.exp(2 * j * k2 * x2) - (
                          ((Z2 ** 2) + 2 * Z1 * Z2 + (Z1 ** 2)) * np.exp(2 * j * k2 * x1)))
    R = Prf1 / Pin
    T = Ptr2 / Pin
    R1 = Prf1Alt
    T1 = Ptr2Alt
    SaveMatrix[0][f] = frequencyRange[f]
    SaveMatrix[1][f] = R
    SaveMatrix[2][f] = T
    SaveMatrix[3][f] = R1
    SaveMatrix[4][f] = T1

f_seq = SaveMatrix[0][:]
R = np.zeros(len(f_seq)+1,dtype=complex)
T = np.zeros(len(f_seq)+1,dtype=complex)


"""Add 0 frequency value"""
Freq = np.zeros(len(f_seq)+1, dtype=complex)
R_freq = SaveMatrix[1][:]
T_freq = SaveMatrix[2][:]
R_freq_Alt = SaveMatrix[3][:]
T_freq_Alt = SaveMatrix[4][:]
R[0] = 0+ 0j
T[0] = 1+ 0j
Freq[0] = 0+0j
R[1::] = R_freq
T[1::] = T_freq
Freq[1::] = f_seq
Freq.tofile("FrequencyRange.csv", sep=",")
R.tofile("ReflectionCoefficient.csv", sep=",")
T.tofile("TransmissionCoefficient.csv", sep=",")
time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)


# """Test and Verification"""
# time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)
# R = np.real(np.fft.ifft(R, (len(R))))
# T = np.real(np.fft.ifft(T, (len(T))))
# F = np.real(Freq)
# R_Alt = np.real(np.fft.ifft(R_freq_Alt, (len(R_freq_Alt))))
# T_Alt = np.real(np.fft.ifft(T_freq_Alt, (len(T_freq_Alt))))
# plt.figure("Diff")
# plt.plot(time[:-1], R[:-1]-R_Alt)
# plt.plot(time[:-1], T[:-1]-T_Alt)
# plt.figure("1")
# plt.plot(time, R, label="Reflection Coefficient")
# plt.plot(time, T, label="Transmission Coefficient")
# plt.legend()
# plt.show()
