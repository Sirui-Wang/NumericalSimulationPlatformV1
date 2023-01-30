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

df = 0.01
maxf = 100

print("Z1 =", Z1)
print("Z2 =", Z2)

frequencyRange = np.arange(df, maxf + df, df)
SaveMatrix = np.zeros((7, len(frequencyRange)), dtype=complex)

t = 1
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
    SaveMatrix[5][f] = R - R1
    SaveMatrix[6][f] = T - T1

f_seq = SaveMatrix[0][:]
R_seq = SaveMatrix[1][:]
T_seq = SaveMatrix[2][:]
R_seq1 = SaveMatrix[3][:]
T_seq1 = SaveMatrix[4][:]
R_seqdiff = SaveMatrix[5][:]
T_seqdiff = SaveMatrix[6][:]
time = np.arange(1 / maxf, (1 / df) + 1 / maxf, 1 / maxf)
fft_R_seq = np.real(np.fft.ifft(R_seq, (len(f_seq))))
fft_T_seq = np.real(np.fft.ifft(T_seq, (len(f_seq))))
fft_R_seq1 = np.real(np.fft.ifft(R_seq1, (len(f_seq))))
fft_T_seq1 = np.real(np.fft.ifft(T_seq1, (len(f_seq))))
fft_R_seqdiff = np.real(np.fft.ifft(R_seqdiff, (len(f_seq))))
fft_T_seqdiff = np.real(np.fft.ifft(T_seqdiff, (len(f_seq))))
np.real(f_seq).tofile("FrequencyRange.csv", sep=",")
R_seq.tofile("ReflectionCoefficient.csv", sep=",")
T_seq.tofile("TransmissionCoefficient.csv", sep=",")
plt.figure("1")
plt.plot(time, fft_R_seq)
plt.plot(time, fft_T_seq)
plt.figure("Alt")
plt.plot(time, fft_R_seq1)
plt.plot(time, fft_T_seq1)
plt.figure("Diff")
plt.plot(time, fft_R_seqdiff)
plt.plot(time, fft_T_seqdiff)
plt.show()
