import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

ro = 1000
""" Main Pipe """
Dm = 0.300
Aream = np.pi * Dm ** 2
am = 1000
Zm = ro * am / Aream
""" Branch Unblocked """
Db1 = 0.300
Areab1 = np.pi * Db1 ** 2
ab1 = 1000
Zb1 = ro * ab1 / Areab1

""" Branch blocked """
Db2 = 0.100
Areab2 = np.pi * Db2 ** 2
ab2 = 1000
Zb2 = ro * ab2 / Areab2

df = 0.1
maxf = 1000

print("Z main =", Zm)
print("Z branch unblocked =", Zb1)
print("Z branch blocked =", Zb2)

frequencyRange = np.arange(0, maxf + df, df)
SaveMatrix = np.zeros((3, len(frequencyRange)), dtype=complex)

t = 0
x1 = 0
B1 = 300
L = 150
B2 = B1 + L
print("x1 = {}, x2 = {}, x3 = {}".format(x1, B1, B2))

for f in tqdm(range(0, len(frequencyRange))):
    omega = 2 * np.pi * frequencyRange[f]
    km = -omega / am
    kb1 = -omega / ab1
    kb2 = -omega / ab2
    j = 1j
    Pin = 1
    """
    P[0] = Prfm
    P[1] = Ptrm
    P[2] = Pb
    P[3] = Pbrf1
    P[4] = Pbtr1
    P[5] = Pbrf2
    P[6] = Pbtr2
    """
    A = [
        [np.exp(-j * km * x1), -np.exp(j * km * x1), 0, 0, 0, 0, 0],
        [0, np.exp(j * km * x1), -np.exp(j * kb1 * x1), -np.exp(-j * kb1 * x1), 0, 0, 0],
        [(np.exp(-j * km * x1) / Zm), np.exp(j * km * x1) / Zm, (np.exp(j * kb1 * x1) / Zb1),
         -np.exp(-j * kb1 * x1) / Zb1, 0, 0, 0],
        [0, 0, np.exp(j * kb1 * B1), np.exp(-j * kb1 * (B1)), -np.exp(j * kb1 * (B1 - B1)),
         -np.exp(-j * kb1 * (B1 - B1)), 0],
        [0, 0, np.exp(j * kb1 * B1) / Zb1, -np.exp(-j * kb1 * (B1)) / Zb1, -np.exp(j * kb1 * (B1 - B1)) / Zb2,
         np.exp(-j * kb1 * (B1 - B1)) / Zb2, 0],
        [0, 0, 0, 0, np.exp(j * kb2 * (B2 - B1)), np.exp(-j * kb2 * (B2 - B1)), -np.exp(j * kb2 * (B2 - (B1 + L)))],
        [0, 0, 0, 0, np.exp(j * kb2 * (B2 - B1)) / Zb2, -np.exp(-j * kb2 * (B2 - B1)) / Zb2,
         -np.exp(j * kb2 * (B2 - (B1 + L))) / Zb1]
    ]

    B = [
        [-Pin * np.exp(j * km * x1)],
        [0],
        [(Pin / Zm) * np.exp(j * km * x1)],
        [0],
        [0],
        [0],
        [0]
    ]

    P = np.linalg.solve(A, B)
    Prfm = P[0][0]
    Ptrm = P[1][0]
    R = Prfm / Pin
    T = Ptrm / Pin
    SaveMatrix[0][f] = frequencyRange[f]
    SaveMatrix[1][f] = R
    SaveMatrix[2][f] = T

f_seq = SaveMatrix[0][:]
R = np.zeros(len(f_seq), dtype=complex)
T = np.zeros(len(f_seq), dtype=complex)

"""Add 0 frequency value"""
Freq = np.zeros(len(f_seq), dtype=complex)
R_freq = SaveMatrix[1][:]
T_freq = SaveMatrix[2][:]
R = R_freq
T = T_freq
Freq = f_seq
Freq.tofile("../Data/BIBFrequencyRangeCONVENTIONAL.csv", sep=",")
R.tofile("../Data/BIBReflectionCoefficientCONVENTIONAL.csv", sep=",")
T.tofile("../Data/BIBTransmissionCoefficientCONVENTIONAL.csv", sep=",")
time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)

"""Test and Verification"""
R = np.real(np.fft.ifft(R, (len(R))))
R.tofile("../Result/BIBReflectionCoefficientCONVENTIONAL.csv", sep=",")
T = np.real(np.fft.ifft(T, (len(T))))
T.tofile("../Result/BIBTransmissionCoefficientCONVENTIONAL.csv", sep=",")
F = np.real(Freq)
plt.figure("1")
plt.plot(time, R, label="Reflection Coefficient", linewidth=3)
plt.plot(time, T, label="Transmission Coefficient")
plt.legend()
plt.show()
