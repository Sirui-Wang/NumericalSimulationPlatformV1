import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

ro = 1000
""" Main Pipe """
Dm = 0.300
Aream = np.pi * Dm ** 2
am = 1000
Zm = ro * am / Aream
""" Branch Pipe """
Db = 0.300
Areab = np.pi * Db ** 2
ab = 1000
Zb = ro * ab / Areab
df = 0.01
maxf = 1000

print("Z main =", Zm)
print("Z branch unblocked =", Zb)

frequencyRange = np.arange(0, maxf + df, df)
SaveMatrix = np.zeros((3, len(frequencyRange)), dtype=complex)

t = 0
x1 = 0

for f in tqdm(range(0, len(frequencyRange))):
    omega = 2 * np.pi * frequencyRange[f]
    k = -omega / am
    j = 1j
    Pin = 1
    """
    P[0] = Prfm
    P[1] = Ptrm
    P[2] = Pb
    """
    A = [
        [np.exp(-j * k * x1), -np.exp(j * k * x1), 0],
        [0, np.exp(j * k * x1), -np.exp(j * k * x1)],
        [(np.exp(-j * k * x1) / Zm), np.exp(j * k * x1) / Zm, (np.exp(j * k * x1) / Zb)],
    ]

    B = [
        [-Pin * np.exp(j * k * x1)],
        [0],
        [(Pin / Zm) * np.exp(j * k * x1)],
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
Freq.tofile("../Data/BranchFrequencyRange.csv", sep=",")
R.tofile("../Data/BranchReflectionCoefficient.csv", sep=",")
T.tofile("../Data/BranchTransmissionCoefficient.csv", sep=",")
time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)

"""Test and Verification"""
R = np.real(np.fft.ifft(R, (len(R))))
R.tofile("../Result/BranchReflectionCoefficient.csv", sep=",")
T = np.real(np.fft.ifft(T, (len(T))))
T.tofile("../Result/BranchTransmissionCoefficient.csv", sep=",")
F = np.real(Freq)
plt.figure("1")
plt.plot(time, R, label="Reflection Coefficient", linewidth=3)
plt.plot(time, T, label="Transmission Coefficient")
plt.legend()
plt.show()
