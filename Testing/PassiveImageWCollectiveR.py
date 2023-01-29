import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

R_Freq = np.array(np.genfromtxt("ReflectionCoefficient.csv", delimiter=",", dtype=complex))
T_Freq = np.array(np.genfromtxt("TransmissionCoefficient.csv", delimiter=",", dtype=complex))

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
h = 300
L = -50
dy = 1
F_Freq = np.zeros(len(omega))
F_Freq[1:] = ((D ** 4) * (a ** 2) * (np.pi ** 2)) / (64 * omega[1:] ** 2)
alpha = 4.49666e-5
Y_SUM = np.zeros(len(F_Freq))
"""Y1, -1550 ~ -300"""

F = np.real(np.fft.fft(F_Freq, (len(F_Freq))) / len(F_Freq))
R = np.real(np.fft.fft(R_Freq, (len(R_Freq))) / len(F_Freq))
T = np.real(np.fft.fft(T_Freq, (len(T_Freq))) / len(F_Freq))
print("Pause")
time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)
for y in tqdm(range(-1550, -300, dy)):
    Direct_Impulse_Freq = T_Freq * np.exp((2 * j * omega * h / a) + (2 * y * alpha))
    Direct_Impulse_Time = np.real(np.fft.ifft(Direct_Impulse_Freq, (len(Direct_Impulse_Freq))))
    Scattered_Impulse_Freq = R_Freq * T_Freq * np.exp(((-(2 * h) + (2 * y - 2 * L)) * alpha) - (2 * j * omega * L / a))
    Scattered_Impulse_Time = np.real(np.fft.ifft(Scattered_Impulse_Freq, (len(Scattered_Impulse_Freq))))
    Total_Impulse_Time = Direct_Impulse_Time + Scattered_Impulse_Time
    Y_SUM += Total_Impulse_Time
    # Y_single = F * T * np.exp((2 * j * omega * h / a) + (2 * y * alpha)) + F * R * T * np.exp(
    #     ((-(2 * h) + (2 * y - 2 * L)) * alpha) - (2 * j * omega * L / a))
"""Y2, -300 ~ -50"""
for y in tqdm(range(-300, -50, dy)):
    Y_single = F * T * np.exp(-(2 * h * a) - (2 * j * y * omega / a)) + F * R * T * np.exp(
        (-2 * h + (2 * y - 2 * L)) * alpha - (2 * j * omega * L / a))
    Y_single = np.real(np.fft.fft(Y_single, (len(Y_single)))) / (len(Y_single))
    Y_SUM += Y_single
"""Y3, -50 ~ 300"""
for y in tqdm(range(-50, 300, dy)):
    Y_single = F * T * R * np.exp(((-2 * h + (2 * L - 2 * y)) * alpha) - (2 * j * omega * L / a)) + F * T * np.exp(
        (-2 * h * alpha) - (2 * j * omega * y / a))
    Y_single = np.real(np.fft.fft(Y_single, (len(Y_single)))) / (len(Y_single))
    Y_SUM += Y_single
"""Y4, 300 ~ 1550"""
for y in tqdm(range(300, 1550 + dy, dy)):
    Y_single = F * T * R * np.exp(((-2 * h + (2 * L - 2 * y)) * alpha) - (2 * j * omega * L / a)) + F * T * np.exp(
        (-2 * y * alpha) - (2 * j * omega * h / a))
    Y_single = np.real(np.fft.fft(Y_single, (len(Y_single)))) / (len(Y_single))
    Y_SUM += Y_single

time = np.arange(0, (1 / df), 1 / maxf)
plt.plot(time, Y_SUM)
plt.show()

# """Y1, -1550 ~ -300"""
# for y in tqdm(range(-1550, -300, dy)):
#     Y_single = F * T * np.exp((2 * j * omega * h / a) + (2 * y * alpha)) + F * R * T * np.exp(
#         ((-(2 * h) + (2 * y - 2 * L)) * alpha) - (2 * j * omega * L / a))
#     Y_SUM += Y_single
# """Y2, -300 ~ -50"""
# for y in tqdm(range(-300, -50, dy)):
#     Y_single = F * T * np.exp(-(2 * h * a) - (2 * j * y * omega / a)) + F * R * T * np.exp(
#         (-2 * h + (2 * y - 2 * L)) * alpha - (2 * j * omega * L / a))
#     Y_SUM += Y_single
# """Y3, -50 ~ 300"""
# for y in tqdm(range(-50, 300, dy)):
#     Y_single = F * T * R * np.exp(((-2 * h + (2 * L - 2 * y)) * alpha) - (2 * j * omega * L / a)) + F * T * np.exp(
#         (-2 * h * alpha) - (2 * j * omega * y / a))
#     Y_SUM += Y_single
# """Y4, 300 ~ 1550"""
# for y in tqdm(range(300, 1550 + dy, dy)):
#     Y_single = F * T * R * np.exp(((-2 * h + (2 * L - 2 * y)) * alpha) - (2 * j * omega * L / a)) + F * T * np.exp(
#         (-2 * y * alpha) - (2 * j * omega * h / a))
