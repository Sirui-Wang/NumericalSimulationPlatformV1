import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from matplotlib.animation import PillowWriter
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
"""Frequency Range Configuration"""
df = 0.01
maxf = 100
frequencyRange = np.arange(df, maxf + df, df)
time = np.arange(1 / maxf, (1 / df) + 1 / maxf, 1 / maxf)
"""Define other parameters"""
t = 1
L = 150


def calculate_R_F(x1_location):
    """Initialize save matrix"""
    SaveMatrix = np.zeros((3, len(frequencyRange)), dtype=complex)
    x1 = x1_location
    x2 = x1 + L
    """Solve for R and T at each frequency"""
    for f in (range(0, len(frequencyRange))):
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
        A = [[np.exp(j * (-k1 * x1 - omega * t)), -np.exp(j * (k2 * x1 - omega * t)),
              -np.exp(j * (-k2 * x1 - omega * t)), 0],
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
        R = Prf1 / Pin
        T = Ptr2 / Pin
        SaveMatrix[0][f] = frequencyRange[f]
        SaveMatrix[1][f] = R
        SaveMatrix[2][f] = T
    f_seq = SaveMatrix[0][:]
    R_seq = SaveMatrix[1][:]
    T_seq = SaveMatrix[2][:]
    fft_R_seq = np.real(np.fft.ifft(R_seq, (len(f_seq))))
    fft_T_seq = np.real(np.fft.ifft(T_seq, (len(f_seq))))
    yR = fft_R_seq
    yT = fft_T_seq
    return yR, yT


def update1(frame_index):
    global RStorageArray, TStorageArray, lineR, lineT, x_text, location_range, lineR1, lineT1, x_text1
    R_frame = RStorageArray[frame_index]
    T_frame = TStorageArray[frame_index]
    print(max(abs(R_frame)))
    lineR.set_data(time, R_frame)
    lineT.set_data(time, T_frame)
    x_text.set_text("x1 ={}m".format(location_range[frame_index]))
    lineR1.set_data(time, R_frame)
    lineT1.set_data(time, T_frame)
    x_text1.set_text("x1 ={}m".format(location_range[frame_index]))
    return lineR, lineT, x_text, lineR1, lineT1, x_text1


def main():
    global RStorageArray, TStorageArray, lineR, lineT, x_text, location_range, lineR1, lineT1, x_text1
    Rstorage = []
    Tstorage = []
    location_range = np.arange(-5000, 5000, 50)
    for x1 in tqdm(location_range):
        yR, yT = calculate_R_F(x1)
        Rstorage.append(yR)
        Tstorage.append(yT)
    RStorageArray = np.array(Rstorage)
    TStorageArray = np.array(Tstorage)
    fig, axs = plt.subplots(2)
    fig.suptitle("Animated R&F with different x1 value")
    ax1 = axs[0]
    ax1.set_xlim([89, 101])
    ax1.set_ylim([-1, 1])
    lineR, = ax1.plot([], [])
    lineT, = ax1.plot([], [])
    x_text = ax1.text(90.5, .5, '', fontsize=15)
    ax2 = axs[1]
    ax2.set_xlim([-1, 11])
    ax2.set_ylim([-1, 1])
    lineR1, = ax2.plot([], [])
    lineT1, = ax2.plot([], [])
    x_text1 = ax2.text(6.5, .5, '', fontsize=15)
    Figureanimation1 = animation.FuncAnimation(fig, update1, interval=100, frames=range(len(RStorageArray)), blit=True)
    Figureanimation1.save("R&F Animated W Dx.gif", dpi=300, writer=PillowWriter(fps=25))
    plt.show()


main()
