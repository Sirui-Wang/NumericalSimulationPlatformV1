import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

global PipeLength, a, D, U, f, A, Q0, dt, MaxT, HA, HB
# HA = 60
# HB = 39.08
PipeLength = 3100
a = 1000  # WaveSpeed
D = 0.3  # Diameter
U = 0.71  # FlowVelocity
f = 0.038  # FrictionFactor
df = 0.0025
MaxF = 500
dt = 1 / MaxF
MaxT = 1 / df
FreqRange = np.arange(0, MaxF + df, df)
A = (np.pi * D ** 2) / 4
Q0 = A * U
S = [[1, 0, 0],
     [0, 1, 1],
     [0, 0, 1]]


def MOC():
    global HA, HB, D, PipeLength, f, a, MaxT, dt
    g = 9.81
    dx = dt * a
    t = np.linspace(0, MaxT, int(MaxT / dt) + 1)
    x = np.linspace(0, PipeLength, int(PipeLength / dx + 1))

    # ------ Head & velocity ------
    U = np.zeros((len(t), len(x)))
    H = np.zeros((len(t), len(x)))
    H_sta = np.zeros((1, len(x)))
    # ------ Steady State ------
    u_sta = np.sqrt((HA - HB) * 2 * g * D / f / PipeLength)
    # u_sta = 1.38
    Q_sta = u_sta * D ** 2 / 4 * np.pi
    U[0, :] = u_sta
    for i in range(len(x)):
        H_sta[0, i] = HA - u_sta ** 2 / 2 / g * f * x[i] / D
    H[0, :] = H_sta

    # ------ Transient Calculation ------
    dH = a * u_sta / g
    U[0, -1] = 0
    H[0, -1] = H[0, -1] + dH

    for j in tqdm(range(len(t))):
        if j == 0:
            continue
        else:
            for i in range(len(x)):
                if i == 0:
                    H[j, i] = HA
                    U[j, i] = (HA - H[j - 1, i + 1] - f * U[j - 1, i + 1] * abs(
                        U[j - 1, i + 1]) / 2 / g / D * dx) * g / a + U[j - 1, i + 1]
                elif i == len(x) - 1:
                    U[j, i] = 0
                    H[j, i] = a / g * U[j - 1, i - 1] - f * U[j - 1, i - 1] * abs(
                        U[j - 1, i - 1]) / 2 / g / D * dx + H[j - 1, i - 1]
                else:
                    coeff = np.array([[a / g, 1], [a / g, -1]])
                    res = np.array([[-f * U[j - 1, i - 1] * abs(U[j - 1, i - 1]) / 2 / g / D * dx + a / g * U[
                        j - 1, i - 1] + H[j - 1, i - 1]],
                                    [-f * U[j - 1, i + 1] * abs(U[j - 1, i + 1]) / 2 / g / D * dx + a / g * U[
                                        j - 1, i + 1] - H[j - 1, i + 1]]])
                    para = np.linalg.solve(coeff, res)
                    U[j, i] = para[0]
                    H[j, i] = para[1]
        # print("t =%-7.6s" % (t[j]))

    # # ------ Linear regression ------
    # index = np.arange(0, len(t), 80)
    # tran_H = H[index, -1]
    # x_lin = t[index] * a
    # y_lin = np.log(tran_H)
    # const = np.ones(len(x_lin))
    # coeff_m = np.column_stack((const, -x_lin))
    # para = np.linalg.inv(coeff_m.T @ coeff_m) @ coeff_m.T @ y_lin
    # Amp = np.exp(para[0])
    # alpha = para[1]
    # tran_m = Amp * np.exp(-alpha * x_lin)
    # print("u =%-7.6s *e-6" % (alpha))

    # ------ Plot figures ------
    plt.close('all')
    plt.figure(1)
    plt.plot(t, H[:, -1], 'b')
    plt.xlabel("time")
    plt.ylabel("Head")
    plt.grid(True)

    # plt.figure(3)
    # plt.plot(t * a, H[:, -1], 'b')
    # plt.plot(x_lin, tran_m, '--r')
    # plt.xlabel("Location")
    # plt.ylabel("Head")
    # plt.grid(True)

def FieldMatrix(freq, L):
    g = 9.81
    n = 2  # empirical value for TM friction term R
    global a, D, U, f, A, Q0
    omega = 2 * np.pi * freq
    R = (n * f * (Q0 ** (n - 1))) / (2 * g * D * (A ** n))
    mu = np.sqrt(-((omega ** 2) / (a ** 2)) + ((1j * g * A * omega * R) / (a ** 2)))
    print(A * R * g / (2 * a))
    Zc = (mu * a ** 2) / (1j * omega * g * A)
    F = np.array([[np.cosh(mu * L), (-1 / Zc) * np.sinh(mu * L), 0],
                  [-Zc * np.sinh(mu * L), np.cosh(mu * L), 0],
                  [0, 0, 1]])
    return F

def DataatSensor(freq, L, q, h):
    F = FieldMatrix(freq, L)
    U = F @ S
    c1, c2, c3, c4, c5, c6, c7, c8, c9 = U.reshape(9)
    """ RPR, Source is delta-h, H and h are fixed at both end"""
    Q = c1 * q + c2
    H = c4 * q + c5
    # """RPV, source is delta-q, H is fixed, q is fixed"""
    # Q = c2*h+c1
    # H = c5*h + c4
    return H, Q


def TM():
    SensorHeadFreq = np.zeros(len(FreqRange), dtype=complex)
    SensorFlowFreq = np.zeros(len(FreqRange), dtype=complex)
    for freq_index in tqdm(np.arange(1, len(FreqRange))):
        freq = FreqRange[freq_index]
        F = FieldMatrix(freq, PipeLength)
        U = F @ S
        c1, c2, c3, c4, c5, c6, c7, c8, c9 = U.reshape(9)
        """ RPR, Source is delta-h, H and h are fixed at both end"""
        q = -c5 / c4
        h = 0
        Q = c1 * q + c2
        H = 0
        # """RPV, source is delta-q, H is fixed, q is fixed"""
        # h = -c4/c5
        # q = 0
        SensorHead, SensorFlow = DataatSensor(freq, 50, q, h)
        SensorHeadFreq[freq_index] = SensorHead
        SensorFlowFreq[freq_index] = SensorFlow
    SensorHeadTime = (np.real(np.fft.ifft(SensorHeadFreq, len(SensorHeadFreq))))
    SensorFlowTime = np.real(np.fft.ifft(SensorFlowFreq, len(SensorFlowFreq)))
    time = np.arange(0, 1 / df + 1 / MaxF, 1 / MaxF)
    plt.figure("Head")
    plt.plot(time, SensorHeadTime)
    # plt.plot(time, np.cumsum(SensorHeadTime))
    # plt.figure("Flow")
    # plt.plot(time, SensorFlowTime)
    """Linear Regression"""
    Indexes = np.argwhere(SensorHeadTime > 0.00325)
    Y = SensorHeadTime[Indexes]
    X = time[Indexes]
    plt.plot(X, Y, ".")
    x_lin = X * a
    y_lin = np.log(Y)
    const = np.ones(len(x_lin))
    coeff_m = np.column_stack((const, -x_lin))
    para = np.linalg.inv(coeff_m.T @ coeff_m) @ coeff_m.T @ y_lin
    Amp = np.exp(para[0])
    alpha = para[1]
    EstY = Amp * np.exp(-alpha * time * a)
    print(alpha)
    plt.plot(time, EstY)


# MOC()
TM()
plt.show()
