import math

import matplotlib.pyplot as plt
import numpy
import numpy as np
from tqdm import tqdm


def ref_within():
    np.random.seed(1107)
    D = 0.3
    A = (np.pi * D ** 2) / 4
    ro = 1000
    a = 1000
    j = 1j
    df = 0.01
    maxf = 1000
    FreqRange = np.arange(0, maxf + df, df)
    omega = 2 * np.pi * FreqRange
    PipeLength = 3100
    h1 = 300
    h2 = 150
    fault_length = 150
    L1 = -50
    dy = 1
    alpha = 4.49666e-5
    R_Freq = np.array(
        np.genfromtxt("../Data/SBReflectionCoefficient.csv", delimiter=",", dtype=complex))  # reflection coefficient
    T_Freq = np.array(np.genfromtxt("../Data/SBTransmissionCoefficient.csv", delimiter=",",
                                    dtype=complex))  # transmission coefficient
    time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)
    time_rotated = time - (1 / df) / 2
    Pert_Sum = np.zeros(len(time))
    """Y1, -1550 ~ -300"""
    for y in tqdm(range(int(-PipeLength / 2), -h1, dy)):
        # for y in tqdm(range(-590, -h1, dy)):
        """Sensor 1"""
        # SgnFunction = (1 / (j * omega))
        Sensor1Direct = ((A * a * np.exp(((alpha * a + j * omega) * (h1 + y)) / a)) / 2)
        Sensor1Scattered = ((A * R_Freq * a * np.exp(-(alpha * a + j * omega) * (h1 + 2 * L1 - y) / a)) / 2)
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
        """Sensor 2 - Complex Conjugate"""
        Sensor2Direct = ((-A * T_Freq[::-1] * a * np.exp((-(omega * j * (y - h2)) / a) + ((y - h2) * alpha))) / 2)
        # Sensor2Direct_NonCongugate = np.exp(-blockage_alpha * 150) * A * T_Freq * a * np.exp(
        #     ((y - h2) * (omega * j + alpha * a)) / a) / 2
        Sensor2Scattered = np.zeros(len(omega), dtype=complex)
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
        # Sensor2SignalTime_NonConjugate = np.real(np.fft.ifft(Sensor2Direct_NonCongugate, (len(Sensor2Direct_NonCongugate))))
        """Frequency domain multiplication - time domain cross correlation"""
        Y1CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y1CC_Time = np.real(np.fft.ifft(Y1CC_Freq, (len(Y1CC_Freq)))) * -1
        # print(math.ceil(len(Y1CC_Time) / 2) - 1)
        Y1CC_Time = numpy.roll(Y1CC_Time, math.ceil(len(Y1CC_Time) / 2) - 1)
        """Debug Verification"""
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
    """Y2, -300 ~ -50"""
    for y in tqdm(range(-h1, L1, dy)):
        """Sensor 1"""
        # SgnFunction = (1 / (j * omega))
        Sensor1Direct = ((A * a * np.exp((-(h1 + y) * (omega * j + alpha * a)) / a)) / 2)
        Sensor1Scattered = ((A * R_Freq * a * np.exp((-(h1 + 2 * L1 - y) * (omega * j + alpha * a)) / a)) / 2)
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
        """Sensor 2 - Complex Conjugate"""
        Sensor2Direct = (-(A * T_Freq[::-1] * a * np.exp(((h2 - y) * (omega * j - alpha * a)) / a)) / 2)
        Sensor2Scattered = np.zeros(len(omega), dtype=complex)
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
        """Frequency domain multiplication - time domain cross correlation"""
        Y2CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y2CC_Time = np.real(np.fft.ifft(Y2CC_Freq, (len(Y2CC_Freq)))) * -1
        # print(math.ceil(len(Y2CC_Time) / 2) - 1)
        Y2CC_Time = numpy.roll(Y2CC_Time, math.ceil(len(Y2CC_Time) / 2) - 1)
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
    """Y3, -50 ~ 300"""
    for y in tqdm(range(L1, h2, dy)):
        """Sensor 1"""
        Sensor1Direct = (A * T_Freq * a * np.exp(-((h1 + y) * (omega * j + alpha * a)) / (a))) / (2)
        Sensor1Scattered = np.zeros(len(omega), dtype=complex)
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
        """Sensor 2 - Complex Conjugate"""
        Sensor2Direct = -(A * a * np.exp((h2 - y) * (omega * j - alpha * a) / a)) / (2)
        Sensor2Scattered = -(A * R_Freq[::-1] * a * np.exp((h2 + y - 2 * L1) * (omega * j - alpha * a) / a)) / (2)
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
        """Frequency domain multiplication - time domain cross correlation"""
        Y3CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y3CC_Time = np.real(np.fft.ifft(Y3CC_Freq, (len(Y3CC_Freq)))) * -1
        # print(math.ceil(len(Y3CC_Time) / 2) - 1)
        Y3CC_Time = numpy.roll(Y3CC_Time, math.ceil(len(Y3CC_Time) / 2) - 1)
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
    """Y4, 300 ~ 1550"""
    for y in tqdm(range(h2, int(PipeLength / 2) - fault_length, dy)):
        """Sensor 1"""
        Sensor1Direct = (A * T_Freq * a * np.exp(-((h1 + y) * (omega * j + alpha * a)) / (a))) / (2)
        Sensor1Scattered = np.zeros(len(omega), dtype=complex)
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
        """Sensor 2 - Complex Conjugate"""
        Sensor2Direct = -(A * a * np.exp(-((h2 - y) * (omega * j - alpha * a)) / a)) / 2
        Sensor2Scattered = -(A * R_Freq[::-1] * a * np.exp((h2 + y - 2 * L1) * (omega * j - alpha * a) / a)) / 2
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
        """Frequency domain multiplication - time domain cross correlation"""
        Y4CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y4CC_Time = np.real(np.fft.ifft(Y4CC_Freq, (len(Y4CC_Freq)))) * -1
        # print(math.ceil(len(Y4CC_Time) / 2) - 1)
        Y4CC_Time = numpy.roll(Y4CC_Time, math.ceil(len(Y4CC_Time) / 2) - 1)
        """Debug Verification"""
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
    Pert_Sum_Normalized = Pert_Sum / np.max(Pert_Sum)
    return Pert_Sum_Normalized


def ref_beyoud():
    np.random.seed(1107)
    D = 0.3
    A = (np.pi * D ** 2) / 4
    ro = 1000
    a = 1000
    Z = ro * a / A
    j = 1j
    df = 0.01
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
    R_Freq = np.array(np.genfromtxt("../Data/SBReflectionCoefficient.csv", delimiter=",", dtype=complex))
    T_Freq = np.array(np.genfromtxt("../Data/SBTransmissionCoefficient.csv", delimiter=",", dtype=complex))
    R = np.real(np.fft.ifft(R_Freq, (len(R_Freq))))
    T = np.real(np.fft.ifft(T_Freq, (len(T_Freq))))
    time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)
    time_rotated = time - (1 / df) / 2

    Pert_Sum = np.zeros(len(time))
    Pert_Sum_withNoise = np.zeros(len(time))
    """Y1, -1550 ~ -300"""
    for y in tqdm(range(int(-PipeLength / 2), -h1, dy)):
        """Sensor 1"""
        Sensor1Direct = ((A * a * np.exp(((alpha * a + j * omega) * (h1 + y)) / a)) / 2)
        Sensor1Scattered = ((A * R_Freq * a * np.exp(-(alpha * a + j * omega) * (h1 + 2 * L - y) / a)) / 2)
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
        """Sensor 2 - Complex Conjugate"""
        Sensor2Direct = ((-A * a * np.exp((h1 - y) * (omega * j - alpha * a) / a)) / 2)
        Sensor2Scattered = -(A * R_Freq[::-1] * a * np.exp(-(h1 + y - 2 * L) * (omega * j - a * alpha) / a)) / 2
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
        """Frequency domain multiplication - time domain cross correlation"""
        Y1CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y1CC_Time = np.real(np.fft.ifft(Y1CC_Freq, (len(Y1CC_Freq)))) * -1
        Y1CC_Time = numpy.roll(Y1CC_Time, math.ceil(len(Y1CC_Time) / 2) - 1)
        # """Debug Verification"""
        # plt.figure(1)
        # plt.plot(time, Sensor1SignalTime, label="sensor1")
        # plt.plot(time, Sensor2SignalTime, label="sensor2")
        # plt.legend()
        # plt.figure(3)
        # plt.plot(time, R)
        # plt.plot(time, T)
        # plt.figure(2)
        # plt.plot(time_rotated, Y1CC_Time)
        # plt.show()
        Pert_Sum += Y1CC_Time
    for y in tqdm(range(-h1, h2, dy)):
        Sensor1Direct = ((A * a * np.exp((-(h1 + y) * (omega * j + alpha * a)) / a)) / 2)
        Sensor1Scattered = ((A * R_Freq * a * np.exp((-(h1 + 2 * L - y) * (omega * j + alpha * a)) / a)) / 2)
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        Sensor2Direct = (-(A * a * np.exp(((h2 - y) * (omega * j - alpha * a)) / a)) / 2)
        Sensor2Scattered = (-(A * R_Freq[::-1] * a * np.exp(-(h2 + y - 2 * L) * (omega * j - alpha * a) / a)) / 2)
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        """Frequency domain multiplication - time domain cross correlation"""
        Y2CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y2CC_Time = np.real(np.fft.ifft(Y2CC_Freq, (len(Y2CC_Freq)))) * -1
        Y2CC_Time = numpy.roll(Y2CC_Time, math.ceil(len(Y2CC_Time) / 2) - 1)
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
    """Y3, -50 ~ 300"""
    for y in tqdm(range(h2, L, dy)):
        """Sensor 1"""
        Sensor1Direct = (A * a * np.exp(-((h1 + y) * (omega * j + alpha * a)) / (a))) / (2)
        Sensor1Scattered = A * R_Freq * a * np.exp(-(h1 + 2 * L - y) * (omega * j + alpha * a) / a) / 2
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        """Sensor 2 - Complex Conjugate"""
        Sensor2Direct = -(A * a * np.exp(-(h2 - y) * (omega * j - alpha * a) / a)) / (2)
        Sensor2Scattered = -(A * R_Freq[::-1] * a * np.exp(-(h2 + y - 2 * L) * (omega * j - alpha * a) / a)) / (2)
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
        """Frequency domain multiplication - time domain cross correlation"""
        Y3CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y3CC_Time = np.real(np.fft.ifft(Y3CC_Freq, (len(Y3CC_Freq)))) * -1
        # print(math.ceil(len(Y3CC_Time) / 2) - 1)
        Y3CC_Time = numpy.roll(Y3CC_Time, math.ceil(len(Y3CC_Time) / 2) - 1)
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
    """Y4, 300 ~ 1550"""
    for y in tqdm(range(L, int(PipeLength / 2) - fault_length, dy)):
        """Sensor 1"""
        Sensor1Direct = (A * T_Freq * a * np.exp(-((h1 + y) * (omega * j + alpha * a)) / (a))) / (2)
        Sensor1Scattered = np.zeros(len(omega), dtype=complex)
        Sensor1SignalFreq = Sensor1Direct + Sensor1Scattered
        Sensor1SignalTime = np.real(np.fft.ifft(Sensor1SignalFreq, (len(Sensor1SignalFreq))))
        """Sensor 2 - Complex Conjugate"""
        Sensor2Direct = -(A * T_Freq[::-1] * a * np.exp(-((h2 - y) * (omega * j - alpha * a)) / a)) / 2
        Sensor2Scattered = np.zeros(len(omega), dtype=complex)
        Sensor2SignalFreq = Sensor2Direct + Sensor2Scattered
        Sensor2SignalTime = np.real(np.fft.ifft(Sensor2SignalFreq, (len(Sensor2SignalFreq))))
        """Frequency domain multiplication - time domain cross correlation"""
        Y4CC_Freq = Sensor1SignalFreq * Sensor2SignalFreq
        Y4CC_Time = np.real(np.fft.ifft(Y4CC_Freq, (len(Y4CC_Freq)))) * -1
        # print(math.ceil(len(Y4CC_Time) / 2) - 1)
        Y4CC_Time = numpy.roll(Y4CC_Time, math.ceil(len(Y4CC_Time) / 2) - 1)
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
    Pert_Sum_Normalized = Pert_Sum / np.max(Pert_Sum)
    return Pert_Sum_Normalized


def PieceWise_integrated_within(a, A, alpha, h1, h2, L, PipeLength, R, T, omega, j):
    As = (A ** 2) * (a ** 2)
    Tconj = T[::-1]  # conjugate of the transmission coefficient
    Rconj = R[::-1]  # conjugate of the reflection coefficient
    """ The following equations are all in the frequency domain"""
    """Piece Wise integration of Function 1: Noise Perturbation on the LHS of both sensor and fault"""
    Y1_Pri = -(Tconj * As * ((np.exp((omega * j * (h2 + h1)) / a - (h2 + h1) * alpha)) - (
        np.exp((omega * j * (h2 + h1)) / a + ((h1 - h2) * alpha) - (alpha * PipeLength))))) / (32 * alpha)
    Y1_Sec = -(Tconj * As * R * (
                (np.exp(-(omega * j * (2 * L + h1 - h2)) / a - (2 * alpha * L + (h2 + 3 * h1) * alpha))) - (np.exp(
            -(omega * j * (2 * L + h1 - h2)) / a - (alpha * PipeLength + 2 * alpha * L + (h2 + h1) * alpha))))) / (
                         32 * alpha)
    """Piece Wise integration of Function 2: Noise Perturbation on the RHS of sensor1 and LHS of fault"""
    Y2_Pri = (Tconj * As * a * j * (np.exp((omega * j * (h2 + h1)) / a - ((h2 + h1) * alpha)) - np.exp(
        -(omega * j * (2 * L + h1 - h2)) / a - ((h2 + h1) * alpha)))) / (32 * omega)
    Y2_Sec = -(Tconj * As * R * (np.exp(-(omega * j * (2 * L + (h1 - h2)) / a) - (h2 + h1) * alpha) - np.exp(
        -(omega * j * (2 * L + h1 - h2) / a) - (2 * alpha * L + (h2 + 3 * h1) * alpha)))) / (32 * alpha)
    """Piece Wise integration of Function 3: Noise Perturbation on the RHS of fault and LHS of Sensor 2"""
    Y3_Pri = -(As * a * T * j * (np.exp(-(omega * j * (h2 + h1) / a) - (h2 + h1) * alpha) - np.exp(
        -(omega * j * (2 * L + (h1 - h2))) / a - (h2 + h1) * alpha))) / (32 * omega)
    Y3_Sec = (Rconj * As * T * (
                np.exp(-(omega * j * (2 * L + (h1 - h2)) / a) + (2 * alpha * L - (3 * h2 + h1) * alpha)) - np.exp(
            -(omega * j * (2 * L + (h1 - h2))) / a - (h2 + h1) * alpha))) / (32 * alpha)
    """Piece Wise integration of Function 4: Noise Perturbation on the RHS of both Sensor 2 and fault"""
    Y4_Pri = (-As * T * (np.exp(-(omega * j * (h2 + h1)) / a - (h2 + h1) * alpha) - np.exp(
        -(omega * j * (h2 + h1)) / a + (2 * alpha * fault_length + (h2 - h1) * alpha - alpha * PipeLength)))) / (
                         32 * alpha)
    Y4_Sec = -(Rconj * As * T * (
                np.exp(-(omega * j * (2 * L + h1 - h2)) / a + (2 * alpha * L - (3 * h2 + h1) * alpha)) - np.exp(
            -(omega * j * (2 * L + h1 - h2)) / a + (
                        2 * alpha * L + 2 * alpha * fault_length - (h2 + h1) * alpha - alpha * PipeLength)))) / (
                         32 * alpha)
    Y = Y1_Pri + Y1_Sec + Y2_Pri + Y2_Sec + Y3_Pri + Y3_Sec + Y4_Pri + Y4_Sec  # Sum of all the piece wise integration
    TimeDomainY = np.fft.ifft(Y, len(Y))  # Inverse Fourier Transform
    TimeDomainY = TimeDomainY.real * -1  # Extract the real part of the time domain signal
    TimeDomainY = np.roll(TimeDomainY,
                          int(len(TimeDomainY) / 2))  # Rotate the time domain signal by 1/2 of the time range
    TimeDomainY = TimeDomainY / max(TimeDomainY)  # Normalization
    return TimeDomainY


def PieceWise_integrated_beyound(a, A, alpha, h1, h2, L, PipeLength, R, T, omega, j):
    As = (A ** 2) * (a ** 2)
    Tconj = T[::-1]  # conjugate of the transmission coefficient
    Rconj = R[::-1]  # conjugate of the reflection coefficient
    """ The following equations are all in the frequency domain"""
    """Piece Wise integration of Function 1: Noise Perturbation on the LHS of both sensor and fault"""
    Y1_P1 = -(As * (np.exp((omega * j * (h2 + h1) / a) - (h2 + h1) * alpha) - np.exp(
        (omega * j * (h2 + h1)) / a + (h1 - h2) * alpha - alpha * PipeLength))) / (32 * alpha)
    Y1_P2 = -(As * Rconj * (np.exp((omega * j * (2 * L + h1 - h2)) / a + (-2 * alpha * L + (h2 - h1) * alpha)) - np.exp(
        (omega * j * (2 * L + h1 - h2)) / a + (-alpha * PipeLength - 2 * alpha * L + (h2 + h1) * alpha)))) / (
                        32 * alpha)
    Y1_P3 = -(As * R * (np.exp(-(omega * j * (2 * L + h1 - h2)) / a - (2 * alpha * L + (h2 + 3 * h1) * alpha)) - np.exp(
        -(omega * j * (2 * L + h1 - h2)) / a - (alpha * PipeLength + 2 * alpha * L + (h2 + h1) * alpha)))) / (
                        32 * alpha)
    Y1_P4 = -(As * R * Rconj * (np.exp(-(omega * j * (h2 + h1)) / a + ((h2 - 3 * h1) * alpha - 4 * alpha * L)) - np.exp(
        -(omega * j * (h2 + h1)) / a + ((h2 - h1) * alpha - 4 * alpha * L - alpha * PipeLength)))) / (32 * alpha)
    """Piece Wise integration of Function 2: Noise Perturbation between sensors and on the LHS of fault"""
    Y2_P1 = As * a * j * (np.exp((h2 + h1) * (omega * j - alpha * a) / a) - np.exp(
        -((h2 + h1) * (omega * j + alpha * a)) / a)) / (32 * omega)
    Y2_P2 = Rconj * As * a * j * (
                -np.exp(omega * j * (2 * L - (3 * h2 + h1)) / a + ((h2 - h1) * alpha - 2 * alpha * L)) + np.exp(
            (2 * L + h1 - h2) * (omega * j - alpha * a) / a)) / (32 * omega)
    Y2_P3 = -(As * R * (np.exp(-(omega * j * (2 * L + (h1 - h2)) / a) + (h2 - h1) * alpha - 2 * alpha * L) - np.exp(
        -(omega * j * (2 * L + h1 - h2) / a) - (2 * alpha * L + (h2 + 3 * h1) * alpha)))) / (32 * alpha)
    Y2_P4 = -(Rconj * As * R * (np.exp(-(omega * j * (h1 + h2)) / a + (3 * h2 - h1) * alpha - 4 * alpha * L) - np.exp(
        -(omega * j * (h1 + h2)) / a + (-4 * alpha * L + (h2 - 3 * h1) * alpha)))) / (32 * alpha)
    """Piece Wise integration of Function 3: Noise Perturbation on the RHS of Sensor 2 and LHS of fault"""
    Y3_P1 = -(As * (np.exp(-(omega * j * (h2 + h1)) / a - (h2 + h1) * alpha) - np.exp(
        -(omega * j * (h2 + h1)) / a + ((h2 - h1) * alpha - 2 * alpha * L)))) / (32 * alpha)
    Y3_P2 = Rconj * As * a * j * (
                np.exp((omega * j * (2 * L - (3 * h2 + h1))) / a + ((h2 - h1) * alpha - 2 * alpha * L)) - np.exp(
            -(omega * j * (h1 + h2)) / a + ((h2 - h1) * alpha - 2 * alpha * L))) / (32 * omega)
    Y3_P3 = -As * R * a * j * (np.exp(-((2 * L + (h1 - h2)) * (omega * j + alpha * a)) / a) - np.exp(
        -(omega * j * (h1 + h2)) / a + ((h2 - h1) * alpha - 2 * alpha * L))) / (32 * omega)
    Y3_P4 = -(Rconj * As * R * (np.exp(-(omega * j * (h1 + h2)) / a + ((h2 - h1) * alpha - 2 * alpha * L)) - np.exp(
        -(omega * j * (h1 + h2)) / a + ((3 * h2 - h1) * alpha - 4 * alpha * L)))) / (32 * alpha)
    """Piece Wise integration of Function 4: Noise Perturbation on the RHS of both Sensor 2 and fault"""
    Y4_P1 = -Tconj * As * T * (np.exp(-(omega * j * (h2 + h1)) / a + ((h2 - h1) * alpha - 2 * alpha * L)) - np.exp(
        -(omega * j * (h2 + h1)) / a + (2 * alpha * fault_length + (h2 - h1) * alpha - alpha * PipeLength))) / (
                        32 * alpha)
    Y = Y1_P1 + Y1_P2 + Y1_P3 + Y1_P4 + Y2_P1 + Y2_P2 + Y2_P3 + Y2_P4 + Y3_P1 + Y3_P2 + Y3_P3 + Y3_P4 + Y4_P1  # Sum of all the piece wise integration
    TimeDomainY = np.fft.ifft(Y, len(Y))  # Inverse Fourier Transform
    TimeDomainY = TimeDomainY.real * -1  # Extract the real part of the time domain signal
    TimeDomainY = np.roll(TimeDomainY,
                          int(len(TimeDomainY) / 2))  # Rotate the time domain signal by 1/2 of the time range
    TimeDomainY = TimeDomainY / max(TimeDomainY)  # Normalization
    return TimeDomainY


def isFaultBeyoundSensor(L, h1, h2):
    if L < -h1 or L > h2:
        print("Fault is beyond the two sensors")
        return True
    else:
        print("Fault is between the two sensors")
        return False


if __name__ == "__main__":
    D = 0.3
    A = (np.pi * D ** 2) / 4
    ro = 1000
    a = 1000
    Z = ro * a / A
    j = 1j
    df = 0.01
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
    R = np.array(
        np.genfromtxt("../Data/SBReflectionCoefficient.csv", delimiter=",", dtype=complex))  # reflection coefficient
    T = np.array(np.genfromtxt("../Data/SBTransmissionCoefficient.csv", delimiter=",",
                               dtype=complex))  # transmission coefficient
    time = np.arange(0, (1 / df) + 1 / maxf, 1 / maxf)  # time range
    time_rotated = time - (1 / df) / 2  # time range rotated by 1/2 of the time range
    omega = np.where(omega == 0, 1e-10, omega)  # replace 0 with 1e-10 to avoid division by zero
    if isFaultBeyoundSensor(L, h1, h2):
        CrossCorrelateFunction = PieceWise_integrated_beyound(a, A, alpha, h1, h2, L, PipeLength, R, T, omega, j)
        ref_result = ref_beyoud()
    else:
        CrossCorrelateFunction = PieceWise_integrated_within(a, A, alpha, h1, h2, L, PipeLength, R, T, omega, j)
        ref_result = ref_within()
    plt.figure("Cross Correlation")
    plt.plot(time_rotated, CrossCorrelateFunction, label="New")  # Plot the time domain signal
    plt.plot(time_rotated, ref_result, label="old", linewidth=3, alpha=0.3)  # Plot the time domain signal
    plt.legend()
    plt.figure("diff")
    plt.plot(time_rotated, CrossCorrelateFunction - ref_result, label="diff")
    plt.legend()
    plt.show()
