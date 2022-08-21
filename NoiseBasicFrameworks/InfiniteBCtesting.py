import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def FieldMatrix(BasicPipeProperties, L, freq):
    D, a, f, U = BasicPipeProperties
    n = 2  # empirical value for TM friction term R
    g = 9.81  # gravitational acceleration
    A = (np.pi * D ** 2) / 4  # Area
    Q0 = A * U  # initial flow rate
    omega = 2 * np.pi * freq  # T = Theoretical period
    R = (n * f * (Q0 ** (n - 1))) / (2 * g * D * (A ** n))
    "Field Matrix for single conduit"
    mu = np.sqrt(-((omega ** 2) / (a ** 2)) + ((1j * g * A * omega * R) / (a ** 2)))
    Zc = (mu * a ** 2) / (1j * omega * g * A)
    F = np.array([[np.cosh(mu * L), (-1 / Zc) * np.sinh(mu * L), 0],
                  [-Zc * np.sinh(mu * L), np.cosh(mu * L), 0],
                  [0, 0, 1]])
    return F, Zc

def main():
    # Define pipe geometry
    Segments = [250, 250, 250]  # Assume the sensor locates at 250 m from the upstream reservoir and the source locates 750 m  from the upstream reservoir.
    Total_Pipe_length = 750
    feature_mapping = ["Normal", "Source", "Normal"] # Mapping of features as we move along the pipe from the downstream end to the upstream
    Pipe_diameter = [0.3, 0.3, 0.3]
    Pipe_friction_factor = [0.038, 0.038, 0.038]
    SteadyStateVelocity = [1.41,1.41,1.41]
    a = 1000
    g = 9.81

    dfreq = 0.025
    maxfreq = 1000
    freq_range = np.arange(0, maxfreq+dfreq, dfreq)
    time = np.arange(0, (1 / dfreq) + (1 / maxfreq), 1 / maxfreq)
    S = [[1, 0, 0],
         [0, 1, 1],
         [0, 0, 1]]

    H_freq_response = np.zeros((len(freq_range), len(Segments)+1), dtype=complex)

    for i in tqdm(range(1, len(freq_range))):
        freq = freq_range[i]
        F = np.identity(3, dtype=complex)
        for segmentID in reversed(range(len(Segments))):
            D = Pipe_diameter[segmentID]
            a = 1000
            f = Pipe_friction_factor[segmentID]
            U = SteadyStateVelocity[segmentID]
            length = Segments[segmentID]
            BasicPipeProperties = D, a, f, U
            if feature_mapping[segmentID] == "Normal":
                F = F
            else:
                F = S@F
            Field, Zc = FieldMatrix(BasicPipeProperties, length, freq)
            F = Field@F
        c1 = F[0][0]
        c2 = F[0][1]
        c3 = F[0][2]
        c4 = F[1][0]
        c5 = F[1][1]
        c6 = F[1][2]
        A = (np.pi * D ** 2) / 4  # Area
        inifiniteBC = True
        if not inifiniteBC:
            EQA = [[-1, c2],
                   [0, c5]]
            EQB = [[-c3],
                   [-c6]]
            Solution = np.linalg.solve(EQA, EQB)
            Q = Solution[0]
            H = [0]
            q = [0]
            h = Solution[1]
        else:
            EQA = [[-1, c1, 0, c2],
                   [0, c4, -1, c5],
                   [Zc, 0, -1, 0],
                   [0, -Zc, 0, -1]]
            EQB = [[-c3],
                   [-c6],
                   [0],
                   [0]]
            Solution = np.linalg.solve(EQA, EQB)
            Q = Solution[0]
            q = Solution[1]
            H = Solution[2]
            h = Solution[3]
        # else:
        #     EQA = [[A*g/a, -(c1*A*g/a)-c2],
        #            [1, -(c4*A*g/a)-c5]]
        #     EQB = [[c3],
        #            [c6],]
        #     Solution = np.linalg.solve(EQA, EQB)
        #     H = Solution[0]
        #     h = Solution[1]
        #     Q = (A*g/a)*H
        #     q = (A*g/a)*h
        downstreamQ = q
        downstreamH = h
        H_freq_response[i][-1]=downstreamH
        for segmentID in reversed(range(len(Segments))):
            F = np.identity(3, dtype=complex)
            D = Pipe_diameter[segmentID]
            a = 1000
            f = Pipe_friction_factor[segmentID]
            U = SteadyStateVelocity[segmentID]
            length = Segments[segmentID]
            BasicPipeProperties = D, a, f, U
            if feature_mapping[segmentID] == "Normal":
                F = F
            else:
                F = S@F
            F = FieldMatrix(BasicPipeProperties, length, freq)[0]@F
            Result = F@[downstreamQ,downstreamH,[1]]
            downstreamQ = Result[0]
            downstreamH = Result[1]
            H_freq_response[i][segmentID] = downstreamH

    H_time_response = np.zeros((len(freq_range), len(Segments)+1))

    for col in tqdm(range(len(H_freq_response[0]))):
        H_in_freq = H_freq_response[:,col]
        H_in_time = np.real(np.fft.ifft(H_in_freq, len(H_in_freq)))
        H_time_response[:,col] = H_in_time
        plt.figure(col)
        plt.plot(time, H_time_response[:,col])


    plt.show()
    return H_freq_response, H_time_response
