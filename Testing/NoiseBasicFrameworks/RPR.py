import os
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pyexcel
from tqdm import tqdm

S = [[1, 0, 0],
     [0, 1, 1],
     [0, 0, 1]]


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
    return F


def PipeSegmentation(Perts, L):
    Perts.sort()
    PipeSegmentsLength = {0: Perts[0]}
    SegmentID = 1
    for i in range(1, len(Perts)):
        PreviousPert = Perts[i - 1]
        CurrentPert = Perts[i]
        LengthBtwnPerts = CurrentPert - PreviousPert
        PipeSegmentsLength.update({SegmentID: LengthBtwnPerts})
        SegmentID += 1
    PipeSegmentsLength.update({SegmentID: L - Perts[-1]})
    return PipeSegmentsLength


def HeadAtSensor(BasicPipeProperties, q, h, sensor, freq, PipeSegments):
    Segments = list(reversed(sorted(PipeSegments.keys())))
    # Length2Down = PipeSegments[Segments[0]]
    Length2Down = 0
    global S
    U = np.identity(3, dtype=complex)
    for Segment in Segments:
        length = PipeSegments[Segment]
        Length2Down += length
        if sensor > Length2Down:
            U = S @ FieldMatrix(BasicPipeProperties, length, freq) @ U
        else:
            if Segment == Segments[0]:
                length = sensor
            else:
                length = sensor - Length2Down + length
            U = FieldMatrix(BasicPipeProperties, length, freq) @ U
            break
    Result = U @ [[q], h, [1]]
    head = Result[1]
    return head


def analysis(BasicPipeProperties, PipeSegmentsLength, Freq_range, Sensor1, Sensor2):
    global S
    Segments = list(reversed(sorted(PipeSegmentsLength.keys())))
    Sensor1_Head_fResponse = np.zeros((len(Freq_range)), dtype=complex)
    Sensor2_Head_fResponse = np.zeros((len(Freq_range)), dtype=complex)
    for i in tqdm(range(1, len(Freq_range))):
        freq = Freq_range[i]
        D, a, f, U = BasicPipeProperties
        g = 9.81
        A = (np.pi * D ** 2) / 4  # Area
        U = FieldMatrix(BasicPipeProperties, PipeSegmentsLength[Segments[0]], freq)
        for SegmentID in Segments[1::]:
            SegmentLength = PipeSegmentsLength[SegmentID]
            F = FieldMatrix(BasicPipeProperties, SegmentLength, freq)
            U = F @ S @ U
        c1 = U[0][0]
        c2 = U[0][1]
        c3 = U[0][2]
        c4 = U[1][0]
        c5 = U[1][1]
        c6 = U[1][2]
        inifiniteBC = False
        if inifiniteBC:
            EQA = [[-1, c2],
                 [0, c5]]
            EQB = [[-c3],
                 [-c6]]
            Solution = np.linalg.solve(EQA, EQB)
            Q = Solution[0]
            H = 0
            q = 0
            h = Solution[1]
        else:
            EQA = [[-1, c1, 0, c2],
                 [0, c4, -1, c5],
                 [1,0,-(g/a)*A, 0],
                 [0, 1, 0, -(g/a)*A]]
            EQB = [[-c3],
                 [-c6],
                 [0],
                 [0]]
            Solution = np.linalg.solve(EQA, EQB)
            Q = Solution[0]
            q = Solution[1]
            H = Solution[2]
            h = Solution[3]
        Sensor1_Head_fResponse[i] = HeadAtSensor(BasicPipeProperties, q=q, h=h, sensor=Sensor1, freq=freq,
                                                 PipeSegments=PipeSegmentsLength)
        Sensor2_Head_fResponse[i] = HeadAtSensor(BasicPipeProperties, q=q, h=h, sensor=Sensor2, freq=freq,
                                                 PipeSegments=PipeSegmentsLength)
    return Sensor1_Head_fResponse, Sensor2_Head_fResponse


def get_Pert_Name(Perts):
    SheetName = "-"
    for PertLoc in Perts:
        SheetName += str(PertLoc) + "-"
    return SheetName


import math


def round_nearest2(x, a):
    rounded_x = np.zeros(np.shape(x))
    for row_index in range(len(x)):
        for col_index in range(len(x[row_index])):
            rounded_x[row_index][col_index] = round(round(x[row_index][col_index] / a) * a,
                                                    -int(math.floor(math.log10(a))))
    return rounded_x


def main():
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    """ Define pipe geometry and properties for RPR system"""
    # must be long and heavily damped system, such that a impulse introduced within the middle section does not(have
    # minimum) reflection from the either ends
    # all units are SI unit
    L = 1000  # Full pipe length
    D = 0.1  # diameter
    a = 1000  # wavespeed
    f = 0.02  # friction factor
    U = 1.27  # flow velocity
    df = 0.025
    MaxF = 1000
    Freq_range = np.arange(0, MaxF + df, df)
    BasicPipeProperties = (D, a, f, U)
    # Critical component location (distance from upstream)
    Sensor1_dist = 450  # distance from upstream
    Sensor2_dist = 550  # distance from upstream
    NPerts = 1
    SimulationSize = 5
    grid_resolution = L / MaxF
    MonteCarloPerts = np.random.uniform(Sensor1_dist, Sensor2_dist, (SimulationSize, NPerts))
    RoundedMonteCarloPerts = round_nearest2(MonteCarloPerts, grid_resolution)
    # RoundedMonteCarloPerts = [[5555], [5200], [5700]]
    time = np.arange(0, (1 / df) * 2 + (1 / MaxF), 1 / MaxF)
    SaveSensor1Time = time
    SaveSensor2Time = time
    SaveSensor2Freq = Freq_range
    SaveSensor1Freq = Freq_range
    TimeHeading = ["Time"]
    FreqHeading = ["Freq"]
    SimulationCount = 1
    Sensor1Superpositions = np.zeros(np.shape(time))
    Sensor2Superpositions = np.zeros(np.shape(time))
    for Perts in RoundedMonteCarloPerts:
        PipeSegmentsLength = PipeSegmentation(Perts, L)
        print("Simulation ({}/{})".format(SimulationCount, SimulationSize))
        Sensor1, Sensor2 = analysis(BasicPipeProperties, PipeSegmentsLength, Freq_range, L - Sensor1_dist,
                                    L - Sensor2_dist)
        Noise = np.random.normal(0, 0.1, np.shape(Sensor1))
        Sensor1Time = np.convolve(np.real(np.fft.ifft(Sensor1, len(Sensor1))), Noise)
        Sensor2Time = np.convolve(np.real(np.fft.ifft(Sensor2, len(Sensor2))), Noise)
        Name = get_Pert_Name(Perts)
        TimeHeading.append(Name)
        FreqHeading.append(Name)
        SaveSensor1Time = np.column_stack((SaveSensor1Time, Sensor1Time))
        SaveSensor2Time = np.column_stack((SaveSensor2Time, Sensor2Time))
        SaveSensor1Freq = np.column_stack((SaveSensor1Freq, abs(Sensor1 * Noise)))
        SaveSensor2Freq = np.column_stack((SaveSensor2Freq, abs(Sensor2 * Noise)))
        SimulationCount += 1
    for i in range(len(SaveSensor1Time)):
        Sensor1Superpositions[i] = sum(SaveSensor1Time[i][1:])
        Sensor2Superpositions[i] = sum(SaveSensor2Time[i][1:])
    SaveSensor1Time = np.row_stack((np.array(TimeHeading), SaveSensor1Time))
    SaveSensor1Freq = np.row_stack((np.array(FreqHeading), SaveSensor1Freq))
    SaveSensor2Time = np.row_stack((np.array(TimeHeading), SaveSensor2Time))
    SaveSensor2Freq = np.row_stack((np.array(FreqHeading), SaveSensor2Freq))
    CCSuperpositioned = np.correlate(Sensor1Superpositions, Sensor2Superpositions, mode="full")
    plt.plot(CCSuperpositioned)
    SaveSuperpositioned = np.column_stack((Sensor1Superpositions, Sensor2Superpositions))
    SaveSuperpositionedName = np.array(["Sensor1","Sensor2"])
    SaveSuperpositioned = np.row_stack((SaveSuperpositionedName, SaveSuperpositioned))
    SaveDict = {"Time-Sensor1": SaveSensor1Time.tolist(),
                "Time-Sensor2": SaveSensor2Time.tolist(),
                "Freq-Sensor1": SaveSensor1Freq.tolist(),
                "Freq-Sensor2": SaveSensor2Freq.tolist(), }
    SuperpositionedData = {"SaveSuperpositioned":SaveSuperpositioned.tolist()}
    pyexcel.isave_book_as(bookdict=SuperpositionedData, dest_file_name=Output_loc[:-5] + "Superpositioned.xlsx")
    plt.show()
    pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    pyexcel.free_resources()
    print("File Saved")


main()
