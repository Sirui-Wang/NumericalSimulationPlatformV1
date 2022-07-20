import os
from copy import deepcopy
from tkinter import filedialog
import networkx as nx
import numpy as np
import pyexcel

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
    Length2Down = PipeSegments[Segments[0]]
    global S
    for Segment in Segments:
        U = np.identity(3, dtype=complex)
        if sensor > Length2Down:
            length = PipeSegments[Segment]
            U = S @ FieldMatrix(BasicPipeProperties, length, freq) @ U
            Length2Down += length
        else:
            if Segment == Segments[0]:
                length = sensor
            else:
                length = sensor - Length2Down
            U = FieldMatrix(BasicPipeProperties, length, freq) @ U
            break
    Result = U @ [[q], [h]]
    head = Result[1]
    return head


def analysis(BasicPipeProperties, PipeSegmentsLength, Freq_range, Sensor1, Sensor2):
    global S
    Segments = list(reversed(sorted(PipeSegmentsLength.keys())))
    Sensor1_Head_fResponse = np.zeros((len(Freq_range), 1))
    Sensor2_Head_fResponse = np.zeros((len(Freq_range), 1))
    for i in range(1, len(Freq_range)):
        freq = Freq_range[i]
        U = FieldMatrix(BasicPipeProperties, PipeSegmentsLength[Segments[0]], freq)
        for SegmentID in Segments[1::]:
            SegmentLength = PipeSegmentsLength[SegmentID]
            F = FieldMatrix(BasicPipeProperties, SegmentLength, freq)
            U = F @ S @ U
        A = [[-1, U[0][1]], [0, U[1][1]]]
        B = [[-U[0][2]], [-U[1][2]]]
        Solution = np.linalg.solve(A, B)
        Sensor1_Head_fResponse[i] = HeadAtSensor(BasicPipeProperties, q=0, h=Solution[1], sensor=Sensor1, freq=freq,
                                                 PipeSegments=PipeSegmentsLength)
        Sensor2_Head_fResponse[i] = HeadAtSensor(BasicPipeProperties, q=0, h=Solution[1], sensor=Sensor2, freq=freq,
                                                 PipeSegments=PipeSegmentsLength)
    return Sensor1_Head_fResponse, Sensor2_Head_fResponse


def main():
    """ Define pipe geometry and properties for RPR system"""
    # must be long and heavily damped system, such that a impulse introduced within the middle section does not(have
    # minimum) reflection from the either ends
    # all units are SI unit
    L = 10000  # Full pipe length
    D = 0.3  # diameter
    a = 1000  # wavespeed
    f = 0.2  # friction factor
    U = 3  # flow velocity
    df = 0.0025
    MaxF = 100
    Freq_range = np.arange(0, MaxF + df, df)
    BasicPipeProperties = (D, a, f, U)
    Perts = [4500, 5000, 5500]  # distance from upstream
    PipeSegmentsLength = PipeSegmentation(Perts, L)
    # Critical component location (distance from upstream)
    Sensor1 = L - 4000  # distance from downstream
    Sensor2 = L - 6000  # distance from downstream
    Sensor1, Sensor2 = analysis(BasicPipeProperties, PipeSegmentsLength, Freq_range, Sensor1, Sensor2)


main()
