import multiprocessing as mp
import os
import time
from tkinter import filedialog

import pyexcel
from matplotlib import pyplot as plt
from tqdm import tqdm

import TransferMatrix.NoiseAnalysis as NoiseAnalysis
import TransferMatrix.NormalAnalysis as NormalAnalysis
from TransferMatrix.TMTools import *


def round_nearest2(x, a):
    rounded_x = round(x / a) * a
    return rounded_x


def editPipeLength(CopyG, Key, SourceLoc):
    source, target = Key
    PipeLength = float(CopyG.edges[(source, "Pert")]["length"])
    LengthBF = SourceLoc
    LengthAFT = PipeLength - SourceLoc
    CopyG.edges[(source, "Pert")]["length"] = LengthBF
    CopyG.edges[("Pert", target)]["length"] = LengthAFT
    return CopyG


def plotImpulseResponse(time, Sensor1, Sensor2, Title):
    # ZeroPeak1 = np.argwhere(abs(Sensor1) < 10)
    # ZeroPeak2 = np.argwhere(abs(Sensor2) < 60)
    # Sensor1[ZeroPeak1] = 0
    # Sensor2[ZeroPeak2] = 0
    plt.figure(Title)
    plt.plot(time, Sensor1, label="Sensor1")
    plt.plot(time, Sensor2, label="Sensor2")
    plt.legend()


def plotCorrelation(time, Sensor1, Sensor2, Title):
    # ZeroPeak1 = np.argwhere(abs(Sensor1) < 10)
    # ZeroPeak2 = np.argwhere(abs(Sensor2) < 60)
    # Sensor1[ZeroPeak1] = 0
    # Sensor2[ZeroPeak2] = 0
    CrossCorrelation = np.correlate(Sensor1, Sensor2, mode="same")
    CCIndexese = np.argwhere(abs(CrossCorrelation) > 100000000)
    for i in CCIndexese:
        print("CC", time[i] - int(time[-1] / 2), CrossCorrelation[i])
    print("end-------------------------------")
    plt.figure(Title)
    plt.plot(time - int(time[-1] / 2), CrossCorrelation)


def worker(key_index):
    global dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1, timeArray2
    key = list(PossibleCaseDict.keys())[key_index]
    # print(key_index, key)
    Graph, PipeLength, wavespeed = PossibleCaseDict[key]
    SimulationSize = max(int(SimulationSizePerMeter * PipeLength), 1)
    np.random.seed(key_index)
    PertLocations = np.random.uniform(0, PipeLength, SimulationSize)
    SumSensor1 = np.zeros(len(timeArray2))
    SumSensor2 = np.zeros(len(timeArray2))
    GridSize = wavespeed / MaxFreq
    np.random.seed(key_index + 10)
    # print(PertLocations)
    print(key)
    for simulation in tqdm(range(len(PertLocations))):
        Pertlocation = PertLocations[simulation]
        Pertlocation = round_nearest2(Pertlocation, GridSize)
        NewGraph = editPipeLength(copy.deepcopy(Graph), key, Pertlocation)
        SensorResult = NoiseAnalysis.main(NewGraph, Envir, freq_range)
        HFreqResultS1 = SensorResult[Sensor1]["hfreq"]
        HFreqResultS2 = SensorResult[Sensor2]["hfreq"]
        Noise = np.random.normal(0, 0.1, len(timeArray1))
        NPower = np.sum((abs(Noise)) ** 2) / len(Noise)
        Noise = Noise / np.sqrt(NPower)
        Sensor1Time = np.real(np.fft.ifft(HFreqResultS1, len(HFreqResultS1)))
        Sensor2Time = np.real(np.fft.ifft(HFreqResultS2, len(HFreqResultS2)))
        Sensor1WNoise = np.convolve(Sensor1Time, Noise, mode="full")
        Sensor2WNoise = np.convolve(Sensor2Time, Noise, mode="full")
        SumSensor1 = np.add(SumSensor1, Sensor1WNoise)
        SumSensor2 = np.add(SumSensor2, Sensor2WNoise)
        # plotCorrelation(timeArray2, Sensor1WNoise, Sensor2WNoise,
        #                 "Sensor Correlation With Noise of {}, source {}".format(key, Pertlocation))
        print(Pertlocation)
        # plotCorrelation(timeArray1, Sensor1Time, Sensor2Time,
        #                 "Sensor Correlation W/O Noise of {}, source {}".format(key, Pertlocation))
        plotImpulseResponse(timeArray1, Sensor1Time, Sensor2Time,
                            "ImpulseResponse of {} with Source {}".format(key, Pertlocation))
        plt.show()
    return SumSensor1, SumSensor2, PertLocations, key


def init_worker(shared_data):
    global dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1, timeArray2
    dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1, timeArray2 = shared_data


def PossibleSourceLocations(G):
    CaseDict = {}
    for possibleCase in list(G.edges):
        CopyG = copy.deepcopy(G)
        CopyG = SplitEdge(CopyG, possibleCase)
        CaseDict[possibleCase] = (CopyG, G.edges[possibleCase]["length"], G.edges[possibleCase]["wave_velocity"])
    return CaseDict


def CorrelationAnalysis(G, Envir, freq_range, dFreq, MaxFreq, timeArray1, timeArray2):
    SaveDict = {}
    start_time = time.time()
    SimulationSizePerMeter = float(Envir["SimSize"])
    """Start MultiProcessing by start multiple processs"""
    Sensor1 = Envir["Sensor1"]
    Sensor2 = Envir["Sensor2"]
    SumSensor1 = np.zeros(len(timeArray2))
    SumSensor2 = np.zeros(len(timeArray2))
    PossibleCaseDict = PossibleSourceLocations(G)
    SourceLocationRecord = {}
    # CoreCount = min(len(G.edges), 12)
    CoreCount = 1
    with mp.Pool(processes=CoreCount, initializer=init_worker, initargs=(
            (dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1,
             timeArray2),)) as pool:
        for result in pool.map(worker, range(len(list(PossibleCaseDict.keys())))):
            Sensor1Result, Sensor2Result, SourceDist, key = result
            SumSensor1 = np.add(SumSensor1, Sensor1Result)
            SumSensor2 = np.add(SumSensor2, Sensor2Result)
            SourceLocationRecord[key] = SourceDist
    CrossCorrelatedResult = np.correlate(SumSensor1, SumSensor2, mode="same")
    AutoCorrelatedResultSensor1 = np.correlate(SumSensor1, SumSensor1, mode="same")
    AutoCorrelatedResultSensor2 = np.correlate(SumSensor2, SumSensor2, mode="same")
    CorrelatedTime = timeArray2 - ((1 / dFreq))
    SaveTime = np.column_stack(
        (CorrelatedTime, CrossCorrelatedResult, AutoCorrelatedResultSensor1, AutoCorrelatedResultSensor2))
    plt.figure("Correlated Result from Sensor1 and Sensor2")
    plt.plot(CorrelatedTime, CrossCorrelatedResult)
    plt.figure("AutoCorrelated result from Sensor1 and Sensor2")
    plt.plot(CorrelatedTime, AutoCorrelatedResultSensor1, label=Sensor1)
    plt.plot(CorrelatedTime, AutoCorrelatedResultSensor2, label=Sensor2)
    plt.legend()
    SaveDict["Result"] = SaveTime.tolist()
    print("--- %s seconds ---" % (time.time() - start_time))
    return SaveDict


def NormalAnalysis(G, Envir, SubProgressBar, MainProgressBar, ProgressPage, dFreq, freq_range, timeArray1):
    SaveDict = {}
    result_dict, all_result_dict = NormalAnalysis.main(G, Envir, SubProgressBar, MainProgressBar, ProgressPage, dFreq,
                                                       freq_range)
    for edge in G.edges:
        source, target = edge
        SensorHfreq = result_dict[edge]["hfreq"]
        SensorQfreq = result_dict[edge]["qfreq"]
        SourceHfreq = all_result_dict[edge]["source"]["hfreq"]
        SourceQfreq = all_result_dict[edge]["source"]["qfreq"]
        TargetHfreq = all_result_dict[edge]["target"]["hfreq"]
        TargetQfreq = all_result_dict[edge]["target"]["qfreq"]
        SensorHTime = np.real(np.fft.ifft(SensorHfreq, len(SensorHfreq)))
        SensorQTime = np.real(np.fft.ifft(SensorQfreq, len(SensorQfreq)))
        SourceHTime = np.real(np.fft.ifft(SourceHfreq, len(SourceHfreq)))
        SourceQTime = np.real(np.fft.ifft(SourceQfreq, len(SourceQfreq)))
        TargetHTime = np.real(np.fft.ifft(TargetHfreq, len(TargetHfreq)))
        TargetQTime = np.real(np.fft.ifft(TargetQfreq, len(TargetQfreq)))
        plt.figure("({},{}) {}".format(source, target, "Source"))
        plt.plot(timeArray1, SourceHTime)
        plt.figure("({},{}) {}".format(source, target, "Target"))
        plt.plot(timeArray1, TargetHTime)
        time_name = "Pipe {0}-{1} Time".format(source, target)
        freq_name = "Pipe {0}-{1} Freq".format(source, target)
        timeHeading = np.array(
            ["Time", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor",
             "Flow Sensor",
             "Head Target({})".format(target), "Flow Target({})".format(target)])
        freqHeading = np.array(
            ["Freq", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor",
             "Flow Sensor",
             "Head Target({})".format(target), "Flow Target({})".format(target)])
        SaveTime = np.column_stack(
            (timeArray1, SourceHTime, SourceQTime, SensorHTime, SensorQTime, TargetHTime, TargetQTime))
        SaveFreq = np.column_stack(
            (freq_range, abs(SourceHfreq), abs(SourceQfreq), abs(SensorHfreq), abs(SensorQfreq), abs(TargetHfreq),
             abs(TargetQfreq)))
        SaveTime = np.row_stack((timeHeading, SaveTime))
        SaveFreq = np.row_stack((freqHeading, SaveFreq))
        SaveDict[time_name] = SaveTime.tolist()
        SaveDict[freq_name] = SaveFreq.tolist()
    return SaveDict


def main(Graph, Envir, SubProgressBar, MainProgressBar, ProgressPage):
    G = Graph
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    dFreq = float(Envir["df"])
    MaxFreq = int(Envir["MaxFreq"])
    freq_range = np.arange(0, MaxFreq + dFreq, dFreq)
    timeArray1 = np.arange(0, (1 / dFreq) + (1 / MaxFreq), 1 / MaxFreq)
    timeArray2 = np.arange(0, (1 / dFreq) * 2 + (1 / MaxFreq), 1 / MaxFreq)
    timeArray4 = np.arange(0, (1 / dFreq) * 4 + (1 / MaxFreq), 1 / MaxFreq)
    if Envir["FreqMode"] == "Randomized Noise":
        SaveDict = CorrelationAnalysis(G, Envir, freq_range, dFreq, MaxFreq, timeArray1, timeArray2)
    else:
        SaveDict = NormalAnalysis(G, Envir, SubProgressBar, MainProgressBar, ProgressPage, dFreq, freq_range,
                                  timeArray1)
    pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    ProgressPage.destroy()
    pyexcel.free_resources()
    print("File Saved")
    plt.show()
