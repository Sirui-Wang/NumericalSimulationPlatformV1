import math
import multiprocessing as mp
import os
import threading
import time
from multiprocessing.pool import ThreadPool
from tkinter import filedialog

import pyexcel
from matplotlib import pyplot as plt
from tqdm import tqdm

import TransferMatrix.TMCoreEngine as NoiseAnalysisEngine
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
    plt.figure(Title)
    plt.plot(time, Sensor1, label="Sensor1")
    plt.plot(time, Sensor2, label="Sensor2")
    plt.legend()


def plotCorrelation(time, Sensor1, Sensor2, Title):
    CrossCorrelation = np.correlate(Sensor1, Sensor2, mode="same")
    CCIndexese = np.argwhere(abs(CrossCorrelation) > 100000000)
    for i in CCIndexese:
        print("CC", time[i] - int(time[-1] / 2), CrossCorrelation[i])
    print("end-------------------------------")
    plt.figure(Title)
    plt.plot(time - int(time[-1] / 2), CrossCorrelation)


def init_worker2(shared_data2):
    global GridSize, Graph, key, Envir, freq_range, ThreadID, key_index
    GridSize, Graph, key, Envir, freq_range, ThreadID, key_index = shared_data2


def worker2(Pertlocation):
    global GridSize, Graph, key, Envir, freq_range, ThreadID, key_index
    Pertlocation = round_nearest2(Pertlocation, GridSize)
    NewGraph = editPipeLength(copy.deepcopy(Graph), key, Pertlocation)
    Sensors_list = [Envir["Sensor1"], Envir["Sensor2"]]
    SensorResult = NoiseAnalysisEngine.main(NewGraph, Envir, freq_range, Sensors_list)
    return SensorResult, Pertlocation


def worker(key_index):
    global dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1, timeArray2, ref_length
    ThreadID = threading.get_ident()
    key = list(PossibleCaseDict.keys())[key_index]
    Graph, PipeLength, wavespeed = PossibleCaseDict[key]
    SimulationSize = max(int(SimulationSizePerMeter * PipeLength), 1)
    np.random.seed(key_index + 10)
    PertLocations = np.random.uniform(0, PipeLength, SimulationSize)
    SumSensor1 = np.zeros(len(timeArray1))
    SumSensor2 = np.zeros(len(timeArray1))
    GridSize = wavespeed / MaxFreq
    try:
        SubCoreCount = min(math.ceil(SimulationSize / ref_length), 12)
    except ZeroDivisionError:
        SubCoreCount = 1
    with mp.Pool(processes=SubCoreCount, initializer=init_worker2, initargs=(
            (GridSize, Graph, key, Envir, freq_range, ThreadID, key_index),)) as pool:
        for result in tqdm(pool.imap(worker2, PertLocations), total=SimulationSize):
            SensorResult, Pertlocation = result
            np.random.seed(np.int(100000 * key_index) + np.int(Pertlocation))
            HFreqResultS1 = SensorResult[Sensor1]["hfreq"]
            HFreqResultS2 = SensorResult[Sensor2]["hfreq"]
            Sensor1Time = np.real(np.fft.ifft(HFreqResultS1, (len(HFreqResultS1))))
            Sensor2Time = np.real(np.fft.ifft(HFreqResultS2, (len(HFreqResultS2))))
            SumSensor1 = np.add(SumSensor1, Sensor1Time)
            SumSensor2 = np.add(SumSensor2, Sensor2Time)
    return SumSensor1, SumSensor2, PertLocations, key


def init_worker(shared_data):
    global dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1, timeArray2, ref_length
    dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1, timeArray2, ref_length = shared_data


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
    """Find GCD for subprocess core count"""
    length_list = []
    for edge in list(G.edges.keys()):
        length_list.append(int(G.edges[edge]["length"]))
    length_median = np.median(length_list)
    ref_length = int(length_median * SimulationSizePerMeter)
    print(ref_length)
    """Start MultiProcessing by start multiple processs"""
    Sensor1 = Envir["Sensor1"]
    Sensor2 = Envir["Sensor2"]
    SumSensor1 = np.zeros(len(timeArray1))
    SumSensor2 = np.zeros(len(timeArray1))
    PossibleCaseDict = PossibleSourceLocations(G)
    SourceLocationRecord = {}
    CoreCount = min(len(G.edges), 12)
    # CoreCount = 1
    with ThreadPool(processes=CoreCount, initializer=init_worker, initargs=(
            (dFreq, MaxFreq, freq_range, Envir, PossibleCaseDict, SimulationSizePerMeter, Sensor1, Sensor2, timeArray1,
             timeArray2, ref_length),)) as pool:
        for result in pool.map(worker, range(len(list(PossibleCaseDict.keys())))):
            Sensor1Result, Sensor2Result, SourceDist, key = result
            SumSensor1 = np.add(SumSensor1, Sensor1Result)
            SumSensor2 = np.add(SumSensor2, Sensor2Result)
            SourceLocationRecord[key] = SourceDist
    CrossCorrelatedResult = np.correlate(SumSensor1, SumSensor2, mode="full")
    AutoCorrelatedResultSensor1 = np.correlate(SumSensor1, SumSensor1, mode="full")
    AutoCorrelatedResultSensor2 = np.correlate(SumSensor2, SumSensor2, mode="full")
    CorrelatedTime = timeArray2 - np.median(timeArray2)
    # SaveTime = np.column_stack(
    #     (CorrelatedTime, CrossCorrelatedResult, AutoCorrelatedResultSensor1))
    SaveTime = np.column_stack(
        (timeArray1, SumSensor1, SumSensor2))
    plt.figure("Cross Correlated Result from Sensor1 and Sensor2")
    plt.plot(CorrelatedTime, CrossCorrelatedResult)
    plt.figure("AutoCorrelated result from Sensor1 and Sensor2")
    plt.plot(CorrelatedTime, AutoCorrelatedResultSensor1, label=Sensor1)
    plt.plot(CorrelatedTime, AutoCorrelatedResultSensor2, label=Sensor2)
    plt.legend()
    SaveDict["Result"] = SaveTime.tolist()
    print("--- %s seconds ---" % (time.time() - start_time))
    return SaveDict


def NormalAnalysis(G, Envir, freq_range, timeArray1):
    SaveDict = {}
    Sensors_list = []
    for edge in list(G.edges()):
        if G.edges[edge]["HasSensor"]:
            Sensors_list.append(G.edges[edge]["SensorLocation"])
    SaveTime = timeArray1
    for i in range(len(Sensors_list)):
        Sensor = Sensors_list[i]
        SensorResult = NoiseAnalysisEngine.main(G, Envir, freq_range, Sensors_list)
        HFreqResultS1 = SensorResult[Sensor]["hfreq"]
        Sensor1Time = np.real(np.fft.ifft(HFreqResultS1, (len(HFreqResultS1))))
        SaveTime = np.column_stack((SaveTime, Sensor1Time))
        plt.figure("IRF at {}".format(Sensor))
        plt.plot(timeArray1, Sensor1Time)
    SaveDict["Result"] = SaveTime.tolist()
    return SaveDict


def main(Graph, Envir, SubProgressBar, MainProgressBar, ProgressPage):
    print("yes")
    G = Graph
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    dFreq = float(Envir["df"])
    MaxFreq = int(Envir["MaxFreq"])
    freq_range = np.arange(0, MaxFreq + dFreq, dFreq)
    timeArray1 = np.arange(0, (1 / dFreq) + (1 / MaxFreq), (1 / MaxFreq))
    timeArray2 = np.arange(0, (1 / dFreq) * 2 + (1 / MaxFreq), 1 / MaxFreq)
    if Envir["FreqMode"] == "Randomized Noise":
        SaveDict = CorrelationAnalysis(G, Envir, freq_range, dFreq, MaxFreq, timeArray1, timeArray2)
    else:
        SaveDict = NormalAnalysis(G, Envir, freq_range, timeArray1)
    pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    ProgressPage.destroy()
    pyexcel.free_resources()
    print("File Saved")
    plt.show()
    active = mp.active_children()
    for child in active:
        child.terminate()
