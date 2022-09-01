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


def worker(Simulations):
    global gG, gMaxFreq, gfreq_range, gNumberOfEdges, gEnvir, gSortedEdges
    SplitedG, PertEdge, PertLocation = RandomPertsinPipes(gG, gMaxFreq, gSortedEdges, gNumberOfEdges)
    SensorResult = NoiseAnalysis.main(SplitedG, gEnvir, gfreq_range)
    SensorA = gEnvir["Sensor1"]
    SensorB = gEnvir["Sensor2"]
    HFreqResultS1 = SensorResult[SensorA]["hfreq"]
    HFreqResultS2 = SensorResult[SensorB]["hfreq"]
    return HFreqResultS1, HFreqResultS2


def init_worker(shared_data):
    global gG, gMaxFreq, gfreq_range, gNumberOfEdges, gEnvir, gSortedEdges
    gG, gMaxFreq, gfreq_range, gNumberOfEdges, gEnvir, gSortedEdges = shared_data


def track_job(job, update_interval=3):
    while job._number_left > 0:
        print("Tasks remaining = {0}".format(
            job._number_left * job._chunksize))
        time.sleep(update_interval)

def main(Graph, Envir, SubProgressBar, MainProgressBar, ProgressPage):
    G = Graph
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    start_time = time.time()
    dFreq = float(Envir["df"])
    maxFreq = int(Envir["MaxFreq"])
    freq_range = np.arange(0, maxFreq + dFreq, dFreq)
    time_range = np.arange(0, (1 / dFreq) + (1 / maxFreq), 1 / maxFreq)
    if Envir["FreqMode"] == "Randomized Noise":
        SaveDict = {}
        MaxFreq = float(Envir["MaxFreq"])
        NumberOfEdges = len(list(G.edges))
        SortedEdges = sorted(list(G.edges))
        SimulationSize = int(Envir["SimSize"])
        Simulations = range(0, SimulationSize, 1)
        Sensor1Superpositioned = np.zeros(len(freq_range))  # * 2 - 1)
        Sensor2Superpositioned = np.zeros(len(freq_range))  # * 2 - 1)
        """Start MultiProcessing by start multiple processs"""
        CoreCount = 3
        chunksize = 10  # round(SimulationSize/CoreCount/20)
        Sensors = [Envir["Sensor1"], Envir["Sensor2"]]
        SensorA = Sensors[0]
        SensorB = Sensors[1]
        with mp.Pool(processes=CoreCount, initializer=init_worker,
                     initargs=((G, MaxFreq, freq_range, NumberOfEdges, Envir, SortedEdges),)) as pool:
            for result in tqdm(pool.imap_unordered(worker, Simulations, chunksize=chunksize), total=SimulationSize):
                Noise = np.random.normal(0, 0.1, np.shape(freq_range))
                HFreqResultS1, HFreqResultS2 = result
                Sensor1Time = np.convolve(np.real(np.fft.ifft(HFreqResultS1, len(HFreqResultS1))), Noise, mode="same")
                Sensor2Time = np.convolve(np.real(np.fft.ifft(HFreqResultS2, len(HFreqResultS2))), Noise, mode="same")
                Sensor1Superpositioned = np.add(Sensor1Superpositioned, Sensor1Time)
                Sensor2Superpositioned = np.add(Sensor2Superpositioned, Sensor2Time)
                MainProgressBar["value"] += 100 / SimulationSize
                ProgressPage.update()
        print("--- %s seconds ---" % (time.time() - start_time))
        CrossCorrelatedResult = np.correlate(Sensor1Superpositioned, Sensor2Superpositioned, mode="same")
        CorrelatedTime = np.arange(0, (1 / dFreq) + (1 / MaxFreq), 1 / MaxFreq)
        SaveTime = np.column_stack((CorrelatedTime, CrossCorrelatedResult))
        plt.figure("Random Source")
        plt.plot(time_range, np.real(np.fft.ifft(HFreqResultS1, len(HFreqResultS1))), label=SensorA)
        plt.plot(time_range, np.real(np.fft.ifft(HFreqResultS2, len(HFreqResultS2))), label=SensorB)
        plt.legend()
        plt.figure("Correlated Result from Sensor1 and Sensor2")
        plt.plot(CorrelatedTime, CrossCorrelatedResult)
        SaveDict["Result"] = SaveTime.tolist()
        pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    else:
        result_dict, all_result_dict = NormalAnalysis.main(Graph, Envir, SubProgressBar, MainProgressBar, ProgressPage, dFreq, freq_range)
        SaveDict = {}
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
            plt.plot(time_range, SourceHTime)
            plt.figure("({},{}) {}".format(source, target, "Target"))
            plt.plot(time_range, TargetHTime)
            time_name = "Pipe {0}-{1} Time".format(source, target)
            freq_name = "Pipe {0}-{1} Freq".format(source, target)
            timeHeading = np.array(
                ["Time", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor", "Flow Sensor",
                 "Head Target({})".format(target), "Flow Target({})".format(target)])
            freqHeading = np.array(
                ["Freq", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor",
                 "Flow Sensor",
                 "Head Target({})".format(target), "Flow Target({})".format(target)])
            SaveTime = np.column_stack(
                (time_range, SourceHTime, SourceQTime, SensorHTime, SensorQTime, TargetHTime, TargetQTime))
            SaveFreq = np.column_stack(
                (freq_range, abs(SourceHfreq), abs(SourceQfreq), abs(SensorHfreq), abs(SensorQfreq), abs(TargetHfreq),
                 abs(TargetQfreq)))
            SaveTime = np.row_stack((timeHeading, SaveTime))
            SaveFreq = np.row_stack((freqHeading, SaveFreq))
            SaveDict[time_name] = SaveTime.tolist()
            SaveDict[freq_name] = SaveFreq.tolist()
        pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    ProgressPage.destroy()
    pyexcel.free_resources()
    print("File Saved")
    plt.show()

