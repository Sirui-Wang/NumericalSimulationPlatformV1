import os
import time as timeMod
from tkinter import filedialog
import multiprocessing as mp
import pyexcel
from matplotlib import pyplot as plt
from tqdm import tqdm
import TransferMatrix.NoiseAnalysis as NoiseAnalysis
import TransferMatrix.NormalAnalysis as NormalAnalysis
from TransferMatrix.TMTools import *


def worker(Simulations, G, MaxFreq, SortedEdges, NumberOfEdges, Envir, freq_range, Sensor1Superpositioned, Sensor2Superpositioned):
    print(f"Process ID: {os.getpid()}")
    for Simulation in tqdm(Simulations):
        # SubProgressBar["value"] = 0
        SplitedG, PertEdge, PertLocation = RandomPertsinPipes(G, MaxFreq, SortedEdges, NumberOfEdges)
        # print(is_picklable((SplitedG, Envir, freq_range)))
        SensorResult = NoiseAnalysis.main(SplitedG, Envir, freq_range)# , SubProgressBar, ProgressPage)
        # MainProgressBar["value"] += 100 / Dist_Sim_Size
        # ProgressPage.update()
        Sensors = [Envir["Sensor1"], Envir["Sensor2"]]
        HFreqResultS1 = SensorResult[Sensors[0]]["hfreq"]
        HFreqResultS2 = SensorResult[Sensors[1]]["hfreq"]
        Noise = np.random.normal(0, 0.1, np.shape(HFreqResultS1))
        Sensor1Time = np.convolve(np.real(np.fft.ifft(HFreqResultS1, len(HFreqResultS1))), Noise)
        Sensor2Time = np.convolve(np.real(np.fft.ifft(HFreqResultS2, len(HFreqResultS2))), Noise)
        Sensor1Superpositioned = np.add(Sensor1Superpositioned, Sensor1Time)
        Sensor2Superpositioned = np.add(Sensor2Superpositioned, Sensor2Time)
        del (SplitedG)
    return (Sensor1Superpositioned, Sensor2Superpositioned)



def main(Graph, Envir, SubProgressBar, MainProgressBar, ProgressPage):
    G = Graph
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    dFreq = float(Envir["df"])
    maxFreq = int(Envir["MaxFreq"])
    freq_range = np.arange(0, maxFreq + dFreq, dFreq)
    time = np.arange(0, (1 / dFreq) + (1 / maxFreq), 1 / maxFreq)
    if Envir["FreqMode"] == "Randomized Noise":
        SaveDict = {}
        MaxFreq = float(Envir["MaxFreq"])
        NumberOfEdges = len(list(G.edges))
        SortedEdges = sorted(list(G.edges))
        SimulationSize = int(Envir["SimSize"])
        Simulations = range(0, SimulationSize, 1)
        Sensor1Superpositioned = np.zeros(len(freq_range)*2-1)
        Sensor2Superpositioned = np.zeros(len(freq_range) * 2 - 1)
        start_time = timeMod.time()
        """Start MultiProcessing by start multiple processs"""
        CoreCount = 2
        # mp.set_start_method("fork")
        pool = mp.Pool(processes=CoreCount)
        ThreadSize = round(SimulationSize/CoreCount)
        NThreads = range(0, SimulationSize, ThreadSize)
        # for i in range(0, CoreCount-1):
        #     p = pool.apply_async(worker, (Simulations[NThreads[i]:NThreads[i+1]], G, MaxFreq, SortedEdges, NumberOfEdges, Envir, freq_range, Sensor1Superpositioned, Sensor2Superpositioned))
        #     Sensor1MPResult, Sensor2MPResult = p.get()
        #     Sensor1Superpositioned = np.add(Sensor1Superpositioned, Sensor1MPResult)
        #     Sensor2Superpositioned = np.add(Sensor2Superpositioned, Sensor2MPResult)
        # p = pool.apply_async(worker, (Simulations[NThreads[i+1]:], G, MaxFreq, SortedEdges, NumberOfEdges, Envir, freq_range,Sensor1Superpositioned, Sensor2Superpositioned))
        # Sensor1MPResult, Sensor2MPResult = p.get()
        input = []
        for i in range(0, CoreCount - 1):
            input.append((Simulations[NThreads[i]:NThreads[i + 1]], G, MaxFreq, SortedEdges, NumberOfEdges, Envir,
                          freq_range, Sensor1Superpositioned, Sensor2Superpositioned))
        # pool.starmap(worker, input)
        results = pool.starmap(worker, input)
        Sensor1MPResult = np.zeros(len(freq_range)*2-1)
        Sensor2MPResult = np.zeros(len(freq_range)*2-1)
        for Sensor1MPResults, Sensor2MPResults in results:
            print(Sensor1MPResults)
            print(Sensor2MPResults)
            Sensor1MPResult = np.add(Sensor1MPResult, Sensor1MPResults)
            Sensor2MPResult = np.add(Sensor2MPResult, Sensor2MPResults)

        Sensor1Superpositioned = np.add(Sensor1Superpositioned, Sensor1MPResult)
        Sensor2Superpositioned = np.add(Sensor2Superpositioned, Sensor2MPResult)
        print("--- %s seconds ---" % (timeMod.time() - start_time))
        # plt.figure("Pert on {}, {}m from {}".format(PertEdge, PertLocation, PertEdge[0]))
        # plt.plot(time, np.real(np.fft.ifft(HFreqResultS1, len(HFreqResultS1))), label="Sensor at {}".format(Envir["Sensor1"]))
        # plt.plot(time, np.real(np.fft.ifft(HFreqResultS2, len(HFreqResultS2))), label="Sensor at {}".format(Envir["Sensor2"]))
        # plt.legend()
        CrossCorrelatedResult = np.correlate(Sensor1Superpositioned, Sensor2Superpositioned, mode="full")
        CorrelatedTime = np.arange(0, (1 / dFreq) * 4 + (1 / MaxFreq), 1 / MaxFreq)
        SaveTime = np.column_stack((CorrelatedTime, CrossCorrelatedResult))
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
            plt.plot(time, SourceHTime)
            plt.figure("({},{}) {}".format(source, target, "Target"))
            plt.plot(time, TargetHTime)
            time_name = "Pipe {0}-{1} Time".format(source, target)
            freq_name = "Pipe {0}-{1} Freq".format(source, target)
            timeHeading = np.array(
                ["Time", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor", "Flow Sensor",
                 "Head Target({})".format(target), "Flow Target({})".format(target)])
            freqHeading = np.array(
                ["Freq", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor", "Flow Sensor",
                 "Head Target({})".format(target), "Flow Target({})".format(target)])
            SaveTime = np.column_stack((time, SourceHTime, SourceQTime, SensorHTime, SensorQTime, TargetHTime, TargetQTime))
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

