import os
from tkinter import filedialog

import pyexcel
from matplotlib import pyplot as plt

import TransferMatrix.NoiseAnalysis as NoiseAnalysis
import TransferMatrix.NormalAnalysis as NormalAnalysis
from TransferMatrix.TMTools import *


def main(Graph, Envir, progress_bar, ProgressPage):
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
        MaxFreq = float(Envir["MaxFreq"])
        NumberOfEdges = len(list(G.edges))
        SortedEdges = sorted(list(G.edges))
        SimulationSize = int(Envir["SimSize"])
        Sensor1Superpositioned = np.zeros(np.shape(freq_range))
        Sensor2Superpositioned = np.zeros(np.shape(freq_range))
        for Simulation in range(SimulationSize):
            SplitedG = RandomPertsinPipes(G, MaxFreq, SortedEdges, NumberOfEdges)
            SensorResult = NoiseAnalysis.main(SplitedG, Envir, freq_range)
            progress_bar["value"] += 100 / SimulationSize
            ProgressPage.update()
            Sensors = [Envir["Sensor1"], Envir["Sensor2"]]
            HFreqResultS1 = SensorResult[Sensors[0]]["hfreq"]
            HFreqResultS2 = SensorResult[Sensors[1]]["hfreq"]
            Noise = np.random.normal(0, 0.1, np.shape(HFreqResultS1))
            Sensor1Time = np.convolve(np.real(np.fft.ifft(HFreqResultS1, len(HFreqResultS1))), Noise)
            Sensor2Time = np.convolve(np.real(np.fft.ifft(HFreqResultS2, len(HFreqResultS2))), Noise)
            print("test")
    else:
        result_dict, all_result_dict = NormalAnalysis.main(Graph, Envir, progress_bar, ProgressPage, dFreq, freq_range)
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
