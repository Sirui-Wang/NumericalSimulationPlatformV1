import os
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyexcel

S = [[1, 0, 0],
     [0, 1, 1],
     [0, 0, 1]]


def insertSource(SourceLoc, PipePD):
    index = 0
    SourceLoc = 650
    while sum(PipePD["L"][0:index]) < SourceLoc:
        index += 1
        print(index, sum(PipePD["L"][0:index]))
    DistBF = SourceLoc - sum(PipePD["L"][0:index - 1])
    DistAFT = sum(PipePD["L"][0:index]) - SourceLoc
    PipePD["L"][index] = DistBF
    NewRow = np.array(PipePD.iloc[index])
    InsertPD = pd.DataFrame({"L": NewRow[0], "D": NewRow[1], "a": NewRow[2], "f": NewRow[3], "U": NewRow[4]},
                            index=[index])
    PipePD = pd.concat([PipePD.iloc[:index], InsertPD, PipePD.iloc[index:]]).reset_index(drop=True)
    PipePD["L"][index + 1] = DistAFT
    # DistBf =
    # DistAft
    # Insert New Row
    print()
    print("pause")


def main():
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    """ Define pipe geometry and properties for RPR system"""
    Sections = np.array([0, 1, 2, 3, 4])
    Sections_L = np.array([250, 250, 100, 250, 250])
    Sections_D = np.array([0.3, 0.3, 0.2, 0.3, 0.3])
    Sections_a = np.full(5, 1000)
    Sections_f = np.array([0.038, 0.038, 0.044, 0.038, 0.038])
    Sections_U = np.array([1.41, 1.41, 3.18, 1.41, 1.41])
    Data = np.column_stack((Sections_L, Sections_D, Sections_a, Sections_f, Sections_U))
    PipePD = pd.DataFrame(data=Data, index=Sections, columns=["L", "D", "a", "f", "U"])
    df = 0.025
    MaxF = 1000
    Freq_range = np.arange(0, MaxF + df, df)
    # Critical component location (distance from upstream)
    Sensor1_Seg = 0
    Sensor2_Seg = 4
    SimulationSizeBYDSensors = 1500
    SimulationSizeBTWSensors = 3000
    SourcesBYDSensor1 = np.random.uniform(0, 250, SimulationSizeBYDSensors)
    SourcesBYDSensor2 = np.random.uniform(850, 1100, SimulationSizeBYDSensors)
    SourcesBTWSensors = np.random.uniform(250, 850, SimulationSizeBTWSensors)
    CombinedSourceLocation = np.concatenate((SourcesBYDSensor1, SourcesBYDSensor2, SourcesBTWSensors), axis=None)
    for SourceLoc in CombinedSourceLocation:
        NewPipePD = insertSource(SourceLoc, PipePD)
    print("pause")

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
    SaveSuperpositionedName = np.array(["Sensor1", "Sensor2"])
    SaveSuperpositioned = np.row_stack((SaveSuperpositionedName, SaveSuperpositioned))
    SaveDict = {"Time-Sensor1": SaveSensor1Time.tolist(),
                "Time-Sensor2": SaveSensor2Time.tolist(),
                "Freq-Sensor1": SaveSensor1Freq.tolist(),
                "Freq-Sensor2": SaveSensor2Freq.tolist(), }
    SuperpositionedData = {"SaveSuperpositioned": SaveSuperpositioned.tolist()}
    pyexcel.isave_book_as(bookdict=SuperpositionedData, dest_file_name=Output_loc[:-5] + "Superpositioned.xlsx")
    plt.show()
    pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    pyexcel.free_resources()
    print("File Saved")


main()
