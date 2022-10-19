import os
from tkinter import filedialog
from tqdm import tqdm
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyexcel
import time as TimeMod

S = [[1, 0, 0],
     [0, 1, 1],
     [0, 0, 1]]
isinifiniteBC = True

def insertSource(SourceLoc, PipePD,Sensor1_Loc, Sensor2_Loc):
    Sensor1Seg, Loc1 = Sensor1_Loc
    Sensor2Seg, Loc2 = Sensor2_Loc
    index = 0
    while sum(PipePD["L"][0:index]) < SourceLoc:
        index += 1
    DistBF = SourceLoc - sum(PipePD["L"][0:index - 1])
    DistAFT = sum(PipePD["L"][0:index]) - SourceLoc
    PipePD["L"][index-1] = DistBF
    NewRow = np.array(PipePD.iloc[index-1])
    InsertPD = pd.DataFrame({"L": NewRow[0], "D": NewRow[1], "A": NewRow[2],"a": NewRow[3], "f": NewRow[4], "U": NewRow[5],"Q0":NewRow[6], "S":1},
                            index=[index-1])
    PipePD = pd.concat([PipePD.iloc[:index-1], InsertPD, PipePD.iloc[index-1:]]).reset_index(drop=True)
    PipePD["L"][index] = DistAFT
    if index-1 <= Sensor1Seg:
        Sensor1Seg +=1
        Sensor1_Loc = (Sensor1Seg, Loc1)
    if index <= Sensor2Seg:
        Sensor2Seg +=1
        Sensor2_Loc = (Sensor2Seg, Loc2)
    return PipePD, Sensor1_Loc, Sensor2_Loc


def FieldMatrix(Pipe_Data, freq):
    g = 9.81
    n = 2  # empirical value for TM friction term R
    omega = 2 * np.pi * freq  # T = Theoretical period
    L,D,A,a,f,U,Q0, S= Pipe_Data
    R = (n * f * (Q0 ** (n - 1))) / (2 * g * D * (A ** n))
    "Field Matrix for single conduit"
    mu = np.sqrt(-((omega ** 2) / (a ** 2)) + ((1j * g * A * omega * R) / (a ** 2)))
    Zc = (mu * a ** 2) / (1j * omega * g * A)
    F = np.array([[np.cosh(mu * L), (-1 / Zc) * np.sinh(mu * L), 0],
                  [-Zc * np.sinh(mu * L), np.cosh(mu * L), 0],
                  [0, 0, 1]])
    return F, Zc


def HeadAtSensor(PipePD, q, h, Sensor_Loc, freq):
    global S
    SensorSeg, SensorLoc = Sensor_Loc
    U = np.identity(3, dtype=complex)
    PipeIndex = PipePD.index[-1]
    while SensorSeg != PipeIndex:
        F, Zc = FieldMatrix(PipePD.iloc[PipeIndex], freq)
        isSource = PipePD["S"][PipeIndex]
        if isSource==1:
            U = F @ S @ U
        else:
            U = F @ U
        PipeIndex -=1
    F, Zc = FieldMatrix(PipePD.iloc[PipeIndex], freq)
    isSource = PipePD["S"][PipeIndex]
    if SensorLoc == 0:
        if isSource==1:
            U = F @ S @ U
        else:
            U = F @ U
    else:
        pass
    Result = U @ [q, h, [1]]
    head = Result[1]
    return head[0]


def analysis(PipePD, Freq_range,Sensor1_Loc, Sensor2_Loc):
    global S, isinifiniteBC
    Sensor1_Head_fResponse = np.zeros((len(Freq_range)), dtype=complex)
    Sensor2_Head_fResponse = np.zeros((len(Freq_range)), dtype=complex)
    for i in range(1, len(Freq_range[1::])):
        freq = Freq_range[i]
        U, Zc = FieldMatrix(PipePD.iloc[-1], freq)
        Zcs = [Zc]
        for PipeIndex in np.array(PipePD.index)[::-1][1::]:
            F, Zc = FieldMatrix(PipePD.iloc[PipeIndex], freq)
            isSource = PipePD["S"][PipeIndex]
            Zcs.append(Zc)
            if isSource == 1:
                U = F @ S @ U
            else:
                U = F @ U
        c1 = U[0][0]
        c2 = U[0][1]
        c3 = U[0][2]
        c4 = U[1][0]
        c5 = U[1][1]
        c6 = U[1][2]
        if isinifiniteBC == False:
            EQA = [[-1, c2],
                 [0, c5]]
            EQB = [[-c3],
                 [-c6]]
            Solution = np.linalg.solve(EQA, EQB)
            Q = Solution[0]
            H = 0
            q = [0]
            h = Solution[1]
        else:
            Zcup = Zcs[-1]
            Zcdown = Zcs[0]
            EQA = np.array([[-1, c1, 0, c2],
                 [0, c4, -1, c5],
                      [Zcup, 0, -1, 0],
                      [0, -Zcdown, 0, -1]])
            EQB = np.array([[-c3],
                 [-c6],
                 [0],
                 [0]])
            Solution = np.linalg.solve(EQA, EQB)
            Q = Solution[0]
            q = Solution[1]
            H = Solution[2]
            h = Solution[3]
        Sensor1_Head_fResponse[i] = HeadAtSensor(PipePD, q=q, h=h, Sensor_Loc=Sensor1_Loc, freq=freq)
        Sensor2_Head_fResponse[i] = HeadAtSensor(PipePD, q=q, h=h, Sensor_Loc=Sensor2_Loc, freq=freq)
    return Sensor1_Head_fResponse, Sensor2_Head_fResponse


def worker(SourceLoc):
    global PipePD, Freq_range,time, Sensor1_Loc, Sensor2_Loc, Sensor1Superpositions, Sensor2Superpositions
    # SourceLoc = 110
    NewPipePD, Sensor1_Loc, Sensor2_Loc = insertSource(SourceLoc, PipePD.copy(deep=True), Sensor1_Loc, Sensor2_Loc)
    Sensor1Head, Sensor2Head = analysis(NewPipePD, Freq_range,Sensor1_Loc,Sensor2_Loc)
    Noise = np.random.normal(0, 1, np.shape(Sensor1Head))
    Sensor1Time = np.convolve(np.real(np.fft.ifft(Sensor1Head, len(Sensor1Head))), Noise, mode="same")
    Sensor2Time = np.convolve(np.real(np.fft.ifft(Sensor2Head, len(Sensor2Head))), Noise, mode="same")
    Sensor1Superpositions = np.add(Sensor1Superpositions,Sensor1Time)
    Sensor2Superpositions = np.add(Sensor2Superpositions,Sensor2Time)
    return Sensor1Superpositions, Sensor2Superpositions, Sensor1Head, Sensor2Head


def init_worker(shared_data):
    global PipePD, Freq_range,time, Sensor1_Loc, Sensor2_Loc, Sensor1Superpositions, Sensor2Superpositions
    PipePD, Freq_range,time, Sensor1_Loc, Sensor2_Loc, Sensor1Superpositions, Sensor2Superpositions = shared_data



def main():
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    """ Define pipe geometry and properties for RPR system"""
    start_time = TimeMod.time()
    Sections = np.array([0, 1, 2, 3, 4])
    Sections_L = np.array([250, 250, 100, 250, 250])
    Sections_D = np.array([0.3, 0.3, 0.2, 0.3, 0.3])
    Sections_a = np.full(len(Sections), 1000)
    Sections_f = np.array([0.038, 0.038, 0.044, 0.038, 0.038])
    Sections_U = np.array([1.41, 1.41, 3.18, 1.41, 1.41])
    Sections_A = (np.pi * Sections_D ** 2) / 4
    Sections_Q0 = Sections_A * Sections_U  # initial flow rate
    # Sections_ZC = np.zeros(len(Sections),dtype=complex)
    Sections_S = np.zeros(len(Sections))
    Data = np.column_stack((Sections_L, Sections_D, Sections_A, Sections_a, Sections_f, Sections_U, Sections_Q0, Sections_S))
    PipePD = pd.DataFrame(data=Data, index=Sections, columns=["L", "D", "A", "a", "f", "U", "Q0", "S"])
    df = 0.0025
    MaxF = 100
    Freq_range = np.arange(0, MaxF + df, df)
    # Critical component location (distance from upstream)
    BaseSensor1_Loc = (0,-1)
    BaseSensor2_Loc = (4,0)
    SimulationSizeBYDSensors = 1000
    SimulationSizeBTWSensors = 1000
    SourcesBYDSensor1 = np.random.uniform(0, 250, SimulationSizeBYDSensors)
    SourcesBYDSensor2 = np.random.uniform(850, 1100, SimulationSizeBYDSensors)
    SourcesBTWSensors = np.random.uniform(250, 850, SimulationSizeBTWSensors)
    CombinedSourceLocation = np.concatenate((SourcesBYDSensor1, SourcesBYDSensor2, SourcesBTWSensors), axis=None)
    global time
    time = np.arange(0, (1 / df)  + (1 / MaxF), 1 / MaxF)
    Sensor1Superpositions = np.zeros(np.shape(time))
    Sensor2Superpositions = np.zeros(np.shape(time))
    CoreCount = 6
    chunksize = 10
    with mp.Pool(processes=CoreCount, initializer=init_worker, initargs=((PipePD.copy(deep=True), Freq_range, time, BaseSensor1_Loc, BaseSensor2_Loc, Sensor1Superpositions, Sensor2Superpositions),)) as pool:
        for result in tqdm(pool.imap(worker, CombinedSourceLocation, chunksize=chunksize), total=len(CombinedSourceLocation)):
            Sensor1Superpositions ,Sensor2Superpositions, Sensor1Head, Sensor2Head = result
    # for SourceLoc in CombinedSourceLocation:
    #     NewPipePD, Sensor1_Loc, Sensor2_Loc = insertSource(SourceLoc, PipePD.copy(deep=True), BaseSensor1_Loc, BaseSensor2_Loc)

    print("--- %s seconds ---" % (TimeMod.time() - start_time))
    CrossCorrelation = np.correlate(Sensor1Superpositions, Sensor2Superpositions, mode="same")
    Sensor1AutoCorrelation = np.correlate(Sensor1Superpositions, Sensor1Superpositions, mode="same")
    Sensor2AutoCorrelation = np.correlate(Sensor2Superpositions, Sensor2Superpositions, mode="same")
    plt.figure("CrossCorrelation")
    plt.plot(time, CrossCorrelation)
    plt.figure("AutoCorrelation")
    plt.plot(time, Sensor1AutoCorrelation, label="Sensor1")
    plt.plot(time, Sensor2AutoCorrelation, label="Sensor2")
    plt.legend()
    plt.figure("TimeDomainImpulseResponse")
    Sensor1Head_time = np.real(np.fft.ifft(Sensor1Head, len(Sensor1Head)))
    Sensor2Head_time = np.real(np.fft.ifft(Sensor2Head, len(Sensor2Head)))
    plt.plot(time, np.real(np.fft.ifft(Sensor1Head, len(Sensor1Head))), label="Sensor1")
    plt.plot(time, np.real(np.fft.ifft(Sensor2Head, len(Sensor2Head))), label="Sensor2")
    SaveData = np.column_stack((time, Sensor1Head_time, Sensor2Head_time, CrossCorrelation, Sensor1AutoCorrelation, Sensor2AutoCorrelation))
    SaveDict = {}
    SaveDict["Result"] = SaveData.tolist()
    pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    pyexcel.free_resources()
    print("File Saved")
    plt.show()

if __name__ == '__main__':
    main()
