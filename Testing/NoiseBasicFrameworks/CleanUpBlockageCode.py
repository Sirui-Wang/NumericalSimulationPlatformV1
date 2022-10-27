import copy
import multiprocessing as mp
import os
import time
from tkinter import filedialog

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import pyexcel
from tqdm import tqdm

S = [[1, 0, 0],
     [0, 1, 1],
     [0, 0, 1]]

isinifiniteBC = True


def definePipe():
    Sections = np.array([0, 1, 2, 3, 4])
    Sources = np.array(["A", "B", "C", "D", "E"])
    Targets = np.array(["B", "C", "D", "E", "F"])
    Sections_L = np.array([250, 250, 100, 250, 250])
    Sections_a = np.full(len(Sections), 1000)
    # Sections_D = np.array([0.3, 0.3, 0.2, 0.3, 0.3]) # OG
    # Sections_U = np.array([1.41, 1.41, 3.18, 1.41, 1.41]) # OG
    # Sections_f = np.array([0.038, 0.038, 0.041, 0.038, 0.038]) # OG
    """0.9R"""
    # Sections_D = np.array([0.435, 0.435, 0.1, 0.435, 0.435])
    # Sections_f = np.array([0.035, 0.035, 0.057, 0.035, 0.035])
    # Sections_U = np.array([0.13, 0.13, 2.55, 0.13, 0.13])
    """0.85R"""
    # Sections_D = np.array([0.351, 0.351, 0.1, 0.351, 0.351])
    # Sections_f = np.array([0.037, 0.037, 0.057, 0.037, 0.037])
    # Sections_U = np.array([0.21, 0.21, 2.55, 0.21, 0.21])
    """0.8R"""
    # Sections_D = np.array([0.3, 0.3, 0.1, 0.3, 0.3])
    # Sections_f = np.array([0.039, 0.039, 0.057, 0.039, 0.039])
    # Sections_U = np.array([0.28, 0.28, 2.55, 0.28, 0.28])
    """0.7R"""
    # Sections_D = np.array([0.238, 0.238, 0.1, 0.238, 0.238])
    # Sections_f = np.array([0.042, 0.042, 0.057, 0.042, 0.042])
    # Sections_U = np.array([0.45, 0.2451, 2.55, 0.45, 0.45])
    """0.5R"""
    Sections_D = np.array([0.173, 0.173, 0.1, 0.173, 0.173])
    # Sections_f = np.array([0.047, 0.047, 0.057, 0.047, 0.047])
    Sections_f = np.array([0.047, 0.047, 0.057, 0.047, 0.047])
    Sections_U = np.array([0.85, 0.85, 2.55, 0.85, 0.85])
    Sections_A = (np.pi * Sections_D ** 2) / 4
    Sections_Q0 = Sections_A * Sections_U
    Sections_S = np.zeros(len(Sections))
    Data = np.column_stack(
        (Sections_L, Sections_D, Sections_A, Sections_a, Sections_f, Sections_U, Sections_Q0, Sections_S, Sources,
         Targets))
    PipePD = pd.DataFrame(data=Data, index=Sections,
                          columns=["L", "D", "A", "a", "f", "U", "Q0", "isSource", "sources", "targets"])
    G = nx.from_pandas_edgelist(PipePD, source="sources", target="targets", edge_attr=True, create_using=nx.DiGraph())
    return G


def SpliteEdge(G, SelectedEdge):
    source, target = SelectedEdge
    SelectedEdgeProperty = G.edges[SelectedEdge]
    G.remove_edge(source, target)
    G.add_node("Pert")
    G.add_edge(source, "Pert",
               L=SelectedEdgeProperty["L"],
               D=SelectedEdgeProperty["D"],
               A=SelectedEdgeProperty["A"],
               a=SelectedEdgeProperty["a"],
               f=SelectedEdgeProperty["f"],
               U=SelectedEdgeProperty["U"],
               Q0=SelectedEdgeProperty["Q0"],
               isSource=1, )
    G.add_edge("Pert", target,
               L=SelectedEdgeProperty["L"],
               D=SelectedEdgeProperty["D"],
               A=SelectedEdgeProperty["A"],
               a=SelectedEdgeProperty["a"],
               f=SelectedEdgeProperty["f"],
               U=SelectedEdgeProperty["U"],
               Q0=SelectedEdgeProperty["Q0"],
               isSource=0)
    return G


def PossiblePerts(G):
    GraphDict = {}
    for possibleSecnario in list(G.edges):
        CopyG = copy.deepcopy(G)
        CopyG = SpliteEdge(CopyG, possibleSecnario)
        GraphDict[possibleSecnario] = (CopyG, G.edges[possibleSecnario]["L"])
    return GraphDict


def editPipeLength(CopyG, Key, SourceLoc):
    source, target = Key
    PipeLength = float(CopyG.edges[(source, "Pert")]["L"])
    LengthBF = SourceLoc
    LengthAFT = PipeLength - SourceLoc
    CopyG.edges[(source, "Pert")]["L"] = LengthBF
    CopyG.edges[("Pert", target)]["L"] = LengthAFT
    if LengthAFT < 0:
        print("Wrong")
    return CopyG


def round_nearest2(x, a):
    rounded_x = round(x / a) * a
    return rounded_x


def drawGraph(G, Title):
    plt.figure(Title)
    edge_labels = {i[0:2]: '{}m'.format(i[2]['L']) for i in G.edges(data=True)}
    pos = nx.spring_layout(G, weight="L")
    nx.draw(G, pos, with_labels=True)
    nx.draw_networkx_edge_labels(G, pos, edge_labels)


def drawDistribution(DistArray, Title, bins=10):
    plt.figure(Title)
    plt.hist(DistArray, bins=bins)


def plotImpulseResponse(time, Sensor1, Sensor2, Title):
    plt.figure(Title)
    plt.plot(time, Sensor1, label="Sensor1")
    plt.plot(time, Sensor2, label="Sensor2")
    plt.legend()


def plotCorrelation(time, Sensor1, Sensor2, Title):
    # CCIndexese = np.argwhere(abs(CrossCorrelation)>0.08)
    ZeroPeak1 = np.argwhere(abs(Sensor1) < 0.08)
    ZeroPeak2 = np.argwhere(abs(Sensor2) < 0.08)
    Sensor1[ZeroPeak1] = 0
    Sensor2[ZeroPeak2] = 0
    ZeroPeak1 = np.argwhere(abs(Sensor1) < 0.08)
    CrossCorrelation = np.correlate(Sensor1, Sensor2, mode="same")
    CCIndexese = np.argwhere(abs(CrossCorrelation) > 0.0025)
    for i in CCIndexese:
        print(time[i] - 20, CrossCorrelation[i])
    print("end-------------------------------")
    plt.figure(Title)
    plt.plot(time, CrossCorrelation)


def FieldMatrix(EdgeData, freq):
    g = 9.81
    n = 2  # empirical value for TM friction term R
    omega = 2 * np.pi * freq  # T = Theoretical period
    L, D, A, a, f, U, Q0, S = [float(i) for i in EdgeData.values()]
    R = (n * f * (Q0 ** (n - 1))) / (2 * g * D * (A ** n))
    "Field Matrix for single conduit"
    mu = np.sqrt(-((omega ** 2) / (a ** 2)) + ((1j * g * A * omega * R) / (a ** 2)))
    Zc = (mu * a ** 2) / (1j * omega * g * A)
    F = np.array([[np.cosh(mu * L), (-1 / Zc) * np.sinh(mu * L), 0],
                  [-Zc * np.sinh(mu * L), np.cosh(mu * L), 0],
                  [0, 0, 1]])
    return F, Zc


def analysis(G, Freq_range, Sensor1, Sensor2):
    global S, isinifiniteBC
    Sensor1_Head_fResponse = np.zeros((len(Freq_range)), dtype=complex)
    Sensor2_Head_fResponse = np.zeros((len(Freq_range)), dtype=complex)
    AnalysisPath = nx.dag_longest_path(G)
    Impedance = np.zeros(len(AnalysisPath) - 1, dtype=complex)
    for i in (range(1, len(Freq_range))):
        U = np.identity(3, dtype=complex)
        for j in range(len(AnalysisPath) - 1):
            source = AnalysisPath[-j - 2]
            target = AnalysisPath[-j - 1]
            F, Zc = FieldMatrix(G.edges[source, target], Freq_range[i])
            Impedance[j] = Zc
            if G.edges[source, target]["isSource"] == 1:
                U = F @ S @ U
            else:
                U = F @ U
            if source == Sensor1:
                Sensor1F = U
            elif source == Sensor2:
                Sensor2F = U
        c1, c2, c3, c4, c5, c6, c7, c8, c9 = U.reshape(9)
        if isinifiniteBC == True:
            Zcup = Impedance[0]
            Zcdown = Impedance[-1]
            EQA = np.array([[-1, c1, 0, c2],
                            [0, c4, -1, c5],
                            [Zcup, 0, -1, 0],
                            [0, -Zcdown, 0, -1]])
            EQB = np.array([[-c3],
                            [-c6],
                            [0],
                            [0]])
            Solution = np.linalg.solve(EQA, EQB)
            Q, q, H, h = Solution
        else:
            EQA = [[-1, c2],
                   [0, c5]]
            EQB = [[-c3],
                   [-c6]]
            Solution = np.linalg.solve(EQA, EQB)
            Q, h = Solution
            H, q = [[0], [0]]
        Sensor1Result = Sensor1F @ [q, h, [1]]
        Sensor2Result = Sensor2F @ [q, h, [1]]
        Sensor1_Head_fResponse[i] = Sensor1Result[1][0]
        Sensor2_Head_fResponse[i] = Sensor2Result[1][0]
    return Sensor1_Head_fResponse, Sensor2_Head_fResponse


def init_worker(shared_data):
    global GraphDict, GridSize, SimulationSizePerPipe, Freq_range, Sensor1, Sensor2, timeArray, timeArray1
    GraphDict, GridSize, SimulationSizePerPipe, Freq_range, Sensor1, Sensor2, timeArray, timeArray1 = shared_data


def worker(key_index):
    global GraphDict, GridSize, SimulationSizePerPipe, Freq_range, Sensor1, Sensor2, timeArray, timeArray1
    key = list(GraphDict.keys())[key_index]
    Graph, PipeLength = GraphDict[key]
    if key == ("C", "D"):
        np.random.seed(key_index)
        SimulationSizePerPipe = max(int(0.4 * SimulationSizePerPipe), 1)
        PertLocations = np.random.uniform(0, int(PipeLength), SimulationSizePerPipe)
    elif key == ("A", "B"):
        # np.random.seed(key_index)
        np.random.seed(11)
        SimulationSizePerPipe = SimulationSizePerPipe
        PertLocations = np.random.uniform(0, int(PipeLength), SimulationSizePerPipe)
    elif key == ("E", "F"):
        # np.random.seed(key_index)
        np.random.seed(10)
        SimulationSizePerPipe = SimulationSizePerPipe
        PertLocations = np.random.uniform(0, int(PipeLength), SimulationSizePerPipe)
    elif key == ("B", "C"):
        # np.random.seed(key_index)
        SimulationSizePerPipe = SimulationSizePerPipe
        PertLocations = np.random.uniform(0, int(PipeLength), SimulationSizePerPipe)
    elif key == ("D", "E"):
        # np.random.seed(key_index)
        SimulationSizePerPipe = SimulationSizePerPipe
        PertLocations = np.random.uniform(0, int(PipeLength), SimulationSizePerPipe)
    SumSensor1 = np.zeros(len(timeArray1))
    SumSensor2 = np.zeros(len(timeArray1))
    RoundedPertLocations = np.zeros(len(PertLocations))
    index = 0
    for simulation in tqdm(PertLocations):
        simulation = round_nearest2(simulation, GridSize)
        RoundedPertLocations[index] = simulation
        index += 1
        SplitedG = editPipeLength(copy.deepcopy(Graph), key, simulation)
        Sensor1FreqHead, Sensor2FreqHead = analysis(SplitedG, Freq_range, Sensor1, Sensor2)
        Noise = np.random.normal(0, 0.1, len(Sensor1FreqHead))
        NPower = np.sum((abs(Noise)) ** 2) / len(Noise)
        Noise = Noise / np.sqrt(NPower)
        Sensor1Time = np.real(np.fft.ifft(Sensor1FreqHead, len(Sensor1FreqHead)))
        Sensor2Time = np.real(np.fft.ifft(Sensor2FreqHead, len(Sensor2FreqHead)))
        Sensor1TimeHeadWNoise = np.convolve(np.real(np.fft.ifft(Sensor1FreqHead, len(Sensor1FreqHead))), Noise,
                                            mode="full")
        Sensor2TimeHeadWNoise = np.convolve(np.real(np.fft.ifft(Sensor2FreqHead, len(Sensor2FreqHead))), Noise,
                                            mode="full")
        plotImpulseResponse(timeArray, np.real(np.fft.ifft(Sensor1FreqHead, len(Sensor1FreqHead))),
                            np.real(np.fft.ifft(Sensor2FreqHead, len(Sensor2FreqHead))),
                            "ImpulseResponse {}, {}".format(key, simulation))
        # Sensor1ImpulsesIndexes = np.argwhere(abs(Sensor1Time)>0.02)
        # print("Sensor1", "\n","Head", Sensor1Time[Sensor1ImpulsesIndexes], "\n", "time", timeArray[Sensor1ImpulsesIndexes])
        # Sensor2ImpulsesIndexes = np.argwhere(abs(Sensor2Time)>0.02)
        # print("Sensor2", "\n", "Head", Sensor2Time[Sensor2ImpulsesIndexes], "\n", "time",timeArray[Sensor2ImpulsesIndexes])
        # drawGraph(SplitedG, "Graph {}, {}".format(key, PertLocation))
        # plotCorrelation(timeArray1, Sensor1TimeHeadWNoise, Sensor2TimeHeadWNoise,"CrossCorrelation {}, {}".format(key, simulation))
        plotCorrelation(timeArray, Sensor1Time, Sensor2Time,
                        "CrossCorrelation without noise {}, {}".format(key, simulation))
        plt.show()
        S1Power = np.sum((abs(Sensor1TimeHeadWNoise)) ** 2) / len(Sensor1TimeHeadWNoise)
        S2Power = np.sum((abs(Sensor2TimeHeadWNoise)) ** 2) / len(Sensor2TimeHeadWNoise)
        Sensor1TimeHeadWNoise = Sensor1TimeHeadWNoise / np.sqrt(S1Power)
        Sensor2TimeHeadWNoise = Sensor2TimeHeadWNoise / np.sqrt(S2Power)
        SumSensor1 = np.add(SumSensor1, Sensor1TimeHeadWNoise)
        SumSensor2 = np.add(SumSensor2, Sensor2TimeHeadWNoise)
    return SumSensor1, SumSensor2, RoundedPertLocations, key


def main():
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname
                                                         + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    start_time = time.time()
    G = definePipe()
    df = 0.025
    MaxF = 100
    Freq_range = np.arange(0, MaxF + df, df)
    Sensor1 = "B"
    Sensor2 = "E"
    GraphDict = PossiblePerts(G)
    timeArray = np.arange(0, (1 / df) + (1 / MaxF), 1 / MaxF)
    timeArray1 = np.arange(0, (1 / df) * 2 + (1 / MaxF), 1 / MaxF)
    GridSize = 1000 / MaxF
    SimulationSizePerPipe = 1
    SumSensor1 = np.zeros(len(timeArray1))
    SumSensor2 = np.zeros(len(timeArray1))
    # CoreCount = len(GraphDict.keys())
    CoreCount = 1
    SourceLocationDist = {}
    with mp.Pool(processes=CoreCount, initializer=init_worker, initargs=(
    (GraphDict, GridSize, SimulationSizePerPipe, Freq_range, Sensor1, Sensor2, timeArray, timeArray1),)) as pool:
        for result in pool.map(worker, range(len(list(GraphDict.keys())))):
            Sensor1Result, Sensor2Result, SourceDist, key = result
            SumSensor1 = np.add(SumSensor1, Sensor1Result)
            SumSensor2 = np.add(SumSensor2, Sensor2Result)
            SourceLocationDist[key] = SourceDist
    """ Single Process Code Start """
    # for Key, item in tqdm(GraphDict.items()):
    #     Graph, PipeLength = item
    #     # if Key == ("C","D"):
    #     #     SimulationSizePerPipe = 100
    #     # elif Key == ("A", "B"):
    #     #     SimulationSizePerPipe = 100
    #     # elif Key == ("E", "F"):
    #     #     SimulationSizePerPipe = 100
    #     # else:
    #     #     SimulationSizePerPipe = 100
    #     PertLocations = np.random.uniform(GridSize, float(PipeLength)-GridSize, SimulationSizePerPipe)
    #     # PertLocations=[180]
    #     # drawDistribution(PertLocations, str(Key))
    #     for simulation in PertLocations:
    #         G = copy.deepcopy(Graph)
    #         PertLocation = round_nearest2(simulation, GridSize)
    #         SplitedG = editPipeLength(G, Key, PertLocation)
    #         # drawGraph(SplitedG, "{}, {}".format(Key, PertLocation))
    #         Sensor1FreqHead, Sensor2FreqHead = analysis(SplitedG, Freq_range, Sensor1, Sensor2)
    #         Noise = np.random.normal(0, 0.1, np.shape(Sensor1FreqHead))
    #         NPower = np.sum(abs((Noise)) ** 2) / len(Noise)
    #         Noise = Noise / NPower
    #         Sensor1TimeHeadWNoise = np.convolve(np.real(np.fft.ifft(Sensor1FreqHead, len(Sensor1FreqHead))), Noise,
    #                                             mode="full")
    #         Sensor2TimeHeadWNoise = np.convolve(np.real(np.fft.ifft(Sensor2FreqHead, len(Sensor2FreqHead))), Noise,
    #                                             mode="full")
    #         SumSensor1 = np.add(SumSensor1, Sensor1TimeHeadWNoise)
    #         SumSensor2 = np.add(SumSensor2, Sensor2TimeHeadWNoise)
    #         plotImpulseResponse(timeArray, np.real(np.fft.ifft(Sensor1FreqHead, len(Sensor1FreqHead))),np.real(np.fft.ifft(Sensor2FreqHead, len(Sensor2FreqHead))),"ImpulseResponse {}, {}".format(Key, PertLocation))
    #         drawGraph(SplitedG,"Graph {}, {}".format(Key, PertLocation))
    #         plotCorrelation(timeArray1, Sensor1TimeHeadWNoise, Sensor2TimeHeadWNoise, "CrossCorrelation {}, {}".format(Key, PertLocation))
    """ Single Process Code End """
    # print(np.sum((abs(SumSensor1)) ** 2) / len(SumSensor1), np.sum((abs(SumSensor2)) ** 2) / len(SumSensor2))
    CrossCorrelation = np.correlate(SumSensor1, SumSensor2, mode="same")
    Sensor1AutoCorrelation = np.correlate(SumSensor1, SumSensor1, mode="same")
    Sensor2AutoCorrelation = np.correlate(SumSensor2, SumSensor2, mode="same")
    AnalysisPath = nx.dag_longest_path(G)
    CummulativeLength = 0
    TotalDist = np.array([])
    for j in range(len(AnalysisPath) - 1):
        source = AnalysisPath[j]
        target = AnalysisPath[j + 1]
        key = (source, target)
        TotalDist = np.concatenate((TotalDist, SourceLocationDist[key] + CummulativeLength))
        # print(TotalDist)
        CummulativeLength += int(G.edges[key]["L"])
    drawDistribution(TotalDist, "Total Distribution", int(CummulativeLength))
    # drawDistribution(TotalDist, "Total Distribution")
    plt.figure("CrossCorrelation")
    plt.plot(timeArray1, CrossCorrelation)
    plt.figure("AutoCorrelation")
    plt.plot(timeArray1, Sensor1AutoCorrelation, label="Sensor1{}".format(Sensor1))
    plt.plot(timeArray1, Sensor2AutoCorrelation, label="Sensor2{}".format(Sensor2))
    plt.legend()
    SaveData = np.column_stack((timeArray1, CrossCorrelation, Sensor1AutoCorrelation, Sensor2AutoCorrelation))
    SaveDict = {}
    SaveDict["Result"] = SaveData.tolist()
    pyexcel.isave_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    pyexcel.free_resources()
    print("File Saved")
    print("--------- %s seconds ---------" % (time.time() - start_time))


if __name__ == '__main__':
    main()
    plt.show()
