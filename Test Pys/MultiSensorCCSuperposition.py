from multiprocessing.pool import ThreadPool
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

"""Prompt user to select two csv files"""
file1 = filedialog.askopenfilename()
file2 = filedialog.askopenfilename()
"""Read the csv files into numpy arrays"""
data1 = np.genfromtxt(file1, delimiter=',')
data2 = np.genfromtxt(file2, delimiter=',')
T1, S2, S3 = data1[:, 0], data1[:, 1], data1[:, 2]
T2, S1, S4 = data2[:, 0], data2[:, 1], data2[:, 2]
InputArray = np.array([S1, S2, S3, S4])
InputLabel = np.array([0, 1, 2, 3])
DistancetoSensor0 = np.array([0, 500, 1100, 1400])
dt = T1[1] - T1[0]
wavespeed = 1000


def pairs(InputLabel):
    """Generate all possible pairs of sensors"""
    for i in range(len(InputLabel)):
        for j in range(i + 1, len(InputLabel)):
            yield InputLabel[i], InputLabel[j], DistancetoSensor0[i] - DistancetoSensor0[0], DistancetoSensor0[j] - \
                  DistancetoSensor0[i]


def init_worker(SharedData):
    """Initialize worker"""
    global InputArray
    InputArray = SharedData


def worker(pair):
    """Compute cross-correlation"""
    i, j, disToSensor0, distanceBetweenSensor = pair
    print("Computing cross-correlation between sensor {} and sensor {}".format(i, j))
    corr = np.correlate(InputArray[i], InputArray[j], mode='full')
    time = np.linspace(-T1[-1], T1[-1], len(corr))
    return corr, time, i, j


def plotsignal(InputArray):
    """Plot signals"""
    for i in range(len(InputArray)):
        plt.figure(figsize=(6, 3))  # Set the figure size to be small and suitable for a flowchart
        plt.plot(T1, InputArray[i], linewidth=2)  # Set the line width to make it clearly visible
        plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False,
                        labelleft=False)  # Remove all ticks and tick labels
        """limit x axis to 20s """
        for spine in plt.gca().spines.values():
            spine.set_visible(False)  # Remove all spines
        plt.xlim(0, 20)
        plt.gca().autoscale(enable=True, axis='y', tight=True)  # Set the y-axis to be automatically scaled
        plt.xticks(fontsize=8)  # Set x-axis tick labels font size
        plt.yticks(fontsize=8)  # Set y-axis tick labels font size
        if i == 2:
            plt.title("Sensor {}".format("..."), fontsize=30)
        elif i == 3:
            plt.title("Sensor {}".format("n"), fontsize=30)
        else:
            plt.title("Sensor {}".format(i + 1), fontsize=30)
        plt.xlabel("", fontsize=10)
        plt.ylabel("", fontsize=10)
        plt.tight_layout()  # Adjust the layout to minimize whitespace


def computeCrossCorrelation(InputArray):
    """apply multi-threading using map.pool"""
    CoreCount = 6
    ResultDict = {}
    with ThreadPool(processes=CoreCount, initializer=init_worker, initargs=(InputArray,)) as pool:
        for result in tqdm(pool.map(worker, pairs(InputLabel)), total=len(list(pairs(InputLabel)))):
            ResultDict[result[2], result[3]] = result[0], result[1]
    ResultDict = processCCs(ResultDict)
    ResultDict = shiftOrigin(ResultDict)
    return ResultDict


def processCCs(ResultDict):
    """Process cross-correlation"""
    for pair in pairs(InputLabel):
        i, j, disToSensor0, distanceBetweenSensor = pair
        time = ResultDict[i, j][1]
        distance = time * wavespeed / 2
        ResultDict[i, j] = ResultDict[i, j][0], distance
        print(distance)
    return ResultDict


def shiftOrigin(ResultDict):
    """Shift origin to sensor 0"""
    for pair in pairs(InputLabel):
        i, j, disToSensor0, distanceBetweenSensor = pair
        CC = ResultDict[i, j][0]
        distance = ResultDict[i, j][1]
        dDistance = distance[1] - distance[0]
        index = int((distanceBetweenSensor / 2) / dDistance)
        CC = np.roll(CC, index)
        index = int((disToSensor0 / dDistance))
        CC = np.roll(CC, index)
        ResultDict[i, j] = CC, distance
    return ResultDict


def plotCrossCorrelation(ResultDict):
    """Plot cross-correlations"""
    for pair in pairs(InputLabel):
        i, j, disToSensor0, distanceBetweenSensor = pair
        plt.figure(figsize=(4, 3))  # Set the figure size to be small and suitable for a flowchart
        plt.plot(ResultDict[i, j][1], ResultDict[i, j][0], linewidth=2)  # Set the line width to make it clearly visible
        plt.tick_params(axis='both', which='both', length=0)  # Remove tick marks
        plt.xticks(fontsize=8)  # Set x-axis tick labels font size
        plt.yticks(fontsize=8)  # Set y-axis tick labels font size
        plt.title("Sensor {} & {}".format(i, j), fontsize=12)
        plt.xlabel("Distance (m)", fontsize=10)
        plt.ylabel("Correlation", fontsize=10)
        plt.tight_layout()  # Adjust the layout to minimize whitespace


def superpositionAllCorrelations(ResultDict):
    """Superposition all cross-correlations"""
    for pair in pairs(InputLabel):
        i, j, disToSensor0, distanceBetweenSensor = pair
        if i == 0:
            CC = ResultDict[i, j][0]
            distance = ResultDict[i, j][1]
        else:
            CC += ResultDict[i, j][0]

    plt.figure(figsize=(4, 3))  # Set the figure size to be small and suitable for a flowchart
    plt.plot(distance, CC, linewidth=2)  # Set the line width to make it clearly visible
    plt.tick_params(axis='both', which='both', length=0)  # Remove tick marks
    plt.xticks(fontsize=8)  # Set x-axis tick labels font size
    plt.yticks(fontsize=8)  # Set y-axis tick labels font size
    plt.title("Processed", fontsize=12)
    plt.xlabel("Distance (m)", fontsize=10)
    plt.ylabel("Correlation", fontsize=10)
    plt.tight_layout()  # Adjust the layout to minimize whitespace

    plt.show()


if __name__ == '__main__':
    plotsignal(InputArray)
    plt.show()
    ResultDict = computeCrossCorrelation(InputArray)
    plotCrossCorrelation(ResultDict)
    superpositionAllCorrelations(ResultDict)
