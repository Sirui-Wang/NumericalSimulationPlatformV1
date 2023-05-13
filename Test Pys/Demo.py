from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np

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


# Demo for correlation function superposition technique (post correlation processing)
def plotsignal(InputArray):
    """Plot signals"""
    for i in range(len(InputArray)):
        plt.figure(figsize=(9, 3))  # Set the figure size to be small and suitable for a flowchart
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


if __name__ == '__main__':
    # plotsignal(InputArray)
    """Compute and plot one example correlation function"""
    corr = np.correlate(InputArray[0], InputArray[2], mode='full')
    time = np.linspace(-T1[-1], T1[-1], len(corr))
    time = time * wavespeed / 2 + 550 + 0
    plt.figure(1, figsize=(9, 3))
    plt.plot(time, corr)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)  # Remove all spines
    plt.title("Correlation function", fontsize=20)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False,
                    labelleft=False)  # Remove all ticks and tick labels
    plt.xlim(-1000, 1000)
    # plt.xlim(-1.5, 1.5)
    plt.xticks(fontsize=12)
    # plt.xlabel("Time (s)", fontsize=20)
    # tick_positions, tick_labels = plt.xticks()
    # new_tick_labels = [f"{label.get_text()} × α" for label in tick_labels]
    # plt.xticks(tick_positions, new_tick_labels)
    plt.xlabel("Distance (m)", fontsize=20)
    """Compute and plot auto correlation function"""
    autocorr = np.correlate(InputArray[1], InputArray[1], mode='full')
    time = np.linspace(-T1[-1], T1[-1], len(autocorr))
    time = time * wavespeed / 2 + 500 + 0
    plt.figure(2, figsize=(9, 3))
    plt.plot(time, autocorr)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.title("Auto correlation function", fontsize=20)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False,
                    labelleft=False)
    plt.xlim(-1000, 1000)
    plt.xticks(fontsize=12)
    plt.xlabel("Distance (m)", fontsize=20)
    """Compute and plot auto correlation function"""
    autocorr = np.correlate(InputArray[0], InputArray[0], mode='full')
    time = np.linspace(-T1[-1], T1[-1], len(autocorr))
    time = time * wavespeed / 2 + 0 + 0
    plt.figure(3, figsize=(9, 3))
    plt.plot(time, autocorr)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.title("Auto correlation function", fontsize=20)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False,
                    labelleft=False)
    plt.xlim(-1000, 1000)
    plt.xticks(fontsize=12)
    plt.xlabel("Distance (m)", fontsize=20)
    plt.show()
