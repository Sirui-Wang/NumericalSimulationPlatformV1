import copy
import json
import os
from tkinter import filedialog, messagebox

import numpy as np
import pyexcel

import MOC.CalculateSteadyState as toSteadyState

# from tqdm import tqdm

# Auther: Sirui Wang
# Updated date: 13/July/2022
"""Module Comment"""


# Line Comment

def create_domain(dt, run_time):
    t_range = np.arange(-dt, run_time + dt, dt)  # define time axis for the computational square
    if len(t_range) % 2 == 1:
        t_range = np.linspace(-dt, run_time, int(run_time / dt) + 2)
    G.graph["time_seq"] = t_range  # append this time axis to the graph for future plotting
    for source, target, data in G.edges.data(data=True):
        dx = dt * data["wave_velocity"]  # specific to each pipe as it depends on the wave speed of each pipe
        size = int(np.round(data["length"] / dx) + 1)
        x_range = np.zeros(size)
        # create computational square
        U_mat = np.zeros((len(t_range), len(x_range)))
        H_mat = np.zeros((len(t_range), len(x_range)))
        """Append computational squares to each edge"""
        G.edges[source, target]["H_mat"] = H_mat
        G.edges[source, target]["U_mat"] = U_mat
        """Append spatial axis to each edge in case it is needed in the future"""
        G.edges[source, target]["dx"] = dx
        G.edges[source, target]["x_range"] = x_range


def ApplyExistingIC(ICfile):
    with open(ICfile, "r") as file:
        recovered_dict = json.load(file)
    for edge in G.edges:
        source, target = edge
        G.edges[edge]["H_mat"][-1] = np.array(recovered_dict[source + "," + target]["H"])
        G.edges[edge]["U_mat"][-1] = np.array(recovered_dict[source + "," + target]["U"])


def node_classification(G):
    for node in G.nodes():
        DoF = (len(G.in_edges(node)) + len(G.out_edges(node)))
        if DoF == 1:
            if len(G.out_edges(node)) == 1:
                classification = "Source"
            else:
                classification = "Sink"
        elif DoF == 2:
            classification = "SeriesConnection"
        elif DoF >= 3:
            classification = "BranchConnection"
        print(node, classification)
        G.nodes[node]["classification"] = classification
    return G


def network_decomposition(G):
    decomposed_network = dict.fromkeys(list(G.nodes.keys()), {})
    for node in list(decomposed_network.keys()):
        inedges = G.in_edges(node)
        outedges = G.out_edges(node)
        temp = {}
        for edge in inedges:
            temp[edge] = ["in"]
        for edge in outedges:
            temp[edge] = ["out"]
        decomposed_network[node] = temp
    return decomposed_network


def analysis(t, decomposed_network, G):
    global count
    for node, branch in decomposed_network.items():
        node_class = G.nodes[node]["classification"]
        if node_class == "Sink":
            NodeBCType = G.nodes[node]["BCType"]
            edge = list(branch.keys())[0]
            HA = G.edges[edge]["H_mat"][-t, -2]
            UA = G.edges[edge]["U_mat"][-t, -2]
            a = G.edges[edge]["wave_velocity"]
            f = G.edges[edge]["friction_factor"]
            r = (G.edges[edge]["diameter"]) / 2
            dx = G.edges[edge]["dx"]
            if G.edges[edge]["PertType"] == "":
                """ if there is no specified perturbation."""
                if NodeBCType == "Constant flow":
                    "Constant flow therefore UP is the same all the time, HP is calculated based on this UP and it can vary"
                    # Constant Flow boundary condition, known flow velocity, unknow head
                    UP = G.edges[edge]["U_mat"][-1, -1]
                    HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                    G.edges[edge]["H_mat"][-t - 1, -1] = HP
                    G.edges[edge]["U_mat"][-t - 1, -1] = UP
                else:
                    "Constant head therefore HP is the same all the time, UP is calculated based on this HP"
                    # Constant Flow boundary condition, known head, unknow flow velocity
                    HP = G.edges[edge]["H_mat"][-1, -1]
                    UP = UA + (9.81 / a) * (-HP + HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx)
                    G.edges[edge]["H_mat"][-t - 1, -1] = HP
                    G.edges[edge]["U_mat"][-t - 1, -1] = UP
            elif G.edges[edge]["PertType"] == "Impulse" and G.edges[edge]["PertLocation"] == node:
                """A impulse perturbation is when flow/head changes for 1 time step and resume to no perturbation status"""
                if G.graph["time_seq"][t] >= G.edges[edge]["Start Time"] and count == 0:
                    """to avoid python float inaccuracy, used if the start time is larger than run time to speicify when the perturbation is introduced
                    the count variable is to ensure that the impulse is only introduced one time."""
                    if NodeBCType == "Constant flow":
                        # if node is constant Flow BC, then the impulse is introduced as suddently set flow to 0
                        UP = 0
                        HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
                    else:
                        # if node is constant head BC, then the impulse is introduced as suddently set head to 0
                        HP = 0
                        UP = UA + (9.81 / a) * (-HP + HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
                    count += 1
                else:
                    """Before and after the impulse is introduced, this should act the same as if there are no perturbation in troduced"""
                    if NodeBCType == "Constant flow":
                        # Constant Flow boundary condition, known flow velocity, unknow head
                        UP = G.edges[edge]["U_mat"][-1, -1]
                        HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
                    else:
                        # Constant Flow boundary condition, known head, unknow flow velocity
                        HP = G.edges[edge]["H_mat"][-1, -1]
                        UP = UA + (9.81 / a) * (-HP + HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
            elif G.edges[edge]["PertType"] == "Controlled Flow" and G.edges[edge]["PertLocation"] == node:
                if G.graph["time_seq"][t] >= G.edges[edge]["Start Time"]:
                    """ if past start time, change flow velocity to specified flow velocity, regardless previous node boundary condition"""
                    # Constant Flow boundary condition, known flow velocity, unknow head
                    UP = G.edges[edge]["FlowRate"]
                    HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                    G.edges[edge]["H_mat"][-t - 1, -1] = HP
                    G.edges[edge]["U_mat"][-t - 1, -1] = UP
                else:
                    """if not past start time, act like if there is no perturbation introduced, calculated based on given node boundary condition"""
                    if NodeBCType == "Constant flow":
                        "Constant flow therefore UP is the same all the time, HP is calculated based on this UP and it can vary"
                        # Constant Flow boundary condition, known flow velocity, unknow head
                        UP = G.edges[edge]["U_mat"][-1, -1]
                        HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
                    else:
                        "Constant head therefore HP is the same all the time, UP is calculated based on this HP"
                        # Constant Flow boundary condition, known head, unknow flow velocity
                        HP = G.edges[edge]["H_mat"][-1, -1]
                        UP = UA + (9.81 / a) * (-HP + HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
            elif G.edges[edge]["PertType"] == "Full Closure" and G.edges[edge]["PertLocation"] == node:
                if G.graph["time_seq"][t] >= G.edges[edge]["Start Time"]:
                    """if past specified start time, apply 0 flow velocity condition, regardless what the previous condition is"""
                    # Constant Flow boundary condition, known flow velocity, unknow head
                    UP = 0
                    HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                    G.edges[edge]["H_mat"][-t - 1, -1] = HP
                    G.edges[edge]["U_mat"][-t - 1, -1] = UP
                else:
                    """if not past the specified start time, assume there are no perturbation introduced"""
                    if NodeBCType == "Constant flow":
                        "Constant flow therefore UP is the same all the time, HP is calculated based on this UP and it can vary"
                        # Constant Flow boundary condition, known flow velocity, unknow head
                        UP = G.edges[edge]["U_mat"][-1, -1]
                        HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
                    else:
                        "Constant head therefore HP is the same all the time, UP is calculated based on this HP"
                        # Constant Flow boundary condition, known head, unknow flow velocity
                        HP = G.edges[edge]["H_mat"][-1, -1]
                        UP = UA + (9.81 / a) * (-HP + HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, -1] = HP
                        G.edges[edge]["U_mat"][-t - 1, -1] = UP
            elif G.edges[edge]["PertType"] == "Sinusoidal Head" and G.edges[edge]["PertLocation"] == node:
                """When introducing a Sinusoidal perturbation, it need to start from the very beginning
                since a sinusoidal head perturbation is applied, head is known and flow is unknown
                head vary sinusoidally with a given frequency and amplitude, centered at the steady state"""
                SSReservoir = G.edges[edge]["H_mat"][-1, -1]
                omega = 2 * np.pi * G.edges[edge]["Freq"]
                amp = G.edges[edge]["Amp"]
                HP = SSReservoir + amp * np.sin(omega * t_range[t])
                UP = UA + (9.81 / a) * (-HP + HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx)
                G.edges[edge]["H_mat"][-t - 1, -1] = HP
                G.edges[edge]["U_mat"][-t - 1, -1] = UP
            elif G.edges[edge]["PertType"] == "Sinusoidal Flow" and G.edges[edge]["PertLocation"] == node:
                """Same with Sinusoidal head perturbation but just with flow"""
                SSFlowSpeed = G.edges[edge]["U_mat"][-1, -1]
                omega = 2 * np.pi * G.edges[edge]["Freq"]
                amp = G.edges[edge]["Amp"]
                UP = SSFlowSpeed + amp * np.sin(omega * t_range[t])
                HP = HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx - (a / 9.81) * (UP - UA)
                G.edges[edge]["H_mat"][-t - 1, -1] = HP
                G.edges[edge]["U_mat"][-t - 1, -1] = UP
        elif node_class == "Source":
            NodeBCType = G.nodes[node]["BCType"]
            edge = list(branch.keys())[0]
            HB = G.edges[edge]["H_mat"][-t, 1]
            UB = G.edges[edge]["U_mat"][-t, 1]
            a = G.edges[edge]["wave_velocity"]
            f = G.edges[edge]["friction_factor"]
            r = (G.edges[edge]["diameter"]) / 2
            dx = G.edges[edge]["dx"]
            if G.edges[edge]["PertType"] == "":
                if NodeBCType == "Constant flow":
                    UP = G.edges[edge]["U_mat"][-1, 0]
                    HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                    G.edges[edge]["H_mat"][-t - 1, 0] = HP
                    G.edges[edge]["U_mat"][-t - 1, 0] = UP
                else:
                    # Default is constant reservoir head BC
                    HP = G.edges[edge]["H_mat"][-1, 0]
                    UP = UB + (9.81 / a) * (HP - HB - ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx)
                    G.edges[edge]["H_mat"][-t - 1, 0] = HP
                    G.edges[edge]["U_mat"][-t - 1, 0] = UP
            elif G.edges[edge]["PertType"] == "Impulse" and G.edges[edge]["PertLocation"] == node:
                if G.graph["time_seq"][t] >= G.edges[edge]["Start Time"] and count == 0:
                    if NodeBCType == "Constant flow":
                        UP = 0
                        HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
                    else:
                        # Default is constant reservoir head BC
                        HP = 0
                        UP = UB + (9.81 / a) * (HP - HB - ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
                    count += 1
                else:
                    if NodeBCType == "Constant flow":
                        UP = G.edges[edge]["U_mat"][-1, 0]
                        HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
                    else:
                        # Default is constant reservoir head BC
                        HP = G.edges[edge]["H_mat"][-1, 0]
                        UP = UB + (9.81 / a) * (HP - HB - ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
            elif G.edges[edge]["PertType"] == "Controlled Flow" and G.edges[edge]["PertLocation"] == node:
                if G.graph["time_seq"][t] >= G.edges[edge]["Start Time"]:
                    # Constant Flow boundary condition, known flow velocity, unknow head
                    UP = G.edges[edge]["FlowRate"]
                    HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                    G.edges[edge]["H_mat"][-t - 1, 0] = HP
                    G.edges[edge]["U_mat"][-t - 1, 0] = UP
                else:
                    if NodeBCType == "Constant flow":
                        UP = G.edges[edge]["U_mat"][-1, 0]
                        HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
                    else:
                        # Default is constant reservoir head BC
                        HP = G.edges[edge]["H_mat"][-1, 0]
                        UP = UB + (9.81 / a) * (HP - HB - ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
            elif G.edges[edge]["PertType"] == "Full Closure" and G.edges[edge]["PertLocation"] == node:
                if G.graph["time_seq"][t] >= G.edges[edge]["Start Time"]:
                    # Constant Flow boundary condition, known flow velocity, unknow head
                    UP = 0
                    HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                    G.edges[edge]["H_mat"][-t - 1, 0] = HP
                    G.edges[edge]["U_mat"][-t - 1, 0] = UP
                else:
                    if NodeBCType == "Constant flow":
                        UP = G.edges[edge]["U_mat"][-1, 0]
                        HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
                    else:
                        # Default is constant reservoir head BC
                        HP = G.edges[edge]["H_mat"][-1, 0]
                        UP = UB + (9.81 / a) * (HP - HB - ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx)
                        G.edges[edge]["H_mat"][-t - 1, 0] = HP
                        G.edges[edge]["U_mat"][-t - 1, 0] = UP
            elif G.edges[edge]["PertType"] == "Sinusoidal Head" and G.edges[edge]["PertLocation"] == node:
                SSReservoir = G.edges[edge]["H_mat"][-1, 0]
                omega = 2 * np.pi * G.edges[edge]["Freq"]
                amp = G.edges[edge]["Amp"]
                HP = SSReservoir + amp * np.sin(omega * t_range[t])
                UP = UB + (9.81 / a) * (HP - HB - ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx)
                G.edges[edge]["H_mat"][-t - 1, 0] = HP
                G.edges[edge]["U_mat"][-t - 1, 0] = UP
            elif G.edges[edge]["PertType"] == "Sinusoidal Flow" and G.edges[edge]["PertLocation"] == node:
                SSFlowSpeed = G.edges[edge]["U_mat"][-1, 0]
                omega = 2 * np.pi * G.edges[edge]["Freq"]
                amp = G.edges[edge]["Amp"]
                UP = SSFlowSpeed + amp * np.sin(omega * t_range[t])
                HP = (a / 9.81) * (UP - UB) + HB + ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                G.edges[edge]["H_mat"][-t - 1, 0] = HP
                G.edges[edge]["U_mat"][-t - 1, 0] = UP
        else:
            # doesnt matter what is inputed as BC type, if it is not classified as a edge node, ignore BCtype and proceed with the assumption of it is not a BC node
            size = len(branch) + 1
            A = np.zeros((size, size))
            B = np.zeros((size, 1))
            i = 1
            mapping = dict.fromkeys(list(branch.keys()))
            for edge, direction in branch.items():
                HA = G.edges[edge]["H_mat"][-t, -2]
                UA = G.edges[edge]["U_mat"][-t, -2]
                HB = G.edges[edge]["H_mat"][-t, 1]
                UB = G.edges[edge]["U_mat"][-t, 1]
                a = G.edges[edge]["wave_velocity"]
                f = G.edges[edge]["friction_factor"]
                r = (G.edges[edge]["diameter"]) / 2
                dx = G.edges[edge]["dx"]
                Area = np.pi * r ** 2
                if direction == ["in"]:
                    A[0, i] = Area
                    A[i, 0] = 1
                    A[i, i] = a / 9.81
                    B[i, 0] = (a / 9.81) * UA + HA - ((f * UA * abs(UA)) / (4 * 9.81 * r)) * dx
                    mapping[edge] = i
                    i += 1
                else:
                    A[0, i] = -Area
                    A[i, 0] = -1
                    A[i, i] = a / 9.81
                    B[i, 0] = (a / 9.81) * UB - HB - ((f * UB * abs(UB)) / (4 * 9.81 * r)) * dx
                    mapping[edge] = i
                    i += 1
            C = np.linalg.solve(A, B)
            for edge, direction in branch.items():
                if direction == ["in"]:
                    G.edges[edge]["H_mat"][-t - 1, -1] = C[0][0]
                    G.edges[edge]["U_mat"][-t - 1, -1] = C[mapping[edge]][0]
                else:
                    G.edges[edge]["H_mat"][-t - 1, 0] = C[0][0]
                    G.edges[edge]["U_mat"][-t - 1, 0] = C[mapping[edge]][0]
    for edge in list(G.edges):
        source, target = edge
        for x in np.arange(1, len(G.edges[edge]["x_range"]) - 1):
            a = G.edges[source, target]["wave_velocity"]
            f = G.edges[source, target]["friction_factor"]
            d = G.edges[source, target]["diameter"]
            dx = G.edges[source, target]["dx"]
            UA = G.edges[source, target]["U_mat"][-t, x - 1]
            HA = G.edges[source, target]["H_mat"][-t, x - 1]
            UB = G.edges[source, target]["U_mat"][-t, x + 1]
            HB = G.edges[source, target]["H_mat"][-t, x + 1]
            forward = UA * a / 9.81 + HA - ((f * UA * abs(UA)) / (4 * 9.81 * d / 2)) * dx
            backward = UB * a / 9.81 - HB - ((f * UB * abs(UB)) / (4 * 9.81 * d / 2)) * dx
            B = np.array([forward, backward])
            A = np.array([[a / 9.81, 1], [a / 9.81, -1]])
            UP, HP = np.linalg.solve(A, B)
            G.edges[source, target]["U_mat"][-t - 1, x] = UP
            G.edges[source, target]["H_mat"][-t - 1, x] = HP
    return G


def main(Graph, Envir, progress_bar, ProgressPage):
    global G, t_range, ReservoirHead, dt, count
    G = copy.deepcopy(
        Graph)  # isolated copy of the original graph, used for analysis. This is set because the original graph may be used for steady state calculation
    create_domain(float(Envir["dt"]), int(Envir["TotalTime"]))
    t_range = G.graph["time_seq"]
    haveICFile = messagebox.askyesno(title="Load Initial Condition from Existing File",
                                     message="Click Yes to select initial condition file, no to select location to save a new initial condition file")
    dirname = os.getcwd()
    extensions = [("Json File", ".json")]
    dt = float(Envir["dt"])
    if haveICFile:
        ICFile = filedialog.askopenfilename(initialdir=dirname + "/MOC_IC", title="Open File",
                                            defaultextension=".json",
                                            filetypes=extensions)
    else:
        ICFile = filedialog.asksaveasfilename(initialdir=dirname + "/MOC_IC", title="Save File",
                                              defaultextension=".json",
                                              filetypes=extensions)
        toSteadyState.main(Graph, dt, ICFile)
    ApplyExistingIC(ICFile)
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
    G = node_classification(G)
    decomposed_network = network_decomposition(G)
    count = 0
    for t in range(1, len(t_range)):
        progress_bar["value"] += 100 / len(t_range[1::])
        ProgressPage.update()
        G = analysis(t, decomposed_network, G)
    t_range = G.graph["time_seq"]
    time = t_range[1::]  # to get rid of the -dt term used for steady state
    SaveDict = {}
    if isinstance(Envir["RecordStart"], str):
        saveStartIndex = 0
    else:
        saveStartIndex = int((Envir["RecordStart"]) / dt)
    if isinstance(Envir["RecordEnd"], str):
        saveEndIndex = -1
    else:
        saveEndIndex = int((Envir["RecordEnd"]) / dt)
    for edge in list(G.edges):
        H_mat = G.edges[edge]["H_mat"][1::]  # to get rid of the -dt term used for Steady state
        H_mat = np.flip(H_mat, axis=0)
        source, target = edge
        time_name = "Pipe {0}-{1} Time".format(source, target)
        freq_name = "Pipe {0}-{1} Freq".format(source, target)
        trimed_time = time[saveStartIndex:saveEndIndex]
        trimed_data = H_mat[saveStartIndex:saveEndIndex]
        Nrow, Ncol = trimed_data.shape
        fft_H_mat = np.zeros_like(trimed_data)
        fft_range = np.linspace(0, 1 / dt, len(trimed_time))
        for col in range(Ncol):
            data = trimed_data[:, col]
            data = data - np.mean(data)
            fft_result = abs(np.fft.fft(data)) / (len(fft_range) / 2)
            fft_H_mat[:, col] = fft_result
        SaveTime = np.column_stack((trimed_time, trimed_data))
        SaveFreq = np.column_stack((fft_range, fft_H_mat))
        SaveDict[time_name] = SaveTime.tolist()
        SaveDict[freq_name] = SaveFreq.tolist()
    pyexcel.save_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    ProgressPage.destroy()
    print("File Saved")
