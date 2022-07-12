import json

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


class InconsistantInputFile(Exception):
    pass


def create_domain(G, dt, run_time):
    t_range = np.arange(0, run_time + dt, dt)  # define time axis for the computational square
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
    return G


def IC(G):
    for source, target, data in G.edges.data(data=True):
        v = data["flow_velocity"]
        H_mat = data["H_mat"]
        U_mat = data["U_mat"]
        Ea = G.nodes[source]["head"]
        Eb = G.nodes[target]["head"]
        H_mat[-1, 0] = Ea
        H_mat[-1, -1] = Eb
        dH = (Ea - Eb) / (len(H_mat[-1]) - 1)
        for i in range(1, len(H_mat[-1])):
            H_mat[-1, i] = H_mat[-1, i - 1] - dH
        U_mat[-1, :] = v
        G.edges[source, target]["H_mat"] = H_mat
        G.edges[source, target]["U_mat"] = U_mat
    return G


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


def analysis(t, decomposed_network, G):
    print("pause")
    for node, branch in decomposed_network.items():
        node_class = G.nodes[node]["classification"]
        """Apply a special boundary condition at the sink/source depends on input. Default is Constant Head
        if the node is classified as a sink/source, it must be a boundary condition, even if the input is none. 
        No need to have isBC boolean check because BCtype is created in every nodes' property """
        if node_class == "Sink":
            NodeBCType = G.nodes[node]["BCType"]
            edge = list(branch.keys())[0]
            HA = G.edges[edge]["H_mat"][-t, -2]
            UA = G.edges[edge]["U_mat"][-t, -2]
            a = G.edges[edge]["wave_velocity"]
            f = G.edges[edge]["friction_factor"]
            r = (G.edges[edge]["diameter"]) / 2
            dx = G.edges[edge]["dx"]
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
        elif node_class == "Source":
            NodeBCType = G.nodes[node]["BCType"]
            edge = list(branch.keys())[0]
            HB = G.edges[edge]["H_mat"][-t, 1]
            UB = G.edges[edge]["U_mat"][-t, 1]
            a = G.edges[edge]["wave_velocity"]
            f = G.edges[edge]["friction_factor"]
            r = (G.edges[edge]["diameter"]) / 2
            dx = G.edges[edge]["dx"]
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


def main(G, dt, ICfile):
    """ Analyses """
    G = create_domain(G, dt, 1000)
    G = IC(G)
    G = node_classification(G)
    decomposed_network = network_decomposition(G)
    t_range = G.graph["time_seq"]
    for t in tqdm(range(1, len(t_range))):
        G = analysis(t, decomposed_network, G)
    steadystate_dict = {}
    for edge in G.edges():
        source, target = edge
        steadystate_dict[source + "," + target] = {"H": list(G.edges[edge]["H_mat"][0]),
                                                   "U": list(G.edges[edge]["U_mat"][0])}
    with open(ICfile, "w") as file:
        json.dump(steadystate_dict, file)
    file.close()
    for edge in G.edges():
        plt.figure(str(edge) + "IC")
        plt.plot(t_range, np.flip(G.edges[edge]["H_mat"][:, -2]))
    plt.show()
