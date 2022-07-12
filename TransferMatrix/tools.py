import numpy as np
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import os


def ReadInput(filename):
    """
    This function read the Input file.
    :param filename: location of the excel file describe the adjacency matrix
    :return: the three page (adjacency matrix in pd form, node attribute and edge attribute in dictionary form.)
    """
    adjpd = pd.read_excel(filename, sheet_name="PipeNetwork", index_col=0, header=0)
    NodeAttribute = pd.read_excel(filename, sheet_name="NodeAttribute", index_col=0, header=0)
    NodeAttribute_dict = NodeAttribute.to_dict(orient="index")
    EdgeAttribute = pd.read_excel(filename, sheet_name="EdgeAttribute", header=0)
    EdgeAttribute_dict = EdgeAttribute.to_dict(orient="index")
    return adjpd, NodeAttribute_dict, EdgeAttribute_dict


def CreateGraph(adjpd, NodeAttribute, EdgeAttribuite):
    """
    Take the Input parameters and create a graph with all field filled.
    :param adjpd: adjancency matrix in pandas dataframe form
    :param NodeAttribute: Node Attribute in Dictionary form
    :param EdgeAttribuite: Edge Attribute in Diectionary form
    :return: G, networkx graph with all field filled.
    """
    G = nx.convert_matrix.from_pandas_adjacency(adjpd, create_using=nx.DiGraph)
    nx.set_node_attributes(G, NodeAttribute)
    for i, item in EdgeAttribuite.items():
        G[item["source"]][item["target"]]["length"] = item["length"]
        G[item["source"]][item["target"]]["diameter"] = item["diameter"]
        G[item["source"]][item["target"]]["material"] = item["material"]
        G[item["source"]][item["target"]]["wave_velocity"] = item["wave_velocity"]
        G[item["source"]][item["target"]]["friction_factor"] = item["friction_factor"]
        G[item["source"]][item["target"]]["flow_velocity"] = item["flow_velocity"]
    return G


def main(filename):
    adjpd, NodeAttribute, EdgeAttribuite = ReadInput(filename)
    G = CreateGraph(adjpd, NodeAttribute, EdgeAttribuite)
    return G


def get_G(file):
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, file)
    G = main(filename)
    return G


def field_matrix_single(G, source, target, freq):
    L = G.edges[source, target]["length"]
    a = G.edges[source, target]["wave_velocity"]
    D = G.edges[source, target]["diameter"]
    fc = G.edges[source, target]["friction_factor"]
    U = G.edges[source, target]["flow_velocity"]
    n = 2  # empirical number
    j = 1j
    g = 9.81
    A = (np.pi * D ** 2) / 4  # D = diameter of the pipe
    Q0 = A * U
    # Q0 = 0.01
    omega = 2 * np.pi * freq  # T = Theoretical period
    R = (n * fc * (Q0 ** (n - 1))) / (2 * g * D * (A ** n))
    "Field Matrix for single conduit"
    mu = np.sqrt(-((omega ** 2) / (a ** 2)) + ((j * g * A * omega * R) / (a ** 2)))
    Zc = (mu * a ** 2) / (j * omega * g * A)
    F = np.array([[np.cosh(mu * L), (-1 / Zc) * np.sinh(mu * L), 0],
                  [-Zc * np.sinh(mu * L), np.cosh(mu * L), 0],
                  [0, 0, 1]])
    return F


def Reverse_field_matrix_single(G, source, target, freq, length):
    L = length
    a = G.edges[source, target]["wave_velocity"]
    D = G.edges[source, target]["diameter"]
    fc = G.edges[source, target]["friction_factor"]
    U = G.edges[source, target]["flow_velocity"]
    n = 2  # empirical number
    j = 1j
    g = 9.81
    A = (np.pi * D ** 2) / 4  # D = diameter of the pipe
    Q0 = A * U
    # Q0 = 0.01
    omega = 2 * np.pi * freq  # T = Theoretical period
    R = (n * fc * (Q0 ** (n - 1))) / (2 * g * D * (A ** n))
    "Field Matrix for single conduit"
    mu = np.sqrt(-((omega ** 2) / (a ** 2)) + ((j * g * A * omega * R) / (a ** 2)))
    Zc = (mu * a ** 2) / (j * omega * g * A)
    F = np.array([[np.cosh(mu * L), (-1 / Zc) * np.sinh(mu * L), 0],
                  [-Zc * np.sinh(mu * L), np.cosh(mu * L), 0],
                  [0, 0, 1]])
    return F


def source_matrix(isPertFlow):
    if isPertFlow == True:
        PMat = [[1, 0, 1],
                [0, 1, 0],
                [0, 0, 1]]
    else:
        PMat = [[1, 0, 0],
                [0, 1, 1],
                [0, 0, 1]]
    return PMat
