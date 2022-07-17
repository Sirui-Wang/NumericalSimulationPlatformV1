import networkx as nx

import MOC.MOCEngine as MOCEngine
import TransferMatrix.TransferMatrixEngine as TMEngine

# Auther: Sirui Wang
# Updated date: 13/July/2022
"""Module Comment"""
# Line Comment
global G
G = nx.DiGraph()


def main(data_dict, progress_bar, ProgressPage):
    """Match GUI backend save dictionary to networkx G"""
    NodesDictionary = data_dict["Nodes"]
    LinksDictionary = data_dict["Links"]
    EnvirDictionary = data_dict["Environment"]
    AnalysisMode = data_dict["Mode"]
    for node, NodeAttributes in NodesDictionary.items():
        G.add_node(node)
        G.nodes[node]["head"] = float(NodeAttributes["head"])
        G.nodes[node]["isBC"] = NodeAttributes["isBC"]
        G.nodes[node]["BCType"] = NodeAttributes["BCType"]
    for edge, LinkAttributes in LinksDictionary.items():
        edge = edge.split(",")
        G.add_edge(edge[0], edge[1])
        G.edges[edge[0], edge[1]]["length"] = float(LinkAttributes["Length"])
        G.edges[edge[0], edge[1]]["diameter"] = float(LinkAttributes["Diameter"])
        G.edges[edge[0], edge[1]]["wave_velocity"] = float(LinkAttributes["Wave speed"])
        G.edges[edge[0], edge[1]]["friction_factor"] = float(LinkAttributes["Friction factor"])
        G.edges[edge[0], edge[1]]["flow_velocity"] = float(LinkAttributes["Flow Velocity"])
        G.edges[edge[0], edge[1]]["PertType"] = LinkAttributes["PertType"]
        if AnalysisMode:
            if LinkAttributes["PertType"] == "":
                pass
            else:
                G.edges[edge[0], edge[1]]["PertLocation"] = LinkAttributes["Location"]
            G.edges[edge[0], edge[1]]["HasSensor"] = LinkAttributes["HasSensor"]
            if LinkAttributes["HasSensor"]:
                G.edges[edge[0], edge[1]]["SensorLocation"] = LinkAttributes["SensorLocation"]
                G.edges[edge[0], edge[1]]["SensorDist"] = float(LinkAttributes["SensorDist"])
        else:
            if LinkAttributes["PertType"] == "None":
                pass
            elif LinkAttributes["PertType"] == "Impulse":
                G.edges[edge[0], edge[1]]["Start Time"] = float(LinkAttributes["Time"])
                G.edges[edge[0], edge[1]]["PertLocation"] = LinkAttributes["Location"]
            elif LinkAttributes["PertType"] == "Full Closure":
                G.edges[edge[0], edge[1]]["Start Time"] = float(LinkAttributes["Time"])
                G.edges[edge[0], edge[1]]["PertLocation"] = LinkAttributes["Location"]
            elif LinkAttributes["PertType"] == "Sinusoidal Head":
                G.edges[edge[0], edge[1]]["Freq"] = float(LinkAttributes["Freq"])
                G.edges[edge[0], edge[1]]["Amp"] = float(LinkAttributes["Amp"])
                G.edges[edge[0], edge[1]]["PertLocation"] = LinkAttributes["Location"]
            elif LinkAttributes["PertType"] == "Sinusoidal Flow":
                G.edges[edge[0], edge[1]]["Freq"] = float(LinkAttributes["Freq"])
                G.edges[edge[0], edge[1]]["Amp"] = float(LinkAttributes["Amp"])
                G.edges[edge[0], edge[1]]["PertLocation"] = LinkAttributes["Location"]
            elif LinkAttributes["PertType"] == "Controlled Flow":
                G.edges[edge[0], edge[1]]["Start Time"] = float(LinkAttributes["Time"])
                G.edges[edge[0], edge[1]]["PertLocation"] = LinkAttributes["Location"]
                G.edges[edge[0], edge[1]]["FlowRate"] = float(LinkAttributes["FlowRate"])
    if AnalysisMode:
        TMEngine.main(G, EnvirDictionary, progress_bar, ProgressPage)
    else:
        MOCEngine.main(G, EnvirDictionary, progress_bar, ProgressPage)
