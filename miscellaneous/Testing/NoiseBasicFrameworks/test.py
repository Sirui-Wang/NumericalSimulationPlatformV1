import networkx as nx
import numpy as np
import pandas as pd

Sections = np.array([0, 1, 2, 3, 4])
Sources = np.array(["A", "B", "B", "D", "E"])
Targets = np.array(["B", "C", "D", "E", "F"])
Sections_L = np.array([250, 250, 100, 250, 250])
Sections_a = np.full(len(Sections), 1000)
Sections_D = np.array([0.173, 0.173, 0.1, 0.173, 0.173])
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

for edge in G.edges:
    length = int(G.edges[edge]["L"])
    lengthDist = range(0, length, 1)
    dist = np.random.uniform(0, 10, len(lengthDist))
    PressureDist = np.column_stack((lengthDist, dist))
    G.edges[edge]["PressureDist"] = PressureDist

# plt.figure(1)
# edge_labels = {i[0:2]: '{}m'.format(i[2]['L']) for i in G.edges(data=True)}
# pos = nx.spring_layout(G, weight="L")
# nx.draw(G, pos, with_labels=True)
# nx.draw_networkx_edge_labels(G, pos, edge_labels)
# plt.show()
