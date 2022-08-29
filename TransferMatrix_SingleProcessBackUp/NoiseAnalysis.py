from TransferMatrix.TMTools import *


def TMCalculation(data_pack):
    global SensorResult, i
    freq, branches_dict, sub_matrixes = data_pack
    juncs = list(sub_matrixes.keys())
    A = sub_matrixes[juncs[0]]["aa"]
    for junc in juncs[1::]:
        A = block_diag_selfmade(A, sub_matrixes[junc]["aa"])
    A = np.array(A, dtype=complex)
    B = np.zeros((len(A), 1), dtype=complex)
    IndexMap = np.empty((np.shape(A)[1]), dtype=object)
    RowIterable = 0
    ColIterable = 0
    column2del = []
    h_placed = dict.fromkeys(list(branches_dict.keys()), [False, 0])
    junc_index = 1
    JuncIterable = dict.fromkeys(list(branches_dict.keys()))
    for key in list(JuncIterable.keys()):
        JuncIterable[key] = junc_index
        junc_index += 1
    for junc, branches in branches_dict.items():
        records = (RowIterable, ColIterable, JuncIterable, column2del, h_placed, IndexMap)
        for node, data in branches.items():
            """ node here is essentially the end/start node opposite to the junction
            """
            paths, direction = data
            for path in paths:
                Impedances = np.zeros(2, dtype=complex)
                U = np.identity(3)
                for edge in path[-1::-1]:
                    source, target = edge
                    FieldMatrix, Zc = field_matrix_single(G, source, target, freq)
                    if direction == "upstream" and len(path) == 1:
                        Impedances[0] = Zc
                    elif direction == "downstream" and len(path) == 1:
                        Impedances[-1] = -Zc
                    elif direction == "simple" and len(path) == 1:
                        Impedances[0] = Zc
                        Impedances[-1] = -Zc
                    elif edge == path[-1]:
                        Impedances[-1] = -Zc
                    elif edge == path[0]:
                        Impedances[0] = Zc
                    if G.edges[edge]["PertType"] != "":
                        PertLocation = G.edges[edge]["PertLocation"]
                        PertType = G.edges[edge]["PertType"]
                        S = source_matrix(PertType == "Flow")
                        if PertLocation == source:
                            U = S @ FieldMatrix @ U
                        else:
                            U = FieldMatrix @ S @ U
                    else:
                        U = FieldMatrix @ U
                records, A, B = creatingMatrix(G, U, direction, path, junc, node, records, A, B, Impedances)
        RowIterable, ColIterable, JuncIterable, column2del, h_placed, IndexMap = records
    rows2del = np.argwhere(np.all(np.array(A) == 0, axis=1))
    A = np.delete(A, column2del, 1)
    A = np.delete(A, rows2del, 0)
    B = np.delete(B, rows2del, 0)
    IndexMap = np.delete(IndexMap, column2del, 0)
    Solution = np.linalg.solve(A, B)
    SensorPath = []
    SensorResultSaved = dict.fromkeys(Sensors_list, False)
    for sensor in Sensors_list:
        inedges = list(G.in_edges(sensor))
        outedges = list(G.out_edges(sensor))
        EdgesFromSensor = inedges + outedges
        for junc, branches in branches_dict.items():
            for node, path in branches.items():
                path, direction = path
                for edges in path:
                    if bool(set(edges) & set(EdgesFromSensor)) and not SensorResultSaved[sensor]:
                        SensorPath.append((sensor, edges))
                        SensorResultSaved[sensor]=True
    for path in SensorPath:
        sensor, edges = path
        downstreamNode = edges[-1][-1]
        try:
            flow_index = np.where(IndexMap == "{} flow, in {}".format(downstreamNode, edges))[0][0]
            end_node_flow = Solution[flow_index][0]
        except IndexError:
            end_node_flow = 0
        try:
            head_index = np.where(IndexMap == "{} head".format(downstreamNode))[0][0]
            end_node_head = Solution[head_index][0]
        except IndexError:
            end_node_head = 0
        U = np.identity(3)
        for edge in edges[-1::-1]:
            source, target = edge
            FieldMatrix, Zc = field_matrix_single(G, source, target, freq)
            if G.edges[edge]["PertType"] != "":
                PertLocation = G.edges[edge]["PertLocation"]
                PertType = G.edges[edge]["PertType"]
                S = source_matrix(PertType == "Flow")
                if PertLocation == source:
                    U = S @ FieldMatrix @ U
                elif PertLocation == target:
                    U = FieldMatrix @ S @ U
            else:
                U = FieldMatrix @ U
            if source == sensor:
                C = U @ [[end_node_flow], [end_node_head], [1]]
                SensorResult[sensor]["hfreq"][i] = C[1][0]
                SensorResult[sensor]["qfreq"][i] = C[0][0]
                break
            elif target == sensor:
                SensorResult[sensor]["hfreq"][i] = end_node_head
                SensorResult[sensor]["qfreq"][i] = end_node_flow
                break


def main(Graph, Envir, freq_range, SubProgressBar, ProgressPage):
    global G, Sensors_list, SensorResult, i
    G = Graph
    G = node_classification(G)
    branches_dict = all_paths(G)
    branches_dict = assign_direction(G, branches_dict)
    branches_dict = SortByJunc(G, branches_dict)
    sub_matrixes = InitializeSubMatrixes(branches_dict, G)
    Sensors_list = [Envir["Sensor1"], Envir["Sensor2"]]
    SensorResult = dict.fromkeys(Sensors_list)
    for sensor in Sensors_list:
        SensorResult[sensor] = {"hfreq": np.zeros(len(freq_range), dtype=complex),
                                "qfreq": np.zeros(len(freq_range), dtype=complex)}
    for i in range(1, len(freq_range[1::])):
        SubProgressBar["value"] += 100 / len(freq_range[1::])
        ProgressPage.update()
        freq = freq_range[i]
        data = (freq, branches_dict, sub_matrixes)
        TMCalculation(data)
    return SensorResult
