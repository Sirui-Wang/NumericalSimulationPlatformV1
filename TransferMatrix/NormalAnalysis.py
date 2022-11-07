from TransferMatrix.TMTools import *


def TMCalculation(data_pack, SubProgressBar, ProgressPage):
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
    Result_dict = dict.fromkeys(list(G.edges), {})
    for junc, branches in branches_dict.items():
        SubProgressBar["value"] += 50 / len(list(branches_dict.keys()))
        ProgressPage.update()
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
    for junc, branches in branches_dict.items():
        SubProgressBar["value"] += 50 / len(list(branches_dict.keys()))
        ProgressPage.update()
        for node, path in branches.items():
            path, direction = path
            for edges in path:
                end_node = edges[-1][-1]
                try:
                    flow_index = np.where(IndexMap == "{} flow, in {}".format(end_node, edges))[0][0]
                    end_node_flow = Solution[flow_index][0]
                except IndexError:
                    end_node_flow = 0
                try:
                    head_index = np.where(IndexMap == "{} head".format(end_node))[0][0]
                    end_node_head = Solution[head_index][0]
                except IndexError:
                    end_node_head = 0
                for edge in edges[-1::-1]:
                    U = np.identity(3)
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
                    C = U @ [[end_node_flow], [end_node_head], [1]]
                    temp_dict = {source: {"flow": C[0][0], "head": C[1][0]},
                                 target: {"flow": end_node_flow, "head": end_node_head}}
                    Result_dict[edge] = temp_dict
                    end_node_flow = C[0][0]
                    end_node_head = C[1][0]
    return Result_dict


def main(Graph, Envir, SubProgressBar, MainProgressBar, ProgressPage, dFreq, freq_range):
    global G
    G = Graph
    G = node_classification(G)
    branches_dict = all_paths(G)
    branches_dict = assign_direction(G, branches_dict)
    branches_dict = SortByJunc(G, branches_dict)
    sub_matrixes = InitializeSubMatrixes(branches_dict, G)
    result_dict = dict.fromkeys(list(G.edges))
    all_result_dict = dict.fromkeys(list(G.edges))
    for key in result_dict.keys():
        result_dict[key] = {"hfreq": np.zeros(len(freq_range), dtype=complex),
                            "qfreq": np.zeros(len(freq_range), dtype=complex)}
        all_result_dict[key] = {"source": {"hfreq": np.zeros(len(freq_range), dtype=complex),
                                           "qfreq": np.zeros(len(freq_range), dtype=complex)},
                                "target": {"hfreq": np.zeros(len(freq_range), dtype=complex),
                                           "qfreq": np.zeros(len(freq_range), dtype=complex)}}
    if Envir["FreqMode"] == "MultiFrequency":
        for i in range(1, len(freq_range[1::])):
            SubProgressBar["value"] = 0
            freq = freq_range[i]
            MainProgressBar["value"] += 100 / len(freq_range[1::])
            ProgressPage.update()
            data = (freq, branches_dict, sub_matrixes)
            Solutions = TMCalculation(data, SubProgressBar, ProgressPage)
            CalculateAtSensor(G, Solutions, i, freq, result_dict)
            CalculateAllNode(G, Solutions, i, all_result_dict)
    else:
        freq = np.array([float(Envir["ExcitationFreq"])])
        data = (freq, branches_dict, sub_matrixes)
        Solutions = TMCalculation(data)
        MainProgressBar["value"] += 100
        ProgressPage.update()
        CalculateAtSensor(Solutions, int(freq / dFreq), freq, result_dict)
        CalculateAllNode(Solutions, int(freq / dFreq), all_result_dict)
    return result_dict, all_result_dict
