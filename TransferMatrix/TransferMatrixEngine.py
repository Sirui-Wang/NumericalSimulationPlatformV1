from copy import deepcopy

import networkx as nx
import numpy as np
import scipy

import TransferMatrix.tools as tools


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
        else:
            classification = "DoF Error"
        print(node, classification)
        G.nodes[node]["classification"] = classification
    return G


def contained_in(lst, sub):
    return ','.join(map(str, sub)) in ','.join(map(str, lst))


def IdentJunction(G):
    JunctionList = []
    for node in G.nodes():
        if G.nodes[node]["classification"] == "BranchConnection":
            JunctionList.append(node)
    return JunctionList


def path2edge_simple(node_list):
    path = []
    for i in range(len(node_list) - 1):
        path.append((node_list[i], node_list[i + 1]))
    return path


def check_for_loop(G, JunctionList, junc):
    temp = JunctionList.copy()
    temp.remove(junc)
    loopdict = {}
    for i in temp:
        loops = list(nx.algorithms.simple_paths.all_simple_paths(G, junc, i)) + list(
            nx.algorithms.simple_paths.all_simple_paths(G, i, junc))
        temp_dict = {}
        loopCount = 0
        for loop in loops:
            temp_dict.update({"Loop{}".format(loopCount): path2edge_simple(loop)})
            loopCount += 1
        loopdict[tuple(set((junc, i)))] = temp_dict
    return loopdict


def all_paths(G):
    nodes = list(G.nodes)
    all_dict = {}
    for node in nodes:
        temp_list = nodes.copy()
        temp_list.remove(node)
        node_dict = {}
        for i in temp_list:
            simplePath = list(nx.algorithms.simple_paths.all_simple_edge_paths(G, node, i))
            node_dict[i] = simplePath
        all_dict[node] = node_dict
    copy_dict = deepcopy(all_dict)
    trimed = trim_all_path(G, copy_dict)
    return trimed


def JunctionNOTInBetween(G, source, target):
    simplePath = list(nx.algorithms.simple_paths.all_simple_paths(G, source, target))
    JunctionList = IdentJunction(G)
    status = []
    junction = []
    for index in range(len(simplePath)):
        if set(simplePath[index][1:-1]).isdisjoint((JunctionList)):
            status.append(True)
        else:
            status.append(False)
            junction.append(list(set(simplePath[index][1:-1]).intersection((JunctionList))))
    return status, junction


def trim_all_path(G, copy_dict):
    for node in list(copy_dict.keys()):
        for i in list(copy_dict[node].keys()):
            if G.nodes[i]["classification"] == "SeriesConnection":
                del copy_dict[node][i]
            elif len(copy_dict[node][i]) == 0:
                del copy_dict[node][i]
        if len(copy_dict[node]) == 0:
            del copy_dict[node]
        elif G.nodes[node]["classification"] == "SeriesConnection":
            del copy_dict[node]
    for source in list(copy_dict.keys()):
        for target in list(copy_dict[source].keys()):
            status, JunctionInBetween = JunctionNOTInBetween(G, source, target)
            paths = copy_dict[source][target]
            ref_path = deepcopy(paths)
            for index in range(len(paths)):
                if not status[index]:
                    copy_dict[source][target].remove(ref_path[index])
            if len(copy_dict[source][target]) == 0:
                del copy_dict[source][target]
    return copy_dict


def assign_direction(G, branches_dict):
    for source, branches in branches_dict.items():
        for target, path in branches.items():
            if G.nodes[source]["classification"] == "BranchConnection" and G.nodes[target][
                "classification"] != "BranchConnection":
                direction = "downstream"
            elif G.nodes[source]["classification"] != "BranchConnection" and G.nodes[target][
                "classification"] == "BranchConnection":
                direction = "upstream"
            elif G.nodes[source]["classification"] == "Source" and G.nodes[target]["classification"] == "Sink":
                direction = "Simple"
            else:
                direction = "midstream"
            branches_dict[source][target] = (path, direction)
    return branches_dict


def SortByJunc(G, branches_dict):
    Sorted_dict = dict.fromkeys(IdentJunction(G), {})
    JunctionList = IdentJunction(G)
    if len(JunctionList) == 0 and len(list(branches_dict.keys())) == 1:
        Sorted_dict = branches_dict
    else:
        included = set()
        for junc in JunctionList:
            temp = {}
            for source, branches in branches_dict.items():
                for target, data in branches.items():
                    if not tuple(sorted((source, target))) in included:
                        if junc == source:
                            temp.update({target: data})
                        elif junc == target:
                            temp.update({source: data})
                        included.add(tuple(sorted((source, target))))
            Sorted_dict[junc] = temp
    return Sorted_dict


class UnknownConditionDetected(Exception):
    pass


def InitializeSubMatrixes(branches_dict, G):
    SubMatrixes = dict.fromkeys(list(branches_dict.keys()))
    for junc, branches in branches_dict.items():
        subbranches = len(branches)
        for node, data in branches.items():
            path, direction = data
            if len(path) > 1:
                subbranches += len(path) - 1
        nCol = subbranches * 4
        nRow = subbranches * 2 + 1  # this would work for a junction, but not for a simple/ series pipe. for simple/series pipe, nRow should be 2
        aa = np.zeros((nRow, nCol), dtype=complex)
        row = 0
        col = 0
        SubMatrixes[junc] = {"nCol": nCol, "nRow": nRow, "aa": aa, "row": row, "col": col}
    return SubMatrixes


def TMCalculation(data_pack):
    freq, branches_dict, sub_matrixes = data_pack
    juncs = list(sub_matrixes.keys())
    if len(juncs) == 1 and juncs[0] == "series":
        A = sub_matrixes["series"]["aa"]
    else:
        A = sub_matrixes[juncs[0]]["aa"]
        for junc in juncs[1::]:
            A = scipy.linalg.block_diag(A, sub_matrixes[junc]["aa"])
    A = np.array(A, dtype=complex)
    B = np.zeros((len(A), 1), dtype=complex)
    global IndexMap
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
        records = (RowIterable, ColIterable, JuncIterable, column2del, h_placed)
        for node, data in branches.items():
            paths, direction = data
            for path in paths:
                U = np.identity(3)
                for edge in path[-1::-1]:
                    source, target = edge
                    if G.edges[edge]["PertType"] != "":
                        PertLocation = G.edges[edge]["PertLocation"]
                        PertType = G.edges[edge]["PertType"]
                        if PertType == "Head":
                            S = tools.source_matrix(False)
                        else:
                            S = tools.source_matrix(True)
                        if PertLocation == source:
                            U = S @ tools.field_matrix_single(G, source, target, freq) @ U
                        else:
                            U = tools.field_matrix_single(G, source, target, freq) @ S @ U
                    else:
                        U = tools.field_matrix_single(G, source, target, freq) @ U
                records, A, B = creatingMatrix(U, direction, path, junc, node, records, A, B)
        RowIterable, ColIterable, JuncIterable, column2del, h_placed = records
    A = np.delete(A, column2del, 1)
    IndexMap = np.delete(IndexMap, column2del, 0)
    Solution = np.linalg.solve(A, B)
    for junc, branches in branches_dict.items():
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
                    if G.edges[edge]["PertType"] != "":
                        PertLocation = G.edges[edge]["PertLocation"]
                        PertType = G.edges[edge]["PertType"]
                        if PertType == "Head":
                            S = tools.source_matrix(False)
                        else:
                            S = tools.source_matrix(True)
                        FieldMatrix = tools.field_matrix_single(G, source, target, freq)
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


def creatingMatrix(U, direction, path, junc, node, records, A, B):
    global IndexMap
    RowIterable, ColIterable, JuncIterable, column2del, h_placed = records
    a = U[0][0]
    b = U[0][1]
    c = U[0][2]
    d = U[1][0]
    e = U[1][1]
    f = U[1][2]
    if direction == "upstream":
        # upstream relative to the junction, ie q(1) and h(3) are must have and Q(0) and H(2) are optional depends on boundary condition
        startNode = path[0][0]
        NodeBCType = G.nodes[startNode]["BCType"]
        if NodeBCType == "Constant Flow":
            IndexMap[ColIterable] = "{} flow, in {}".format(path[0][0], path)
            IndexMap[ColIterable + 1] = "{} flow, in {}".format(path[-1][-1], path)
            IndexMap[ColIterable + 2] = "{} head".format(path[0][0])
            IndexMap[ColIterable + 3] = "{} head".format(path[-1][-1])
            a_temp = [[0, a, 0, b],
                      [0, d, -1, e]]
            b_temp = [[-c], [-f]]
            A[RowIterable:RowIterable + 2, ColIterable:ColIterable + 4] = a_temp
            B[RowIterable:RowIterable + 2] = b_temp
            A[-JuncIterable[junc]][ColIterable + 1] = 1
            col_del_index = ColIterable + 0
            column2del.append(col_del_index)
            if h_placed[junc][0]:
                col_del_index = ColIterable + 3  # remove h (index 3) as this is the upstream of the junc
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[junc][1]] = A[RowIterable:RowIterable + 2, ColIterable + 3]
            else:
                h_placed[junc] = [True, ColIterable + 3]
            ColIterable += 4
            ColIterable += 2
        else:
            IndexMap[ColIterable] = "{} flow, in {}".format(path[0][0], path)
            IndexMap[ColIterable + 1] = "{} flow, in {}".format(path[-1][-1], path)
            IndexMap[ColIterable + 2] = "{} head".format(path[0][0])
            IndexMap[ColIterable + 3] = "{} head".format(path[-1][-1])
            a_temp = [[-1, a, 0, b],
                      [0, d, 0, e]]  # index 2 need to be deleted as H at the opposite of junction is reservoir hence 0
            b_temp = [[-c], [-f]]
            A[RowIterable:RowIterable + 2, ColIterable:ColIterable + 4] = a_temp
            B[RowIterable:RowIterable + 2] = b_temp
            A[-JuncIterable[junc]][ColIterable + 1] = 1
            col_del_index = ColIterable + 2
            column2del.append(col_del_index)
            if h_placed[junc][0]:
                col_del_index = ColIterable + 3
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[junc][1]] = A[RowIterable:RowIterable + 2, ColIterable + 3]
            else:
                h_placed[junc] = [True, ColIterable + 3]
            ColIterable += 4
            RowIterable += 2
    elif direction == "midstream":
        if junc == path[0][0]:  # the branch analyzed is downstream of the junction
            IndexMap[ColIterable] = "{} flow, in {}".format(path[0][0], path)
            IndexMap[ColIterable + 1] = "{} flow, in {}".format(path[-1][-1], path)
            IndexMap[ColIterable + 2] = "{} head".format(path[0][0])
            IndexMap[ColIterable + 3] = "{} head".format(path[-1][-1])
            a_temp = [[-1, a, 0, b], [0, d, -1, e]]
            b_temp = [[-c], [-f]]
            A[RowIterable:RowIterable + 2, ColIterable:ColIterable + 4] = a_temp
            B[RowIterable:RowIterable + 2] = b_temp
            A[-JuncIterable[junc]][ColIterable + 0] = -1
            A[-JuncIterable[node]][ColIterable + 1] = 1
            if h_placed[junc][0]:
                col_del_index = ColIterable + 2
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[junc][1]] = A[RowIterable:RowIterable + 2, ColIterable + 2]
            else:
                h_placed[junc] = [True, ColIterable + 2]
            if h_placed[node][0]:
                col_del_index = ColIterable + 3
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[node][1]] = A[RowIterable:RowIterable + 2, ColIterable + 3]
            else:
                h_placed[node] = [True, ColIterable + 3]
            ColIterable += 4
            RowIterable += 2
        elif junc == path[-1][-1]:
            """
            path is upstream of the junction
            """
            """Assume we cant model valve closure at the middle"""
            IndexMap[ColIterable] = "{} flow, in {}".format(path[0][0], path)
            IndexMap[ColIterable + 1] = "{} flow, in {}".format(path[-1][-1], path)
            IndexMap[ColIterable + 2] = "{} head".format(path[0][0])
            IndexMap[ColIterable + 3] = "{} head".format(path[-1][-1])
            a_temp = [[-1, a, 0, b], [0, d, 0, e]]
            b_temp = [[-c], [-f]]
            A[RowIterable:RowIterable + 2, ColIterable:ColIterable + 4] = a_temp
            B[RowIterable:RowIterable + 2] = b_temp
            A[-JuncIterable[junc]][ColIterable + 1] = 1
            A[-JuncIterable[node]][ColIterable + 0] = -1
            if h_placed[junc][0]:
                col_del_index = ColIterable + 3
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[junc][1]] = A[RowIterable:RowIterable + 2, ColIterable + 3]
            else:
                h_placed[junc] = [True, ColIterable + 3]
            if h_placed[node][0]:
                col_del_index = ColIterable + 2
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[node][1]] = A[RowIterable:RowIterable + 2, ColIterable + 2]
            else:
                h_placed[node] = [True, ColIterable + 2]
            ColIterable += 4
            RowIterable += 2
    elif direction == "downstream":
        endNode = path[-1][-1]
        NodeBCType = G.nodes[endNode]["BCType"]
        if NodeBCType == "Constant Flow":
            IndexMap[ColIterable] = "{} flow, in {}".format(path[0][0], path)
            IndexMap[ColIterable + 1] = "{} flow, in {}".format(path[-1][-1], path)
            IndexMap[ColIterable + 2] = "{} head".format(path[0][0])
            IndexMap[ColIterable + 3] = "{} head".format(path[-1][-1])
            a_temp = [[-1, 0, 0, b], [0, 0, -1, e]]
            b_temp = [[-c], [-f]]
            A[RowIterable:RowIterable + 2, ColIterable:ColIterable + 4] = a_temp
            B[RowIterable:RowIterable + 2] = b_temp
            A[-JuncIterable[junc]][ColIterable + 0] = -1
            col_del_index = ColIterable + 1
            column2del.append(col_del_index)
            if h_placed[junc][0]:
                col_del_index = ColIterable + 2
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[junc][1]] = A[RowIterable:RowIterable + 2, ColIterable + 2]
            else:
                h_placed[junc] = [True, ColIterable + 2]
            ColIterable += 4
            RowIterable += 2
            """
            remember to change all indexes
            """
        else:
            IndexMap[ColIterable] = "{} flow, in {}".format(path[0][0], path)
            IndexMap[ColIterable + 1] = "{} flow, in {}".format(path[-1][-1], path)
            IndexMap[ColIterable + 2] = "{} head".format(path[0][0])
            IndexMap[ColIterable + 3] = "{} head".format(path[-1][-1])
            a_temp = [[-1, a, 0, 0], [0, d, -1, 0]]
            b_temp = [[-c], [-f]]
            A[RowIterable:RowIterable + 2, ColIterable:ColIterable + 4] = a_temp
            B[RowIterable:RowIterable + 2] = b_temp
            A[-JuncIterable[junc]][ColIterable + 0] = -1
            col_del_index = ColIterable + 3
            column2del.append(col_del_index)
            if h_placed[junc][0]:
                col_del_index = ColIterable + 2
                column2del.append(col_del_index)
                A[RowIterable:RowIterable + 2, h_placed[junc][1]] = A[RowIterable:RowIterable + 2, ColIterable + 2]
            else:
                h_placed[junc] = [True, ColIterable + 2]
            ColIterable += 4
            RowIterable += 2
    elif direction == "Simple":
        """ Because the code is set up in such a way that the junction node will always be the upstream node.
        """
        print("placeholder")
    records = (RowIterable, ColIterable, JuncIterable, column2del, h_placed)
    return records, A, B


def main(Graph, Envir, progress_bar, ProgressPage):
    global G
    G = Graph
    G = node_classification(G)
    branches_dict = all_paths(G)
    branches_dict = assign_direction(G, branches_dict)
    branches_dict = SortByJunc(G, branches_dict)
    sub_matrixes = InitializeSubMatrixes(branches_dict, G)
    dFreq = float(Envir["df"])
    maxFreq = int(Envir["MaxFreq"])
    freq_range = np.arange(0, maxFreq + dFreq, dFreq)
    hfreq_array = np.zeros(len(freq_range), dtype=complex)
    qfreq_array = np.zeros(len(freq_range), dtype=complex)
    result_dict = dict.fromkeys(G.edges, {"hfreq": hfreq_array, "qfreq": qfreq_array, "SensorLocation": ""})
    if Envir["FreqMode"] == "MultiFrequency":
        for i in range(1, len(freq_range[1::])):
            freq = freq_range[i]
            progress_bar["value"] += 100 / len(freq_range[1::])
            ProgressPage.update()
            data = (freq, branches_dict, sub_matrixes)
            Solutions = TMCalculation(data)
    else:
        freq = np.array([Envir["ExcitationFreq"]])
        data = (freq, branches_dict, sub_matrixes)
        Solutions = TMCalculation(data)

    #     for edge in G.edges:
    #         HasSensor = G.edges[edge]["HasSensor"]
    #         if HasSensor == True or HasSensor == "normal":
    #             Sensor_loc = G.edges[edge]["SensorLocation"]
    #             dx = G.edges[edge]["SensorDist"]
    #             HasValve = G.edges[edge]["HasValve"]
    #             valveLocation = G.edges[edge]["ValveLocation"]
    #             freqHead = Solutions[edge][Sensor_loc]["head"]
    #             freqFlow = Solutions[edge][Sensor_loc]["flow"]
    #             if HasValve:
    #                 S = tools.source_matrix(G.edges[edge]["ValveMovement"] == "constant")
    #                 if Sensor_loc == edge[0]:
    #                     length = G.edges[edge]["length"] - dx
    #                 else:
    #                     length = dx
    #                 if valveLocation == edge[0]:
    #                     FieldMatrix = tools.Reverse_field_matrix_single(G, edge[0], edge[1], freq, length)
    #                     # backCalculatedResult = S @ FieldMatrix @ [[freqFlow], [freqHead], [1]]
    #                     backCalculatedResult = FieldMatrix @ [[freqFlow], [freqHead], [1]]
    #                     # ASK!!! if we back calculate from the downstream, if the valve/perturbation is at the
    #                     # upstream, should we multiply by source matrix again?
    #                 else:
    #                     FieldMatrix = tools.Reverse_field_matrix_single(G, edge[0], edge[1], freq, length)
    #                     backCalculatedResult = FieldMatrix @ S @ [[freqFlow], [freqHead], [1]]
    #             else:
    #                 if Sensor_loc == edge[0]:
    #                     length = G.edges[edge]["length"] - dx
    #                 else:
    #                     length = dx
    #                 FieldMatrix = tools.Reverse_field_matrix_single(G, edge[0], edge[1], freq, length)
    #                 backCalculatedResult = FieldMatrix @ [[freqFlow], [freqHead], [1]]
    #             result_dict[edge]["hfreq"][i] = backCalculatedResult[1]
    #             result_dict[edge]["qfreq"][i] = backCalculatedResult[0]
    #             # result_dict[edge]["hfreq"][i] = freqHead
    #             # result_dict[edge]["qfreq"][i] = freqFlow
    # time = np.arange(0, (1 / dFreq) + 1 / (maxFreq), 1 / (maxFreq))
    # df1 = pd.DataFrame()
    # df1["Time"] = time
    # for key, item in result_dict.items():
    #     hfreq = item["hfreq"]
    #     qfreq = item["qfreq"]
    #     sensor_loc = item["SensorLocation"]
    #     hTime_default = np.fft.ifft(hfreq, len(hfreq))
    #     qTime_default = np.fft.ifft(qfreq, len(qfreq))
    #     result = itertools.accumulate(hTime_default)
    #     result = np.real(list(result))
    #     plt.figure("H at {} in Frequency Domain".format(sensor_loc))
    #     plt.plot(freq_range, abs(hfreq))
    #     plt.figure("Q at {} in Frequency Domain".format(sensor_loc))
    #     plt.plot(freq_range, abs(qfreq))
    #     plt.figure("H at {}".format(sensor_loc))
    #     plt.plot(time, result)  # accumulated Result
    #     plt.plot(time, np.real(hTime_default))  # step change result
    #     plt.figure("Q at {}".format(sensor_loc))
    #     plt.plot(time, np.real(list(qTime_default)))
    #     df1["delta_H at {}".format(sensor_loc)] = np.real(hTime_default)
    #     df1["accu_Hat at {}".format(sensor_loc)] = result
    #     df1["delta_Qat at {}".format(sensor_loc)] = np.real(qTime_default)
    #     df1["HFreq at {}".format(sensor_loc)] = abs(hfreq)
    #     df1["QFreq at {}".format(sensor_loc)] = abs(qfreq)
    # dirname = os.getcwd()
    # extensions = [("Excel File", ".xlsx")]
    # Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", title="Save File",
    #                                           defaultextension=".xlsx", filetypes=extensions)
    # df1.to_excel(Output_loc, sheet_name="TM")
    # plt.show()
    # ProgressPage.destroy()
