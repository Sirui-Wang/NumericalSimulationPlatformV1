import os
from copy import deepcopy
from tkinter import filedialog

import networkx as nx
import numpy as np
import pyexcel


def block_diag_selfmade(*arrs):
    """Create a block diagonal matrix from the provided arrays.

    Given the inputs `A`, `B` and `C`, the output will have these
    arrays arranged on the diagonal::

        [[A, 0, 0],
         [0, B, 0],
         [0, 0, C]]

    If all the input arrays are square, the output is known as a
    block diagonal matrix.

    Parameters
    ----------
    A, B, C, ... : array-like, up to 2D
        Input arrays.  A 1D array or array-like sequence with length n is
        treated as a 2D array with shape (1,n).

    Returns
    -------
    D : ndarray
        Array with `A`, `B`, `C`, ... on the diagonal.  `D` has the
        same dtype as `A`.

    References
    ----------
    .. [1] Wikipedia, "Block matrix",
           http://en.wikipedia.org/wiki/Block_diagonal_matrix
    """
    if arrs == ():
        arrs = ([],)
    arrs = [np.atleast_2d(a) for a in arrs]

    bad_args = [k for k in range(len(arrs)) if arrs[k].ndim > 2]
    if bad_args:
        raise ValueError("arguments in the following positions have dimension "
                         "greater than 2: %s" % bad_args)

    shapes = np.array([a.shape for a in arrs])
    out = np.zeros(np.sum(shapes, axis=0), dtype=arrs[0].dtype)

    r, c = 0, 0
    for i, (rr, cc) in enumerate(shapes):
        out[r:r + rr, c:c + cc] = arrs[i]
        r += rr
        c += cc
    return out


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
    A = sub_matrixes[juncs[0]]["aa"]
    for junc in juncs[1::]:
        A = block_diag_selfmade(A, sub_matrixes[junc]["aa"])
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
            """ node here is essentially the end/start node opposite to the junction
            """
            paths, direction = data
            for path in paths:
                U = np.identity(3)
                for edge in path[-1::-1]:
                    source, target = edge
                    if G.edges[edge]["PertType"] != "":
                        PertLocation = G.edges[edge]["PertLocation"]
                        PertType = G.edges[edge]["PertType"]
                        S = source_matrix(PertType == "Flow")
                        if PertLocation == source:
                            U = S @ field_matrix_single(G, source, target, freq) @ U
                        else:
                            U = field_matrix_single(G, source, target, freq) @ S @ U
                    else:
                        U = field_matrix_single(G, source, target, freq) @ U
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
                        S = source_matrix(PertType == "Flow")
                        FieldMatrix = field_matrix_single(G, source, target, freq)
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
        # startNode = path[0][0]
        NodeBCType = G.nodes[node]["BCType"]
        if NodeBCType == "Constant flow":
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
        # endNode = path[-1][-1]
        NodeBCType = G.nodes[node]["BCType"]
        if NodeBCType == "Constant flow":
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
        IndexMap[ColIterable] = "{} flow, in {}".format(path[0][0], path)
        IndexMap[ColIterable + 1] = "{} flow, in {}".format(path[-1][-1], path)
        IndexMap[ColIterable + 2] = "{} head".format(path[0][0])
        IndexMap[ColIterable + 3] = "{} head".format(path[-1][-1])
        """ Because the code is set up in such a way that the junction node will always be the upstream node."""
        upBCType = G.nodes[junc]["BCType"]
        downBCType = G.nodes[node]["BCType"]
        if upBCType == "Constant head" and downBCType == "Constant head":
            # since this is a simple path, there are and should only be two unknow for 2 eqn to solve.
            # where previously it all connect to a junction, hence there will be 3 or 4 unknown needed solving
            # H = 0 and h = 0, Q and q need solving resulting in index 2 and 3 being empty
            a_temp = [[-1, a, 0, 0], [0, d, 0, 0]]
            b_temp = [[-c], [-f]]
            column2del.append(ColIterable + 2)
            column2del.append(ColIterable + 3)
        elif upBCType == "Constant head" and downBCType == "Constant flow":
            a_temp = [[-1, 0, 0, b], [0, 0, 0, e]]
            b_temp = [[-c], [-f]]
            column2del.append(ColIterable + 1)
            column2del.append(ColIterable + 2)
        elif upBCType == "Constant flow" and downBCType == "Constant head":
            a_temp = [[0, a, 0, 0], [0, d, -1, 0]]
            b_temp = [[-c], [-f]]
            column2del.append(ColIterable + 0)
            column2del.append(ColIterable + 3)
        elif upBCType == "Constant flow" and downBCType == "Constant flow":
            a_temp = [[0, 0, 0, b], [0, 0, -1, e]]
            b_temp = [[-c], [-f]]
            column2del.append(ColIterable + 0)
            column2del.append(ColIterable + 1)
        A[RowIterable:RowIterable + 2, ColIterable:ColIterable + 4] = a_temp
        B[RowIterable:RowIterable + 2] = b_temp
        A = np.delete(A, 2, 0)
        B = np.delete(B, 2, 0)
    records = (RowIterable, ColIterable, JuncIterable, column2del, h_placed)
    return records, A, B


def CalculateAtSensor(Solutions, i, freq):
    global result_dict
    for edge in G.edges:
        HasSensor = G.edges[edge]["HasSensor"]
        if HasSensor:
            SensorLocation = G.edges[edge]["SensorLocation"]
            SensorDist = G.edges[edge]["SensorDist"]
            PertType = G.edges[edge]["PertType"]
            """Always take the value from the downstream end to forward calculate value at the sensor"""
            freqHead = Solutions[edge][edge[1]]["head"]
            freqFlow = Solutions[edge][edge[1]]["flow"]
            if SensorLocation == edge[0]:
                """always distance from the downstream end, if sensor is at the upstream, then the distance to 
                downstream node is length - SensorDist """
                length = G.edges[edge]["length"] - SensorDist
            else:
                """if sensor is at the downstream, then the length is just SensorDist since we are always calculate 
                from the downstream end """
                length = SensorDist
            FieldMatrix = Reverse_field_matrix_single(G, edge[0], edge[1], freq, length)
            if PertType != "":
                PertLocation = G.edges[edge]["PertLocation"]
                S = source_matrix(PertType == "Flow")
                if PertLocation == edge[1]:
                    backCalculatedResult = FieldMatrix @ S @ [[freqFlow], [freqHead], [1]]
                else:
                    backCalculatedResult = FieldMatrix @ [[freqFlow], [freqHead], [1]]
            else:
                backCalculatedResult = FieldMatrix @ [[freqFlow], [freqHead], [1]]
            result_dict[edge]["hfreq"][i] = backCalculatedResult[1]
            result_dict[edge]["qfreq"][i] = backCalculatedResult[0]


def CalculateAllNode(Solutions, i):
    global all_result_dict
    for edge in G.edges:
        source, target = edge
        sourceFlow = Solutions[edge][source]["flow"]
        sourceHead = Solutions[edge][source]["head"]
        targetFlow = Solutions[edge][target]["flow"]
        targetHead = Solutions[edge][target]["head"]
        all_result_dict[edge]["source"]["hfreq"][i] = sourceHead
        all_result_dict[edge]["source"]["qfreq"][i] = sourceFlow
        all_result_dict[edge]["target"]["hfreq"][i] = targetHead
        all_result_dict[edge]["target"]["qfreq"][i] = targetFlow


def main(Graph, Envir, progress_bar, ProgressPage):
    dirname = os.getcwd()
    extensions = [("Excel File", ".xlsx")]
    Output_loc = filedialog.asksaveasfilename(initialdir=dirname + "/Results", initialfile="Result", title="Save File",
                                              defaultextension=".xlsx",
                                              filetypes=extensions)  # Save File
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
    global result_dict, all_result_dict
    result_dict = dict.fromkeys(list(G.edges))
    all_result_dict = dict.fromkeys(list(G.edges))
    for key in result_dict.keys():
        result_dict[key] = {"hfreq": np.zeros(len(freq_range), dtype=complex), "qfreq": np.zeros(len(freq_range), dtype=complex)}
        all_result_dict[key] = {"source": {"hfreq": np.zeros(len(freq_range), dtype=complex), "qfreq": np.zeros(len(freq_range), dtype=complex)},
                                              "target": {"hfreq": np.zeros(len(freq_range), dtype=complex), "qfreq": np.zeros(len(freq_range), dtype=complex)}}
    if Envir["FreqMode"] == "MultiFrequency":
        for i in range(1, len(freq_range[1::])):
            freq = freq_range[i]
            progress_bar["value"] += 100 / len(freq_range[1::])
            ProgressPage.update()
            data = (freq, branches_dict, sub_matrixes)
            Solutions = TMCalculation(data)
            CalculateAtSensor(Solutions, i, freq)
            CalculateAllNode(Solutions, i)
    else:
        freq = np.array([float(Envir["ExcitationFreq"])])
        data = (freq, branches_dict, sub_matrixes)
        Solutions = TMCalculation(data)
        progress_bar["value"] += 100
        ProgressPage.update()
        CalculateAtSensor(Solutions, int(freq / dFreq), freq)
        CalculateAllNode(Solutions, int(freq / dFreq))
    time = np.arange(0, (1 / dFreq) + (1 / maxFreq), 1 / maxFreq)
    SaveDict = {}
    for edge in G.edges:
        source, target = edge
        SensorHfreq = result_dict[edge]["hfreq"]
        SensorQfreq = result_dict[edge]["qfreq"]
        SourceHfreq = all_result_dict[edge]["source"]["hfreq"]
        SourceQfreq = all_result_dict[edge]["source"]["qfreq"]
        TargetHfreq = all_result_dict[edge]["target"]["hfreq"]
        TargetQfreq = all_result_dict[edge]["target"]["qfreq"]
        SensorHTime = np.real(np.fft.ifft(SensorHfreq, len(SensorHfreq)))
        SensorQTime = np.real(np.fft.ifft(SensorQfreq, len(SensorQfreq)))
        SourceHTime = np.real(np.fft.ifft(SourceHfreq, len(SourceHfreq)))
        SourceQTime = np.real(np.fft.ifft(SourceQfreq, len(SourceQfreq)))
        TargetHTime = np.real(np.fft.ifft(TargetHfreq, len(TargetHfreq)))
        TargetQTime = np.real(np.fft.ifft(TargetQfreq, len(TargetQfreq)))
        time_name = "Pipe {0}-{1} Time".format(source, target)
        freq_name = "Pipe {0}-{1} Freq".format(source, target)
        timeHeading = np.array(
            ["Time", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor", "Flow Sensor",
             "Head Target({})".format(target), "Flow Target({})".format(target)])
        freqHeading = np.array(
            ["Freq", "Head Source({})".format(source), "Flow Source({})".format(source), "Head Sensor", "Flow Sensor",
             "Head Target({})".format(target), "Flow Target({})".format(target)])
        SaveTime = np.column_stack((time, SourceHTime, SourceQTime, SensorHTime, SensorQTime, TargetHTime, TargetQTime))
        SaveFreq = np.column_stack(
            (freq_range, abs(SourceHfreq), abs(SourceQfreq), abs(SensorHfreq), abs(SensorQfreq), abs(TargetHfreq), abs(TargetQfreq)))
        SaveTime = np.row_stack((timeHeading, SaveTime))
        SaveFreq = np.row_stack((freqHeading, SaveFreq))
        SaveDict[time_name] = SaveTime.tolist()
        SaveDict[freq_name] = SaveFreq.tolist()
    pyexcel.save_book_as(bookdict=SaveDict, dest_file_name=Output_loc)
    ProgressPage.destroy()
    print("File Saved")
