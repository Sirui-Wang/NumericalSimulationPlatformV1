# used for solving orifice equation, is commented out as in this version, we do not need the orifice equation (steady state)
"""def func(x, t, node, edges, G):
    tau = tau_array[t]
    eqn = []
    xIterable = 0
    h_placed = {"valve": [False, 0], "P": [False, 0]}
    index_recorder = dict.fromkeys(list(edges.keys()), [0, 0])
    global mapping
    mapping = dict.fromkeys(list(edges.keys()), {})
    for edge, status in edges.items():
        source, target = edge
        direction = status[0]
        a = G.edges[source, target]["wave_velocity"]
        f = G.edges[source, target]["friction_factor"]
        d = G.edges[source, target]["diameter"]
        area = (np.pi * d ** 2) / 4
        if G.edges[edge]["HasValve"]:
            Cv = G.edges[edge]["Cv"]
            ReservoirHead = G.edges[edge]["ReservoirHead"]
        dx = G.edges[source, target]["dx"]
        UA = G.edges[source, target]["U_mat"][-t, -2]
        HA = G.edges[source, target]["H_mat"][-t, -2]
        UB = G.edges[source, target]["U_mat"][-t, 1]
        HB = G.edges[source, target]["H_mat"][-t, 1]
        # print(tau)
        if len(edges) == 1:
            if direction == "in":
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    MOC = (a / 9.81) * (x[xIterable] - UA) + (x[xIterable + 1] - HA) + (
                            ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    variableBC = x[xIterable] * area - Cv * tau * np.sqrt(abs(x[xIterable + 1] - ReservoirHead))
                    eqn.append(variableBC)
                    mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                else:
                    G.edges[source, target]["H_mat"][-t - 1, -1] = G.edges[source, target]["H_mat"][-t, -1]
                    HP = G.edges[source, target]["H_mat"][-t, -1]
                    MOC = (a / 9.81) * (x[xIterable] - UA) + (HP - HA) + (
                            ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    mapping[edge] = {"UP": xIterable}
            else:
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    MOC = (a / 9.81) * (x[xIterable] - UB) - (x[xIterable + 1] - HB) + (
                            ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    variableBC = x[xIterable] * area - Cv * tau * np.sqrt(abs(ReservoirHead - x[xIterable + 1]))
                    eqn.append(variableBC)
                    mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                else:
                    G.edges[source, target]["H_mat"][-t - 1, 0] = G.edges[source, target]["H_mat"][-t, 0]
                    HP = G.edges[source, target]["H_mat"][-t, 0]
                    MOC = (a / 9.81) * (x[xIterable] - UB) - (HP - HB) + (
                            ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    mapping[edge] = {"UP": xIterable}
        else:
            if direction == "in":
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    index_recorder[edge] = [xIterable, area]
                    if h_placed["valve"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[h_placed["valve"][1]] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HU": h_placed["valve"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[xIterable + 1] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["valve"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HU": xIterable + 1}
                        xIterable += 2
                    if h_placed["P"][0]:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[h_placed["valve"][1]] - x[h_placed["P"][1]]))
                    else:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[h_placed["valve"][1]] - x[xIterable]))
                        h_placed["P"] = [True, xIterable]
                        xIterable += 1
                    eqn.append(MOC)
                    eqn.append(variableBC)
                else:
                    index_recorder[edge] = [xIterable, area]
                    if h_placed["P"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[h_placed["P"][1]] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HP": h_placed["P"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[xIterable + 1] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["P"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                        xIterable += 2
                    eqn.append(MOC)
            else:
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    index_recorder[edge] = [xIterable, -area]
                    if h_placed["valve"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[h_placed["valve"][1]] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HP": h_placed["valve"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[xIterable + 1] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["valve"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                        xIterable += 2
                    if h_placed["P"][0]:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[h_placed["P"][1]] - x[h_placed["valve"][1]]))
                    else:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[xIterable] - x[h_placed["valve"][1]]))
                        h_placed["P"] = [True, xIterable]
                        xIterable += 1
                    eqn.append(MOC)
                    eqn.append(variableBC)
                else:
                    index_recorder[edge] = [xIterable, -area]
                    if h_placed["P"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[h_placed["P"][1]] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HP": h_placed["P"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[xIterable + 1] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["P"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                        xIterable += 2
                    eqn.append(MOC)
    indexes = list(index_recorder.keys())
    if len(index_recorder) == 2:
        equiv = x[index_recorder[indexes[0]][0]] * index_recorder[indexes[0]][1] + x[index_recorder[indexes[1]][0]] * \
                index_recorder[indexes[1]][1]
        eqn.append(equiv)
    elif len(index_recorder) == 3:
        equiv = x[index_recorder[indexes[0]][0]] * index_recorder[indexes[0]][1] + x[index_recorder[indexes[1]][0]] * \
                index_recorder[indexes[1]][1] + x[index_recorder[indexes[2]][0]] * index_recorder[indexes[2]][1]
        eqn.append(equiv)
    elif len(index_recorder) == 4:
        equiv = x[index_recorder[indexes[0]][0]] * index_recorder[indexes[0]][1] + x[index_recorder[indexes[1]][0]] * \
                index_recorder[indexes[1]][1] + x[index_recorder[indexes[2]][0]] * index_recorder[indexes[2]][1] + x[
                    index_recorder[indexes[3]][0]] * index_recorder[indexes[3]][1]
        eqn.append(equiv)
    elif len(index_recorder) == 1:
        pass
    else:
        raise TooManyBranches
    return eqn
"""
"""
class TooManyBranches(Exception):
    pass


class WrongRootError:
    pass


class NoSolution():
    pass


def rootVerification(G, t, node, root, POI):
    ReservoirHead = G.edges[POI]["ReservoirHead"]
    valve_map_indexes = mapping[POI]
    valve_loc = POI.index(node)
    UP = root[valve_map_indexes["UP"]]
    HP = root[valve_map_indexes["HP"]]
    status = True
    if valve_loc == 1:
        if HP > ReservoirHead and UP < 0:
            status = False
        elif HP < ReservoirHead and UP > 0:
            status = False
    else:
        if HP > ReservoirHead and UP > 0:
            status = False
        elif HP < ReservoirHead and UP < 0:
            status = False
    return status


def rootVerification_series(root):
    status = True
    HA = root[1]
    try:
        HB = root[3]
    except IndexError:
        print("pause")
    UA = root[0]
    UB = root[2]
    if HA > HB:
        if UA < 0 or UB < 0:
            status = False
    elif HA<HB:
        if UA>0 or UB>0:
            status = False
    return status

def quadratic_root(G, edge, t):
    print("Analytical Solution at t = {}".format(t))
    UA = G.edges[edge]["U_mat"][-t, -2]
    HA = G.edges[edge]["H_mat"][-t, -2]
    ReservoirHead = G.edges[edge]["ReservoirHead"]
    l = G.edges[edge]["length"]
    d = G.edges[edge]["diameter"]
    f = G.edges[edge]["friction_factor"]
    dx = G.edges[edge]["dx"]
    a = G.edges[edge]["wave_velocity"]
    f_term = (f*UA*abs(UA))/(4*9.81*0.5*d)*dx
    Area = (np.pi*d**2)/4
    Cv = G.edges[edge]["Cv"]
    tau = G.edges[edge]["tau"]
    K = G.edges[edge]["K"]
    a_factor = (Area/(Cv*tau))**2
    b_factor = a/9.81
    c_factor = ((-a/9.81)*UA)+ReservoirHead-HA+f_term
    roots = np.zeros((4,2))
    try:
        roots[0,0] = (-b_factor+np.sqrt((b_factor**2)-4*a_factor*c_factor))/(2*a_factor)
        roots[0,1] = (Area/(tau*Cv))**2*roots[0,0]**2+ReservoirHead
    except TypeError:
        roots[0,0] = None
        roots[0,1] = None
    try:
        roots[1,0] =(-b_factor-np.sqrt((b_factor**2)-4*a_factor*c_factor))/(2*a_factor)
        roots[1,1] = (Area/(tau*Cv))**2*roots[1,0]**2+ReservoirHead
    except TypeError:
        roots[1,0] = None
        roots[1,1] = None
    try:
        roots[2,0] = (-b_factor+np.sqrt((b_factor**2)-4*(-1*a_factor)*c_factor))/(-2*a_factor)
        roots[2,1] = ReservoirHead - (Area/(tau*Cv))**2*roots[2,0]**2
    except TypeError:
        roots[2,0] = None
        roots[2,1] = None
    try:
        roots[3,0] = (-b_factor-np.sqrt((b_factor**2)-4*(-1*a_factor)*c_factor))/(-2*a_factor)
        roots[3,1] = ReservoirHead - (Area/(tau*Cv))**2*roots[3,0]**2
    except  TypeError:
        roots[3,0] = None
        roots[3,1] = None
    isCorrectRoot = False
    for root in roots:
        U = root[0]
        H = root[1]
        if U == U:
            if H > ReservoirHead and U > 0:
                correct_root = np.array([U, H])
                isCorrectRoot = True
            elif H < ReservoirHead and U < 0:
                correct_root = np.array([U, H])
                isCorrectRoot = True
    return correct_root, isCorrectRoot
"""
"""def analysis(t, decomposed_network, G):
    for node, branch in decomposed_network.items():
        hasValve = False
        for edge in branch:
            source, target = edge
            if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                hasValve = True
        if len(branch) == 1:
            if hasValve:
                size = len(branch) + 1
            else:
                size = len(branch)
        elif G.nodes[node]["classification"] == "SeriesConnection" and hasValve:
            size = len(branch) + 2
        else:
            size = len(branch) + 1
        InitGuess = np.full(size, 1)
        root = fsolve(func, InitGuess, args=(t, node, branch, G))
        global POIs
        for POI, LOI in POIs.items():
            if POI in mapping.keys() and LOI == node:
                # Verification
                ReservoirHead = G.edges[edge]["ReservoirHead"]
                try:
                    check_list = [-1, -0.1, -0.01, -0.001, -0.0001, 0, 0.0001, 0.001, 0.01, 0.1, 1, ReservoirHead, ReservoirHead - 1, ReservoirHead + 1, ReservoirHead - 10,
                                  ReservoirHead + 10, ReservoirHead - 5, ReservoirHead + 5, ReservoirHead - 50,
                                  ReservoirHead + 50, ReservoirHead - 100, ReservoirHead + 100]
                    count = 0
                    isCorrectRoot = rootVerification(G, t, node, root, POI)
                    while not isCorrectRoot:
                        if count == len(check_list):
                            root, isCorrectRoot = quadratic_root(G, edge, t)
                        else:
                            InitGuess = np.full(size, check_list[count])
                            root = fsolve(func, InitGuess, args=(t, node, branch, G))
                            isCorrectRoot = rootVerification(G, t, node, root, POI)
                            count += 1
                except TypeError:
                    trails = np.arange(0, 1000, 0.01)
                    isCorrectRoot = rootVerification_series(root)
                    count = 0
                    while not isCorrectRoot:
                        InitGuess = np.full(size, trails[count])
                        root = fsolve(func, InitGuess, args=(t, node, branch, G))
                        isCorrectRoot = rootVerification_series(root)
                        count+= 1
        for key, item in mapping.items():
            source, target = key
            loc = key.index(node)
            for key, index in item.items():
                if key == "UP":
                    G.edges[source, target]["U_mat"][-t - 1, -loc] = root[index]
                elif key == "HP":
                    G.edges[source, target]["H_mat"][-t - 1, -loc] = root[index]
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
"""

"""
def func(x, t, node, edges, G):
    eqn = []
    xIterable = 0
    h_placed = {"valve": [False, 0], "P": [False, 0]}
    index_recorder = dict.fromkeys(list(edges.keys()), [0, 0])
    global mapping
    mapping = dict.fromkeys(list(edges.keys()), {})
    for edge, status in edges.items():
        source, target = edge
        base_tau = G.edges[edge]["tau"]
        Cv = G.edges[edge]["Cv"]
        ReservoirHead = G.edges[edge]["ReservoirHead"]
        if G.edges[edge]["ValveMovement"] == "sine":
            global sineFreq, amp
            omega = 2 * np.pi * sineFreq
            tau = float(base_tau) - amp * np.sin(omega * t_range[t])
        elif G.edges[edge]["ValveMovement"] == "Instant close open":
            if t == 1:
                tau = base_tau
            else:
                tau = 1
        elif G.edges[edge]["ValveMovement"] == "constant":
            tau = base_tau
        direction = status[0]
        a = G.edges[source, target]["wave_velocity"]
        f = G.edges[source, target]["friction_factor"]
        d = G.edges[source, target]["diameter"]
        area = (np.pi * d ** 2) / 4
        dx = G.edges[source, target]["dx"]
        UA = G.edges[source, target]["U_mat"][-t, -2]
        HA = G.edges[source, target]["H_mat"][-t, -2]
        UB = G.edges[source, target]["U_mat"][-t, 1]
        HB = G.edges[source, target]["H_mat"][-t, 1]
        # print(tau)
        if len(edges) == 1:
            if direction == "in":
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    MOC = (a / 9.81) * (x[xIterable] - UA) + (x[xIterable + 1] - HA) + (
                            ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    variableBC = x[xIterable] * area - Cv * tau * np.sqrt(abs(x[xIterable + 1] - ReservoirHead))
                    eqn.append(variableBC)
                    mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                else:
                    G.edges[source, target]["H_mat"][-t - 1, -1] = G.edges[source, target]["H_mat"][-t, -1]
                    HP = G.edges[source, target]["H_mat"][-t, -1]
                    MOC = (a / 9.81) * (x[xIterable] - UA) + (HP - HA) + (
                            ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    mapping[edge] = {"UP": xIterable}
            else:
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    MOC = (a / 9.81) * (x[xIterable] - UB) - (x[xIterable + 1] - HB) + (
                            ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    variableBC = x[xIterable] * area - Cv * tau * np.sqrt(abs(ReservoirHead - x[xIterable + 1]))
                    eqn.append(variableBC)
                    mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                else:
                    G.edges[source, target]["H_mat"][-t - 1, 0] = G.edges[source, target]["H_mat"][-t, 0]
                    HP = G.edges[source, target]["H_mat"][-t, 0]
                    MOC = (a / 9.81) * (x[xIterable] - UB) - (HP - HB) + (
                            ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                    eqn.append(MOC)
                    mapping[edge] = {"UP": xIterable}
        else:
            if direction == "in":
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    index_recorder[edge] = [xIterable, area]
                    if h_placed["valve"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[h_placed["valve"][1]] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HU": h_placed["valve"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[xIterable + 1] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["valve"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HU": xIterable + 1}
                        xIterable += 2
                    if h_placed["P"][0]:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[h_placed["valve"][1]] - x[h_placed["P"][1]]))
                    else:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[h_placed["valve"][1]] - x[xIterable]))
                        h_placed["P"] = [True, xIterable]
                        xIterable += 1
                    eqn.append(MOC)
                    eqn.append(variableBC)
                else:
                    index_recorder[edge] = [xIterable, area]
                    if h_placed["P"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[h_placed["P"][1]] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HP": h_placed["P"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UA) + (x[xIterable + 1] - HA) + (
                                ((f * UA * abs(UA)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["P"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                        xIterable += 2
                    eqn.append(MOC)
            else:
                if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                    index_recorder[edge] = [xIterable, -area]
                    if h_placed["valve"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[h_placed["valve"][1]] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HP": h_placed["valve"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[xIterable + 1] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["valve"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                        xIterable += 2
                    if h_placed["P"][0]:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[h_placed["P"][1]] - x[h_placed["valve"][1]]))
                    else:
                        variableBC = x[index_recorder[edge][0]] * area - Cv * tau * np.sqrt(
                            abs(x[xIterable] - x[h_placed["valve"][1]]))
                        h_placed["P"] = [True, xIterable]
                        xIterable += 1
                    eqn.append(MOC)
                    eqn.append(variableBC)
                else:
                    index_recorder[edge] = [xIterable, -area]
                    if h_placed["P"][0]:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[h_placed["P"][1]] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        mapping[edge] = {"UP": xIterable, "HP": h_placed["P"][1]}
                        xIterable += 1
                    else:
                        MOC = (a / 9.81) * (x[xIterable] - UB) - (x[xIterable + 1] - HB) + (
                                ((f * UB * abs(UB)) / (4 * 9.81 * (d / 2))) * dx)
                        h_placed["P"] = [True, xIterable + 1]
                        mapping[edge] = {"UP": xIterable, "HP": xIterable + 1}
                        xIterable += 2
                    eqn.append(MOC)
    indexes = list(index_recorder.keys())
    if len(index_recorder) == 2:
        equiv = x[index_recorder[indexes[0]][0]] * index_recorder[indexes[0]][1] + x[index_recorder[indexes[1]][0]] * \
                index_recorder[indexes[1]][1]
        eqn.append(equiv)
    elif len(index_recorder) == 3:
        equiv = x[index_recorder[indexes[0]][0]] * index_recorder[indexes[0]][1] + x[index_recorder[indexes[1]][0]] * \
                index_recorder[indexes[1]][1] + x[index_recorder[indexes[2]][0]] * index_recorder[indexes[2]][1]
        eqn.append(equiv)
    elif len(index_recorder) == 4:
        equiv = x[index_recorder[indexes[0]][0]] * index_recorder[indexes[0]][1] + x[index_recorder[indexes[1]][0]] * \
                index_recorder[indexes[1]][1] + x[index_recorder[indexes[2]][0]] * index_recorder[indexes[2]][1] + x[
                    index_recorder[indexes[3]][0]] * index_recorder[indexes[3]][1]
        eqn.append(equiv)
    elif len(index_recorder) == 1:
        pass
    else:
        raise TooManyBranches
    return eqn


class TooManyBranches(Exception):
    pass


class NoSolution():
    pass


def rootVerification(G, t, node, root, POI):
    ReservoirHead = G.edges[POI]["ReservoirHead"]
    valve_map_indexes = mapping[POI]
    valve_loc = POI.index(node)
    UP = root[valve_map_indexes["UP"]]
    HP = root[valve_map_indexes["HP"]]
    status = True
    if valve_loc == 1:
        if HP > ReservoirHead and UP < 0:
            status = False
        elif HP < ReservoirHead and UP > 0:
            status = False
    else:
        if HP > ReservoirHead and UP > 0:
            status = False
        elif HP < ReservoirHead and UP < 0:
            status = False
    return status


def rootVerification_series(root):
    status = True
    HA = root[1]
    HB = root[3]
    UA = root[0]
    UB = root[2]
    if HA > HB:
        if UA < 0 or UB < 0:
            status = False
    elif HA < HB:
        if UA > 0 or UB > 0:
            status = False
    return status


def analysis(t, decomposed_network, G):
    for node, branch in decomposed_network.items():
        hasValve = False
        for edge in branch:
            source, target = edge
            if G.edges[source, target]["HasValve"] and G.edges[source, target]["ValveLocation"] == node:
                hasValve = True
        if len(branch) == 1:
            if hasValve:
                size = len(branch) + 1
            else:
                size = len(branch)
        elif G.nodes[node]["classification"] == "SeriesConnection" and hasValve:
            size = len(branch) + 2
        else:
            size = len(branch) + 1
        InitGuess = np.full(size, 1)
        root = fsolve(func, InitGuess, args=(t, node, branch, G))
        global POIs
        for POI, LOI in POIs.items():
            if POI in mapping.keys() and LOI == node:
                # Verification
                ReservoirHead = G.edges[edge]["ReservoirHead"]
                try:
                    check_list = [-1, -0.1, -0.01, -0.001, -0.0001, 0, 0.0001, 0.001, 0.01, 0.1, 1, ReservoirHead,
                                  ReservoirHead - 1,
                                  ReservoirHead + 1, ReservoirHead - 10,
                                  ReservoirHead + 10, ReservoirHead - 5, ReservoirHead + 5, ReservoirHead - 50,
                                  ReservoirHead + 50, ReservoirHead - 100, ReservoirHead + 100]
                    count = 0
                    isCorrectRoot = rootVerification(G, t, node, root, POI)
                    while not isCorrectRoot:
                        if count == len(check_list):
                            raise NoSolution
                        InitGuess = np.full(size, check_list[count])
                        root = fsolve(func, InitGuess, args=(t, node, branch, G))
                        isCorrectRoot = rootVerification(G, t, node, root, POI)
                        count += 1
                except TypeError:
                    trails = np.arange(0, 1000, 0.01)
                    isCorrectRoot = rootVerification_series(root)
                    count = 0
                    while not isCorrectRoot:
                        InitGuess = np.full(size, trails[count])
                        root = fsolve(func, InitGuess, args=(t, node, branch, G))
                        isCorrectRoot = rootVerification_series(root)
                        count += 1
        for key, item in mapping.items():
            source, target = key
            loc = key.index(node)
            for key, index in item.items():
                if key == "UP":
                    G.edges[source, target]["U_mat"][-t - 1, -loc] = root[index]
                elif key == "HP":
                    G.edges[source, target]["H_mat"][-t - 1, -loc] = root[index]
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
"""  # Orifice Eqn (transient analysis)
