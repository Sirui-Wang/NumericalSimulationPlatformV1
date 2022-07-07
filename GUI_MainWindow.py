import json
import os
import tkinter
from tkinter import *
from tkinter import filedialog, ttk, messagebox

import numpy as np

import GUI_Engine_Interface as Interface

# Auther: Sirui Wang
# Updated date: 07/July/2022
"""Module Comment"""
# Line Comment

global NodesDictionary, LinksDictionary, EnvirDictionary
NodesDictionary = {}  # Main storage format for Nodes(BC/Junctions) & its property
LinksDictionary = {}  # Main storage format for Links(Pipes) & its property
EnvirDictionary = {}  # Main storage format for Run Configurations

""" Sub Functions - Open, Close, Save, Refresh"""


def on_closing():
    """When closing the window, pop up message box ask if want to save file"""
    if messagebox.askyesno("Quit", message="Do you want to save before quit?"):
        SaveFile()
        root.destroy()
    else:
        root.destroy()


def OpenFile():
    """open json file and map to its corresponding entry box and canvas drawing"""
    print("Placeholder")


def SaveFile():
    """Save inputs as json file"""
    current_dir = os.getcwd()
    saveConfig()
    Save_Dict = {"Nodes": NodesDictionary, "Links": LinksDictionary, "Environment": EnvirDictionary,
                 "Mode": AnalysisMode.get()}
    extensions = [("Json Files", ".json")]
    savedFileName = filedialog.asksaveasfilename(initialdir=current_dir + "/Save", title="Save File",
                                                 defaultextension=".json",
                                                 filetypes=extensions)
    with open(savedFileName, "w") as file:
        json.dump(Save_Dict, file)
    file.close()


def Refresh_Display(*args):
    """
    This function is used to refresh main panel display on changing transient analysis mode
    switching between input field for MOC modelling and Transfer matrix modelling
    :return:
    """
    MOCFrame.pack_forget()
    TMFrame.pack_forget()
    MultiFrame.grid_forget()
    SingleFrame.grid_forget()
    if AnalysisMode.get():
        TMFrame.pack(fill=BOTH, expand=1)
        if frequencyModeVar.get() == "MultiFrequency":
            MultiFrame.grid(row=2, column=0, columnspan=4, sticky="ew", padx=20)
        else:
            SingleFrame.grid(row=2, column=0, columnspan=4, sticky="ew", padx=20)
    else:
        MOCFrame.pack(fill=BOTH, expand=1)


def Refresh_Link_Window(*args):
    """
    This function is used to refresh link pop up window on changing transient analysis mode as well as other parameters such as HasPerturbation, HasSensor, etc
    :param args:
    :return:
    """
    pass


def saveConfig():
    """Save Run configuration to dictionary for later use or overall file save"""
    global EnvirDictionary
    if AnalysisMode.get():
        if frequencyModeVar.get() == "MultiFrequency":
            EnvirDictionary = {"df": freqStepEntry.get(), "MaxFreq": FreqRangeEntry.get(),
                               "FreqMode": frequencyModeVar.get()}
        else:
            EnvirDictionary = {"ExcitationFreq": FreqEntry.get(), "FreqMode": frequencyModeVar.get()}
    else:
        EnvirDictionary = {"dt": timestepEntry.get(), "TotalTime": runtimeEntry.get(),
                           "RecordStart": record_startEntry.get(), "RecordEnd": record_endEntry.get()}


def startAnalysis():
    """pass all information to a interface between the GUI and engine to start transient analysis"""
    saveConfig()
    global Progresses
    Progresses = Toplevel()
    Progresses.title("Progresses")
    my_progress = ttk.Progressbar(Progresses, orient=HORIZONTAL, length=400, mode="determinate")
    my_progress.pack(fill=BOTH, expand=1)
    compacted_data = {"Nodes": NodesDictionary, "Links": LinksDictionary, "Environment": EnvirDictionary,
                      "Mode": AnalysisMode.get()}
    Interface.main(compacted_data, my_progress, Progresses)


"""Sub Functions - Create, edit, delete network items"""


def node_BC_field(displayWindow, occurance, IS=True):
    global BCFrame, BCVar  # declared for use in other function such as drawNodeNSave
    if occurance == "New":
        if isBCVar.get():
            BCFrame = LabelFrame(displayWindow, text="Select BC type")
            Label(BCFrame, text="BC Type").grid(row=0, column=0, sticky="w")
            BCVar = StringVar()
            BCVar.set("Constant head")
            BCSelection = OptionMenu(BCFrame, BCVar, "Constant head", "Constant flow")
            BCSelection.grid(row=0, column=1, sticky="ew")
            BCFrame.grid(row=3, column=0, columnspan=2, sticky="ew")
            BCFrame.grid_columnconfigure(0, weight=1)
        else:
            BCFrame.grid_forget()
    else:
        NodeName = occurance
        if not IS:
            pass
        elif not isBCVar.get():
            BCFrame.grid_forget()
        else:
            BCFrame = LabelFrame(displayWindow, text="Select BC type")
            Label(BCFrame, text="BC Type").grid(row=0, column=0, sticky="w")
            BCVar = StringVar()
            BCVar.set(NodesDictionary[NodeName]["BCType"])
            BCSelection = OptionMenu(BCFrame, BCVar, "Constant head", "Constant flow")
            BCSelection.grid(row=0, column=1, sticky="ew")
            BCFrame.grid(row=3, column=0, columnspan=2, sticky="ew")
            BCFrame.grid_columnconfigure(0, weight=1)


def new_node():
    global add_node
    add_node = Toplevel()
    add_node.resizable(0, 0)
    add_node.geometry(f"+{root.winfo_x() + 200}+{root.winfo_y() + 200}")
    add_node.title("Create New Node")
    Label(add_node, text="Specify node name: ").grid(row=0, column=0)
    Label(add_node, text="Specify node head (m): ").grid(row=1, column=0)
    node_name_entry = Entry(add_node)
    node_name_entry.grid(row=0, column=1)
    node_head_entry = Entry(add_node)
    node_head_entry.grid(row=1, column=1)
    Label(add_node, text="is a BC").grid(row=2, column=0)
    global isBCVar
    isBCVar = BooleanVar()
    isBCVar.set(False)
    isBCCheck = Checkbutton(add_node, variable=isBCVar, command=lambda: node_BC_field(add_node, "New"))
    isBCCheck.grid(row=2, column=1, sticky="ew")
    saveNlocation = Label(add_node, text="Save & Select node location on screen")
    saveNlocation_btn = Button(add_node, text="Save & Select",
                               command=lambda: getorigin(node_name_entry.get(), node_head_entry.get(), ))
    saveNlocation.grid(row=4, column=0, columnspan=2)
    saveNlocation_btn.grid(row=5, column=0, columnspan=2)


def getorigin(NodeName, NodeHead):
    add_node.destroy()
    drawArea.bind("<Button 1>", lambda self: drawNodeNSave(self, NodeName, NodeHead))


def drawNodeNSave(eventorigin, NodeName, NodeHead):
    x = eventorigin.x
    y = eventorigin.y
    NodesDictionary[NodeName] = {"head": NodeHead}
    NodesDictionary[NodeName].update({"Coord": (x, y)})
    if isBCVar.get():
        NodesDictionary[NodeName].update({"isBC": True, "BCType": BCVar.get()})
    else:
        NodesDictionary[NodeName].update({"isBC": False, "BCType": None})
    drawArea.create_oval(x - 10, y - 10, x + 10, y + 10, fill="blue", tag=NodeName)
    drawArea.create_text(x, y + 30, fill="blue", text=NodeName, tag=NodeName + "text")
    drawArea.unbind("<Button 1>")
    edit_node_menu.add_command(label=NodeName, command=lambda: editNode(NodeName))


def editNode(key):
    global editNodePage
    editNodePage = Toplevel()
    editNodePage.resizable(0, 0)
    editNodePage.geometry(f"+{root.winfo_x() + 200}+{root.winfo_y() + 200}")
    editNodePage.title(key)
    node_name_label = Label(editNodePage, text="Node Name: " + key)
    node_name_label.grid(row=0, column=0, columnspan=2)
    node_head_label = Label(editNodePage, text="Node Head: ")
    node_head_label.grid(row=1, column=0)
    node_head_entry = Entry(editNodePage)
    node_head_entry.grid(row=1, column=1)
    node_head_entry.insert(0, NodesDictionary[key]["head"])
    Label(editNodePage, text="is a BC").grid(row=2, column=0)
    global isBCVar
    node_BC_field(editNodePage, key, IS=NodesDictionary[key]["isBC"])
    isBCVar.set(bool(NodesDictionary[key]["isBC"]))
    isBCCheck = Checkbutton(editNodePage, variable=isBCVar, command=lambda: node_BC_field(editNodePage, key))
    isBCCheck.grid(row=2, column=1, sticky="ew")
    saveNodePropertyBtn = Button(editNodePage, text="Save Changes",
                                 command=lambda: saveNodeEdit(key, node_head_entry.get()))
    saveNodePropertyBtn.grid(row=6, column=0, columnspan=2)


def saveNodeEdit(key, node_head_entry):
    NodesDictionary[key]["head"] = node_head_entry
    if isBCVar.get():
        NodesDictionary[key].update({"isBC": True, "BCType": BCVar.get()})
    else:
        NodesDictionary[key].update({"isBC": False, "BCType": None})
    editNodePage.destroy()


def remove_node():
    global delete_node
    delete_node = Toplevel()
    delete_node.resizable(0, 0)
    delete_node.geometry(f"+{root.winfo_x() + 200}+{root.winfo_y() + 200}")
    delete_node.title("Delete specified node")
    node_name_label = Label(delete_node, text="Specify node name to delete: ")
    node_name_entry = Entry(delete_node)
    node_name_label.grid(row=0, column=0)
    node_name_entry.grid(row=0, column=1)
    delete_btn = Button(delete_node, text="Delete", command=lambda: destroy_node(node_name_entry.get()))
    delete_btn.grid(row=1, column=0, columnspan=2)


def destroy_node(NodeName):
    drawArea.delete(NodeName)
    drawArea.delete(NodeName + "text")
    edit_node_menu.delete(NodeName)
    del NodesDictionary[NodeName]
    delete_node.destroy()


def link_MOC_property_field(displayWindow):
    global MOCAdditionFrame
    MOCAdditionFrame = LabelFrame(displayWindow, text="Additional Boundary Condition", highlightthickness=2)
    Label(MOCAdditionFrame, text="Has perturbation").grid(row=0, column=0, sticky="w")
    Label(MOCAdditionFrame, text="NOTE: This only applies to boundary nodes").grid(row=1, column=0, columnspan=2,
                                                                                   sticky="ew")
    global MOCPerturbationTypeVar
    MOCPerturbationTypeVar = StringVar()
    MOCPerturbationTypeVar.set("None")
    PertTypeList = ["None", "Impulse", "Full Closure", "Sinusoidal Head", "Sinusoidal Flow", "Controlled Flow"]
    PerturbationSelection = OptionMenu(MOCAdditionFrame, MOCPerturbationTypeVar, *PertTypeList,
                                       command=lambda v: MOC_perturbation_field((displayWindow, v)))
    PerturbationSelection.grid(row=0, column=1, sticky="ew")
    MOCAdditionFrame.grid(row=2, sticky="nsew")
    MOCAdditionFrame.grid_columnconfigure(0, weight=1)


def MOC_perturbation_field(passed):
    displayWindow, PerturbationType = passed
    global PertFrame, Location
    try:
        PertFrame.destroy()
    except NameError:
        pass
    PertFrame = LabelFrame(displayWindow, text=PerturbationType, highlightthickness=0)
    Label(PertFrame, text="Location").grid(row=0, column=0, sticky="w")  # perturbation location label
    Location = Entry(PertFrame)
    Location.grid(row=0, column=1, sticky="e")
    if PerturbationType == "Impulse":
        ImpulseFrame = LabelFrame(PertFrame, highlightthickness=0, borderwidth=0)
        Label(ImpulseFrame, text="Time").grid(row=0, column=0, sticky="w")  # impulse time
        StartTime = Entry(ImpulseFrame)
        StartTime.grid(row=0, column=1, sticky="w")  # impulse time entry
        ImpulseFrame.grid(row=1, column=0, columnspan=2, sticky="nsew")
        ImpulseFrame.grid_columnconfigure(0, weight=1)
        PertFrame.grid(row=10, column=0, sticky="ew")
        PertFrame.grid_columnconfigure(0, weight=1)
    elif PerturbationType == "Full Closure":
        global StartTime
        FCFrame = LabelFrame(PertFrame, highlightthickness=0, borderwidth=0)
        Label(FCFrame, text="Time").grid(row=0, column=0, sticky="w")  # impulse time
        StartTime = Entry(FCFrame)
        StartTime.grid(row=0, column=1, sticky="w")  # impulse time entry
        FCFrame.grid(row=1, column=0, columnspan=2, sticky="ew")
        FCFrame.grid_columnconfigure(0, weight=1)
        PertFrame.grid(row=10, column=0, sticky="ew")
        PertFrame.grid_columnconfigure(0, weight=1)
    elif PerturbationType == "Controlled Flow":
        global StartTime, FlowRate
        CFFrame = LabelFrame(PertFrame, highlightthickness=0, borderwidth=0)
        Label(CFFrame, text="Time").grid(row=0, column=0, sticky="w")  # impulse time
        StartTime = Entry(CFFrame)
        StartTime.grid(row=0, column=1, sticky="e")  # impulse time entry
        Label(CFFrame, text="Flow rate").grid(row=1, column=0, sticky="w")
        FlowRate = Entry(CFFrame)
        FlowRate.grid(row=1, column=1, sticky="e")
        CFFrame.grid(row=1, column=0, columnspan=2, sticky="ew")
        CFFrame.grid_columnconfigure(0, weight=1)
        PertFrame.grid(row=10, column=0, sticky="ew")
        PertFrame.grid_columnconfigure(0, weight=1)
    elif PerturbationType == "None":
        pass
    else:
        global Frequency, Amplitude
        SineFrame = LabelFrame(PertFrame, highlightthickness=2, borderwidth=0)
        Label(SineFrame, text="Frequency").grid(row=0, column=0, sticky="w")
        Label(SineFrame, text="Amplitude").grid(row=1, column=0, sticky="w")
        Frequency = Entry(SineFrame)
        Amplitude = Entry(SineFrame)
        Frequency.grid(row=0, column=1, sticky="e")
        Amplitude.grid(row=1, column=1, sticky="e")
        SineFrame.grid(row=2, column=0, columnspan=2, sticky="ew")
        SineFrame.grid_columnconfigure(0, weight=1)
        PertFrame.grid(row=10, column=0, sticky="ew")
        PertFrame.grid_columnconfigure(0, weight=1)


def link_TM_property_field(displayWindow):
    global TMAdditionFrame
    TMAdditionFrame = LabelFrame(displayWindow, text="Additional Boundary Condition", highlightthickness=2)
    Label(TMAdditionFrame, text="Has perturbation").grid(row=0, column=0, sticky="w")
    Label(TMAdditionFrame, text="NOTE: This only applies to boundary nodes").grid(row=1, column=0, columnspan=2,
                                                                                  sticky="ew")
    Label(TMAdditionFrame, text="Has Sensor").grid(row=2, column=0, sticky="w")
    # AdditionFrame.pack(fill=BOTH)
    global TMPerturbationTypeVar
    TMPerturbationTypeVar = StringVar()
    TMPerturbationTypeVar.set("None")
    PerturbationSelection = OptionMenu(TMAdditionFrame, TMPerturbationTypeVar, "None", "Flow", "Head",
                                       command=lambda v: TM_perturbation_field((displayWindow, v)))
    PerturbationSelection.grid(row=0, column=1, sticky="ew")
    SensorVar = BooleanVar()
    SensorVar.set(False)
    HasSensor = Checkbutton(TMAdditionFrame, variable=SensorVar,
                            command=lambda: TM_sensor_field(displayWindow, SensorVar.get()))
    HasSensor.grid(row=2, column=1, sticky="ew")
    TMAdditionFrame.grid(row=2, sticky="nsew")
    TMAdditionFrame.grid_columnconfigure(0, weight=1)


def TM_perturbation_field(passed):
    displayWindow, PerturbationType = passed
    global PertFrame
    try:
        PertFrame.destroy()
    except NameError:
        pass
    PertFrame = LabelFrame(displayWindow, text=PerturbationType, highlightthickness=0)
    Label(PertFrame, text="Location").grid(row=0, column=0, sticky="w")  # perturbation location label
    Location = Entry(PertFrame)
    Location.grid(row=0, column=1, sticky="e")
    if PerturbationType == "Head":
        PertFrame.grid(row=10, column=0, sticky="ew")
        PertFrame.grid_columnconfigure(0, weight=1)
    elif PerturbationType == "Flow":
        PertFrame.grid(row=10, column=0, sticky="ew")
        PertFrame.grid_columnconfigure(0, weight=1)
    else:
        pass


def TM_sensor_field(displayWindow, SensorVar):
    global SensorFrame, SensorLocation
    try:
        SensorFrame.destroy()
    except NameError:
        pass
    SensorFrame = LabelFrame(displayWindow, highlightthickness=0)
    Label(SensorFrame, text="Sensor Location").grid(row=0, column=0, sticky="w")  # Sensor location label
    SensorLocation = Entry(SensorFrame)
    SensorLocation.grid(row=0, column=1, sticky="e")
    if SensorVar:
        SensorFrame.grid(row=15, sticky="ew")
        SensorFrame.grid_columnconfigure(0, weight=1)
    else:
        pass


# noinspection PyGlobalUndefined
def link_Basic_property_field(displayWindow):
    BasicFrame = LabelFrame(displayWindow, text="Basic Properties", highlightthickness=2)
    Label(BasicFrame, text="Specify Pipe Length in m").grid(row=0, column=0, sticky="w")  # pipe length label
    Label(BasicFrame, text="Specify Pipe Diameter in m").grid(row=1, column=0, sticky="w")  # pipe diameter label
    Label(BasicFrame, text="Specify Wave Velocity in m/s").grid(row=2, column=0, sticky="w")  # wave speed label
    Label(BasicFrame, text="Specify Friction Factor").grid(row=3, column=0, sticky="w")  # friction factor label
    Label(BasicFrame, text="Specify Flow Velocity").grid(row=4, column=0, sticky="w")  # flow velocity label
    # BasicFrame.pack(fill=BOTH,expand=True)
    BasicFrame.grid(row=1, sticky="nsew")
    BasicFrame.grid_columnconfigure(0, weight=1)
    global l, d, a, f, u
    l = Entry(BasicFrame, justify="right")
    d = Entry(BasicFrame)
    a = Entry(BasicFrame)
    f = Entry(BasicFrame)
    u = Entry(BasicFrame)
    l.grid(row=0, column=1, sticky="ew")
    d.grid(row=1, column=1, sticky="e")
    a.grid(row=2, column=1, sticky="e")
    f.grid(row=3, column=1, sticky="e")
    u.grid(row=4, column=1, sticky="e")


def new_link():
    global add_link
    add_link = Toplevel()
    add_link.resizable(0, 0)
    add_link.geometry(f"+{root.winfo_x() + 200}+{root.winfo_y() + 200}")
    add_link.title("Create New Link")
    ConnectionFrame = LabelFrame(add_link, text="Connection", highlightthickness=2)
    Label(ConnectionFrame, text="Specify Starting Node (upstream node)").grid(row=0, column=0,
                                                                              sticky="w")  # staring node label
    Label(ConnectionFrame, text="Specify Ending Node (downstream node)").grid(row=1, column=0,
                                                                              sticky="w")  # ending node label
    global startNode, endNode
    startNode = Entry(ConnectionFrame)
    endNode = Entry(ConnectionFrame)
    startNode.grid(row=0, column=1, sticky="w")
    endNode.grid(row=1, column=1, sticky="w")
    # ConnectionFrame.pack(fill=BOTH,expand=True)
    ConnectionFrame.grid(row=0, sticky="nsew")
    ConnectionFrame.grid_columnconfigure(0, weight=1)
    link_Basic_property_field(add_link)
    if AnalysisMode.get():
        link_TM_property_field(add_link)
    else:
        link_MOC_property_field(add_link)
    save_btn = Button(add_link, text="Save",
                      command=lambda: drawLink())
    save_btn.grid(row=20, column=0, columnspan=2)


def drawLink():
    x0, y0 = NodesDictionary[startNode]["Coord"]
    x1, y1 = NodesDictionary[endNode]["Coord"]
    drawArea.create_line(x0, y0, x1, y1, tag="{},{}".format(startNode, endNode))
    LinksDictionary["{},{}".format(startNode, endNode)] = {"l": l, "d": d, "a": a, "f": f, "U": u}
    if AnalysisMode.get():
        LinksDictionary["{},{}".format(startNode, endNode)].update({"PertType":TMPerturbationTypeVar.get(), "SensorLocation":SensorLocation.get()})
    else:
        LinksDictionary["{},{}".format(startNode, endNode)].update({"Location":Location.get()})
        if MOCPerturbationTypeVar.get()=="Impulse":
            LinksDictionary["{},{}".format(startNode, endNode)].update({"Location": Location.get()})
        elif MOCPerturbationTypeVar.get()=="Full Closure":
            LinksDictionary["{},{}".format(startNode, endNode)].update({"Location": Location.get()})
        elif MOCPerturbationTypeVar.get()=="Controlled Flow":
            LinksDictionary["{},{}".format(startNode, endNode)].update({"Location": Location.get()})
        elif MOCPerturbationTypeVar.get()=="None":
            LinksDictionary["{},{}".format(startNode, endNode)].update({"Location": Location.get()})
        else:
            LinksDictionary["{},{}".format(startNode, endNode)].update({"Location": Location.get()})

    x = (x0 + x1) / 2
    y = (y0 + y1) / 2
    drawArea.create_text(x, y + 30, text=l + "m", tag="{},{}".format(startNode, endNode) + "text")
    edit_link_menu.add_command(label=("{},{}".format(startNode, endNode)), command=lambda: editLink("{},{}".format(startNode, endNode)))
    add_link.destroy()


def remove_link():
    print("Placeholder")


def ViewLink():
    print("Placeholder")


def ViewNode():
    print("Placeholder")


def clear():
    global NodesDictionary, LinksDictionary, EnvirDictionary
    for key in NodesDictionary.keys():
        edit_node_menu.delete(key)
    for key in LinksDictionary.keys():
        key = key.split(",")
        start, end = key
        edit_link_menu.delete("{},{}".format(start, end))
    drawArea.delete("all")
    NodesDictionary = {}
    LinksDictionary = {}
    EnvirDictionary = {}


""" Main Window """
root = Tk()
root.title("Transient Analysis")
root.geometry("1080x768")  # Define initial window size
root.protocol("WM_DELETE_WINDOW", on_closing)
mainMenu = Menu(root)
root.config(menu=mainMenu)

""" Define main menu """
mainMenu = Menu(root)
root.config(menu=mainMenu)

""" Create a menu item """
# file menu
file_menu = Menu(mainMenu, tearoff=False)
mainMenu.add_cascade(label="File", menu=file_menu)
file_menu.add_command(label="Open...", command=OpenFile)
file_menu.add_command(label="Save", command=SaveFile)
# view menu
view_menu = Menu(mainMenu, tearoff=False)
mainMenu.add_cascade(label="View", menu=view_menu)
view_menu.add_command(label="Links", command=ViewLink)
view_menu.add_command(label="Nodes", command=ViewNode)
# option menu
option_menu = Menu(mainMenu, tearoff=False)
mainMenu.add_cascade(label="Options", menu=option_menu)
method_menu = Menu(option_menu, tearoff=False)
option_menu.add_cascade(label="Numerical Method", menu=method_menu)
global AnalysisMode  # declare global because it will be used in Function: RefreshDisplay, SaveFile, OpenFile, saveConfig, startAnalysis
AnalysisMode = BooleanVar()  # Define boolean variable to determine which mode is active, since there is only two mode, boolean variable is used
AnalysisMode.set(True)  # Default mode Transfer Matrix Analysis
method_menu.add_checkbutton(label="TransferMatrix", onvalue=1, offvalue=0, variable=AnalysisMode,
                            command=Refresh_Display)
method_menu.add_checkbutton(label="Method of Characteristic", onvalue=0, offvalue=1, variable=AnalysisMode,
                            command=Refresh_Display)
# define geometry menu
geometry_menu = Menu(mainMenu, tearoff=False)
mainMenu.add_cascade(label="Geometry", menu=geometry_menu)
geometry_menu.add_command(label="add Node", command=new_node)
geometry_menu.add_command(label="add Link", command=new_link)
geometry_menu.add_command(label="Remove Node", command=remove_node)
geometry_menu.add_command(label="Remove Link", command=remove_link)
geometry_menu.add_command(label="Clear all", command=clear)
# edit geometry menu
global edit_node_menu, edit_link_menu, edit_menu  # declare global because it will be used in Function: add/destroy node/link, refreshDisplay
edit_menu = Menu(mainMenu, tearoff=False)
mainMenu.add_cascade(label="Edit", menu=edit_menu)
edit_node_menu = Menu(edit_menu, tearoff=0)
edit_link_menu = Menu(edit_menu, tearoff=0)
edit_menu.add_cascade(label="Edit Node", menu=edit_node_menu)
edit_menu.add_cascade(label="Edit Link", menu=edit_link_menu)

""" Define panel space """
# Define Paned Windows
pw = PanedWindow(bd=2, relief="flat", orient=tkinter.HORIZONTAL)
LeftFrame = LabelFrame(root, text="Draw network here")
RightFrame = LabelFrame(root, text="Change run configuration here")
pw.add(LeftFrame, stretch="always")
pw.add(RightFrame, width=400, stretch="never")
pw.pack(fill=BOTH, expand=True)
global MOCFrame, TMFrame  # declared global because it will be used in Function: Refresh
MOCFrame = LabelFrame(RightFrame, borderwidth=0, highlightthickness=0)
TMFrame = LabelFrame(RightFrame, borderwidth=0, highlightthickness=0)

"""Define MOC Frame layout"""
# Define entry boxes
Label(MOCFrame, text="MOC Run Configuration", font=("Helvetica", "24")).grid(row=0, column=0, columnspan=4, pady=10,
                                                                             sticky="W")  # Title
Label(MOCFrame, text="dt = ", anchor="e").grid(row=1, column=0, pady=10, sticky="W")  # dt Label
Label(MOCFrame, text="second").grid(row=1, column=3, pady=10)  # unit Label
Label(MOCFrame, text="Total Run Time").grid(row=2, column=0, pady=10, sticky="W")  # total time Label
Label(MOCFrame, text="second").grid(row=2, column=3, pady=10)  # unit Label
Label(MOCFrame, text="Record data from").grid(row=3, column=0, pady=10, sticky="W")  # record range Label
Label(MOCFrame, text="to").grid(row=3, column=2, pady=10)  # "to" Label
Label(MOCFrame, text="Start/second").grid(row=4, column=1)  # unit Label
Label(MOCFrame, text="End/second").grid(row=4, column=3)  # unit Label
timestepEntry = Entry(MOCFrame, width=10)
timestepEntry.grid(row=1, column=1, columnspan=2, pady=10, padx=0, sticky="WE")
runtimeEntry = Entry(MOCFrame, width=10)
runtimeEntry.grid(row=2, column=1, columnspan=2, pady=10, ipadx=0, sticky="WE")
record_startEntry = Entry(MOCFrame, width=10)
record_endEntry = Entry(MOCFrame, width=10)
record_startEntry.grid(row=3, column=1, pady=(10, 0), ipadx=0, sticky="WE")  # pady can take a tuple (top,bottom)
record_endEntry.grid(row=3, column=3, pady=(10, 0), ipadx=0, sticky="WE")

EmptyLabel = Label(MOCFrame, text="")
EmptyLabel.grid(row=5, column=0,
                columnspan=4)  # Left sufficient place to add future required entry field, row starting from 5 ends at 19.

# MOC save and start analysis
saveconfigBtn = Button(MOCFrame, text="Save Configuration", command=saveConfig)
saveconfigBtn.grid(row=20, column=1, columnspan=2, pady=(50, 0))
StartBtn = Button(MOCFrame, text="Start Analysis", command=startAnalysis)
StartBtn.grid(row=21, column=1, columnspan=2, pady=25)

"""Define TransferMatrx Frame layout"""
Label(TMFrame, text="TM Run Configuration", font=("Helvetica", "24")).grid(row=0, column=0, columnspan=4, pady=10,
                                                                           sticky="ew")  # Title
SelectionFrame = LabelFrame(TMFrame, borderwidth=0, highlightthickness=0)
Label(SelectionFrame, text="FrequencyMode").grid(row=0, column=0)
frequencyModeVar = StringVar()
frequencyModeVar.set("MultiFrequency")
FrequencyModeSelection = OptionMenu(SelectionFrame, frequencyModeVar, "MultiFrequency", "SingleFrequency")
FrequencyModeSelection.grid(row=0, column=1, pady=10, padx=10, sticky="ew")
SelectionFrame.grid(row=1, column=0, columnspan=4, padx=20, sticky="ew")

# Create Frame for multi-frequency analysis
MultiFrame = LabelFrame(TMFrame, borderwidth=0, highlightthickness=0)
Label(MultiFrame, text="dfreq = ", anchor="e").grid(row=0, column=0, pady=10, sticky="W")
Label(MultiFrame, text="Max Freq").grid(row=1, column=0, pady=10, sticky="W")
freqStepEntry = Entry(MultiFrame, width=40)
freqStepEntry.grid(row=0, column=1, pady=10, padx=10, columnspan=2, sticky="WE")
FreqRangeEntry = Entry(MultiFrame, width=40)
FreqRangeEntry.grid(row=1, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")

# Create Frame for single frequency analysis
SingleFrame = LabelFrame(TMFrame, borderwidth=0, highlightthickness=0)
Label(SingleFrame, text="Excitation Frequency").grid(row=0, column=0, pady=10)
FreqEntry = Entry(SingleFrame, width=30)
FreqEntry.grid(row=0, column=1, pady=10, columnspan=2, padx=10, sticky="WE")
frequencyModeVar.trace("w", Refresh_Display)

# TM save and start analysis
saveconfigBtn = Button(TMFrame, text="Save Configuration", command=saveConfig)
saveconfigBtn.grid(row=10, column=1, columnspan=2, pady=(50, 0))
StartBtn = Button(TMFrame, text="Start Analysis", command=startAnalysis)
StartBtn.grid(row=11, column=1, columnspan=2, pady=50)

""" Display selected analysis mode entry boxes on screen """
AnalysisMode.set(False)
if AnalysisMode.get():
    TMFrame.pack(fill=BOTH, expand=1)
    if frequencyModeVar.get() == "MultiFrequency":
        MultiFrame.grid(row=2, column=0, columnspan=4, sticky="ew", padx=20)
    else:
        SingleFrame.grid(row=2, column=0, columnspan=4, sticky="ew", padx=20)
else:
    MOCFrame.pack(fill=BOTH, expand=1)

# define canvas (draw area)
drawArea = Canvas(LeftFrame, bg="grey")
drawArea.pack(fill=BOTH, expand=1)

root.mainloop()
