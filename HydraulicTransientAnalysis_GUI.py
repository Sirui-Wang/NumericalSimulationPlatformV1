import json
import os
import tkinter
from ctypes import windll
from tkinter import *
from tkinter import filedialog, ttk, messagebox

import convertion.GUI_Engine_Interface as Interface

# Auther: Sirui Wang
# Updated date: 13/July/2022
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
    if messagebox.askyesno("Open New", message="Do you want to save current layout before open"):
        SaveFile()
    clear(Warning=False)
    current_dir = os.getcwd()
    filename = filedialog.askopenfilename(initialdir=current_dir + "/Save", title="Select a File",
                                          filetypes=(("Json Files", "*.json"), ("All Files", "*.*")))
    with open(filename, "r") as file:
        recovered_dict = json.load(file)
    file.close()
    global NodesDictionary, LinksDictionary, EnvirDictionary
    NodesDictionary = recovered_dict["Nodes"]
    LinksDictionary = recovered_dict["Links"]
    EnvirDictionary = recovered_dict["Environment"]
    for node in NodesDictionary.keys():
        """Read saved file and create node"""
        x, y = NodesDictionary[node]["Coord"]
        drawArea.create_oval(x - 10, y - 10, x + 10, y + 10, fill="blue", tag=node)
        drawArea.create_text(x, y + 30, fill="blue", text=node, tag=node + "text")
        global isBCVar
        isBCVar = BooleanVar()
        isBCVar.set(NodesDictionary[node]["isBC"])
        # print(isBCVar.get())
        edit_node_menu.add_command(label=node, command=lambda n=node: editNode(n))
    for edge in LinksDictionary.keys():
        """Read saved file and create link"""
        edgestr = edge.split(",")
        startNode = edgestr[0]
        endNode = edgestr[1]
        l = LinksDictionary[edge]["Length"]
        x0, y0 = NodesDictionary[startNode]["Coord"]
        x1, y1 = NodesDictionary[endNode]["Coord"]
        drawArea.create_line(x0, y0, x1, y1, tag="{},{}".format(startNode, endNode))
        x = (x0 + x1) / 2
        y = (y0 + y1) / 2
        drawArea.create_text(x, y + 30, text=l + "m", tag="{},{}".format(startNode, endNode) + "text")
        edit_link_menu.add_command(label=("{},{}".format(startNode, endNode)),
                                   command=lambda e="{},{}".format(startNode, endNode): editLink(e))
    AnalysisMode.set(recovered_dict["Mode"])
    Refresh_display()
    if not recovered_dict["Mode"]:
        """ (MOC) Read saved file in MOC mode and store relevant information in dictionary"""
        timestepEntry.delete(0, END)
        timestepEntry.insert(END, EnvirDictionary["dt"])
        runtimeEntry.delete(0, END)
        runtimeEntry.insert(END, EnvirDictionary["TotalTime"])
        record_startEntry.delete(0, END)
        record_startEntry.insert(END, EnvirDictionary["RecordStart"])
        record_endEntry.delete(0, END)
        record_endEntry.insert(END, EnvirDictionary["RecordEnd"])
    else:
        """ (TM) Read saved file in TM mode and store relevant information in dictionary"""
        global frequencyModeVar
        frequencyModeVar.set(EnvirDictionary["FreqMode"])
        Refresh_Panel()
        if EnvirDictionary["FreqMode"] == "MultiFrequency":
            """Multi Frequency"""
            MultifreqStepEntry.delete(0, END)
            MultifreqStepEntry.insert(END, EnvirDictionary["df"])
            MultiFreqRangeEntry.delete(0, END)
            MultiFreqRangeEntry.insert(END, EnvirDictionary["MaxFreq"])
        elif EnvirDictionary["FreqMode"] == "SingleFrequency":
            """Single Frequency"""
            FreqEntry.delete(0, END)
            FreqEntry.insert(END, EnvirDictionary["ExcitationFreq"])
            SinglefreqStepEntry.delete(0, END)
            SinglefreqStepEntry.insert(END, EnvirDictionary["df"])
            SingleFreqRangeEntry.delete(0, END)
            SingleFreqRangeEntry.insert(END, EnvirDictionary["MaxFreq"])
        else:
            """Randomized Noise"""
            NoiseFreqStepEntry.delete(0, END)
            NoiseFreqStepEntry.insert(END, EnvirDictionary["df"])
            NoiseFreqRangeEntry.delete(0, END)
            NoiseFreqRangeEntry.insert(END, EnvirDictionary["MaxFreq"])
            SimSizeEntry.delete(0, END)
            SimSizeEntry.insert(END, EnvirDictionary["SimSize"])
            MeanNoiseEntry.delete(0,END)
            MeanNoiseEntry.insert(END, EnvirDictionary["NoiseMean"])
            NoiseSTDEntry.delete(0,END)
            NoiseSTDEntry.insert(END, EnvirDictionary["NoiseSTD"])
            Sensor1Entry.delete(0,END)
            Sensor1Entry.insert(END, EnvirDictionary["Sensor1"])
            Sensor2Entry.delete(0,END)
            Sensor2Entry.insert(END, EnvirDictionary["Sensor2"])


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


def Refresh_Panel(*args):
    """ Display the correct Environment panel based on environment frequency mode (single frequency / multi frequency)"""
    MultiFrame.grid_forget()
    SingleFrame.grid_forget()
    NoiseFrame.grid_forget()
    if frequencyModeVar.get() == "MultiFrequency":
        MultiFrame.grid(row=2, column=0, columnspan=4, sticky="ew", padx=20)
    elif frequencyModeVar.get() == "SingleFrequency":
        SingleFrame.grid(row=2, column=0, columnspan=4, sticky="ew", padx=20)
    else:
        NoiseFrame.grid(row=2, column=0, columnspan=4, sticky="ew", padx=20)


def Refresh_display():
    """ Display the correct Environment panel based on Transient analysis mode ( MOC / TM )"""
    MOCFrame.pack_forget()
    TMFrame.pack_forget()
    Refresh_Panel()
    if AnalysisMode.get():
        TMFrame.pack(fill=BOTH, expand=1)
    else:
        MOCFrame.pack(fill=BOTH, expand=1)


def ChangeMode(Warning=True, *args):
    """ This function is to make sure that every time analysis mode is changed, the backend dictionary is refreshed
    This is because if mode is changed after a network is created, it could cause error because some field do not exist"""
    Answer = clear(Warning, flip=True)
    if Answer:
        Refresh_display()


def saveConfig():
    """Save Run configuration to dictionary for later use or overall file save"""
    global EnvirDictionary
    if AnalysisMode.get():
        if frequencyModeVar.get() == "MultiFrequency":
            EnvirDictionary = {}
            EnvirDictionary = {"df": MultifreqStepEntry.get(), "MaxFreq": MultiFreqRangeEntry.get(),
                               "FreqMode": frequencyModeVar.get()}
        elif frequencyModeVar.get() == "SingleFrequency":
            EnvirDictionary = {}
            EnvirDictionary = {"ExcitationFreq": FreqEntry.get(), "df": SinglefreqStepEntry.get(),
                               "MaxFreq": SingleFreqRangeEntry.get(), "FreqMode": frequencyModeVar.get()}
        else:
            EnvirDictionary = {"df": NoiseFreqStepEntry.get(), "MaxFreq": NoiseFreqRangeEntry.get(),
                               "SimSize":SimSizeEntry.get(), "NoiseMean":MeanNoiseEntry.get(),
                               "NoiseSTD":NoiseSTDEntry.get(), "Sensor1":Sensor1Entry.get(),
                               "Sensor2":Sensor2Entry.get(),
                               "FreqMode": frequencyModeVar.get()}
    else:
        EnvirDictionary = {"dt": timestepEntry.get(), "TotalTime": runtimeEntry.get(),
                           "RecordStart": record_startEntry.get(), "RecordEnd": record_endEntry.get()}


def startAnalysis():
    """pass all information to a convertion between the GUI and engine to start transient analysis"""
    saveConfig()
    global SubProgresses, MainProgresses
    Progresses = LabelFrame(RightFrame, text="Progresses", height=50)
    SubProgresses = LabelFrame(Progresses, text="SubProgress", height=20)
    MainProgresses = LabelFrame(Progresses,  text="MainProgress", height=20)
    SubProgressBar = ttk.Progressbar(SubProgresses, orient=HORIZONTAL, length=400, mode="determinate")
    SubProgressBar.pack(fill=BOTH, expand=1)
    MainProgressBar = ttk.Progressbar(MainProgresses, orient=HORIZONTAL, length=400, mode="determinate")
    MainProgressBar.pack(fill=BOTH, expand=1)
    SubProgresses.pack(fill=BOTH)
    MainProgresses.pack(fill=BOTH)
    Progresses.pack(fill=BOTH, expand=1)
    compacted_data = {"Nodes": NodesDictionary, "Links": LinksDictionary, "Environment": EnvirDictionary,
                      "Mode": AnalysisMode.get()}
    Interface.main(compacted_data, SubProgressBar, MainProgressBar, Progresses)


"""Sub Functions - Create, edit, delete network items"""


def node_BC_field(displayWindow, occurance, IS=True):
    """ Modular function to display node Boundary condition, is called by both new node and edit existing node"""
    global BCFrame, BCVar  # declared for use in other function such as drawNodeNSave
    if occurance == "New":
        if isBCVar.get():
            BCFrame = LabelFrame(displayWindow, text="Select BC type")
            Label(BCFrame, text="BC Type").grid(row=0, column=0, sticky="w")
            BCVar = StringVar()
            BCVar.set("Constant head")
            BCSelection = OptionMenu(BCFrame, BCVar, "Constant head", "Constant flow", "Infinite Non-Reflecting")
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
            BCSelection = OptionMenu(BCFrame, BCVar, "Constant head", "Constant flow","Infinite Non-Reflecting")
            BCSelection.grid(row=0, column=1, sticky="ew")
            BCFrame.grid(row=3, column=0, columnspan=2, sticky="ew")
            BCFrame.grid_columnconfigure(0, weight=1)


def new_node():
    """Gui where new node is created"""
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
    """Function is needed as a mid step between save_node button, and create the node on canvas"""
    add_node.destroy()
    drawArea.bind("<Button 1>", lambda self: drawNodeNSave(self, NodeName, NodeHead))


def drawNodeNSave(eventorigin, NodeName, NodeHead):
    """The final step for a new node is created, actually create the node on canvas, and store all field information to background dictionary"""
    x = eventorigin.x
    y = eventorigin.y
    NodesDictionary[NodeName] = {"head": NodeHead}
    NodesDictionary[NodeName].update({"Coord": (x, y)})
    if isBCVar.get():
        NodesDictionary[NodeName].update({"isBC": isBCVar.get(), "BCType": BCVar.get()})
    else:
        NodesDictionary[NodeName].update({"isBC": isBCVar.get(), "BCType": None})
    drawArea.create_oval(x - 10, y - 10, x + 10, y + 10, fill="blue", tag=NodeName)
    drawArea.create_text(x, y + 30, fill="blue", text=NodeName, tag=NodeName + "text")
    drawArea.unbind("<Button 1>")
    edit_node_menu.add_command(label=NodeName, command=lambda: editNode(NodeName))


def editNode(key):
    """Front end gui to edit node"""
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
    isBCVar.set(bool(NodesDictionary[key]["isBC"]))
    node_BC_field(editNodePage, key, IS=NodesDictionary[key]["isBC"])
    isBCCheck = Checkbutton(editNodePage, variable=isBCVar, command=lambda: node_BC_field(editNodePage, key))
    isBCCheck.grid(row=2, column=1, sticky="ew")
    saveNodePropertyBtn = Button(editNodePage, text="Save Changes",
                                 command=lambda: saveNodeEdit(key, node_head_entry.get()))
    saveNodePropertyBtn.grid(row=6, column=0, columnspan=2)


def saveNodeEdit(key, node_head_entry):
    """end function to save the new edited node"""
    NodesDictionary[key]["head"] = node_head_entry
    if isBCVar.get():
        NodesDictionary[key].update({"isBC": True, "BCType": BCVar.get()})
    else:
        NodesDictionary[key].update({"isBC": False, "BCType": None})
    editNodePage.destroy()


def remove_node():
    """front end gui to delete node"""
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
    """back end function to remove node from canvas, and clear dictionary corresponding to the node"""
    drawArea.delete(NodeName)
    drawArea.delete(NodeName + "text")
    edit_node_menu.delete(NodeName)
    del NodesDictionary[NodeName]
    delete_node.destroy()


def link_MOC_property_field(displayWindow, key="New", HasPert=False):
    """modular function to display front end gui of MOC property"""
    global MOCAdditionFrame
    MOCAdditionFrame = LabelFrame(displayWindow, text="Additional Boundary Condition", highlightthickness=2)
    Label(MOCAdditionFrame, text="Has perturbation").grid(row=0, column=0, sticky="w")
    Label(MOCAdditionFrame, text="NOTE: This only applies to boundary nodes").grid(row=1, column=0, columnspan=2,
                                                                                   sticky="ew")
    global MOCPerturbationTypeVar
    if not HasPert:
        MOCPerturbationTypeVar = StringVar()
        MOCPerturbationTypeVar.set("None")
        PertTypeList = ["None", "Impulse", "Full Closure", "Sinusoidal Head", "Sinusoidal Flow", "Controlled Flow"]
        PerturbationSelection = OptionMenu(MOCAdditionFrame, MOCPerturbationTypeVar, *PertTypeList,
                                           command=lambda v: MOC_perturbation_field((displayWindow, v, "New")))
        PerturbationSelection.grid(row=0, column=1, sticky="ew")
    else:
        MOCPerturbationTypeVar = StringVar()
        MOCPerturbationTypeVar.set(LinksDictionary[key]["PertType"])
        PertTypeList = ["None", "Impulse", "Full Closure", "Sinusoidal Head", "Sinusoidal Flow", "Controlled Flow"]
        passed = (displayWindow, LinksDictionary[key]["PertType"], key)
        MOC_perturbation_field(passed)
        PerturbationSelection = OptionMenu(MOCAdditionFrame, MOCPerturbationTypeVar, *PertTypeList,
                                           command=lambda v: MOC_perturbation_field((displayWindow, v, "New")))
        PerturbationSelection.grid(row=0, column=1, sticky="ew")
    MOCAdditionFrame.grid(row=2, sticky="nsew")
    MOCAdditionFrame.grid_columnconfigure(0, weight=1)


def MOC_perturbation_field(passed):
    """modular function to display front end gui to introduce MOC perturbation"""
    displayWindow, PerturbationType, key = passed
    global PertFrame, Location, FlowRate, Frequency, Amplitude, StartTime, impulseSize
    if key == "New":
        """When a new edge is created, the following code will be used"""
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
            Label(ImpulseFrame, text="Impulse Size").grid(row = 1, column=0, sticky="w")
            StartTime = Entry(ImpulseFrame)
            StartTime.grid(row=0, column=1, sticky="w")  # impulse time entry
            impulseSize = Entry(ImpulseFrame)
            impulseSize.grid(row=1, column=1, sticky="w")
            ImpulseFrame.grid(row=1, column=0, columnspan=2, sticky="nsew")
            ImpulseFrame.grid_columnconfigure(0, weight=1)
            PertFrame.grid(row=10, column=0, sticky="ew")
            PertFrame.grid_columnconfigure(0, weight=1)
        elif PerturbationType == "Full Closure":
            FCFrame = LabelFrame(PertFrame, highlightthickness=0, borderwidth=0)
            Label(FCFrame, text="Time").grid(row=0, column=0, sticky="w")  # impulse time
            StartTime = Entry(FCFrame)
            StartTime.grid(row=0, column=1, sticky="w")  # impulse time entry
            FCFrame.grid(row=1, column=0, columnspan=2, sticky="ew")
            FCFrame.grid_columnconfigure(0, weight=1)
            PertFrame.grid(row=10, column=0, sticky="ew")
            PertFrame.grid_columnconfigure(0, weight=1)
        elif PerturbationType == "Controlled Flow":
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
    else:
        """When editing old links, use the following code, this is because for edited pages, entry needed to be filled with existing info"""
        PertFrame = LabelFrame(displayWindow, text=PerturbationType, highlightthickness=0)
        Label(PertFrame, text="Location").grid(row=0, column=0, sticky="w")  # perturbation location label
        Location = Entry(PertFrame)
        Location.insert(0, LinksDictionary[key]["Location"])
        Location.grid(row=0, column=1, sticky="e")
        if PerturbationType == "Impulse":
            ImpulseFrame = LabelFrame(PertFrame, highlightthickness=0, borderwidth=0)
            Label(ImpulseFrame, text="Time").grid(row=0, column=0, sticky="w")  # impulse time
            Label(ImpulseFrame, text="Impulse Size").grid(row = 1, column=0, sticky="w")
            StartTime = Entry(ImpulseFrame)
            StartTime.insert(0, LinksDictionary[key]["Time"])
            StartTime.grid(row=0, column=1, sticky="w")  # impulse time entry
            impulseSize = Entry(ImpulseFrame)
            impulseSize.grid(row=1, column=1, sticky="w")
            impulseSize.insert(0, LinksDictionary[key]["impulseSize"])
            ImpulseFrame.grid(row=1, column=0, columnspan=2, sticky="nsew")
            ImpulseFrame.grid_columnconfigure(0, weight=1)
            PertFrame.grid(row=10, column=0, sticky="ew")
            PertFrame.grid_columnconfigure(0, weight=1)
        elif PerturbationType == "Full Closure":
            FCFrame = LabelFrame(PertFrame, highlightthickness=0, borderwidth=0)
            Label(FCFrame, text="Time").grid(row=0, column=0, sticky="w")  # impulse time
            StartTime = Entry(FCFrame)
            StartTime.insert(0, LinksDictionary[key]["Time"])
            StartTime.grid(row=0, column=1, sticky="w")  # impulse time entry
            FCFrame.grid(row=1, column=0, columnspan=2, sticky="ew")
            FCFrame.grid_columnconfigure(0, weight=1)
            PertFrame.grid(row=10, column=0, sticky="ew")
            PertFrame.grid_columnconfigure(0, weight=1)
        elif PerturbationType == "Controlled Flow":
            CFFrame = LabelFrame(PertFrame, highlightthickness=0, borderwidth=0)
            Label(CFFrame, text="Time").grid(row=0, column=0, sticky="w")  # impulse time
            StartTime = Entry(CFFrame)
            StartTime.insert(0, LinksDictionary[key]["Time"])
            StartTime.grid(row=0, column=1, sticky="e")  # impulse time entry
            Label(CFFrame, text="Flow rate").grid(row=1, column=0, sticky="w")
            FlowRate = Entry(CFFrame)
            FlowRate.insert(0, LinksDictionary[key]["FlowRate"])
            FlowRate.grid(row=1, column=1, sticky="e")
            CFFrame.grid(row=1, column=0, columnspan=2, sticky="ew")
            CFFrame.grid_columnconfigure(0, weight=1)
            PertFrame.grid(row=10, column=0, sticky="ew")
            PertFrame.grid_columnconfigure(0, weight=1)
        elif PerturbationType == "None":
            try:
                PertFrame.destroy()
            except NameError:
                pass
        else:
            SineFrame = LabelFrame(PertFrame, highlightthickness=2, borderwidth=0)
            Label(SineFrame, text="Frequency").grid(row=0, column=0, sticky="w")
            Label(SineFrame, text="Amplitude").grid(row=1, column=0, sticky="w")
            Frequency = Entry(SineFrame)
            Frequency.insert(0, LinksDictionary[key]["Freq"])
            Amplitude = Entry(SineFrame)
            Amplitude.insert(0, LinksDictionary[key]["Amp"])
            Frequency.grid(row=0, column=1, sticky="e")
            Amplitude.grid(row=1, column=1, sticky="e")
            SineFrame.grid(row=2, column=0, columnspan=2, sticky="ew")
            SineFrame.grid_columnconfigure(0, weight=1)
            PertFrame.grid(row=10, column=0, sticky="ew")
            PertFrame.grid_columnconfigure(0, weight=1)


def link_TM_property_field(displayWindow, key="New", HasPert=False, HasSensor=False):
    """modular function to display front end gui of TM property"""
    global TMAdditionFrame
    TMAdditionFrame = LabelFrame(displayWindow, text="Additional Boundary Condition", highlightthickness=2)
    Label(TMAdditionFrame, text="Has perturbation").grid(row=0, column=0, sticky="w")
    Label(TMAdditionFrame, text="NOTE: This only applies to boundary nodes").grid(row=1, column=0, columnspan=2,
                                                                                  sticky="ew")
    Label(TMAdditionFrame, text="Has Sensor").grid(row=2, column=0, sticky="w")
    # AdditionFrame.pack(fill=BOTH)
    global TMPerturbationTypeVar, SensorVar
    if not HasPert:
        TMPerturbationTypeVar = StringVar()
        TMPerturbationTypeVar.set("None")
        PerturbationSelection = OptionMenu(TMAdditionFrame, TMPerturbationTypeVar, "None", "Flow", "Head",
                                           command=lambda v: TM_perturbation_field((displayWindow, v, "New")))
        PerturbationSelection.grid(row=0, column=1, sticky="ew")
    else:
        TMPerturbationTypeVar = StringVar()
        TMPerturbationTypeVar.set(LinksDictionary[key]["PertType"])
        passed = (displayWindow, LinksDictionary[key]["PertType"], key)
        TM_perturbation_field(passed)
        PerturbationSelection = OptionMenu(TMAdditionFrame, TMPerturbationTypeVar, "None", "Flow", "Head",
                                           command=lambda v: TM_perturbation_field((displayWindow, v, "New")))
        PerturbationSelection.grid(row=0, column=1, sticky="ew")
    if not HasSensor:
        SensorVar = BooleanVar()
        SensorVar.set(False)
        SensorCheck = Checkbutton(TMAdditionFrame, variable=SensorVar,
                                  command=lambda: TM_sensor_field(displayWindow))
        SensorCheck.grid(row=2, column=1, sticky="ew")
    else:
        SensorVar = BooleanVar()
        SensorVar.set(HasSensor)
        TM_sensor_field(displayWindow, key)
        SensorCheck = Checkbutton(TMAdditionFrame, variable=SensorVar,
                                  command=lambda: TM_sensor_field(displayWindow))
        SensorCheck.grid(row=2, column=1, sticky="ew")
    TMAdditionFrame.grid(row=2, sticky="nsew")
    TMAdditionFrame.grid_columnconfigure(0, weight=1)


def TM_perturbation_field(passed):
    """modular function to display front end gui to introduce TM perturbation"""
    displayWindow, PerturbationType, key = passed
    global PertFrame, Location
    if key == "New":
        """When a new link is created, the following code is used"""
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
    else:
        """When a old link is edited, the following code is used"""
        PertFrame = LabelFrame(displayWindow, text=PerturbationType, highlightthickness=0)
        Label(PertFrame, text="Location").grid(row=0, column=0, sticky="w")  # perturbation location label
        Location = Entry(PertFrame)
        Location.insert(0, LinksDictionary[key]["Location"])
        Location.grid(row=0, column=1, sticky="e")
        PertFrame.grid(row=10, column=0, sticky="ew")
        PertFrame.grid_columnconfigure(0, weight=1)


def TM_sensor_field(displayWindow, key="New"):
    """modular function to display front end gui to introduce TM sensor"""
    global SensorFrame, SensorLocation, SensorDist
    if key == "New":
        """When a new link is created, use the following code"""
        try:
            SensorFrame.destroy()
        except NameError:
            pass
        SensorFrame = LabelFrame(displayWindow, highlightthickness=0)
        Label(SensorFrame, text="Sensor Location").grid(row=0, column=0, sticky="w")  # Sensor location label
        SensorLocation = Entry(SensorFrame)
        SensorLocation.grid(row=0, column=1, sticky="e")
        Label(SensorFrame, text="Distance to boundary").grid(row=1, column=0, sticky="w")
        SensorDist = Entry(SensorFrame)
        SensorDist.grid(row=1, column=1, sticky="e")
        if SensorVar.get():
            SensorFrame.grid(row=15, sticky="ew")
            SensorFrame.grid_columnconfigure(0, weight=1)
        else:
            pass
    else:
        """When a old link is edited, use the following code"""
        SensorFrame = LabelFrame(displayWindow, highlightthickness=0)
        Label(SensorFrame, text="Sensor Location").grid(row=0, column=0, sticky="w")  # Sensor location label
        Label(SensorFrame, text="Distance to boundary").grid(row=1, column=0, sticky="w")
        SensorLocation = Entry(SensorFrame)
        SensorLocation.insert(0, LinksDictionary[key]["SensorLocation"])
        SensorDist = Entry(SensorFrame)
        SensorDist.insert(0, LinksDictionary[key]["SensorDist"])
        SensorLocation.grid(row=0, column=1, sticky="e")
        SensorDist.grid(row=1, column=1, sticky="e")
        SensorFrame.grid(row=15, sticky="ew")
        SensorFrame.grid_columnconfigure(0, weight=1)


# noinspection PyGlobalUndefined
def link_Basic_property_field(displayWindow, key="New"):
    """modular function to mount the display for link's basic property"""
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
    l = Entry(BasicFrame)
    d = Entry(BasicFrame)
    a = Entry(BasicFrame)
    f = Entry(BasicFrame)
    u = Entry(BasicFrame)
    if key == "New":
        pass
    else:
        """Fill in info from backend dictionary to entry field when editing links"""
        l.insert(0, LinksDictionary[key]["Length"])
        d.insert(0, LinksDictionary[key]["Diameter"])
        a.insert(0, LinksDictionary[key]["Wave speed"])
        f.insert(0, LinksDictionary[key]["Friction factor"])
        u.insert(0, LinksDictionary[key]["Flow Velocity"])
    l.grid(row=0, column=1, sticky="ew")
    d.grid(row=1, column=1, sticky="e")
    a.grid(row=2, column=1, sticky="e")
    f.grid(row=3, column=1, sticky="e")
    u.grid(row=4, column=1, sticky="e")


def new_link():
    """call all the modular functions to create the front end gui for creating new links"""
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
    """back end function to create link on canvas and save entry field to backend dictionary"""
    start = startNode.get()
    end = endNode.get()
    x0, y0 = NodesDictionary[start]["Coord"]
    x1, y1 = NodesDictionary[end]["Coord"]
    drawArea.create_line(x0, y0, x1, y1, tag="{},{}".format(start, end))
    LinksDictionary["{},{}".format(start, end)] = {"Length": l.get(), "Diameter": d.get(), "Wave speed": a.get(),
                                                   "Friction factor": f.get(), "Flow Velocity": u.get()}
    if not AnalysisMode.get() and MOCPerturbationTypeVar.get() != "None":
        """Save information for MOC"""
        LinksDictionary["{},{}".format(start, end)].update(
            {"PertType": MOCPerturbationTypeVar.get(), "Location": Location.get()})
        if MOCPerturbationTypeVar.get() == "Impulse":
            LinksDictionary["{},{}".format(start, end)].update({"Time": StartTime.get(), "impulseSize":impulseSize.get()})
        elif MOCPerturbationTypeVar.get() == "Full Closure":
            LinksDictionary["{},{}".format(start, end)].update({"Time": StartTime.get()})
        elif MOCPerturbationTypeVar.get() == "Controlled Flow":
            LinksDictionary["{},{}".format(start, end)].update({"Time": StartTime.get(), "FlowRate": FlowRate.get()})
        else:
            LinksDictionary["{},{}".format(start, end)].update({"Freq": Frequency.get(), "Amp": Amplitude.get()})
    elif not AnalysisMode.get() and MOCPerturbationTypeVar.get() == "None":
        """Save information for MOC"""
        LinksDictionary["{},{}".format(start, end)].update({"PertType": ""})
    else:
        """Save information for TM"""
        if TMPerturbationTypeVar.get() != "None":
            LinksDictionary["{},{}".format(start, end)].update(
                {"PertType": TMPerturbationTypeVar.get(), "Location": Location.get()})
        else:
            LinksDictionary["{},{}".format(start, end)].update(
                {"PertType": ""})
        if SensorVar.get():
            LinksDictionary["{},{}".format(start, end)].update(
                {"HasSensor": SensorVar.get(), "SensorLocation": SensorLocation.get(), "SensorDist": SensorDist.get()})
        else:
            LinksDictionary["{},{}".format(start, end)].update({"HasSensor": ""})
    x = (x0 + x1) / 2
    y = (y0 + y1) / 2
    drawArea.create_text(x, y + 30, text=l.get() + "m", tag="{},{}".format(start, end) + "text")
    edit_link_menu.add_command(label=("{},{}".format(start, end)), command=lambda: editLink("{},{}".format(start, end)))
    add_link.destroy()


def editLink(key):
    """call all the modular function to create the front end gui used to edit saved links"""
    global editLinkPage
    editLinkPage = Toplevel()
    editLinkPage.title(key)
    editLinkPage.resizable(0, 0)
    editLinkPage.geometry(f"+{root.winfo_x() + 200}+{root.winfo_y() + 200}")
    start, end = key.split(",")
    Identificationlabel = Label(editLinkPage, text="Source Node = {}, Target Node = {}".format(start, end))
    Identificationlabel.grid(row=0, column=0, columnspan=2)
    HasPert = LinksDictionary[key]["PertType"] != ""
    link_Basic_property_field(editLinkPage, key)
    if AnalysisMode.get():
        link_TM_property_field(editLinkPage, key=key, HasPert=HasPert, HasSensor=LinksDictionary[key]["HasSensor"])
    else:
        link_MOC_property_field(editLinkPage, key=key, HasPert=HasPert)

    save_btn = Button(editLinkPage, text="Save Changes",
                      command=lambda: saveLinkEdit(key))
    save_btn.grid(row=30, column=0, columnspan=2)


def saveLinkEdit(key):
    """Back end function to refresh back end links'dictionary based on entry filled in the entry field in edit link page"""
    start, end = key.split(",")
    LinksDictionary[key] = {}
    LinksDictionary[key]["Length"] = l.get()
    LinksDictionary[key]["Diameter"] = d.get()
    LinksDictionary[key]["Wave speed"] = a.get()
    LinksDictionary[key]["Friction factor"] = f.get()
    LinksDictionary[key]["Flow Velocity"] = u.get()
    if not AnalysisMode.get() and MOCPerturbationTypeVar.get() != "None":
        LinksDictionary["{},{}".format(start, end)].update(
            {"PertType": MOCPerturbationTypeVar.get(), "Location": Location.get()})
        if MOCPerturbationTypeVar.get() == "Impulse":
            LinksDictionary["{},{}".format(start, end)].update({"Time": StartTime.get(), "impulseSize":impulseSize.get()})
        elif MOCPerturbationTypeVar.get() == "Full Closure":
            LinksDictionary["{},{}".format(start, end)].update({"Time": StartTime.get()})
        elif MOCPerturbationTypeVar.get() == "Controlled Flow":
            LinksDictionary["{},{}".format(start, end)].update({"Time": StartTime.get(), "FlowRate": FlowRate.get()})
        else:
            LinksDictionary["{},{}".format(start, end)].update({"Freq": Frequency.get(), "Amp": Amplitude.get()})
    elif not AnalysisMode.get() and MOCPerturbationTypeVar.get() == "None":
        LinksDictionary["{},{}".format(start, end)].update({"PertType": ""})
    else:
        if TMPerturbationTypeVar.get() != "None":
            LinksDictionary["{},{}".format(start, end)].update(
                {"PertType": TMPerturbationTypeVar.get(), "Location": Location.get()})
        else:
            LinksDictionary["{},{}".format(start, end)].update(
                {"PertType": ""})
        if SensorVar.get():
            LinksDictionary["{},{}".format(start, end)].update(
                {"HasSensor": SensorVar.get(), "SensorLocation": SensorLocation.get(), "SensorDist": SensorDist.get()})
        else:
            LinksDictionary["{},{}".format(start, end)].update({"HasSensor": ""})
    drawArea.delete(key)
    drawArea.delete(key + "text")
    start, end = key.split(",")
    x0, y0 = NodesDictionary[start]["Coord"]
    x1, y1 = NodesDictionary[end]["Coord"]
    drawArea.create_line(x0, y0, x1, y1, tag=key)
    x = (x0 + x1) / 2
    y = (y0 + y1) / 2
    drawArea.create_text(x, y + 30, text=l.get() + "m", tag=key + "text")
    editLinkPage.destroy()


def remove_link():
    """front end gui to delete link"""
    global delete_link
    delete_link = Toplevel()
    delete_link.title("Delete specified link")
    Label(delete_link, text="Specify start and end of the link to delete ").grid(row=1, column=0)
    Label(delete_link, text="Start").grid(row=0, column=1)
    Label(delete_link, text="End").grid(row=0, column=2)
    startEntry = Entry(delete_link)
    endEntry = Entry(delete_link)
    startEntry.grid(row=1, column=1)
    endEntry.grid(row=1, column=2)
    delete_btn = Button(delete_link, text="Delete", command=lambda: destroy_link(startEntry.get(), endEntry.get()))
    delete_btn.grid(row=2, column=0, columnspan=3)


def destroy_link(start, end):
    """back end function to remove link from back end dictionary"""
    drawArea.delete("{},{}".format(start, end))
    drawArea.delete("{},{}".format(start, end) + "text")
    edit_link_menu.delete("{},{}".format(start, end))
    del LinksDictionary["{},{}".format(start, end)]
    delete_link.destroy()


def ViewLink():
    """Front end gui to display back end link dictionary"""
    LinkSummary = Toplevel()
    LinkSummary.title("Link Summaries")
    # editLinkPage.resizable(0, 0)
    LinkSummary.geometry(f"+{root.winfo_x() + 200}+{root.winfo_y() + 200}")
    rowindex = 1
    columnindex = 1
    PropertiesMatchingDict = {}
    LinksMatchingDict = {}
    for Link, properties in LinksDictionary.items():
        for property, value in properties.items():
            if not property in PropertiesMatchingDict.keys():
                PropertiesMatchingDict[property] = columnindex
                tempFrame = LabelFrame(LinkSummary, borderwidth=1, highlightthickness=1)
                tempFrame.grid(row=0, column=PropertiesMatchingDict[property], sticky="nsew")
                Label(tempFrame, text=property).pack()
                columnindex += 1
        if not Link in LinksMatchingDict.keys():
            LinksMatchingDict[Link] = rowindex
            nodeFrame = LabelFrame(LinkSummary, borderwidth=1, highlightthickness=1)
            nodeFrame.grid(row=LinksMatchingDict[Link], column=0, sticky="nsew")
            Label(nodeFrame, text=Link).pack()
            rowindex += 1
    for Link, properties in LinksDictionary.items():
        for property, value in properties.items():
            row = LinksMatchingDict[Link]
            column = PropertiesMatchingDict[property]
            valueFrame = LabelFrame(LinkSummary, borderwidth=1, highlightthickness=1)
            valueFrame.grid(row=row, column=column, sticky="nsew")
            Label(valueFrame, text=value).pack()
    LinkSummary.grid_columnconfigure(0, weight=1)


def ViewNode():
    """Front end gui to display back end node dictionary"""
    NodeSummary = Toplevel()
    NodeSummary.title("Node Summaries")
    # editLinkPage.resizable(0, 0)
    NodeSummary.geometry(f"+{root.winfo_x() + 200}+{root.winfo_y() + 200}")
    rowindex = 1
    columnindex = 1
    PropertiesMatchingDict = {}
    nodesMatchingDict = {}
    for node, properties in NodesDictionary.items():
        for property, value in properties.items():
            if not property in PropertiesMatchingDict.keys():
                PropertiesMatchingDict[property] = columnindex
                tempFrame = LabelFrame(NodeSummary, borderwidth=1, highlightthickness=1)
                tempFrame.grid(row=0, column=PropertiesMatchingDict[property], sticky="nsew")
                Label(tempFrame, text=property).pack()
                columnindex += 1
        if not node in nodesMatchingDict.keys():
            nodesMatchingDict[node] = rowindex
            nodeFrame = LabelFrame(NodeSummary, borderwidth=1, highlightthickness=1)
            nodeFrame.grid(row=nodesMatchingDict[node], column=0, sticky="nsew")
            Label(nodeFrame, text=node).pack()
            rowindex += 1
    for node, properties in NodesDictionary.items():
        for property, value in properties.items():
            row = nodesMatchingDict[node]
            column = PropertiesMatchingDict[property]
            valueFrame = LabelFrame(NodeSummary, borderwidth=1, highlightthickness=1)
            valueFrame.grid(row=row, column=column, sticky="nsew")
            Label(valueFrame, text=value).pack()
    NodeSummary.grid_columnconfigure(0, weight=1)


def clear(Warning=True, flip=False):
    """Delete everything on canvas, and clear all back end dictionary"""
    global NodesDictionary, LinksDictionary, EnvirDictionary
    if Warning:
        if messagebox.askyesno("Warning", message="Contiune will delete everything !"):
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
            return True
        else:
            if flip:
                flipedAnalysisMode = not AnalysisMode.get()
                AnalysisMode.set(flipedAnalysisMode)
            else:
                pass
    else:
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
        return True



""" Main Window """
windll.shcore.SetProcessDpiAwareness(1)
root = Tk()
root.title("Hydraulic Transient Analysis")
# Base size
normal_width = 1080
normal_height = 768
# Get screen size
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
# Get percentage of screen size from Base size
percentage_width = screen_width / (normal_width / 100)
percentage_height = screen_height / (normal_height / 100)
# Make a scaling factor, this is bases on average percentage from
# width and height.
scale_factor = 1.1 * (((percentage_width + percentage_height) / 4) / 100)
displayRes = "{}x{}".format(int(normal_width * scale_factor), int(normal_height * scale_factor))
root.geometry(displayRes)  # Define initial window size
# # Set the fontsize based on scale_factor,
# # if the fontsize is less than minimum_size
# # it is set to the minimum size
fontsize = int(10 * scale_factor)
minimum_size = 8
if fontsize < minimum_size:
    fontsize = minimum_size
# Create a style and configure for ttk.Button widget
default_style = ttk.Style()
default_style.configure('.', font=("Helvetica", fontsize))

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
                            command=ChangeMode)
method_menu.add_checkbutton(label="Method of Characteristic", onvalue=0, offvalue=1, variable=AnalysisMode,
                            command=ChangeMode)
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
pw.add(RightFrame, width=400 * scale_factor, stretch="never")
pw.pack(fill=BOTH, expand=True)
global MOCFrame, TMFrame  # declared global because it will be used in Function: Refresh
MOCFrame = LabelFrame(RightFrame, borderwidth=0, highlightthickness=0)
TMFrame = LabelFrame(RightFrame, borderwidth=0, highlightthickness=0)

"""Define MOC Frame layout"""
# Define entry boxes
Label(MOCFrame, text="MOC Run Configuration", font=("Helvetica", "24")).grid(row=0, column=0, columnspan=4, pady=10,
                                                                             sticky="ew")  # Title
Label(MOCFrame, text="dt = ", anchor="e").grid(row=1, column=0, pady=10)  # dt Label
Label(MOCFrame, text="second").grid(row=1, column=3, pady=10)  # unit Label
Label(MOCFrame, text="Total Run Time").grid(row=2, column=0, pady=10)  # total time Label
Label(MOCFrame, text="second").grid(row=2, column=3, pady=10)  # unit Label
Label(MOCFrame, text="Record data from").grid(row=3, column=0, pady=10)  # record range Label
Label(MOCFrame, text="to").grid(row=3, column=2, pady=10)  # "to" Label
Label(MOCFrame, text="Start/second").grid(row=4, column=1)  # unit Label
Label(MOCFrame, text="End/second").grid(row=4, column=3)  # unit Label
timestepEntry = Entry(MOCFrame)
timestepEntry.grid(row=1, column=1, columnspan=2, pady=10, padx=0, sticky="WE")
runtimeEntry = Entry(MOCFrame)
runtimeEntry.grid(row=2, column=1, columnspan=2, pady=10, ipadx=0, sticky="WE")
record_startEntry = Entry(MOCFrame)
record_endEntry = Entry(MOCFrame)
record_startEntry.grid(row=3, column=1, pady=(10, 0))  # pady can take a tuple (top,bottom)
record_endEntry.grid(row=3, column=3, pady=(10, 0))

MOCFrame.grid_columnconfigure(0, weight=1)
# MOC save and start analysis
saveconfigBtn = Button(MOCFrame, text="Save Configuration", command=saveConfig)
saveconfigBtn.grid(row=20, column=0, columnspan=4, pady=(10, 0), sticky="ew")
StartBtn = Button(MOCFrame, text="Start Analysis", command=startAnalysis)
StartBtn.grid(row=21, column=0, columnspan=4, pady=10, sticky="ew")

"""Define TransferMatrx Frame layout"""
Label(TMFrame, text="TM Run Configuration", font=("Helvetica", "24")).grid(row=0, column=0, columnspan=4, pady=10,
                                                                           sticky="ew")  # Title

SelectionFrame = LabelFrame(TMFrame, borderwidth=0, highlightthickness=0)
Label(SelectionFrame, text="FrequencyMode").grid(row=0, column=0, sticky="ew")
frequencyModeVar = StringVar()
frequencyModeVar.set("MultiFrequency")
FrequencyModeSelection = OptionMenu(SelectionFrame, frequencyModeVar, "MultiFrequency", "SingleFrequency", "Randomized Noise")
FrequencyModeSelection.grid(row=0, column=1, pady=10, padx=10, sticky="ew")
SelectionFrame.grid(row=1, column=0, columnspan=4, padx=20, sticky="ew")
SelectionFrame.grid_columnconfigure(0, weight=1)
# Create Frame for multi-frequency analysis
MultiFrame = LabelFrame(TMFrame, borderwidth=0, highlightthickness=0)
Label(MultiFrame, text="dfreq = ").grid(row=0, column=0, pady=10, sticky="ew")
Label(MultiFrame, text="Max Freq").grid(row=1, column=0, pady=10, sticky="ew")
MultifreqStepEntry = Entry(MultiFrame)
MultifreqStepEntry.grid(row=0, column=1, pady=10, padx=10, columnspan=2, sticky="WE")
MultiFreqRangeEntry = Entry(MultiFrame)
MultiFreqRangeEntry.grid(row=1, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
MultiFrame.grid_columnconfigure(0, weight=1)

# Create Frame for single frequency analysis
SingleFrame = LabelFrame(TMFrame, borderwidth=0, highlightthickness=0)
Label(SingleFrame, text="Excitation Frequency").grid(row=0, column=0, pady=10)
FreqEntry = Entry(SingleFrame)
FreqEntry.grid(row=0, column=1, pady=10, columnspan=2, padx=10, sticky="WE")
Label(SingleFrame, text="dfreq = ").grid(row=1, column=0, pady=10)
Label(SingleFrame, text="Max Freq").grid(row=2, column=0, pady=10)
SinglefreqStepEntry = Entry(SingleFrame)
SinglefreqStepEntry.grid(row=1, column=1, pady=10, padx=10, columnspan=2, sticky="WE")
SingleFreqRangeEntry = Entry(SingleFrame)
SingleFreqRangeEntry.grid(row=2, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
frequencyModeVar.trace("w", Refresh_Panel)
SingleFrame.grid_columnconfigure(0, weight=1)


# Create Frame for Randomized Noise analysis
NoiseFrame = LabelFrame(TMFrame, borderwidth=0, highlightthickness=0)
Label(NoiseFrame, text="dfreq = ").grid(row=0, column=0, pady=10, sticky="ew")
Label(NoiseFrame, text="Max Freq").grid(row=1, column=0, pady=10, sticky="ew")
Label(NoiseFrame, text="Simulation Size").grid(row=2, column=0, pady=10, sticky="ew")
Label(NoiseFrame, text="Mean Noise Level").grid(row=3, column=0, pady=10, sticky="ew")
Label(NoiseFrame, text="Noise Standard Deviation").grid(row=4, column=0, pady=10, sticky="ew")
Label(NoiseFrame, text="Sensor 1").grid(row=5, column=0, pady=10, sticky="ew")
Label(NoiseFrame, text="Sensor 2").grid(row=6, column=0, pady=10, sticky="ew")
NoiseFreqStepEntry = Entry(NoiseFrame)
NoiseFreqStepEntry.grid(row=0, column=1, pady=10, padx=10, columnspan=2, sticky="WE")
NoiseFreqRangeEntry = Entry(NoiseFrame)
NoiseFreqRangeEntry.grid(row=1, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
SimSizeEntry = Entry(NoiseFrame)
SimSizeEntry.grid(row=2, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
MeanNoiseEntry = Entry(NoiseFrame)
MeanNoiseEntry.grid(row=3, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
NoiseSTDEntry = Entry(NoiseFrame)
NoiseSTDEntry.grid(row=4, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
Sensor1Entry = Entry(NoiseFrame)
Sensor1Entry.grid(row=5, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
Sensor2Entry = Entry(NoiseFrame)
Sensor2Entry.grid(row=6, column=1, pady=10, padx=10, columnspan=2, ipadx=0, sticky="WE")
frequencyModeVar.trace("w", Refresh_Panel)
NoiseFrame.grid_columnconfigure(0, weight=1)

# TM save and start analysis
saveconfigBtn = Button(TMFrame, text="Save Configuration", command=saveConfig)
saveconfigBtn.grid(row=10, column=0, columnspan=4, pady=(10, 0), sticky="ew")
StartBtn = Button(TMFrame, text="Start Analysis", command=startAnalysis)
StartBtn.grid(row=11, column=0, columnspan=4, pady=10, sticky="ew")
TMFrame.grid_columnconfigure(0, weight=1)

""" Display selected analysis mode entry boxes on screen """
AnalysisMode.set(True)
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
