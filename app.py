import HydraulicTransientAnalysis_GUI
import tkinter
from tkinter import *
from tkinter import filedialog, ttk, messagebox
import json
import os
from HydraulicTransientAnalysis_GUI import *


if __name__ == '__main__':
    """ Main Window """
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
    FrequencyModeSelection = OptionMenu(SelectionFrame, frequencyModeVar, "MultiFrequency", "SingleFrequency",
                                        "Randomized Noise")
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