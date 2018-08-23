﻿#Fixing incompatability issue
import matplotlib
matplotlib.use('TkAgg')
#TKinter Modules
import tkinter as TK 
from tkinter import ttk
from tkinter.filedialog import askopenfilename
#KinMS Modules
from KinMS import KinMS
from makebeam import makebeam
#Importing Numpy, Pandas and Matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#Scipy minimize for fitting 
from scipy.optimize import minimize
#Importing .fits reader
from astropy.io import fits
#Importing threading, webbrowser and time
import threading
import webbrowser 
import time
#Functions required for GUI operation 
import GUIFunctions as GUI


#########################################################################################################################################################
"""
Standard Functions
"""
#########################################################################################################################################################

def gaussian(x,x0,sigma):
    """
    Returns a guasian curve for given inputs
    """
    return np.exp(-((x - x0)/sigma)**2/2.)

def AOI(cube):
    '''
    Function estimates the start and end point of the signal in the fits data.
    A threshold value is determined from the spec mean then the first and last 
    points to cross this threshold are taken as the staring and end points.
    '''
    spec = cube.sum(axis=0).sum(axis=0)
    Threshold = np.mean(spec) * 0.8
    Count = 0
    for i in spec:
        if i > Threshold:
            Start = Count - 5
            if Start < 0:
                Start = 1
            break
        Count += 1
    Count = len(spec)
    for i in spec[::-1]:
        if i > Threshold:
            Stop = Count + 5
            if Stop > len(spec):
                Stop = len(spec) - 1
            break    
        Count -= 1
    Start, Stop = int(Start), int(Stop)
    return Start, Stop  
 
#Extract Values
def ExtractVarFromStr(String):
    '''
    Takes the user input equation and Extracts variables and limits
    '''
    Index = []
    i = String.find('|')
    while i >= 0 :
        Index.append(i)
        i = String.find('|', i+1)
    Values = []
    i = 0
    while i <= len(Index)-1:
        Values.append(String[int(Index[i]+1): int(Index[i+1])])
        i += 2
    Inputs, Bounds = [], []
    for i in Values:
        if ',' in i:
            i = i.replace(' ','')
            i = i.replace('(','')
            i = i.replace(')','')
            segment = i.split(',')
            Inputs.append(eval(segment[0]))
            Bounds.append((eval(segment[1]), eval(segment[2])))
        else:
            Inputs.append(eval(i))
            Bounds.append((0,1e50))  
    return Inputs, Bounds

#Edit String
def EditUserString(String, ValueArray):
    '''
    Takes the user input equation and inputs new variables from ValueArray
    '''
    Index = []
    i = String.find('|')
    while i >= 0 :
        Index.append(i)
        i = String.find('|', i+1)
    NewString = String[:Index[0]+1]
    ValueCount = 0
    IndexCount = 1
    while IndexCount < len(Index) -1 :
        NewString += '%f' + String[Index[int(IndexCount)]:Index[int(IndexCount+1)]]
        NewString = NewString %ValueArray[ValueCount]
        ValueCount += 1
        IndexCount += 2
    NewString += '%f' + String[Index[int(IndexCount)]:]
    NewString = NewString.replace('|', '').replace('\n', '')  %ValueArray[-1]
    return NewString

#########################################################################################################################################################
"""
self for classes Functions
"""
#########################################################################################################################################################

def LaunchWebpage():
    url = 'https://github.com/TristanBurman/KinMSGUI'
    return webbrowser.open(url, new=0, autoraise=True)

#Self functions
def SwitchMode(self, Win):
    if Win == 'Read Mode':
        self.master.withdraw()
        self.newWindow = TK.Toplevel(self.master)
        self.app = ReadWin(self.newWindow)
    if Win == 'Custom Model':
        self.master.withdraw()
        self.newWindow = TK.Toplevel(self.master)
        self.app = CustomWin(self.newWindow)

def KeyWord(self):
    """
    Loads Keyword DataBase and opens new window
    """
    def LookUp(Keywords, Help, option, TextBox):
        #Read selected keyword
        Selected = var.get()
        #Find Index
        Index = Keywords.index(Selected)
        #Update Text Box
        TextBox.delete('1.0', TK.END)
        TextBox.insert(TK.INSERT, Help[Index])
    def Close():
        #Close Window
        KeyWordWin.destroy()
    #Create Window
    KeyWordWin = TK.Tk()
    KeyWordWin.configure(background='black')
    KeyWordWin.title('KinMS GUI - Keyword Lookup')
    KeyWordWin.configure(background='black')
    KeyWordWin.geometry("320x300")
    TK.Label(KeyWordWin, text='Keyword Lookup', bg='black', fg='white', font='none 12 bold').grid(row=0, column=0, columnspan=2, sticky="nesw")
    #load Database
    KeyWordDatabase = pd.read_csv('KeyWordLookUp.csv', delimiter=',')
    Keywords, Help = KeyWordDatabase['Keyword'].tolist(), KeyWordDatabase['Help'].tolist()
    #Create dropdown box
    var = TK.StringVar(KeyWordWin)
    option = TK.OptionMenu(KeyWordWin, var, *Keywords)
    option.grid(row=1, columnspan=2, sticky='nesw')
    #Create Text Box
    TextBox = TK.Text(KeyWordWin, height= 10, width=40, wrap=TK.WORD)
    TextBox.grid(row=2, columnspan=2, sticky='nesw')
    #Add Buttons
    TK.Button(KeyWordWin, text='Look Up', width=10, command=lambda : LookUp(Keywords, Help, option, TextBox)).grid(row=5, column=1, pady=10, sticky = 'w')
    TK.Button(KeyWordWin, text='Close', width=10, command=lambda : Close()).grid(row=5, column=0, pady=10, sticky = 'e')

#Insert Example Equations
def VPExamples(self, Type):
    """
    Fills in text box with example equations
    """
    if Type == 'Arctan':
        Equ='2/np.pi * [AMPLITUDE] * np.arctan(Radius)'
        self.VPEqu.delete(0,TK.END)
        self.VPEqu.insert(0,Equ)
    elif Type == 'Natural Log':
        Equ='np.log(Radius)'
        self.VPEqu.delete(0,TK.END)
        self.VPEqu.insert(0,Equ) 
    
    #Insert Example Equations
def SBExamples(self, Type):
    """
    Fills in text box with example equations
    """
    if Type == 'Exponential':
        Equ='np.exp(-Radius/ [SCALE RADIUS] )'
        self.SBEqu.delete(0,TK.END)
        self.SBEqu.insert(0,Equ)
    elif Type == 'Gaussian':
        Equ='[AMPLITUDE] * gaussian(Radius, [PEAK CENTER], [STANDARD DEVATION])'
        self.SBEqu.delete(0,TK.END)
        self.SBEqu.insert(0,Equ) 
     
#Load SB .csv data
def SBLoadCSV(self):
     """
     Load Surface Brightness .csv
     """
     SB_CSVName = askopenfilename()
     self.SB_CSV = pd.read_csv(SB_CSVName, delimiter=',')
     self.SB_Rad, self.SB_Prof = self.SB_CSV['Radius'].values, self.SB_CSV['Surface Brightness'].values
     #Update Equation Box
     self.SBEqu.delete(0, TK.END)
     self.SBEqu.insert(0, 'File Loaded - %s' %SB_CSVName)     
     #Update StatusBox
     self.StatusBox.delete(0, TK.END)
     self.StatusBox.insert(0, 'Surface Brightness .CSV loaded')
        
#Load VP .csv data
def VPLoadCSV(self):
    """
    Load Velocity Profile .csv
    """
    VP_CSVName = askopenfilename()
    self.VP_CSV = pd.read_csv(VP_CSVName, delimiter=',')
    self.VP_Rad, self.VP_Prof = self.VP_CSV['Radius'].values, self.VP_CSV['Velocity Profile'].values
    #Update Equation Box
    self.VPEqu.delete(0, TK.END)
    self.VPEqu.insert(0, 'File Loaded - %s' %VP_CSVName)     
    #Update StatusBox
    self.StatusBox.delete(0, TK.END)
    self.StatusBox.insert(0, 'Velocity Profile .CSV loaded')

#########################################################################################################################################################
"""
Boot Window Class
"""
#########################################################################################################################################################

#Initial Boot window
class BootWindow:
    #Initialization 
    def __init__(self, master):
        self.master = master
        master.title('KinMS GUI')
        master.configure(background='black')

        #Load Photo
        self.im = TK.PhotoImage(file='Header.gif')
        self.Image = TK.Label(master, image=self.im, bg='black').grid(row=0, column=0, columnspan=2)   

        #Adding Buttons        
        TK.Button(master, text='Model Fitting Mode', width=20, command=lambda : SwitchMode(self, 'Read Mode')).grid(row=1, column=0, sticky="nesw")
        TK.Button(master, text='Custom Model Mode', width=20, command=lambda : SwitchMode(self, 'Custom Model')).grid(row=1, column=1, sticky="nesw")

#########################################################################################################################################################
"""
Model Fitting Class
"""
#########################################################################################################################################################

class ReadWin:
    #Initialization 
    def __init__(self, master):
        #Setting default Optimization settings
        self.method, self.MaxIt, self.Samp = 'Nelder-Mead', 500, 100000
        self.Done = False
        
        #Load window
        self.master = master
        master.title('KinMS GUI Model Fitting')
        master.configure(background='black')
        
        #Create Menu Bar
        menu = TK.Menu(root)
        master.config(menu=menu)
        filemenu = TK.Menu(menu)
        Optionmenu = TK.Menu(menu)
        Helpmenu = TK.Menu(menu)
        menu.add_cascade(label="File", menu=filemenu)
        menu.add_cascade(label="Options", menu=Optionmenu)
        menu.add_cascade(label="Help", menu=Helpmenu)
        
        #Add Options
        filemenu.add_command(label="Load Parameters", command=lambda : self.LoadInputs())
        filemenu.add_command(label="Save Parameters", command=lambda : self.SaveInputs())
        filemenu.add_separator()
        filemenu.add_command(label="View .fits Data", command=lambda : self.FitsInput())
        filemenu.add_separator()
        filemenu.add_command(label="Switch to Custom Model Mode", command=lambda : SwitchMode(self, 'Custom Model'))
        filemenu.delete(0)
        Optionmenu.add_command(label="Simulation Options", command=lambda : self.SimOptions())
        Optionmenu.delete(0)
        Helpmenu.add_command(label="Key Word Lookup", command=lambda : KeyWord(self))
        Helpmenu.add_command(label="Website Documentation", command=lambda : LaunchWebpage())
        Helpmenu.delete(0)
        
        #Load Photo
        self.im = TK.PhotoImage(file='Header.gif')
        self.Image = TK.Label(master, image=self.im, bg='black').grid(row=0, column=0, columnspan=4)   

        #Fits Data
        TK.Label(master, text='File selected:', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, pady=10, sticky="e")
        self.ShowFile = TK.Entry(master, width = 40, font = 'none 12')
        self.ShowFile.grid(row=1, column=1, columnspan = 2, sticky="e")
        TK.Button(master, text='Select .fits file', width=20, command=lambda : self.Loadfits()).grid(row=1, column=3, sticky='w')
          
        #Surface Brightness 
        TK.Label(master, text='Surface Brightness Equation', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, columnspan=4, pady=10, sticky="n")
        SBFrame = TK.Frame(master, bg='black')
        SBFrame.grid(row=4, columnspan=4)
        self.SBEqu = TK.Entry(SBFrame, width=70, bg='white')
        self.SBEqu.grid(row=1, column=0, columnspan=3, sticky='nswe')
        self.SBLoaded = False
        TK.Button(SBFrame, text='Load .csv', width=10, command=lambda : SBLoadCSV(self)).grid(row=1, column=4, sticky="w")
        SBSubFrame = TK.Frame(master, bg='black')
        SBSubFrame.grid(row=5, columnspan=4, sticky='n')
        TK.Button(SBSubFrame, text='Gaussian', width=10, command=lambda : SBExamples(self, Type='Gaussian')).grid(row=2, column=0, sticky="e")
        TK.Button(SBSubFrame, text='Exponential', width=10, command=lambda : SBExamples(self, Type='Exponential')).grid(row=2, column=1, sticky="w")
    
        #Velocity profile
        TK.Label(master, text='Velocity Profile Equation', bg='black', fg='white', font='none 12 bold').grid(row=6, column=0, columnspan=4, pady=10, sticky="n")
        VPFrame = TK.Frame(master, bg='black')
        VPFrame.grid(row=7, columnspan=4)
        self.VPEqu = TK.Entry(VPFrame, width=70, bg='white')
        self.VPEqu.grid(row=1, column=0, columnspan=3, sticky='nwse')
        self.VPLoaded = False
        TK.Button(VPFrame, text='Load .csv', width=10, command=lambda : VPLoadCSV(self)).grid(row=1, column=4, sticky="w")
        VPSubFrame = TK.Frame(master, bg='black')
        VPSubFrame.grid(row=8, columnspan=4, sticky='n')
        TK.Button(VPSubFrame, text='Arctan', width=10, command=lambda : VPExamples(self, Type='Arctan')).grid(row=0, column=0, sticky="ne")
        TK.Button(VPSubFrame, text='Natural Log', width=10, command=lambda : VPExamples(self, Type='Natural Log')).grid(row=0, column=1, sticky="nw")
        
        #Other parameters
        TK.Label(master, text='Other Parameters', bg='black', fg='white', font='none 12 bold').grid(row=9, column=0, columnspan=4, pady=10, sticky="n")
        OtherFrame = TK.Frame(master, bg='black')
        OtherFrame.grid(row=10, columnspan=4)
        TK.Label(OtherFrame, text='Disk Thickness (") =', bg='black', fg='white', font='none 12 bold').grid(row=0, column=0, sticky="ne")
        TK.Label(OtherFrame, text='Position Angle (deg) =', bg='black', fg='white', font='none 12 bold').grid(row=0, column=2, sticky="ne")
        TK.Label(OtherFrame, text='Inclination (deg) =', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, sticky="ne")
        TK.Label(OtherFrame, text='Mask Clip Level =', bg='black', fg='white', font='none 12 bold').grid(row=2, column=2, sticky="ne")
        TK.Label(OtherFrame, text='Velocity Dispersion (km/s) =', bg='black', fg='white', font='none 12 bold').grid(row=2, column=0, sticky="ne")
        TK.Label(OtherFrame, text='Initial Flux (Jy/km/s) =', bg='black', fg='white', font='none 12 bold').grid(row=1, column=2, sticky="ne")

        self.UserDiskThic = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserDiskThic.grid(row=0, column=1, sticky="nwse")
        self.UserDiskAng = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserDiskAng.grid(row=0, column=3, sticky="nwse")
        self.UserInc = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserInc.grid(row=1, column=1, sticky="nwse")
        self.UserRatio = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserRatio.grid(row=2, column=3, sticky="nwse")
        self.UserGasSig = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserGasSig.grid(row=2, column=1, sticky="nwse")
        self.UserFlux = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserFlux.grid(row=1, column=3, sticky="nwse")
                        
        #Lower Buttons 
        Lower = TK.Frame(master, bg='black')
        Lower.grid(row=11, columnspan=4)
        TK.Button(Lower, text='Run', width=15, command=lambda : self.Run()).grid(row=11, column=2, columnspan=2, pady=10, sticky = 'nw')
        TK.Button(Lower, text='Show Results', width=15, command=lambda : self.ShowPlot()).grid(row=11, column=0, columnspan=2, pady=10,sticky="ne")
        
        #Status Message
        self.StatusBox = TK.Entry(master, width=20, bg='white', justify='center')
        self.StatusBox.grid(row=12, columnspan=4, sticky="we")     
        self.StatusBox.insert(0, 'Waiting for User Input')
        
        #Progress Bar
        self.prog_bar = ttk.Progressbar(master, orient="horizontal", length=650, mode="determinate")
        self.prog_bar["value"] = 0
        self.prog_bar["maximum"] = 1
        self.prog_bar.grid(row = 13, columnspan=4)
            
    def SimOptions(self):
        """
        Open Window to edit simulation options
        """
        def Accept(self):
            #Save new values
            self.method = var.get()
            self.MaxIt = int(self.UserMaxIt.get())
            self.Samp = int(self.UserSamp.get())
            #Update StatusBox
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, 'Simulation Options Updated')
            #Close window
            SimOptInputWin.destroy()
        
        #Resets to default values
        def Reset(self):
            self.MaxIt, self.Samp = 500, 100000
            Inputs = [self.UserMaxIt, self.UserSamp]
            Data = [self.MaxIt, self.Samp]
            var.set(self.Methods[0])
            for i in range(len(Inputs)):
                Inputs[i].delete(0, TK.END)
                Inputs[i].insert(0, Data[i])
            
        #Create Window
        SimOptInputWin = TK.Tk()
        SimOptInputWin.configure(background='black')
        SimOptInputWin.title('KinMS GUI - Simulation Options')
        SimOptInputWin.configure(background='black')
        TK.Label(SimOptInputWin, text='Simulation Options', bg='black', fg='white', font='none 12 bold').grid(row=0, column=0, columnspan=2, sticky="n")
        TK.Label(SimOptInputWin, text='Method Of Minimsation = ', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, sticky="e")
        TK.Label(SimOptInputWin, text='Maximum Number or Iterations = ', bg='black', fg='white', font='none 12 bold').grid(row=2, column=0, sticky="e")
        TK.Label(SimOptInputWin, text='Number of KinMS Cloudlet Samples =', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, sticky="e")
        
        #Create dropdown box
        self.Methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP']
        var = TK.StringVar(SimOptInputWin)
        self.UserMethod = TK.OptionMenu(SimOptInputWin, var, *self.Methods)
        self.UserMethod.grid(row=1, column=1, sticky='nesw')
        
        #User Inputs
        self.UserMaxIt = TK.Entry(SimOptInputWin, width=20, bg='white')
        self.UserMaxIt.grid(row=2, column=1, sticky="nwse")
        self.UserSamp = TK.Entry(SimOptInputWin, width=20, bg='white')
        self.UserSamp.grid(row=3, column=1, sticky="nwse")
        
        #Input data
        Inputs = [self.UserMaxIt, self.UserSamp]
        Data = [self.MaxIt, self.Samp]
        var.set(self.method)
        for i in range(len(Inputs)):
            Inputs[i].delete(0, TK.END)
            Inputs[i].insert(0, Data[i])
        
        #Lower Buttons
        OptionBottomButtons = TK.Frame(SimOptInputWin, bg='black')
        OptionBottomButtons.grid(row=5, columnspan=2, sticky="n")
        TK.Button(OptionBottomButtons, text='Accept', width=10, command=lambda : Accept(self)).grid(row=11, column=1, pady=10, sticky = 'w')
        TK.Button(OptionBottomButtons, text='Load Default', width=10, command=lambda : Reset(self)).grid(row=11, column=0, pady=10, sticky = 'e') 
        
    #Edit Fits Values
    def FitsInput(self):
        """
        Open Window to edit .fits Values
        """
        #Update Function
        def UpdateFits(Data, Inputs):
            #Get user inputs
            for i in range(len(Data)):
                Data[i] = float(Inputs[i].get())
            #Update values
            self.xsize = Data[0]
            self.xRef = int(Data[1]) 
            self.xVal = Data[2]
            self.ysize = Data[3] 
            self.yRef = int(Data[4])
            self.yVal = Data[5]
            self.cellsize = Data[6]
            self.MaxRad = Data[7]
            self.vsize = Data[8]
            self.vRef = int(Data[9])
            self.vVal = Data[10]
            self.dv = Data[11]
            self.beamsize = np.array([Data[12], Data[13], Data[14]])
            self.VelStart, self.VelStop = int(Data[15]), int(Data[16])
            self.XStart, self.XStop = int(Data[17]), int(Data[18])
            self.YStart, self.YStop = int(Data[19]), int(Data[20])
            #Update StatusBox
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, '.fits data Updated')
            #Close Window
            FitsInputWin.destroy()  
            
        #Create data
        Data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        
        #Create Window
        FitsInputWin = TK.Tk()
        FitsInputWin.configure(background='black')
        FitsInputWin.title('KinMS GUI - Fits data')
        FitsInputWin.configure(background='black')
        FitsInputCenter = TK.Frame(FitsInputWin, bg='black')
        FitsInputCenter.grid(row=1, sticky="n")
        
        #Section Labels
        TK.Label(FitsInputCenter, text = 'Data from .fits', bg='black', fg='white', font='none 12 bold').grid(row=0, column=0, columnspan=4, sticky='n')
        TK.Label(FitsInputCenter, text = 'Values Derived from .fits', bg='black', fg='white', font='none 12 bold').grid(row=13, column=0, columnspan=4, sticky='n')
        
        #Make Labels 
        TK.Label(FitsInputCenter, text = 'X size (") = ', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'X Reference Cell = ', bg='black', fg='white', font='none 12 bold').grid(row=2, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'X Cell Value = ', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Y size (") = ', bg='black', fg='white', font='none 12 bold').grid(row=4, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Y Reference Cell = ', bg='black', fg='white', font='none 12 bold').grid(row=5, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Y Cell Value = ', bg='black', fg='white', font='none 12 bold').grid(row=6, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Cell Size ("/Pix) = ', bg='black', fg='white', font='none 12 bold').grid(row=7, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Number of Velocity Channels = ', bg='black', fg='white', font='none 12 bold').grid(row=8, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Velocity Reference Cell = ', bg='black', fg='white', font='none 12 bold').grid(row=9, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Velocity Cell Value (km/s) = ', bg='black', fg='white', font='none 12 bold').grid(row=10, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Channel Width (km/s/channel) = ', bg='black', fg='white', font='none 12 bold').grid(row=11, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Beam Size [BMAJ, BMIN, BPA] (deg) = ', bg='black', fg='white', font='none 12 bold').grid(row=12, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Radius (") = ', bg='black', fg='white', font='none 12 bold').grid(row=14, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Channel Start and Stop = ', bg='black', fg='white', font='none 12 bold').grid(row=15, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'X axis Spatial Start and Stop = ', bg='black', fg='white', font='none 12 bold').grid(row=16, column=0, sticky='e')
        TK.Label(FitsInputCenter, text = 'Y axis Spatial Start and Stop = ', bg='black', fg='white', font='none 12 bold').grid(row=17, column=0, sticky='e')
        
        #Creating Entry boxes
        InputXSize = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputXSize.grid(row=1, column=1, columnspan=3, sticky='w')
        InputXRef = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputXRef.grid(row=2, column=1, columnspan=3, sticky='w')
        InputXVal = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputXVal.grid(row=3, column=1, columnspan=3, sticky='w')
        InputYSize = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputYSize.grid(row=4, column=1, columnspan=3, sticky='w')
        InputYRef = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputYRef.grid(row=5, column=1, columnspan=3, sticky='w')
        InputYVal = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputYVal.grid(row=6, column=1, columnspan=3, sticky='w')
        InputCellSize = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputCellSize.grid(row=7, column=1, columnspan=3, sticky='w')
        InputMaxRad = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputMaxRad.grid(row=14, column=1, columnspan=3, sticky='w')
        InputVelChan = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputVelChan.grid(row=8, column=1, columnspan=3, sticky='w')
        InputVelRef = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputVelRef.grid(row=9, column=1, columnspan=3, sticky='w')
        InputVelCell = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputVelCell.grid(row=10, column=1, columnspan=3, sticky='w')
        InputVelWidth = TK.Entry(FitsInputCenter, width=32, bg='white')
        InputVelWidth.grid(row=11, column=1, columnspan=3, sticky='w')
        InputBMAJ = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputBMAJ.grid(row=12, column=1, sticky='w')
        InputBMIN = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputBMIN.grid(row=12, column=2, sticky='w')
        InputBPA = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputBPA.grid(row=12, column=3, sticky='w')
        InputVStart = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputVStart.grid(row=15, column=1, sticky='w')
        InputVStop = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputVStop.grid(row=15, column=2, sticky='w')
        InputXStart = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputXStart.grid(row=16, column=1, sticky='w')
        InputXStop = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputXStop.grid(row=16, column=2, sticky='w')
        InputYStart = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputYStart.grid(row=17, column=1, sticky='w')
        InputYStop = TK.Entry(FitsInputCenter, width=10, bg='white')
        InputYStop.grid(row=17, column=2, sticky='w')
        
        #Consolidating data
        try:
            Data = [self.xsize, self.xRef, self.xVal, self.ysize, self.yRef, self.yVal, \
                    self.cellsize, self.MaxRad, self.vsize, self.vRef, self.vVal, self.dv, \
                    self.beamsize[0], self.beamsize[1], self.beamsize[2], self.VelStart, self.VelStop, \
                    self.XStart, self.XStop, self.YStart, self.YStop]
        except AttributeError:
            Data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        Inputs = [InputXSize, InputXRef, InputXVal, InputYSize, InputYRef, InputYVal, InputCellSize, \
                  InputMaxRad, InputVelChan, InputVelRef, InputVelCell, InputVelWidth, InputBMAJ, \
                  InputBMIN, InputBPA, InputVStart , InputVStop, InputXStart, InputXStop, InputYStart, InputYStop]
        
        #Inputting Fits values to text boxes
        for i in range(len(Inputs)):
            Inputs[i].insert(0, Data[i])
        
        #Adding button to close window and update values
        TK.Button(FitsInputCenter, text='Accept', width=20, command=lambda : UpdateFits(Data, Inputs)).grid(row=18, column=2, columnspan=2, sticky="n")
        
    #Load Fits Data
    def Loadfits(self):
        """
        Load .fits files
        """
        #Select file
        filename = askopenfilename()
        self.ShowFile.delete(0, TK.END) 
        self.ShowFile.insert(0, filename.split('/')[-1]) 
        
        #Try to load .fits
        try:
            hdulist = fits.open(filename)
            #Update StatusBox
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, '.fits Loaded')
        #Produce error code if .fits is not compatible 
        except OSError:
            #Update StatusBox
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, 'File is not compatible')
        
        #Reading data for header
        self.scidata = hdulist[0].data.T
        self.cellsize = abs(hdulist[0].header['CDELT2']) * 3600  #arcseconds/pixel
        self.xsize = hdulist[0].header['NAXIS1'] * self.cellsize
        self.xRef = hdulist[0].header['CRPIX1']
        self.xVal = hdulist[0].header['CRVAL1']
        self.ysize = hdulist[0].header['NAXIS2'] * self.cellsize
        self.yRef = hdulist[0].header['CRPIX2']
        self.yVal = hdulist[0].header['CRVAL2']
        self.vsize = hdulist[0].header['NAXIS3']                        
        self.vRef = hdulist[0].header['CRPIX3']
        self.vVal = hdulist[0].header['CRVAL3'] / 1000
        self.dv = hdulist[0].header['CDELT3'] / 1000  #km/s 
        self.MaxRad = np.sqrt(self.xsize**2 + self.ysize**2)
        try:
            self.beamsize = np.array([hdulist[0].header['BMAJ'] * 3600, hdulist[0].header['BMIN'] * 3600, hdulist[0].header['BPA']])
        except KeyError:
            BeamData = hdulist[1].data
            self.beamsize = np.array([np.mean(BeamData['BMAJ']) * 3600, np.mean(BeamData['BMIN']) * 3600, np.mean(BeamData['BPA'])])
        
        #Removing Stokes axis if applicable 
        if self.scidata.ndim > 3:
            self.scidata = self.scidata[:,:,:,0]
        
        #Determin the relevant velocity channels
        try:
            self.VelStart, self.VelStop = AOI(self.scidata.copy())
        except NameError:
            self.VelStart, self.VelStop = 0, 0
        
        #Determin RMS
        self.RMS = np.sqrt(np.mean(self.scidata.copy()[:,:,0]**2))
        
        #Determin the relevant spacial area
        compressed = np.sum(self.scidata, axis=2)
        x, y = np.sum(compressed, axis=0), np.sum(compressed, axis=1)
        Allowence = np.array([len(x)*0.05, len(y)*0.05])
        x, y = x[int(len(x)*0.1):int(len(x)*0.9)], y[int(len(y)*0.1):int(len(y)*0.9)]  #Excluding outter 10%
        xCutoff, yCutoff = 1.5 * np.mean(x), 1.5 * np.mean(y)
        xAOI, yAOI = np.where(x > xCutoff), np.where(y > yCutoff)
        xBoundry, yBoundry = np.array([np.min(xAOI)-Allowence[0], np.max(xAOI)+Allowence[0]]), np.array([np.min(yAOI)-Allowence[1], np.max(yAOI)+Allowence[1]])
        if xBoundry[0] < 0: xBoundry[0] = 0
        if xBoundry[1] > len(x): xBoundry[1] = int(len(x)-1)
        if yBoundry[0] < 0: yBoundry[0] = 0
        if yBoundry[1] > len(y): yBoundry[1] = int(len(y)-1)
        self.XStart, self.XStop, self.YStart, self.YStop = int(xBoundry[0]+(len(x)*0.1)), int(xBoundry[1]+(len(x)*0.1)), int(yBoundry[0]+(len(y)*0.1)), int(yBoundry[1]+(len(y)*0.1))
        
    #Save paramiters as .txt
    def SaveInputs(self):
        """
        Saves Paramiters file
        """
        #Select location
        Location = TK.filedialog.askdirectory()
        Location = Location + '/KinMSFittingParameters_%i%i%i.txt' %(time.localtime()[3],time.localtime()[4],time.localtime()[5])
        #Wright file
        ParaData = open(Location, "w" )
        ParaData.writelines([' SB Equation¦%s' %self.SBEqu.get(),\
                             '\n VP Equation¦%s' %self.VPEqu.get(),\
                             '\n Disk Thickness¦%s' %self.UserDiskThic.get(), \
                             '\n Position Angle¦%s' %self.UserDiskAng.get(), \
                             '\n Inclination¦%s' %self.UserInc.get(), \
                             '\n Mask Clip Level¦%s' %self.UserRatio.get(),\
                             '\n Velocity Dispersion¦%s' %self.UserGasSig.get(),\
                             '\n Int Flux¦%s' %self.UserFlux.get()])
        ParaData.close()
        #Update StatusBox
        self.StatusBox.delete(0, TK.END)
        self.StatusBox.insert(0, 'Parameters Saved')    
    
    #Load paramiters .txt
    def LoadInputs(self):
        """
        Loads Paramiters file
        """
        #Loads data
        Parafilename = askopenfilename()
        ParaData = open(Parafilename, "r" )
        ParaData = ParaData.readlines()
        inputs = [self.SBEqu,\
                  self.VPEqu,\
                  self.UserDiskThic,\
                  self.UserDiskAng,\
                  self.UserInc,\
                  self.UserRatio,\
                  self.UserGasSig,\
                  self.UserFlux]
        count = 0
        #Inputs data into entry boxes 
        for i in range(len(ParaData)):
            if bool(ParaData[i].strip()) == True:
                inputs[count].delete(0,TK.END)
                inputs[count].insert(0,ParaData[i].split('¦')[-1])
                count += 1
                #Update StatusBox
        self.StatusBox.delete(0, TK.END)
        self.StatusBox.insert(0, 'Parameters loaded')
    
    #Run Fitting        
    def Fitting(self):
        """
        Execute simulation and minimisation
        """
        #Initial Time
        t0 = time.time()  
        #Reset progress bar
        self.prog_bar["value"] = 0
        #Define Fitting Function
        def Chisquare(Values):
            #Seperate Guess
            PosX, PosY, Vsys, Inc, DiskAng, GasSig, Flux, DiskThic, SBArray, VPArray = Values[0], Values[1], Values[2], Values[3], Values[4], Values[5], Values[6], Values[7], Values[8:self.SBlen+8], Values[self.SBlen+8:]
            #Edit User Equation
            Radius = self.Radius
            if self.SBEquString.split('-')[0] == 'File Loaded ':
                SB, SBRadius = self.SB_Prof, self.SB_Rad
            else:
                SB = eval(EditUserString(self.SBEquString, SBArray).replace('|',''))
                SBRadius = Radius
            if self.VPEquString.split('-')[0] == 'File Loaded ':
                VP, VPRadius = self.VP_Prof, self.VP_Rad
            else:
                VP = eval(EditUserString(self.VPEquString, VPArray).replace('|','')) 
                VPRadius = Radius
            #Create Simulation
            self.SimCube = KinMS(self.ClippedXsize, self.ClippedYsize, self.ClippedVsize, self.cellsize, self.dv, self.beamsize, Inc, gasSigma=GasSig,\
                                 sbProf=SB, sbRad=SBRadius, velProf=VP, velRad=VPRadius, diskThick=DiskThic, \
                                 nSamps=self.Samp, posAng=DiskAng, intFlux=Flux, vSys = Vsys, phaseCen=np.array([PosX,PosY]))
            #Determine Chi Square
            Chi = np.sum((self.CubeMask - self.SimCube)**2 / (self.RMS)**2)
            #Update progress bar
            Progress = self.run / self.MaxIt
            self.prog_bar["value"] = Progress
            self.run += 1
            #Save values in case of crash
            self.guess = Values
            return Chi
        #Get user input and bounds
        def GetUserInput(self):
            Userinput = self.UserDiskThic.get(), self.UserDiskAng.get(), self.UserInc.get(), self.UserGasSig.get(), self.UserFlux.get()    
            Ratio = int(self.UserRatio.get())
            defaultbounds = [(0,1e50), (0,360), (0,90), (0,100), (0,1e50)]
            Inputs, Bounds = [], []
            count = 0
            for i in Userinput:
                if ',' in i:
                    i = i.replace(' ','')
                    i = i.replace('(','')
                    i = i.replace(')','')
                    segment = i.split(',')
                    Inputs.append(float(segment[0]))
                    Bounds.append((float(segment[1]), float(segment[2])))
                else:
                    Inputs.append(float(i))
                    Bounds.append(defaultbounds[count])
                count += 1
            return Inputs[0], Inputs[1], Inputs[2], Inputs[3], Inputs[4], Ratio, Bounds
        #Read user inputs
        self.DiskThic, self.DiskAng, Inc, GasSig, Flux, Ratio, InputBounds = GetUserInput(self)
        #Make a copy of the data cube
        Cube = self.scidata.copy()[self.XStart:self.XStop,self.YStart:self.YStop,self.VelStart:self.VelStop]
        #Create cube mask
        Xsize, Ysize = len(Cube[:,0,0])*self.cellsize, len(Cube[0,:,0])*self.cellsize
        psf = makebeam(Xsize, Ysize,[self.beamsize[0]/self.cellsize,self.beamsize[1]/self.cellsize],rot=self.beamsize[2])
        self.CubeMask = GUI.smoothclip(Cube, Ratio, self.beamsize, self.cellsize, self.VelStart, self.VelStop, psf)
        #Clip Cube for simulation
        self.ClippedXsize = len(self.CubeMask[:,0,0]) * self.cellsize
        self.ClippedYsize = len(self.CubeMask[0,:,0]) * self.cellsize
        self.ClippedVsize = len(self.CubeMask[0,0,:]) * self.dv

        #Determin if any .CSV files are loaded
        #Get SB and VP data
        self.SBEquString, self.VPEquString = self.SBEqu.get(), self.VPEqu.get()
        self.Radius = np.arange(0, self.MaxRad, self.cellsize/2)
        #Both SB and VP Loaded
        SBLoaded, VPLoaded = self.SBEquString.split('-')[0], self.VPEquString.split('-')[0] 
        if SBLoaded == 'File Loaded ' and VPLoaded == 'File Loaded ':
            GuessList = [0,0,0, Inc, self.DiskAng, GasSig, Flux, self.DiskThic]
            BNDS = [(-self.xsize, self.xsize), (-self.ysize, self.ysize), (-1e50, 1e50)]
            for i in range(len(InputBounds)):
                BNDS.append(InputBounds[i])
            Guess = np.asarray(GuessList)
            self.SBlen = 0
        #VP Loaded
        elif SBLoaded != 'File Loaded ' and VPLoaded == 'File Loaded ':
            SBArray, SPBounds = ExtractVarFromStr(self.SBEquString)
            GuessList = [0,0,0, Inc, self.DiskAng, GasSig, Flux, self.DiskThic]
            BNDS = [(-self.xsize, self.xsize), (-self.ysize, self.ysize), (-1e50, 1e50)]
            for i in range(len(InputBounds)):
                BNDS.append(InputBounds[i])
            for i in range(len(SBArray)):
                GuessList.append(SBArray[i])
                BNDS.append(SPBounds[i])
            Guess = np.asarray(GuessList)
            self.SBlen = len(SBArray)
        #SB Loaded
        elif SBLoaded == 'File Loaded ' and VPLoaded != 'File Loaded ':
            VPArray, VPBounds = ExtractVarFromStr(self.VPEquString)
            GuessList = [0,0,0, Inc, self.DiskAng, GasSig, Flux, self.DiskThic]
            BNDS = [(-self.xsize, self.xsize), (-self.ysize, self.ysize), (-1e50, 1e50)]
            for i in range(len(InputBounds)):
                BNDS.append(InputBounds[i])
            for i in range(len(VPArray)):
                GuessList.append(VPArray[i])
                BNDS.append(VPBounds[i])
            Guess = np.asarray(GuessList)
            self.SBlen, self.VPlen = 0, len(VPArray)
        #None Loaded
        elif SBLoaded != 'File Loaded ' and VPLoaded != 'File Loaded ':
            #Creating Guess and Bound array
            VPArray, VPBounds = ExtractVarFromStr(self.VPEquString)
            SBArray, SPBounds = ExtractVarFromStr(self.SBEquString)
            GuessList = [0,0,0, Inc, self.DiskAng, GasSig, Flux, self.DiskThic]
            BNDS = [(-self.xsize, self.xsize), (-self.ysize, self.ysize), (-1e50, 1e50)]
            for i in range(len(InputBounds)):
                BNDS.append(InputBounds[i])
            for i in range(len(SBArray)):
                GuessList.append(SBArray[i])
                BNDS.append(SPBounds[i])
            for i in range(len(VPArray)):    
                GuessList.append(VPArray[i])
                BNDS.append(VPBounds[i])
            Guess = np.asarray(GuessList)
            #Save Array Length
            self.SBlen, self.VPlen = len(SBArray), len(VPArray)
                   
        #Minimization
        self.Error = False
        self.run = 0
        #No bounds
        if 'Nelder-Mead' or 'Powell'  == self.method:
            try:
                self.res = minimize(Chisquare, Guess, method = self.method, options={'maxfev': self.MaxIt}) 
            except ValueError:
                self.Error = True
                self.prog_bar["value"] = 0
                #Update StatusBox
                self.StatusBox.delete(0, TK.END)
                self.StatusBox.insert(0, 'Complete - Error in fitting')
        #No Maxfev Option 
        elif 'CG' or 'BFGS' or 'L-BFGS-B' == self.method:
            try:
                self.res = minimize(Chisquare, Guess, method = self.method, bounds = BNDS, options={'maxiter': self.MaxIt}) 
            except ValueError:
                self.Error = True                
                self.prog_bar["value"] = 0
                #Update StatusBox
                self.StatusBox.delete(0, TK.END)
                self.StatusBox.insert(0, 'Complete - Error in fitting')
        #Bounds
        elif 'TNC' or 'SLSQP' or 'COBYLA' == self.method:
            try:
                self.res = minimize(Chisquare, Guess, method = self.method, bounds = BNDS, options={'maxiter': self.MaxIt}) 
            except ValueError:
                self.Error = True
                self.prog_bar["value"] = 0
                #Update StatusBox
                self.StatusBox.delete(0, TK.END)
                self.StatusBox.insert(0, 'Complete - Error in fitting')        
        #Calculate time taken
        t4 = time.time()     
        self.TimeTaken = t4 - t0                                                   
        #State code is complete
        self.prog_bar["value"] = 1
        self.Done = True
        if self.Error == False:
            #Update StatusBox
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, 'Complete')

    #Run code
    def Run(self):
        """
        Creating a seperate thread for Fitting
        """
        #Update StatusBox
        self.StatusBox.delete(0, TK.END)
        self.StatusBox.insert(0, 'Running')
        #Update Done cond
        self.Done = False
        #Run Fitting on new thread
        self.thread = threading.Thread(target=self.Fitting, args=())
        self.thread.daemon = True             
        self.thread.start() 
        
    def ShowPlot(self):
        #Check if thread is complete
        SBLoaded, VPLoaded = self.SBEquString.split('-')[0], self.VPEquString.split('-')[0] 
        if self.Done == True:
            #Create Window
            OutputWindow = TK.Tk()
            OutputWindow.configure(background='black')
            OutputWindow.title('KinMS GUI - Fits data')
            OutputWindow.configure(background='black')
            #Add title
            TK.Label(OutputWindow, text = ('Results using the', self.method, 'Method'), bg='black', fg='white', font='none 12 bold').grid(row=0, column=0, columnspan=2, sticky='n')
            #Create labels
            TK.Label(OutputWindow, text = 'Iterations = ', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Success = ', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Phase Center (X,Y) = ', bg='black', fg='white', font='none 12 bold').grid(row=4, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Vsys = ', bg='black', fg='white', font='none 12 bold').grid(row=5, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Inc = ', bg='black', fg='white', font='none 12 bold').grid(row=6, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Disk Ang = ', bg='black', fg='white', font='none 12 bold').grid(row=7, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Gas Sigma = ', bg='black', fg='white', font='none 12 bold').grid(row=8, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Flux = ', bg='black', fg='white', font='none 12 bold').grid(row=9, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Disk Thickness = ', bg='black', fg='white', font='none 12 bold').grid(row=10, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'SB Equation = ', bg='black', fg='white', font='none 12 bold').grid(row=11, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'VP Equation = ', bg='black', fg='white', font='none 12 bold').grid(row=12, column=0, sticky='ne')
            TK.Label(OutputWindow, text = 'Simulation Time = ', bg='black', fg='white', font='none 12 bold').grid(row=13, column=0, sticky='ne')
            #State Guesses incase of error
            if self.Error == True:
                TK.Label(OutputWindow, text = self.run, bg='black', fg='white', font='none 12 bold').grid(row=1, column=1, sticky='nw')
                TK.Label(OutputWindow, text = 'Failed Due to Fitting error', bg='black', fg='white', font='none 12 bold').grid(row=3, column=1, sticky='nw')
                TK.Label(OutputWindow, text = (self.guess[0], self.guess[1]), bg='black', fg='white', font='none 12 bold').grid(row=4, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.guess[2], bg='black', fg='white', font='none 12 bold').grid(row=5, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.guess[3], bg='black', fg='white', font='none 12 bold').grid(row=6, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.guess[4], bg='black', fg='white', font='none 12 bold').grid(row=7, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.guess[5], bg='black', fg='white', font='none 12 bold').grid(row=8, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.guess[6], bg='black', fg='white', font='none 12 bold').grid(row=9, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.guess[7], bg='black', fg='white', font='none 12 bold').grid(row=10, column=1, sticky='nw')
                if SBLoaded == 'File Loaded ':
                    TK.Label(OutputWindow, text = '.CSV Data Used', bg='black', fg='white', font='none 12 bold').grid(row=11, column=1, sticky='nw')
                else:
                    TK.Label(OutputWindow, text = EditUserString(self.SBEquString, self.guess[8:self.SBlen+8]), bg='black', fg='white', font='none 12 bold').grid(row=11, column=1, sticky='nw')
                if VPLoaded == 'File Loaded ':
                    TK.Label(OutputWindow, text = '.CSV Data Used', bg='black', fg='white', font='none 12 bold').grid(row=11, column=1, sticky='nw')
                else:
                    TK.Label(OutputWindow, text = EditUserString(self.VPEquString, self.guess[self.SBlen+8:]), bg='black', fg='white', font='none 12 bold').grid(row=12, column=1, sticky='nw')
                TK.Label(OutputWindow, text = '%s Seconds' %self.TimeTaken, bg='black', fg='white', font='none 12 bold').grid(row=13, column=1, sticky='nw')
            #State Minimisation Values
            else:
                TK.Label(OutputWindow, text = self.run, bg='black', fg='white', font='none 12 bold').grid(row=1, column=1, sticky='nw')
                TK.Label(OutputWindow, text = '%s' %self.res.success, bg='black', fg='white', font='none 12 bold').grid(row=3, column=1, sticky='nw')
                TK.Label(OutputWindow, text = ('(', self.res.x[0], ',', self.res.x[1], ')'), bg='black', fg='white', font='none 12 bold').grid(row=4, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.res.x[2], bg='black', fg='white', font='none 12 bold').grid(row=5, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.res.x[3], bg='black', fg='white', font='none 12 bold').grid(row=6, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.res.x[4], bg='black', fg='white', font='none 12 bold').grid(row=7, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.res.x[5], bg='black', fg='white', font='none 12 bold').grid(row=8, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.res.x[6], bg='black', fg='white', font='none 12 bold').grid(row=9, column=1, sticky='nw')
                TK.Label(OutputWindow, text = self.res.x[7], bg='black', fg='white', font='none 12 bold').grid(row=10, column=1, sticky='nw')
                if SBLoaded == 'File Loaded ':
                    TK.Label(OutputWindow, text = '.CSV Data Used', bg='black', fg='white', font='none 12 bold').grid(row=11, column=1, sticky='nw')
                else:
                    TK.Label(OutputWindow, text = EditUserString(self.SBEquString, self.guess[8:self.SBlen+8]), bg='black', fg='white', font='none 12 bold').grid(row=11, column=1, sticky='nw')
                if VPLoaded == 'File Loaded ':
                    TK.Label(OutputWindow, text = '.CSV Data Used', bg='black', fg='white', font='none 12 bold').grid(row=11, column=1, sticky='nw')
                else:
                    TK.Label(OutputWindow, text = EditUserString(self.VPEquString, self.guess[self.SBlen+8:]), bg='black', fg='white', font='none 12 bold').grid(row=12, column=1, sticky='nw')
                TK.Label(OutputWindow, text = '%.3f Seconds' %self.TimeTaken, bg='black', fg='white', font='none 12 bold').grid(row=13, column=1, sticky='nw')
        #Plot latest guesses
        try:
            GUI.makeplots(self.SimCube, self.ClippedXsize, self.ClippedYsize, self.ClippedVsize, self.cellsize,\
                          self.dv, self.beamsize ,posang=self.DiskAng, overcube=self.CubeMask, pvdthick=self.DiskThic)      
        #Incase of error
        except AttributeError:
            #Update StatusBox
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, 'Fitting must be run before results displayed')
            
#########################################################################################################################################################
"""
Custom Model Class
"""
#########################################################################################################################################################

class CustomWin:
    #Initialisation
    def __init__(self, master):
       
        self.Done = False
        #Setting default aditional inputs settings
        self.GasSig, self.Flux, self.Vsys, self.DiskThick, self.DiskAng, self.Radius, self.PhaseX, self.PhaseY, self.Samps = 0, 0, 0, 0, 0, 0, 0, 0, int(1e5) 
        #Load window
        self.master = master
        master.title('KinMS GUI Custom Model')
        master.configure(background='black')
        
        #Create Menu Bar
        menu = TK.Menu(root)
        master.config(menu=menu)
        filemenu = TK.Menu(menu)
        Helpmenu = TK.Menu(menu)
        menu.add_cascade(label="File", menu=filemenu)
        menu.add_cascade(label="Help", menu=Helpmenu)
        #Add options
        filemenu.add_command(label="Load Parameters", command=lambda : self.LoadInputs())
        filemenu.add_command(label="Save Parameters", command=lambda : self.SaveInputs())
        filemenu.add_separator()
        filemenu.add_command(label="Switch to Model Fitting Mode", command=lambda : SwitchMode(self, 'Read Mode'))
        Helpmenu.add_command(label="Key Word Lookup", command=lambda : KeyWord(self))
        Helpmenu.add_command(label="Website Documentation", command=lambda : LaunchWebpage())
        
        #Load Photo
        self.im = TK.PhotoImage(file='Header.gif')
        self.Image = TK.Label(master, image=self.im, bg='black').grid(row=0, column=0, columnspan=4)  
        
        #Labels
        TK.Label(master, text = 'X size (") = ', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, sticky='ne')
        TK.Label(master, text = 'Y size (") = ', bg='black', fg='white', font='none 12 bold').grid(row=2, column=0, sticky='ne')
        TK.Label(master, text = 'V size (km/s) = ', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, sticky='ne')
        TK.Label(master, text = 'Cellsize ("/Pix) = ', bg='black', fg='white', font='none 12 bold').grid(row=1, column=2, sticky='ne')
        TK.Label(master, text = 'Inclination (deg) = ', bg='black', fg='white', font='none 12 bold').grid(row=2, column=2, sticky='ne')
        TK.Label(master, text = 'Channel Width (km/s/channel) = ', bg='black', fg='white', font='none 12 bold').grid(row=3, column=2, sticky='ne')
        
        self.UserXsize = TK.Entry(master, width=20, bg='white')
        self.UserXsize.grid(row=1, column=1, sticky="nwse")
        self.UserYsize = TK.Entry(master, width=20, bg='white')
        self.UserYsize.grid(row=2, column=1, sticky="nwse")
        self.UserVsize = TK.Entry(master, width=20, bg='white')
        self.UserVsize.grid(row=3, column=1, sticky="nwse") 
        self.UserCellSize = TK.Entry(master, width=20, bg='white')
        self.UserCellSize.grid(row=1, column=3, sticky="nwse")
        self.UserInc = TK.Entry(master, width=20, bg='white')
        self.UserInc.grid(row=2, column=3, sticky="nwse") 
        self.Userdv = TK.Entry(master, width=20, bg='white')
        self.Userdv.grid(row=3, column=3, sticky="nwse")  
               
        #Surface Brightness 
        TK.Label(master, text='Surface Brightness Equation', bg='black', fg='white', font='none 12 bold').grid(row=4, column=0, columnspan=4, pady=10, sticky="n")
        SBFrame = TK.Frame(master, bg='black', width=600, height=300)
        SBFrame.grid(row=5, columnspan=4)
        self.SBEqu = TK.Entry(SBFrame, width=70, bg='white')
        self.SBEqu.grid(row=1, column=0, columnspan=3, sticky='nswe')
        SBSubFrame = TK.Frame(SBFrame, bg='black', width=600, height=300)
        SBSubFrame.grid(row=6, columnspan=4, sticky='n')
        TK.Button(SBSubFrame, text='Gaussian', width=10, command=lambda : SBExamples(self, Type='Gaussian')).grid(row=2, column=0, sticky="e")
        TK.Button(SBSubFrame, text='Exponential', width=10, command=lambda : SBExamples(self, Type='Exponential')).grid(row=2, column=1, sticky="w")
    
        #Velocity profile
        TK.Label(master, text='Velocity Profile Equation', bg='black', fg='white', font='none 12 bold').grid(row=6, column=0, columnspan=4, pady=10, sticky="n")
        VPFrame = TK.Frame(master, bg='black', width=600, height=300)
        VPFrame.grid(row=7, columnspan=4)
        self.VPEqu = TK.Entry(VPFrame, width=70, bg='white')
        self.VPEqu.grid(row=1, column=0, columnspan=3, sticky='nwse')
        VPSubFrame = TK.Frame(VPFrame, bg='black', width=600, height=300)
        VPSubFrame.grid(row=5, columnspan=4, sticky='n')
        TK.Button(VPSubFrame, text='Arctan', width=10, command=lambda : VPExamples(self, Type='Arctan')).grid(row=0, column=0, sticky="n")
        TK.Button(VPSubFrame, text='Natural Log', width=10, command=lambda : VPExamples(self, Type='Natural Log')).grid(row=0, column=1, sticky="nw")
        
        #Beamsize
        BeamSize = TK.Frame(master, bg='black', width=600, height=300)
        BeamSize.grid(row=8, columnspan=4, pady = 10, sticky='n')
        TK.Label(BeamSize, text = 'Beamsize [BMAJ, BMIN, BPA] (deg) = ', bg='black', fg='white', font='none 12 bold').grid(row=8, column=0, sticky='ne')
        self.UserBMAJ = TK.Entry(BeamSize, width=20, bg='white')
        self.UserBMAJ.grid(row=8, column=1, sticky="n")
        self.UserBMIN = TK.Entry(BeamSize, width=20, bg='white')
        self.UserBMIN.grid(row=8, column=2, sticky="n")
        self.UserBPA = TK.Entry(BeamSize, width=20, bg='white')
        self.UserBPA.grid(row=8, column=3, sticky="n")
        
        #Lower Buttons 
        LowerButton = TK.Frame(master, bg='black', width=600, height=300)
        LowerButton.grid(row=9, columnspan=4, sticky='n')
        TK.Button(LowerButton, text='Run', width=15, command=lambda : self.Run()).grid(row=11, column=2, columnspan=2, pady=10, sticky = 'nw')
        TK.Button(LowerButton, text='Aditional Inputs', width=15, command=lambda : self.AditionalInputs()).grid(row=11, column=0, columnspan=2, pady=10,sticky="ne")
        
        #Status Message
        self.StatusBox = TK.Entry(master, width=20, bg='white', justify='center')
        self.StatusBox.grid(row=12, columnspan=4, sticky="we")     
        self.StatusBox.insert(0, 'Waiting for User Input')
        
        
      #Edit Fits Values
    def AditionalInputs(self):
        """
        Open Window to edit .fits Values
        """
        #Reset to default
        def Default(Inputs):
            Data = [0, 0, 0, 0, 0, 0, 0, 0, int(1e5)]
            for i in range(len(Data)):
                Inputs[i].delete(0, TK.END)
                Inputs[i].insert(0, Data[i])
         
        #Update Function
        def Accept(Data, Inputs):
            #Get user inputs
            for i in range(len(Data)):
                Data[i] = float(Inputs[i].get())
            #Update values
            self.GasSig = Data[0]
            self.Flux = Data[1]
            self.Vsys = Data[2]
            self.DiskThick = Data[3]
            self.DiskAng = Data[4]
            self.Radius = Data[5]
            self.PhaseX, self.PhaseY = Data[6], Data[7]
            self.Samps = int(Data[8])
            #Update StatusBox
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, 'Aditional Inputs Added')
            #Close Window
            AditionalInputs.destroy()  
        
        #Create Window
        AditionalInputs = TK.Tk()
        AditionalInputs.configure(background='black')
        AditionalInputs.title('KinMS GUI - Aditional Inputs')
        AditionalInputs.configure(background='black')
        
        #Creating Labels
        TK.Label(AditionalInputs, text = 'Velocity Dispersion (km/s) = ', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, sticky='ne')
        TK.Label(AditionalInputs, text = 'Initial Flux (Jy/km/s) = ', bg='black', fg='white', font='none 12 bold').grid(row=2, column=0, sticky='ne')
        TK.Label(AditionalInputs, text = 'Systemic velocity (km/s) = ', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, sticky='ne')
        TK.Label(AditionalInputs, text = 'Disk Thickness (") = ', bg='black', fg='white', font='none 12 bold').grid(row=1, column=2, sticky='ne')
        TK.Label(AditionalInputs, text = 'Position Angle (deg) = ', bg='black', fg='white', font='none 12 bold').grid(row=2, column=2, sticky='ne')
        TK.Label(AditionalInputs, text = 'Radius (") = ', bg='black', fg='white', font='none 12 bold').grid(row=3, column=2, sticky='ne')
        
        UserGasSig = TK.Entry(AditionalInputs, width=20, bg='white')
        UserGasSig.grid(row=1, column=1, sticky="n")
        UserFlux = TK.Entry(AditionalInputs, width=20, bg='white')
        UserFlux.grid(row=2, column=1, sticky="n")
        UserVsys = TK.Entry(AditionalInputs, width=20, bg='white')
        UserVsys.grid(row=3, column=1, sticky="n") 
        UserDiskThick = TK.Entry(AditionalInputs, width=20, bg='white')
        UserDiskThick.grid(row=1, column=3, sticky="n")
        UserDiskAng = TK.Entry(AditionalInputs, width=20, bg='white')
        UserDiskAng.grid(row=2, column=3, sticky="n") 
        UserRadius = TK.Entry(AditionalInputs, width=20, bg='white')
        UserRadius.grid(row=3, column=3, sticky="n")          
        
        #Phase Center
        TK.Label(AditionalInputs, text = 'Morphological Centre = ', bg='black', fg='white', font='none 12 bold').grid(row=4, columnspan=2, sticky='ne')
        UserPhaseX = TK.Entry(AditionalInputs, width=20, bg='white')
        UserPhaseX.grid(row=4, column=2, sticky="n")
        UserPhaseY = TK.Entry(AditionalInputs, width=20, bg='white')
        UserPhaseY.grid(row=4, column=3, sticky="n")
        
        #Samps
        TK.Label(AditionalInputs, text = 'Number of KinMS Samples = ', bg='black', fg='white', font='none 12 bold').grid(row=5, columnspan=2, sticky='ne')
        UserSamps = TK.Entry(AditionalInputs, width=20, bg='white')
        UserSamps.grid(row=5, column=2, columnspan=2, sticky="n")
        
        #Consolidating data
        try:
            Data = [self.GasSig, self.Flux, self.Vsys, self.DiskThick, self.DiskAng, self.Radius, self.PhaseX, self.PhaseY, self.Samps]
        except AttributeError:
            Data = [0, 0, 0, 0, 0, 0, 0, 0, int(1e5)]
        Inputs = [UserGasSig, UserFlux, UserVsys, UserDiskThick, UserDiskAng, UserRadius, UserPhaseX, UserPhaseY, UserSamps]
        
        #Inputting Fits values to text boxes
        for i in range(len(Inputs)):
            Inputs[i].insert(0, Data[i])
        
        #Adding button to close window and update values
        TK.Button(AditionalInputs, text='Reset to Default', width=20, command=lambda : Default(Inputs)).grid(row=18, column=0, columnspan=2, sticky="ne")      
        TK.Button(AditionalInputs, text='Accept', width=20, command=lambda : Accept(Data, Inputs)).grid(row=18, column=2, columnspan=2, sticky="nw")      
    
    #Load paramiters .txt
    def LoadInputs(self):
        """
        Loads Paramiters file
        """
        #Loads data
        Parafilename = askopenfilename()
        ParaData = open(Parafilename, "r" )
        ParaData = ParaData.readlines()
        inputs = [self.UserXsize,\
                  self.UserYsize,\
                  self.UserVsize,\
                  self.UserCellSize,\
                  self.UserInc,\
                  self.Userdv,\
                  self.UserBMAJ,\
                  self.UserBMIN,\
                  self.UserBPA,\
                  self.SBEqu,\
                  self.VPEqu]
        Values = [0,0,0,0,0,0,0,0,0]
        count = 0
        i = 0
        #Inputs data into entry boxes
        while count < len(inputs):
            if bool(ParaData[i].strip()) == True:
                inputs[count].delete(0,TK.END)
                inputs[count].insert(0,ParaData[i].split('¦')[-1])
                count += 1
            i += 1
        count = 0
        while count < len(Values):
            if bool(ParaData[i].strip()) == True:
                Values[count] = ParaData[i].split('¦')[-1].replace('\n','')
                count += 1
            i += 1
        self.GasSig, self.Flux, self.Vsys, self.DiskThick, self.DiskAng, self.Radius, self.PhaseX, self.PhaseY, self.Samps = float(Values[0]), float(Values[1]), float(Values[2]), float(Values[3]), float(Values[4]), float(Values[5]), float(Values[6]), float(Values[7]), int(Values[8])
        #Update StatusBox
        self.StatusBox.delete(0, TK.END)
        self.StatusBox.insert(0, 'Parameters loaded')
    
 
    #Save paramiters as .txt
    def SaveInputs(self):
        """
        Saves Paramiters file
        """
        #Select location
        Location = TK.filedialog.askdirectory()
        Location = Location + '/KinMSCustomModel_%i%i%i.txt' %(time.localtime()[3],time.localtime()[4],time.localtime()[5])
        #Wright file
        ParaData = open(Location, "w" )
        ParaData.writelines([' Xsize¦%s' %self.UserXsize.get(),\
                             '\n Ysize¦%s' %self.UserYsize.get(),\
                             '\n Vsize¦%s' %self.UserVsize.get(),\
                             '\n CellSize¦%s' %self.UserCellSize.get(),\
                             '\n Inc¦%s' %self.UserInc.get(),\
                             '\n dv¦%s' %self.Userdv.get(),\
                             '\n BMJ¦%s' %self.UserBMAJ.get(),\
                             '\n BMIN¦%s' %self.UserBMIN.get(),\
                             '\n BPA¦%s' %self.UserBPA.get(),\
                             '\n SB Equation¦%s' %self.SBEqu.get(),\
                             '\n VP Equation¦%s' %self.VPEqu.get(),\
                             
                             '\n Disk Thickness¦%s' %self.DiskThick,\
                             '\n Position Angle¦%s' %self.DiskAng,\
                             '\n Velocity Dispersion¦%s' %self.GasSig,\
                             '\n Int Flux¦%s' %self.Flux,\
                             '\n Vsys¦%s' %self.Vsys,\
                             '\n Radius¦%s' %self.Radius,\
                             '\n PhaseX¦%s' %self.PhaseX,\
                             '\n PhaseY¦%s' %self.PhaseY,\
                             '\n Samps¦%s' %self.Samps])
        ParaData.close()
        #Update StatusBox
        self.StatusBox.delete(0, TK.END)
        self.StatusBox.insert(0, 'Parameters Saved')   

    def Run(self):
        #Read User Inputs
        #Load X and Y Sizes
        Xsize = float(self.UserXsize.get())
        Ysize = float(self.UserYsize.get())
        #Cheack if Radius need to be found
        if self.Radius == 0:
            self.Radius = np.sqrt(Xsize**2 + Ysize**2)
        #Get last of user input
        Vsize = float(self.UserVsize.get())
        CellSize = float(self.UserCellSize.get())
        print(self.Radius)
        print(type(self.Radius))
        Radius = np.arange(0, float(self.Radius), CellSize/2)
        print(self.Radius)
        print(type(self.Radius))
        Inc = float(self.UserInc.get())
        dv = float(self.Userdv.get())
        BMAJ, BMIN, BPA = float(self.UserBMAJ.get()), float(self.UserBMIN.get()), float(self.UserBPA.get())
        Beam = np.array([BMAJ, BMIN, BPA])
        SB = eval(self.SBEqu.get())
        VP = eval(self.VPEqu.get())
        #Preform KinMS Simulation
        SimCube = KinMS(Xsize, Ysize, Vsize, CellSize, dv, Beam, Inc, gasSigma=self.GasSig,\
                             sbProf=SB, sbRad=Radius, velProf=VP, velRad=Radius, diskThick=self.DiskThick, \
                             nSamps=self.Samps, posAng=self.DiskAng, intFlux=self.Flux, vSys = self.Vsys, phaseCen=np.array([self.PhaseX,self.PhaseY]))
        #Plot Results
        try:
            GUI.makeplots(SimCube, Xsize, Ysize, Vsize, CellSize,\
                          dv, Beam ,posang=self.DiskAng, pvdthick=self.DiskThick)  
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, 'Done')  
        except ValueError:
            self.StatusBox.delete(0, TK.END)
            self.StatusBox.insert(0, 'Additional inputs are required for more plots')  

#########################################################################################################################################################
"""
Init GUI
"""
#########################################################################################################################################################

root = TK.Tk()  #Change to TK.Tk() for running ouside of spyder and TK.Toplevel() for inside
my_gui = BootWindow(root)
root.mainloop()


