import os

os.chdir("GUI Resources/")

import tkinter as TK 
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import numpy as np
import pandas as pd
import time
from KinMS import KinMS
from makebeam import makebeam
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import GUIFunctions as GUI

import threading


#from IPython import get_ipython
#ipython = get_ipython()
#ipython.magic("gui tk")

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

def EditUserString(String, ValueArray, ReturnUser = False):
    '''
    Takes the user input equation and adjust values for accurate fitting
    '''
    Index = []
    i = String.find('|')
    while i >= 0 :
        Index.append(i)
        i = String.find('|', i+1)
    if len(Index) == 0:
        return String
    elif len(Index)//2 == 1:
        if ValueArray.size > 1:
            ValueArray = ValueArray[0]
        NewString = String[:Index[0]+1] + '%f' + String[Index[1]:]
        NewString = NewString %ValueArray
        Init = np.array(eval(String[Index[0]+1:Index[1]]))
    elif len(Index)//2 == 2:
        NewString = String[:Index[0]] + '%f' + String[Index[1]+1:Index[2]] +\
                                        '%f' + String[Index[3]+1:]
        NewString =  NewString %(ValueArray[0], ValueArray[1])
        Init = np.array([eval(String[Index[0]+1:Index[1]]), eval(String[Index[2]+1:Index[3]])])  
    elif len(Index)//2 == 3:
        NewString = String[:Index[0]] + '%f' + String[Index[1]+1:Index[2]] +\
                                        '%f' + String[Index[3]+1:Index[4]] +\
                                        '%f' + String[Index[5]+1:]
        NewString =  NewString %(ValueArray[0], ValueArray[1], ValueArray[2])
        Init = np.array([eval(String[Index[0]+1:Index[1]]), eval(String[Index[2]+1:Index[3]]), eval(String[Index[4]+1:Index[5]])])
    elif len(Index)//2 == 4:
        NewString = String[:Index[0]] + '%f' + String[Index[1]+1:Index[2]] +\
                                        '%f' + String[Index[3]+1:Index[4]] +\
                                        '%f' + String[Index[5]+1:Index[6]] +\
                                        '%f' + String[Index[7]+1]
        NewString =  NewString %(ValueArray[0], ValueArray[1], ValueArray[2], ValueArray[3])   
        Init = np.array([eval(String[Index[0]+1:Index[1]]), eval(String[Index[2]+1:Index[3]]), eval(String[Index[4]+1:Index[5]]), eval(String[Index[6]+1:Index[7]])])
    elif len(Index)//2 == 5:
        NewString = String[:Index[0]] + '%f' + String[Index[1]+1:Index[2]] +\
                                        '%f' + String[Index[3]+1:Index[4]] +\
                                        '%f' + String[Index[5]+1:Index[6]] +\
                                        '%f' + String[Index[7]+1:Index[8]] +\
                                        '%f' + String[Index[9]+1]
        NewString =  NewString %(ValueArray[0], ValueArray[1], ValueArray[2], ValueArray[3], ValueArray[4])
        Init = np.array([eval(String[Index[0]+1:Index[1]]), eval(String[Index[2]+1:Index[3]]), eval(String[Index[4]+1:Index[5]]), eval(String[Index[6]+1:Index[7]]), eval(String[Index[8]+1:Index[9]])]) 
    if ReturnUser == False:
        return NewString
    else:
        return Init


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
        TK.Button(master, text='Read .fits Mode', width=20, command=lambda : self.OpenWindow('Read Mode')).grid(row=1, column=0, sticky="nesw")
        TK.Button(master, text='Custom Model Mode', width=20, command=lambda : self.OpenWindow('Custom Model')).grid(row=1, column=1, sticky="nesw")
    #Opening windows
    def OpenWindow(self, Win):
        if Win == 'Read Mode':
            self.master.withdraw()
            self.newWindow = TK.Toplevel(self.master)
            self.app = ReadWin(self.newWindow)
        if Win == 'Custom Model':
            self.master.withdraw()
            self.newWindow = TK.Toplevel(self.master)
            self.app = CustomWin(self.newWindow)


class ReadWin:
    #Initialization 
    def __init__(self, master):
        #Setting default Optimization settings
        self.method, self.MaxIt, self.Samp = 'Nelder-Mead', 500, 100000
        #Load window
        self.master = master
        master.title('KinMS GUI Read fits')
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
        
        filemenu.add_command(label="Load Parameters", command=lambda : self.LoadInputs())
        filemenu.add_command(label="Save Parameters", command=lambda : self.SaveInputs())
        filemenu.add_separator()
        filemenu.add_command(label="View .fits Data", command=lambda : self.FitsInput())
        
        Optionmenu.add_command(label="Simulation Options", command=lambda : self.SimOptions())
        
        Helpmenu.add_command(label="Key Word Lookup", command=lambda : self.KeyWord())
        Helpmenu.add_command(label="Website Documentation", command=lambda : self.SimOptions())
        
        
        #Load Photo
        self.im = TK.PhotoImage(file='Header.gif')
        self.Image = TK.Label(master, image=self.im, bg='black').grid(row=0, column=0, columnspan=4)   

        #Fits Data
        TK.Label(master, text='File selected:', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, pady=10, sticky="e")
        self.ShowFile = TK.Entry(master, width = 30, font = 'none 12')
        self.ShowFile.grid(row=1, column=1, columnspan = 2, sticky="e")
        TK.Button(master, text='Select .fits file', width=20, command=lambda : self.Loadfits()).grid(row=1, column=3, sticky='w')
          
        #Surface Brightness 
        TK.Label(master, text='Surface Brigtness Equation', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, columnspan=4, pady=10, sticky="n")
        SBFrame = TK.Frame(master, bg='black', width=600, height=300)
        SBFrame.grid(row=4, columnspan=4)
        self.SBEqu = TK.Entry(SBFrame, width=70, bg='white')
        self.SBEqu.grid(row=1, column=0, columnspan=3, sticky='nse')
        self.SBLoaded = False
        TK.Button(SBFrame, text='Load .csv', width=10, command=lambda : self.SBLoadCSV()).grid(row=1, column=4, sticky="w")
        SBSubFrame = TK.Frame(master, bg='black', width=600, height=300)
        SBSubFrame.grid(row=5, columnspan=4, sticky='n')
        TK.Button(SBSubFrame, text='Gaussian', width=10, command=lambda : self.SBExamples(Type='Gaussian')).grid(row=2, column=0, sticky="e")
        TK.Button(SBSubFrame, text='Exponential', width=10, command=lambda : self.SBExamples(Type='Exponential')).grid(row=2, column=1, sticky="w")
    
        #Velocity profile
        TK.Label(master, text='Velocity Profile Equation', bg='black', fg='white', font='none 12 bold').grid(row=6, column=0, columnspan=4, pady=10, sticky="n")
        VPFrame = TK.Frame(master, bg='black', width=600, height=300)
        VPFrame.grid(row=7, columnspan=4)
        self.VPEqu = TK.Entry(VPFrame, width=70, bg='white')
        self.VPEqu.grid(row=1, column=0, columnspan=3, sticky='nse')
        self.VPLoaded = False
        TK.Button(VPFrame, text='Load .csv', width=10, command=lambda : self.VPLoadCSV()).grid(row=1, column=4, sticky="w")
        VPSubFrame = TK.Frame(VPFrame, bg='black', width=600, height=300)
        VPSubFrame.grid(row=5, columnspan=4, sticky='n')
        TK.Button(VPSubFrame, text='Arctan', width=10, command=lambda : self.VPExamples(Type='Arctan')).grid(row=0, column=0, sticky="n")
        
        #Other parameters
        TK.Label(master, text='Other Parameters', bg='black', fg='white', font='none 12 bold').grid(row=8, column=0, columnspan=4, pady=10, sticky="n")
        OtherFrame = TK.Frame(master, bg='black', width=600, height=300)
        OtherFrame.grid(row=9, columnspan=4)
        TK.Label(OtherFrame, text='Disk Thickness (") =', bg='black', fg='white', font='none 12 bold').grid(row=0, column=0, sticky="ne")
        TK.Label(OtherFrame, text='Angle of Disk (deg) =', bg='black', fg='white', font='none 12 bold').grid(row=0, column=2, sticky="ne")
        TK.Label(OtherFrame, text='Inclination (deg) =', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, sticky="ne")
        TK.Label(OtherFrame, text='Signal/Noise =', bg='black', fg='white', font='none 12 bold').grid(row=1, column=2, sticky="ne")

        TK.Label(OtherFrame, text='Gas Sigma (Jy/km/s) =', bg='black', fg='white', font='none 12 bold').grid(row=2, column=0, sticky="ne")
        TK.Label(OtherFrame, text='Int Flux =', bg='black', fg='white', font='none 12 bold').grid(row=2, column=2, sticky="ne")
        
        self.UserDiskThic = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserDiskThic.grid(row=0, column=1, sticky="nwse")
        self.UserDiskAng = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserDiskAng.grid(row=0, column=3, sticky="nwse")
        self.UserInc = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserInc.grid(row=1, column=1, sticky="nwse")
        self.UserRatio = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserRatio.grid(row=1, column=3, sticky="nwse")
        self.UserGasSig = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserGasSig.grid(row=2, column=1, sticky="nwse")
        self.UserFlux = TK.Entry(OtherFrame, width=20, bg='white')
        self.UserFlux.grid(row=2, column=3, sticky="nwse")
                
        #Lower Buttons 
        TK.Button(master, text='Run', width=10, command=lambda : self.Run()).grid(row=11, column=0, columnspan=4, pady=10, sticky = 's')
    
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
        KeyWordWin.geometry("350x400")
        TK.Label(KeyWordWin, text='Keyword Lookup', bg='black', fg='white', font='none 12 bold').grid(row=0, column=0, columnspan=2, sticky="nesw")
        #load Database
        KeyWordDatabase = pd.read_csv('KeyWordLookUp.csv', delimiter=',')
        Keywords, Help = KeyWordDatabase['Keyword'].tolist(), KeyWordDatabase['Help'].tolist()
        #Create dropdown box
        var = TK.StringVar(KeyWordWin)
        option = TK.OptionMenu(KeyWordWin, var, *Keywords)
        option.grid(row=1, columnspan=2, sticky='nesw')
        #Create Text Box
        TextBox = TK.Text(KeyWordWin, height= 10, width=40)
        TextBox.grid(row=2, columnspan=2, sticky='nesw')
        #Add Buttons
        TK.Button(KeyWordWin, text='Look Up', width=10, command=lambda : LookUp(Keywords, Help, option, TextBox)).grid(row=5, column=1, pady=10, sticky = 'w')
        TK.Button(KeyWordWin, text='Close', width=10, command=lambda : Close()).grid(row=5, column=0, pady=10, sticky = 'e')
        
    def SimOptions(self):
        """
        Open Window to edit simulation options
        """
        def Accept(self):
            #Save new values
            self.method = self.Usermethod.get()
            self.MaxIt = int(self.UserMaxIt.get())
            self.Samp = int(self.UserSamp.get())
            #Close window
            SimOptInputWin.destroy()
        
        def Reset(self):
            self.method, self.MaxIt, self.Samp = 'Nelder-Mead', 500, 100000
            Inputs = [self.Usermethod, self.UserMaxIt, self.UserSamp]
            Data = [self.method, self.MaxIt, self.Samp]
            
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
        
        self.Usermethod = TK.Entry(SimOptInputWin, width=20, bg='white')
        self.Usermethod.grid(row=1, column=1, sticky="nwse")
        self.UserMaxIt = TK.Entry(SimOptInputWin, width=20, bg='white')
        self.UserMaxIt.grid(row=2, column=1, sticky="nwse")
        self.UserSamp = TK.Entry(SimOptInputWin, width=20, bg='white')
        self.UserSamp.grid(row=3, column=1, sticky="nwse")
        
        Inputs = [self.Usermethod, self.UserMaxIt, self.UserSamp]
        Data = [self.method, self.MaxIt, self.Samp]
        
        for i in range(len(Inputs)):
            Inputs[i].delete(0, TK.END)
            Inputs[i].insert(0, Data[i])
        
        OptionBottomButtons = TK.Frame(SimOptInputWin, bg='black')
        OptionBottomButtons.grid(row=5, columnspan=2, sticky="n")
        TK.Button(OptionBottomButtons, text='Accept', width=10, command=lambda : Accept(self)).grid(row=11, column=1, pady=10, sticky = 'w')
        TK.Button(OptionBottomButtons, text='Load Default', width=10, command=lambda : Reset(self)).grid(row=11, column=0, pady=10, sticky = 'e')      
        return
        
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
            #Close Window
            FitsInputWin.destroy()        

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
        TK.Button(FitsInputCenter, text='Accept', width=20, command=lambda : UpdateFits(Data, Inputs)).grid(row=18, column=0, columnspan=4, sticky="n")
        
        
    #Load Fits Data
    def Loadfits(self):
        """
        Load .fits files
        """
        filename = askopenfilename()
        self.ShowFile.delete(0, TK.END) 
        self.ShowFile.insert(0, filename.split('/')[-1]) 
        hdulist = fits.open(filename)
        self.scidata = hdulist[0].data.T
        #Reading data for header
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
        self.MaxRad = np.sqrt(self.xsize**2 + self.ysize**2) * 2
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
        self.RMS = np.sqrt(np.mean(self.scidata.copy()[:,:,:self.VelStart]**2))
        
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
        
        
    #Insert Example Equations
    def VPExamples(self, Type):
        """
        Fills in text box with example equations
        """
        if Type == 'Arctan':
            Equ='2/np.pi * [AMPLITUDE] * np.arctan(Radius)'
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
        self.SBLoaded = True
        print('SB Loaded')
        
    #Load VP .csv data
    def VPLoadCSV(self):
        """
        Load Velocity Profile .csv
        """
        VP_CSVName = askopenfilename()
        self.VP_CSV = pd.read_csv(VP_CSVName, delimiter=',')
        self.VP_Rad, self.VP_Prof = self.VP_CSV['Radius'].values, self.VP_CSV['Velocity Profile'].values
        self.VPLoaded = True
        print('VP Loaded')
        
        
    #Save paramiters as .txt
    def SaveInputs(self):
        """
        Saves Paramiters file
        """
        Location = TK.filedialog.askdirectory()
        Location = Location + '/KinMSParameters_%i%i%i.txt' %(time.localtime()[3],time.localtime()[4],time.localtime()[5])
        ParaData = open(Location, "w" )
        ParaData.writelines([' SB Equation¦%s' %self.SBEqu.get(),\
                             '\n VP Equation¦%s' %self.VPEqu.get(),\
                             '\n Disk Thickness¦%s' %self.UserDiskThic.get(), \
                             '\n Disk Angle¦%s' %self.UserDiskAng.get(), \
                             '\n Inclination¦%s' %self.UserInc.get(), \
                             '\n Ratio¦%s' %self.UserRatio.get(),\
                             '\n Gas Sigma¦%s' %self.UserGasSig.get(),\
                             '\n Int Flux¦%s' %self.UserFlux.get()])
        ParaData.close()
    
    
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
        return
    
    
    #Run code
    def Run(self):
        """
        Execute simulation and minimisation
        """
##############################################################################
        t0 = time.time()                                                     #
##############################################################################
        #Define Fitting Function
        def Fitting(Values):#Inc, DiskAng, SBArray, VPArray):
            #Seperate Guess
            X, Y, Inc, DiskAng, GasSig, Flux, DiskThic, SBArray, VPArray = Values[0], Values[1], Values[2], Values[3], Values[4], Values[5], Values[6], Values[7:self.SBlen+7], Values[self.SBlen+7:]
            #Edit User Equation
            Radius = self.Radius
            EditedSB = eval(EditUserString(self.SBEquString, SBArray).replace('|','')) 
            EditedVP = eval(EditUserString(self.VPEquString, VPArray).replace('|','')) 
            #Create Simulation
            SimCube = KinMS(self.ClippedXsize, self.ClippedYsize, self.ClippedVsize, self.cellsize, self.dv, self.beamsize, Inc, gasSigma=GasSig,\
                            sbProf=EditedSB, sbRad=self.Radius, velProf=EditedVP, velRad=self.Radius, diskThick=DiskThic, \
                            nSamps=self.Samp, posAng=DiskAng, intFlux=Flux, phaseCen=np.array([X,Y]))
            self.SimCube = SimCube
            #Determine Chi Square
            Chi = np.sum((self.CubeMask - SimCube)**2 / (self.RMS))
            return Chi
        #Read user inputs
        DiskThic, DiskAng, Inc, Ratio, GasSig, Flux = float(self.UserDiskThic.get()), float(self.UserDiskAng.get()), float(self.UserInc.get()), float(self.UserRatio.get()), float(self.UserGasSig.get()), float(self.UserFlux.get())
        #Make a copy of the data cube
        Cube = self.scidata.copy()[self.XStart:self.XStop,self.YStart:self.YStop,self.VelStart:self.VelStop]
        #Create cube mask
        Xsize, Ysize = len(Cube[:,0,0])*self.cellsize, len(Cube[0,:,0])*self.cellsize
        psf = makebeam(Xsize, Ysize,[self.beamsize[0]/self.cellsize,self.beamsize[1]/self.cellsize],rot=self.beamsize[2])
        self.CubeMask = GUI.smoothclip(Cube, Ratio, self.beamsize, self.cellsize, self.VelStart, self.VelStop, psf)

#        TotMask = np.sum(CubeMask, axis = 2)
#        X = np.arange(self.XStart, self.XStop, 1)
#        Y = np.arange(self.YStart, self.YStop, 1)
#        X, Y = np.meshgrid(X, Y)
#        plt.figure('Mask')
#        plt.contourf(X, Y, TotMask.T, linewidth=0, antialiased=False)

        #Clip Cube for simulation
        self.ClippedXsize = len(self.CubeMask[:,0,0]) * self.cellsize
        self.ClippedYsize = len(self.CubeMask[0,:,0]) * self.cellsize
        self.ClippedVsize = len(self.CubeMask[0,0,:]) * self.dv

        #Get SB and VP data
        self.Radius = np.arange(0, self.MaxRad, self.cellsize/2)
        #Both SB and VP Loaded
        if self.SBLoaded == True and self.VPLoaded == True:
            Guess = np.array([0,0, Inc, DiskAng, GasSig, Flux, DiskThic])
        #VP Loaded
        elif self.SBLoaded == False and self.VPLoaded == True:
            self.SBEquString = self.SBEqu.get()
            SBArray = EditUserString(self.SBEqu.get(), np.zeros(10, dtype=float), ReturnUser = True)
            if SBArray.size == 0:
                Guess = np.array([0,0, Inc, DiskAng, GasSig, Flux, DiskThic])
            else:
                Guess = np.zeros(len(SBArray)+7)
                Guess[2], Guess[3], Guess[4], Guess[5], Guess[6], Guess[7:] = Inc, DiskAng, GasSig, Flux, DiskThic, SBArray 
                self.SBlen = SBArray.size
        #SB Loaded
        elif self.SBLoaded == True and self.VPLoaded == False:
            self.VPEquString = self.VPEqu.get()
            VPArray = EditUserString(self.VPEqu.get(), np.zeros(10, dtype=float), ReturnUser = True)
            if VPArray.size == 0:
                Guess = np.array([0,0, Inc, DiskAng, GasSig, Flux, DiskThic])
            else:
                Guess = np.zeros(len(VPArray)+7)
                Guess[2], Guess[3], Guess[4], Guess[5], Guess[6], Guess[7:] = Inc, DiskAng, GasSig, Flux, DiskThic, VPArray  
                self.VPlen = VPArray.size
        #None Loaded
        elif self.SBLoaded == False and self.VPLoaded == False:
            self.VPEquString = self.VPEqu.get()
            VPArray = EditUserString(self.VPEqu.get(), np.zeros(10, dtype=float), ReturnUser = True)
            self.SBEquString = self.SBEqu.get()
            SBArray = EditUserString(self.SBEqu.get(), np.zeros(10, dtype=float), ReturnUser = True)
            if SBArray.size == 0 or VPArray.size == 0:
                if SBArray.size == 0 and VPArray.size == 0:
                    Guess = np.array([0,0, Inc, DiskAng, GasSig, Flux, DiskThic])
                    self.SBlen, self.VPlen = SBArray.size, VPArray.size
                elif SBArray.size == 0:
                    Guess = np.zeros(VPArray.size + 7)
                    Guess[2], Guess[3], Guess[4], Guess[5], Guess[6], Guess[7:] = Inc, DiskAng, GasSig, Flux, DiskThic, VPArray
                    self.SBlen, self.VPlen = SBArray.size, VPArray.size
                elif VPArray.size == 0:
                    Guess = np.zeros(SBArray.size + 2)
                    Guess[2], Guess[3], Guess[4], Guess[5], Guess[6], Guess[7:] = Inc, DiskAng, GasSig, Flux, DiskThic, SBArray  
                    self.SBlen, self.VPlen = SBArray.size, VPArray.size
            else:
                EqVar = np.zeros(SBArray.size + VPArray.size)
                EqVar[:SBArray.size], EqVar[SBArray.size:] = SBArray, VPArray
                Guess = np.zeros(len(EqVar)+7)
                Guess[2], Guess[3], Guess[4], Guess[5], Guess[6], Guess[7:] = Inc, DiskAng, GasSig, Flux, DiskThic, EqVar
                self.SBlen, self.VPlen = SBArray.size, VPArray.size
        
        #Minimizationwv 
        res = minimize(Fitting, Guess, method = self.method, options={'maxfev': self.MaxIt})
#        print(res)
        #Print Results
        print('Success = ', res.success)
        print('Phase Center =', res.x[0], res.x[1])
        print('Inc = ', res.x[2])
        print('Disk Ang = ', res.x[3])
        print('Gas Sigma = ', res.x[4])
        print('Flux = ', res.x[5])
        print('Disk Thickness = ', res.x[6])
        print('SB', EditUserString(self.SBEquString, res.x[7:self.SBlen+7]))
        print('VP', EditUserString(self.VPEquString, res.x[self.SBlen+7:]))
        
        GUI.makeplots(self.SimCube, self.ClippedXsize, self.ClippedYsize, self.ClippedVsize, self.cellsize,\
                      self.dv, self.beamsize ,posang=DiskAng, overcube=self.CubeMask, pvdthick=DiskThic)        
        
##############################################################################
        t4 = time.time()                                                     #  
        print('Time Taken = %.3f s' %(t4-t0))                                #
        print('Time Taken = %.3f min' %((t4-t0)/60))                         #
##############################################################################
        



























class CustomWin:
    #Initialisation
    def __init__(self, master):
        self.master = master
    
        #Surface Brightness 
        TK.Label(master, text='Surface Brigtness Equation', bg='black', fg='white', font='none 12 bold').grid(row=3, column=0, columnspan=4, pady=10, sticky="n")
        SBFrame = TK.Frame(master, bg='black', width=600, height=300)
        SBFrame.grid(row=4, columnspan=4)
        self.SBEqu = TK.Entry(SBFrame, width=70, bg='white')
        self.SBEqu.grid(row=1, column=0, columnspan=3, sticky='nse')
        self.SBLoaded = False
        TK.Button(SBFrame, text='Load .csv', width=10, command=lambda : self.SBLoadCSV()).grid(row=1, column=4, sticky="w")
        SBSubFrame = TK.Frame(master, bg='black', width=600, height=300)
        SBSubFrame.grid(row=5, columnspan=4, sticky='n')
        TK.Button(SBSubFrame, text='Gaussian', width=10, command=lambda : self.SBExamples(Type='Gaussian')).grid(row=2, column=0, sticky="e")
        TK.Button(SBSubFrame, text='Exponential', width=10, command=lambda : self.SBExamples(Type='Exponential')).grid(row=2, column=1, sticky="w")
    
        #Velocity profile
        TK.Label(master, text='Velocity Profile Equation', bg='black', fg='white', font='none 12 bold').grid(row=6, column=0, columnspan=4, pady=10, sticky="n")
        VPFrame = TK.Frame(master, bg='black', width=600, height=300)
        VPFrame.grid(row=7, columnspan=4)
        self.VPEqu = TK.Entry(VPFrame, width=70, bg='white')
        self.VPEqu.grid(row=1, column=0, columnspan=3, sticky='nse')
        self.VPLoaded = False
        TK.Button(VPFrame, text='Load .csv', width=10, command=lambda : self.VPLoadCSV()).grid(row=1, column=4, sticky="w")
        VPSubFrame = TK.Frame(VPFrame, bg='black', width=600, height=300)
        VPSubFrame.grid(row=5, columnspan=4, sticky='n')
        TK.Button(VPSubFrame, text='Arctan', width=10, command=lambda : self.VPExamples(Type='Arctan')).grid(row=0, column=0, sticky="n")
    
    #Insert Example Equations
    def VPExamples(self, Type):
        """
        Fills in text box with example equations
        """
        if Type == 'Arctan':
            Equ='2/np.pi * [AMPLITUDE] * np.arctan(Radius)'
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
            
            
            
root = TK.Toplevel()  #Change to TK.Tk() for running ouside of spyder
my_gui = BootWindow(root)
root.mainloop()


