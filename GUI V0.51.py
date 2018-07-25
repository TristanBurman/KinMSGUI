import tkinter as TK 
from tkinter.filedialog import askopenfilename
from astropy.io import fits
import numpy as np
import pandas as pd
import time
from KinMS import KinMS


#from IPython import get_ipython
#ipython = get_ipython()
#ipython.magic("gui tk")

#Import GUI functions 
import GUIFunctions as GUI

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


#Initial Boot window
class BootWindow:
    #Initialization 
    def __init__(self, master):
        self.master = master
        master.title('KinMS GUI')
        master.configure(background='black')
##############################################################################
        #Load Photo
        self.im = TK.PhotoImage(file='Header.gif')
        self.Image = TK.Label(master, image=self.im, bg='black').grid(row=0, column=0, columnspan=2)   
##############################################################################
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
        self.master = master
        master.title('KinMS GUI Read fits')
        master.configure(background='black')
##############################################################################
        #Load Photo
        self.im = TK.PhotoImage(file='Header.gif')
        self.Image = TK.Label(master, image=self.im, bg='black').grid(row=0, column=0, columnspan=4)   
##############################################################################
        #Fits Data
        TK.Label(master, text='File selected:', bg='black', fg='white', font='none 12 bold').grid(row=1, column=0, columnspan=2, pady=10, sticky="ne")
        TK.Button(master, text='View .fits Inputs', width=20, command=lambda : self.FitsInput()).grid(row=2, column=0, columnspan=2, sticky="ne")
        TK.Button(master, text='Select .fits file', width=20, command=lambda : self.Loadfits()).grid(row=2, column=2, columnspan=2, sticky="nw")
        
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
        self.UserDiskThic = TK.Entry(OtherFrame, width=25, bg='white')
        self.UserDiskThic.grid(row=0, column=1, sticky="nwse")
        self.UserDiskAng = TK.Entry(OtherFrame, width=25, bg='white')
        self.UserDiskAng.grid(row=0, column=3, sticky="nwse")
        self.UserInc = TK.Entry(OtherFrame, width=25, bg='white')
        self.UserInc.grid(row=1, column=1, sticky="nwse")
        self.UserRatio = TK.Entry(OtherFrame, width=25, bg='white')
        self.UserRatio.grid(row=1, column=3, sticky="nwse")
                
        #Lower Buttons 
        TK.Button(master, text='Load', width=10, command=lambda : self.LoadInputs()).grid(row=10, column=0, columnspan=2, sticky='e', pady=10)
        TK.Button(master, text='Save', width=10, command=lambda : self.SaveInputs()).grid(row=10, column=2, columnspan=2, sticky='w', pady=10)       
        TK.Button(master, text='Run', width=10, command=lambda : self.Run()).grid(row=11, column=0, columnspan=4, pady=10)
        

    #Edit Fits Values
    def FitsInput(self):
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
        except NameError:
            Data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
        filename = askopenfilename()
        TK.Label(self.master, text=filename.split('/')[-1], bg='black', fg='white', font='none 12 bold').grid(row=1, column=2, columnspan=2, sticky="w") 
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
        
        #Determin the relevant spacial area
        compressed = np.sum(self.scidata, axis=2)
        x, y = np.sum(compressed, axis=0), np.sum(compressed, axis=1)
        Allowence = np.array([len(x)*0.05, len(y)*0.05])
        x, y = x[int(len(x)*0.1):int(len(x)*0.9)], y[int(len(y)*0.1):int(len(y)*0.9)]  #Excluding outter 10%
        xCutoff, yCutoff = 1.5 * np.mean(x), 1.5 * np.mean(y)
        xAOI, yAOI = np.where(x > xCutoff), np.where(y > yCutoff)
        print(Allowence)
        xBoundry, yBoundry = np.array([np.min(xAOI)-Allowence[0], np.max(xAOI)+Allowence[0]]), np.array([np.min(yAOI)-Allowence[1], np.max(yAOI)+Allowence[1]])
        if xBoundry[0] < 0: xBoundry[0] = 0
        if xBoundry[1] > len(x): xBoundry[1] = int(len(x)-1)
        if yBoundry[0] < 0: yBoundry[0] = 0
        if yBoundry[1] > len(y): yBoundry[1] = int(len(y)-1)
        self.XStart, self.XStop, self.YStart, self.YStop = int(xBoundry[0]+(len(x)*0.1)), int(xBoundry[1]+(len(x)*0.1)), int(yBoundry[0]+(len(y)*0.1)), int(yBoundry[1]+(len(y)*0.1))
        
    #Load SB .csv data
    def SBLoadCSV(self):
        SB_CSVName = askopenfilename()
        self.SB_CSV = pd.read_csv(SB_CSVName, delimiter=',')
        self.SB_Rad, self.SB_Prof = self.SB_CSV['Radius'].values, self.SB_CSV['Surface Brightness'].values
        self.SBLoaded = True
        print('SB Loaded')
    
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
    
    #Load VP .csv data
    def VPLoadCSV(self):
        VP_CSVName = askopenfilename()
        self.VP_CSV = pd.read_csv(VP_CSVName, delimiter=',')
        self.VP_Rad, self.VP_Prof = self.VP_CSV['Radius'].values, self.VP_CSV['Velocity Profile'].values
        self.VPLoaded = True
        print('VP Loaded')
        
    #Insert Example Equations
    def VPExamples(self, Type):
        """
        Fills in text box with example equations
        """
        if Type == 'Arctan':
            Equ='2/np.pi * [AMPLITUDE] * np.arctan(Radius)'
            self.VPEqu.delete(0,TK.END)
            self.VPEqu.insert(0,Equ)

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
                             '\n Ratio¦%s' %self.UserRatio.get()])
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
                  self.UserRatio]
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
##############################################################################
        t0 = time.time()                                                     #
##############################################################################
        #Read user inputs
        DiskThic, DiskAng, Inc, Ratio = float(self.UserDiskThic.get()), float(self.UserDiskAng.get()), float(self.UserInc.get()), float(self.UserRatio.get())
        #Make a copy of the data cube
        Cube = self.scidata.copy()[self.XStart:self.XStop,self.YStart:self.YStop,self.VelStart:self.VelStop]
##############################################################################    
        t1 = time.time()                                                     #
        print('Get user input = ', (t1-t0))                                  #
##############################################################################
        #Create cube mask
        CubeMask = GUI.smoothclip(Cube, Ratio, self.beamsize, self.cellsize, self.VelStart, self.VelStop)
##############################################################################
        t2 = time.time()                                                     #
        print('Create Mask = ', (t2-t1))                                     #
##############################################################################
        #Get SB and VP data
        if self.SBLoaded == True:
            SBRadius = self.SB_Rad
            SB = self.SB_Prof
        else:
            Radius = np.arange(0, self.MaxRad, self.cellsize/2)
            SBRadius, SB = Radius, eval(self.SBEqu.get())
        if self.VPLoaded == True:
            VPRadius = self.VP_Rad
            VP = self.VP_Prof
        else:
            Radius = np.arange(0, self.MaxRad, self.cellsize/2)
            VPRadius, VP = Radius, eval(self.SBEqu.get())
##############################################################################
        t3 = time.time()                                                     #
        print('Load SB and VP = ', t3-t2)                                    #
##############################################################################
        #Create the Sim Cube
        ClippedXsize = len(Cube[:,0,0])*self.cellsize
        ClippedYsize = len(Cube[0,:,0])*self.cellsize
        ClippedVsize = len(Cube[0,0,:])*self.dv
            
        SimCube = KinMS(ClippedXsize, ClippedYsize, ClippedVsize, self.cellsize, self.dv, self.beamsize, Inc,\
                        sbProf=SB, sbRad=SBRadius, \
                        velProf=VP, velRad=VPRadius,\
                        posAng=DiskAng)
##############################################################################
        t4 = time.time()                                                     #
        print('Create Sim = ', t4-t3)                                        #
        print('Total =', t4-t0)                                              #
##############################################################################
        
#        GUI.makeplots(SimCube, ClippedXsize, ClippedYsize, ClippedVsize, self.cellsize,\
#                      self.dv, self.beamsize ,posang=DiskAng, overcube=CubeMask, pvdthick=DiskThic)













class CustomWin:
    #Initialisation
    def __init__(self, master):
        self.master = master
        
        


    
    
    













root = TK.Toplevel()  #Change to TK.Tk() for running ouside of spyder
my_gui = BootWindow(root)
root.mainloop()


