import sys
import os
import configparser
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import messagebox
from PIL import ImageTk, Image
from customWidgets import *
from simCalculationsV2 import *


#GLOBAL VARIABLES:
AppVersionTitle = "SimCalc Beta"
defaultWindowSize = "1450x900"
minWindowSizex = 1450
minWindowSizeY = 900
supportedSimTypes = ["2DOF (Symmetric)",]
cfgPath = r'cfg.py'
faviconPath = r'data\favicon.ico'
dataPath = r'data'

#general Default Params:
unitsDefault = 'metric'

#CTC Default params:
ctcLengthDefault = 5
motorTorqueDefault = 40
motorSpeedDefault = 110
maxAngleDefault = 45
currentAngleDefault = 0

#Sim Geometry Defaults
rodMtxDefault = 8.2
rodMtyDefault = 21.99
rodMtzDefault = 27.68

motorPosxDefault = 22.4
motorPosyDefault = -45.26
motorPoszDefault = -8.0

motorAngleDefault = 0
ctcRestAngleDefault = 45
rodLenDefault = 47.44

#MAGIC CONSTANTS
LARGE_FONT=12
EMPTY = ''
IN_2_MM = 25.4
IN_2_M = 0.0254
MM_2_IN = 0.0393701
M_2_IN = 39.3701
N_2_LBS = 0.224808942443
LBS_2_N = 4.448221628254617
LBF_2_NM = 1.35581794833
NM_2_LBF = 0.7375621492780272
N_2_KG = 0.101971621298
KG_2_N = 9.80665
MM_2_M = 1/1000
M_2_MM = 1000


class Application(tk.Frame):

#CALLBACKS AND FUNCTIONS
########################################################################
    def setUnitsMetric(self):
        #change unit constants
        self.unitLengthShort.set('(mm)')
        self.unitLengthLong.set('(m)')
        self.unitForce.set('(N)')
        self.unitWeight.set('(kg)')
        self.unitTorque.set('(N*m)')
        self.unitSpeed.set('(mm/s)')
        self.units='metric'

        
        print('units changed to metric')
        
    def setUnitsImperial(self):
        #change unit constants
        self.unitLengthShort.set('(in)')
        self.unitLengthLong.set('(in)')
        self.unitForce.set('(lbs)')
        self.unitWeight.set('(lbs)')
        self.unitTorque.set('(lbs*ft)')
        self.unitSpeed.set('(in/s)')
        self.units='imperial'
        print('units changed to imperial')

        #change variable values
    

#ALL CALCULATIONS DONE IN METRIC, CONVERT TO MM IF NEEDED
    def ctcPageCalculate(self):
        print('calculating for CTC page...')

        

        #get necessary user input values from the fields in the GUI
        try:
            ctcLen = float(self.ctcLengthVar.get())
            print('CTC length = ',ctcLen)
            mTorque = float(self.motorTorqueVar.get())
            print('Motor torque = ',mTorque)
            mSpeed = float(self.motorSpeedVar.get())
            print('Motor speed = ',mSpeed)
            maxAngle = float(self.maxAngleVar.get())
            print('Max CTC angle = ',maxAngle)
            currentAngle = float(self.currentAngleVar.get())
            print('Current CTC angle = ',currentAngle)
        except ValueError:
            print('Value Error...');
            messagebox.showerror('Entry Error', 'Correctly fill out all entries before calculating')
            return

#CONVERT UNITS TO STANDARD METRIC FOR CALCULATION
        if (self.units == 'imperial'): #convert to metric
            ctcLen = ctcLen*IN_2_M
            mTorque = mTorque*LBF_2_NM
        if (self.units == 'metric'):
            ctcLen = ctcLen*MM_2_M
            
        #getCTCparams will return the parameters indicated below       
        ctcParams = getCTCParams(ctcLen, maxAngle, mSpeed, mTorque, currentAngle, unit=None)

        travel = ctcParams[0]
        maxVelocity = ctcParams[1]
        maxForce = ctcParams[2]
        currentForce = ctcParams[3]
        currentVelocity = ctcParams[4]

#CONVERT UNITS BACK TO DESIRED VALUES FOR OUTPUT
        if (self.units == 'imperial'): #convert back to imperial
            travel = travel*M_2_IN
            maxVelocity = maxVelocity*M_2_IN
            maxForce = maxForce*N_2_LBS
            currentForce = currentForce*N_2_LBS
            currentVelocity = currentVelocity*M_2_IN
        if (self.units == 'metric'):
            travel = travel*M_2_MM
            maxVelocity = maxVelocity*M_2_MM
            currentVelocity = currentVelocity*M_2_MM

        #change the displayed values in the output containers (the float formatting limits it to 4 decimals)
        self.outputSpeedLinearVar.set('{:f}'.format(currentVelocity))
        self.outputForceLinearVar.set('{:f}'.format(currentForce))
        self.outputLinearTravelVar.set('{:f}'.format(travel*2)) #travel*2 since the function only returns travel in one direction

        return 

#FUNCTION TO DO CALCULATIONS FOR ALL VARIABLES ON THE GEOMETRY PAGE
    def geometryPageCalculate(self):
        print('calculating for geometry page...')

        #Get relevant information from the input fields:
        try:
            ctcLen = float(self.ctcLengthVar.get())
            print('CTC length = ',ctcLen)
            mTorque = float(self.motorTorqueVar.get())
            print('Motor torque = ',mTorque)
            mSpeed = float(self.motorSpeedVar.get())
            print('Motor speed = ',mSpeed)
            maxAngle = float(self.maxAngleVar.get())
            print('Max CTC angle = ',maxAngle)
            currentAngle = float(self.currentAngleVar.get())
            print('Current CTC angle = ',currentAngle)
            restAngle = float(self.ctcRestAngleVar.get())
            print('CTC resting angle is ',restAngle)

            #Geometry Params
            mtrAngle =  float(self.motorAngleVar.get())
            print('Motor angle = ',mtrAngle,' inwards')
            
            rodMtx = float(self.rodMtx.get())
            rodMty = float(self.rodMty.get())
            rodMtz = float(self.rodMtz.get())
            print('Rod mount coords: x=',rodMtx,', y=',rodMty,', z=',rodMtz)
            
            
            mtrCoordsx = float(self.motorPosx.get())
            mtrCoordsy = float(self.motorPosy.get())
            mtrCoordsz = float(self.motorPosz.get())
            print('Motor center coords: x=',mtrCoordsx,', y=',mtrCoordsy,', z=',mtrCoordsz)
        except ValueError:
            print('Value Error...');
            messagebox.showerror('Entry Error', 'Correctly fill out all entries before calculating')
            return

        #temporary bandaid while I rewrite the math code:
        if (mtrCoordsy > 0):
            restAngle = 180-restAngle
            mtrAngle = -mtrAngle

        
        #DO UNIT CONVERSION TO STANDARD METRIC BEFORE INPUT TO FUNCTION
        if (self.units == 'imperial'): #convert back to imperial
            ctcLen = ctcLen*IN_2_M
            mTorque = mTorque*LBF_2_NM

            rodMtx *= IN_2_M
            rodMty *= IN_2_M
            rodMtz *= IN_2_M

            mtrCoordsx *= IN_2_M
            mtrCoordsy *= IN_2_M
            mtrCoordsz *= IN_2_M
            
        if (self.units == 'metric'):
            ctcLen *= MM_2_M
            
        #Create vector objects from the rod mount and motor coordinate data, these will be used in the functions below
        rodMt = Vector(rodMtx, rodMty, rodMtz)
        mtrPosition = Vector(mtrCoordsx, mtrCoordsy, mtrCoordsz) 
        
        #get the actuator geometry:
        ctcParams = getCTCParams(ctcLen, maxAngle, mSpeed, mTorque, currentAngle, unit=None)


        #TO BE USED WHEN DYNAMIC SIMULATION IS INTRODUCED#################################################################
        ctctravel = ctcParams[0]
        ctcMaxVelocity = ctcParams[1]
        ctcMaxForce = ctcParams[2] 
        ctcCurrentForce = ctcParams[3]
        ctcCurrentVelocity = ctcParams[4]
        #TO BE USED WHEN DYNAMIC SIMULATION IS INTRODUCED#################################################################

        #DO UNIT CONVERSION BACK TO DESIRED BEFORE OUTPUTS
        if (self.units == 'imperial'): #convert back to imperial
            ctctravel *= M_2_IN
            ctcMaxVelocity *= M_2_IN
            ctcMaxForce *= N_2_LBS
            ctcCurrentForce *= N_2_LBS
            ctcCurrentVelocity *= M_2_IN
            
        if (self.units == 'metric'):
            ctctravel *= M_2_MM
            ctcMaxVelocity *= M_2_MM
            ctcCurrentVelocity *= M_2_MM
        
        #get the CTC position at rest position
        #getCTCGeometry3d(length, restAngle, ctcAngle, mtrPos, mtrAngle=0)
        ctcRestGeometry = getCTCGeometry3d(ctcLen, restAngle, 0, mtrPosition, mtrAngle) #the 0 is the current angle of the lever, we want it to just be the rest angle
        ctcRestPosition = ctcRestGeometry[0] #position vector of the lever end in space
        restForceVect = ctcRestGeometry[1]
        print('CTC rest position is:')
        #ctcRestPosition.print()

        #get the motor rod length from the resting CTC position and the rot mounting point
        rodLen = vSubtract(ctcRestPosition, rodMt)
        rodLen = rodLen.length
        print('\n\n\nRod length is',rodLen,'\n\n\n')

        #get the CTC position at maximum position
        ctcMaxGeometry = getCTCGeometry3d(ctcLen, restAngle, maxAngle, mtrPosition, mtrAngle)
        ctcMaxPosition = ctcMaxGeometry[0]
        maxForceVect = ctcMaxGeometry[1]
        print('CTC max position is:')
        ctcMaxPosition.print()
        
        #get the CTC position at minimum position
        ctcMinGeometry = getCTCGeometry3d(ctcLen, restAngle, -maxAngle, mtrPosition, mtrAngle)
        ctcMinPosition = ctcMinGeometry[0]
        minForceVect = ctcMinGeometry[1]
        print('CTC min position is:')
        ctcMinPosition.print()

        #get travel parameters for pitch and roll
        #feed in minimum and maximum travel positions as calculated before:
        #simTravel(ctcPosi, ctcPosf, rodMt, rodLen, pivLenUser=None)
        try:
            travParams = simTravel(ctcMinPosition, ctcMaxPosition, rodMt, rodLen)
        except ValueError:
            print('something went wrong in simTravel...')
            messagebox.showerror('Geometry Error', 'You may have entered invalid geometry, please double check your values')
            return
        
        #pitch linear travel and pitch angular travel
        pitchParams = travParams[0]
        pitchLinTravel = pitchParams[0]
        pitchAngle = pitchParams[1]
        print('Pitch Travel is',pitchLinTravel,' or ',pitchAngle,'deg')

        #roll linear travel and roll angular travel
        rollParams = travParams[1]
        rollLinTravel = rollParams[0]
        rollAngle = rollParams[1]
        print('Roll Travel is',rollLinTravel,' or ',rollAngle,'deg')

        #do force calculations now:
        #forceCalcs(ctcPos, rodMt, force, forceDir):
        #return(pitchEff, pitchForce, rollEff, rollForce)
        forceParams = forceCalcs(ctcRestPosition, rodMt, ctcCurrentForce, restForceVect)
        pitchEff = 100*forceParams[0]
        pitchForce = 2*forceParams[1] #return twice the force value since we're counting for 2 motors
        pitchTorque = 2*forceParams[4]
        print('Pitch force of ',pitchForce,' with efficiency of ',pitchEff,'%','eq torque is ',pitchTorque)
        
        rollEff = 100*forceParams[2]
        rollForce = 2*forceParams[3] #return twice the force value since we're counting for 2 motors
        rollTorque = 2*forceParams[5]
        print('Roll force of ',pitchForce,' with efficiency of ',rollEff,'%','eq torque is ',rollTorque)

        #note: Roll efficiency and pitch efficiency denote how much force is being wasted simply applying tension to the frame under loading
        
        #DO UNIT CONVERSION BACK TO DESIRED BEFORE OUTPUTS
        if (self.units == 'imperial'): #convert back to imperial
            pitchForce *= NM_2_LBF
            rollForce *= NM_2_LBF
            pitchLinTravel *= M_2_IN
            rollLinTravel *= M_2_IN
            pitchTorque *= NM_2_LBF
            rollTorque *= NM_2_LBF
            
        if (self.units == 'metric'):
            pitchLinTravel *= M_2_MM
            rollLinTravel *= M_2_MM

        #finally set the output values to the calculated states
        self.pitchEff.set('{:f}'.format(pitchEff))
        self.rollEff.set('{:f}'.format(rollEff))
        self.pitchForce.set('{:f}'.format(pitchForce))
        self.rollForce.set('{:f}'.format(rollForce))
        self.pitchRot.set('{:f}'.format(pitchAngle))
        self.rollRot.set('{:f}'.format(rollAngle))
        self.pitchTrav.set('{:f}'.format(pitchLinTravel))
        self.rollTrav.set('{:f}'.format(rollLinTravel))
        self.pitchTorque.set('{:f}'.format(pitchTorque))
        self.rollTorque.set('{:f}'.format(rollTorque))

        return

    def geometryPageTooltip(self):
        messagebox.showinfo(
            "Help",
            """Some explanation about what these output calculations mean:

Pitch/Roll Efficiency:
    In most designs, the axis of applied force is not perfectly
    perpendicular to the line of pivot. Thus some of the force
    applied to the frame will be wasted in creating tension or
    compression rather than creating movement. For standard
    2DOF designs, pitch efficiency should be near 90% and roll
    efficiency is usually much lower (ranges between 20% and
    60%). Make sure that your CTC rest angle is input
    correctly, otherwise these values will read much lower
    than their maximums. As a general rule, if you want to
    improve the pitch efficiency of a rig, you can usually get
    near 100% by changing the rest angle of the CTC lever to
    the proper position. To get higher roll efficiency, try
    tilting the motors inwards (motor angle), moving the
    actuator mounts further apart, and moving the rod mounts
    further apart.
    

Pitch/Roll Force:
    The applied force in either only the pitch or only the
    roll direction, at the point of attachment between
    the motor arm and the rig.
    

Pitch/Roll Travel:
    The total amount of linear travel in either the pitch
    or roll direction of your rig, measured at the point
    of attachment between the motor arm and the rig.

Pitch/Roll Rotation:
    The total amount of rotation in either the pitch or
    the roll direction

Pitch/Roll Torque:
    The equivalent torque at the universal joint in the
    pitch and roll directions.

            """)




            
        

#saves default values to a .ini file
    def saveDefaults(self):       
        print('saving defaults')

        config = configparser.ConfigParser()
        
        config['GENERAL PARAMETERS'] = {'unitsDefault':self.units}

        config['CTC PARAMETERS'] = {
            'ctcLengthDefault':self.ctcLengthVar.get(),
            'motorTorqueDefault':self.motorTorqueVar.get(),
            'motorSpeedDefault':self.motorSpeedVar.get(),
            'currentAngleDefault':self.currentAngleVar.get(),
            'maxAngleDefault':self.maxAngleVar.get()}

        config['SIM GEOMETRY PARAMETERS'] = {
            'rodMtxDefault':self.rodMtx.get(),
            'rodMtyDefault':self.rodMty.get(),
            'rodMtzDefault':self.rodMtz.get(),
            'motorPosxDefault':self.motorPosx.get(),
            'motorPosyDefault':self.motorPosy.get(),
            'motorPoszDefault':self.motorPosz.get(),
            'motorAngleDefault':self.motorAngleVar.get(),
            'ctcRestAngleDefault':self.ctcRestAngleVar.get(),
            'rodLenDefault':self.rodLen.get()}

        with open('cfg.ini', 'w') as file:
            config.write(file)            
        
        return

#reads back the .ini file and sets the values of the corresponding variables
    def getDefaults(self):
        config = configparser.ConfigParser()

        if not os.path.isfile('cfg.ini'):
            return
        
        config.read('cfg.ini')

        print('setting defaults...')

        self.units = config['GENERAL PARAMETERS']['unitsDefault']
        print(self.units)

        self.ctcLengthVar.set(config['CTC PARAMETERS']['ctclengthdefault'])
        self.motorTorqueVar.set(config['CTC PARAMETERS']['motortorquedefault'])
        self.currentAngleVar.set(config['CTC PARAMETERS']['currentangledefault'])
        self.motorSpeedVar.set(config['CTC PARAMETERS']['motorspeeddefault'])
        self.currentAngleVar.set(config['CTC PARAMETERS']['currentangledefault'])
        self.maxAngleVar.set(config['CTC PARAMETERS']['maxangledefault'])

        self.rodMtx.set(config['SIM GEOMETRY PARAMETERS']['rodmtxdefault'])
        self.rodMty.set(config['SIM GEOMETRY PARAMETERS']['rodmtydefault'])
        self.rodMtz.set(config['SIM GEOMETRY PARAMETERS']['rodmtzdefault'])

        self.motorPosx.set(config['SIM GEOMETRY PARAMETERS']['motorposxdefault'])
        self.motorPosy.set(config['SIM GEOMETRY PARAMETERS']['motorposydefault'])
        self.motorPosz.set(config['SIM GEOMETRY PARAMETERS']['motorposzdefault'])

        self.motorAngleVar.set(config['SIM GEOMETRY PARAMETERS']['motorangledefault'])
        self.ctcRestAngleVar.set(config['SIM GEOMETRY PARAMETERS']['ctcrestangledefault'])
        self.rodLen.set(config['SIM GEOMETRY PARAMETERS']['rodlendefault'])
        
        return
        

#GUI Stuff: ROOT FRAME SETUP
########################################################################
    def __init__(self, master=None, *args, **kwargs):

        tk.Frame.__init__(self, master=None, *args, **kwargs)
        self.master = master
   
        #initializes the main window
        self.master.title(AppVersionTitle)
        self.grid(row=0, column=0)

        #root frame grid configuriation, let it stretch
        self.master.columnconfigure(0, weight=1)
        self.master.columnconfigure(1, weight=1) 
        self.master.rowconfigure(0, weight=1)
        self.master.rowconfigure(1, weight=1)

        topCanvas = tk.Canvas(self)
        topCanvas.configure(scrollregion=topCanvas.bbox('all'))

        topScrollx = tk.Scrollbar(self.master, orient=tk.HORIZONTAL, command=topCanvas.xview)
        topScrollx.grid(row=5, column=0, columnspan=5)
        
        topScrolly = tk.Scrollbar(self.master, orient=tk.VERTICAL, command=topCanvas.yview)
        topScrolly.grid(row=0, column=3)

        topCanvas.configure(yscrollcommand = topScrolly.set)
        topCanvas.configure(xscrollcommand = topScrollx.set)

        topCanvas.grid(row=0, column=0, sticky='nsew')
        

        
        
#GLOBAL PARAMETERS, CALCULATION INPUTS
####################################################################

        #unit stringvars for changing calculator units
        self.unitLengthShort = tk.StringVar()
        self.unitLengthShort.set('(mm)')
        self.unitLengthLong = tk.StringVar()
        self.unitLengthLong.set('(m)')
        self.unitForce = tk.StringVar()
        self.unitForce.set('N')
        self.unitWeight = tk.StringVar()
        self.unitWeight.set('kg')
        self.unitTorque = tk.StringVar()
        self.unitTorque.set('N*m')
        self.unitSpeed = tk.StringVar()
        self.unitSpeed.set('(mm/s)')

        self.units = 'metric'

        #CTC Page variables
        #inputs
        self.ctcLengthVar = tk.StringVar()
        self.motorTorqueVar = tk.StringVar()
        self.motorSpeedVar = tk.StringVar()
        self.maxAngleVar = tk.StringVar()
        self.currentAngleVar = tk.StringVar()
        
        #outputs
        self.outputSpeedLinearVar = tk.StringVar()
        self.outputForceLinearVar = tk.StringVar()
        self.outputLinearTravelVar = tk.StringVar()

        #GEOMETRY Page variables
        #inputs
        self.rodMtx = tk.StringVar()
        self.rodMty = tk.StringVar()
        self.rodMtz = tk.StringVar()
        
        self.motorPosx = tk.StringVar()
        self.motorPosy = tk.StringVar()
        self.motorPosz = tk.StringVar()
        
        self.motorAngleVar = tk.StringVar()
        self.ctcRestAngleVar = tk.StringVar()
        self.rodLen = tk.StringVar()

        

        #outputs
        self.pitchEff = tk.StringVar()
        self.rollEff = tk.StringVar()
        self.pitchForce = tk.StringVar()
        self.rollForce = tk.StringVar()
        self.pitchRot = tk.StringVar()
        self.rollRot = tk.StringVar()
        self.pitchTrav = tk.StringVar()
        self.rollTrav = tk.StringVar()
        self.pitchTorque = tk.StringVar()
        self.rollTorque = tk.StringVar()

        #SET DEFAULTS

        self.getDefaults()
        
        if self.units == 'metric':
            print('default units are is metric')
            self.setUnitsMetric() #change the units          
            
        if self.units == 'imperial':
            print('default units are imperial')
            self.setUnitsImperial()

        if (self.units!= 'metric') and (self.units!='imperial'):
            self.units = 'metric'

        
            
        
        #rest of GUI goes here:
####################################################################
        #MENU BAR SECTION
        menuBar = tk.Menu(self.master) #init the top level menu bar

        #file menu and cascade
        file = tk.Menu(menuBar, tearoff=False)
        file.add_command(label = "Save Values to Default", command = self.saveDefaults)        
        menuBar.add_cascade(label = "File", menu = file)

        #Edit menu and cascade
        #edit = tk.Menu(menuBar, tearoff=False)
        #edit.add_command(label = "test")
        #menuBar.add_cascade(label = "Edit", menu = edit)

        #options menu and cascade
        #options = tk.Menu(menuBar, tearoff=False)
        #options.add_command(label = "test")
        #menuBar.add_cascade(label = "Options", menu = options)
        

        #show the menubar
        self.master.config(menu = menuBar)
#####################################################################



#####################################################################
        #TOPLEVEL OPTIONS FRAME, RESIDES IN ROOT(0,0)
        opsFrame = ttk.Frame(topCanvas, height=100, relief='groove')
        opsFrame.grid(row=0, column=0, sticky = "nsew", ipadx=3)
        

        #unit label, resides in opsFrame(0,0)
        unitsLabel = ttk.Label(opsFrame, text="Units:")
        unitsLabel.grid(row=0, column=0, padx=2, pady=2)
        
        #unit selection
        self.unitsVar = tk.StringVar()
        
        #frame to hold radiobuttons, resides in opsFrame(0,1)
        ubFrame = ttk.Frame(opsFrame)
        ubFrame.grid(row=0, column=1, pady=3, sticky="ns")
        
        unitsMetric = ttk.Radiobutton(ubFrame, text=" Metric", variable=self.unitsVar, value="metric", command=self.setUnitsMetric)
        unitsImperial = ttk.Radiobutton(ubFrame, text=" Imperial", variable=self.unitsVar, value="imperial", command=self.setUnitsImperial)        
        unitsMetric.grid(row=0, column=1, sticky="w", padx=2, pady=2)
        unitsImperial.grid(row=1, column=1, sticky="w", padx=2, pady=2)

        #ensure that one of the radiobuttons is active on startup
        if self.units == 'metric':
            unitsMetric.invoke()                   
        if self.units == 'imperial':
            unitsImperial.invoke()
        

        #spacer at opsFrame(0,2)
        spacer1 = ttk.Frame(opsFrame, width=20)
        spacer1.grid(row=0, column=2)

        #simulator type label, resides in opFrame(0,3)
        simLabel = ttk.Label(opsFrame, text="Select Simulator Type:")
        simLabel.grid(row=0, column=3, padx=2, pady=2)

        
        #simulator type combobox, resides in opFrame(0,4)
        self.simVar = tk.StringVar()
        simType = ttk.Combobox(opsFrame, textvariable=self.simVar, state="readonly")
        simType['values'] = supportedSimTypes
        simType.current(0)
        simType.grid(row=0, column=4, padx=2, pady=2)

        #Spacer at opsFrame(0,5)
        spacer2 = ttk.Frame(opsFrame, width=20)
        spacer2.grid(row=0, column=5)
        
        #actuator type selection label, resided in opsFrame(0,6)
        acLabel = ttk.Label(opsFrame, text="Select Actuator Type:")
        acLabel.grid(row=0, column=6, padx=2)

        #actuator radiobutton frame
        acFrame = ttk.Frame(opsFrame)
        acFrame.grid(row=0, column=7, pady=3, sticky="ns")

        #Actuator Type
        self.actuatorVar = tk.StringVar()
        
        actuatorCTC = ttk.Radiobutton(acFrame, text=" Motor with CTC Lever", variable=self.actuatorVar, value="ctc")
        actuatorLinear = ttk.Radiobutton(acFrame, text=" Linear Actuator", variable=self.actuatorVar, value="linear", state='disabled')
        actuatorCTC.grid(row=0, column=0, sticky='w', padx=2, pady=2)
        actuatorLinear.grid(row=1, column=0, sticky='w', padx=2, pady=2)
        actuatorCTC.invoke()
    
#########################################################################


##CALCULATIONS FRAME, HOLDS NOTEBOOK THAT CALCULATOR INTERFACE IS STORED ON, RESIDES IN ROOT(1,0)
#########################################################################
        calcFrame = ttk.Frame(topCanvas, relief="groove")
        calcFrame.grid(row=1, column=0, padx=5, pady=5, sticky="NSEW")

        calcTabs = ttk.Notebook(calcFrame, padding=10,)
        calcTabs.grid(row=0, column=0, padx=10, pady=10)

#CTC CALCULATION FRAME RESIDES IN CALCFRAME
#########################################################################
        ctcCalcFrameContainer = ttk.Frame(calcFrame)
        ctcCalcFrame = CalculationFrame(ctcCalcFrameContainer)

        #create all objects necessary to put an image in the frame
        ctcImgVar = ImageTk.PhotoImage(Image.open('diagrams\ctcDiag.png'))
        
        ctcImg = ttk.Label(ctcCalcFrame.imgFrame, image=ctcImgVar)
        ctcImg.image = ctcImgVar
        ctcImg.grid(row=0, column=0, padx=2, pady=2)

        #caption to go underneath the image
        imgCapt1="""
Use this tab to calculate your CTC lever parameters.
Some important things to note:

1)  When your lever is at the zero position, always make
     sure that the rod makes a 90 degree angle with the CTC lever.
    
2)  The further away from the rest position your CTC is,
     the less speed and torque you get.
    
3)  Your motor cannot output both maximum torque and maximum RPM at
     the same time, as one goes up, the other goes down."""

        ctcImgCaption = ttk.Label(ctcCalcFrame.captionFrame, text=imgCapt1)
        ctcImgCaption.grid(row=2, column=0, padx=4, pady=2, sticky='n')


        #create entry frames for the various inputs:
        capt1 = 'The length of your CTC lever \n(measured from the motor shaft centre to the rose joint centre)'
        self.ctcLength = EntryFrame(ctcCalcFrame.inputsFrame, title='CTC Length', entryvar=self.ctcLengthVar, unitvar=self.unitLengthShort, caption=capt1)

        capt2 = 'The torque output of your motor'
        motorTorque = EntryFrame(ctcCalcFrame.inputsFrame, title='Motor Torque', entryvar=self.motorTorqueVar, unitvar=self.unitTorque, caption=capt2)

        capt3= 'The speed of your motor'
        motorSpeed = EntryFrame(ctcCalcFrame.inputsFrame, title='Motor Speed', entryvar=self.motorSpeedVar, unit='(rpm)', caption=capt3)

        capt4= 'The maximum angle of travel of your CTC lever \n(measured in degrees away from the zero position, either clockwise or counter-clockwise)'
        maxAngle = EntryFrame(ctcCalcFrame.inputsFrame, title='Max. Angle', entryvar=self.maxAngleVar, unit='(Deg)', caption=capt4)

        capt4a= 'The current angle of your CTC lever, measured in degrees away from the zero position'
        currAngle = EntryFrame(ctcCalcFrame.inputsFrame, title='Current Angle', entryvar = self.currentAngleVar, unit='(Deg)', caption=capt4a)

        #grid the entries into the frame
        self.ctcLength.grid(row=1, column=0, padx=10, pady=3, sticky='w')
        motorTorque.grid(row=2, column=0, padx=10, pady=3, sticky='w')
        motorSpeed.grid(row=3, column=0, padx=10, pady=3, sticky='w')
        maxAngle.grid(row=4, column=0, padx=10, pady=3, sticky='w')
        currAngle.grid(row=5, column=0, padx=10, pady=3, sticky='w')

        #button to do the calculation
        ctcCalcButt = ttk.Button(ctcCalcFrame.outputsFrame, text='Calculate', command=self.ctcPageCalculate)
        ctcCalcButt.grid(row=1, column=2, padx=10, pady=15)

        #create various output frames to hold the calculation results
        capt5 = 'The equivalent "Straight-Line" speed of your CTC lever'
        linearSpeed = OutputFrame(ctcCalcFrame.outputsFrame, title='Linear Speed', variable=self.outputSpeedLinearVar, unitvar=self.unitSpeed, caption=capt5)

        capt6 = 'The equivalent force of your actuator'
        outputForce = OutputFrame(ctcCalcFrame.outputsFrame, title='Output Force', variable=self.outputForceLinearVar, unitvar=self.unitForce, caption=capt6)

        capt7 = 'The equivalent "Straight-Line" Travel of your CTC lever'
        maxLinTravel = OutputFrame(ctcCalcFrame.outputsFrame, title='Linear Travel', variable=self.outputLinearTravelVar, unitvar=self.unitLengthShort, caption=capt7)
       
        #grid the output frames into the output container
        linearSpeed.grid(row=1, column=0, padx=10, pady=3, sticky='w')
        outputForce.grid(row=2, column=0, padx=10, pady=3, sticky='w')
        maxLinTravel.grid(row=3, column=0, padx=10, pady=3, sticky='w')        


#GEOMETRY CALCULATION FRAME
#########################################################################
        geometryCalcFrameContainer = ttk.Frame(calcFrame)
        geometryCalcFrame = CalculationFrame(geometryCalcFrameContainer)

        #image stuff
        geoImgVar = ImageTk.PhotoImage(Image.open('diagrams/geoImg.PNG'))
        geoImg = ctcImg = ttk.Label(geometryCalcFrame.imgFrame, image=geoImgVar)
        geoImg.Image = geoImgVar
        geoImg.grid(row=0, column=0, padx=2, pady=2)

        imgCapt2="""
Use this tab to calculate the travel geometry of your sim.
Some important things to note:
1) Measure everything using the pivot of your sim as the reference point.
2) The X axis of the sim runs horizontally from side-to-side
3) The Y axis of the sim runs horizontally from front-to-back
4) The Z axis of the sim runs vertically
5) When measuring your motor position, measure from the center of the CTC lever.
6) Note that your equivalent linear force and travel numbers are calculated for the point where the rod attaches to the frame, so these numbers will differ
   depending on where you've attached yours"""

        #self.rodMtx = tk.StringVar()
        #self.rodMty = tk.StringVar()
        #self.rodMtz = tk.StringVar()
        #self.motorPosx = tk.StringVar()
        #self.motorPosy = tk.StringVar()
        #self.motorPosz = tk.StringVar()
        
        #self.motorAngleVar = tk.StringVar()
        #self.ctcRestAngleVar = tk.StringVar()

        geoImgCaption = ttk.Label(geometryCalcFrame.captionFrame, text=imgCapt2)
        geoImgCaption.grid(row=0, column=0)

        #create subcontainer to hold the inputs
        inputsSubContainer = ttk.Frame(geometryCalcFrame.inputsFrame)
        inputsSubContainer.grid(row=1, column=0, padx=5)

        #input frame for geometry parameters
        geoIpFrame = ttk.Frame(inputsSubContainer, relief='groove')
        geoIpFrame.grid(row=0, column=0, ipadx=7, ipady=7, padx=7)

        geoIpTitle = ttk.Label(geoIpFrame, text='Geometry Info')
        geoIpTitle.grid(row=0, column=0, pady=5)

        #inputs for different geometry positions
        rodMtCoords = TripleEntryFrame(geoIpFrame, title='Rod Mount Position', xvariable=self.rodMtx, yvariable=self.rodMty, zvariable=self.rodMtz, unitvar = self.unitLengthLong)
        rodMtCoords.grid(row=1, column=0, padx=7, pady=7)
        
        mtrCoords = TripleEntryFrame(geoIpFrame, title='Motor Position', xvariable=self.motorPosx, yvariable=self.motorPosy, zvariable=self.motorPosz,unitvar = self.unitLengthLong)
        mtrCoords.grid(row=2, column=0, padx=7, pady=7)

        mtrAngle = EntryFrame(geoIpFrame, title='Motor Angle', entryvar=self.motorAngleVar, unit='(deg)', caption='')
        mtrAngle.grid(row=3, column=0, padx=6, pady=7, sticky='e')
        
      
        #create subframe for CTC related inputs
        ctcParamsFrame = ttk.Frame(inputsSubContainer, relief='groove')
        ctcParamsFrame.grid(row=0, column=1, ipadx=7, ipady=7, padx=7)
        ctcParamslb = ttk.Label(ctcParamsFrame, text='CTC Parameters')
        ctcParamslb.grid(row=0, column=0, pady=5)

        #NOT USED
        useOldCTCParams = tk.IntVar()
        #checkbutton sends a 1 if the user wants to use the parameters from the other page DISABLED FOR BETA BUILD
        usePageParams = ttk.Checkbutton(ctcParamsFrame, variable=useOldCTCParams)
        #usePageParams.grid(row=0, column=0, sticky='e')
        chButtlb = tk.Label(ctcParamsFrame,text='Use Parameters from Page 1')
        #chButtlb.grid(row=0, column=1, sticky='w', pady=3)
        #NOT USED

        #define input entries
        ctcLenGeo = EntryFrame(ctcParamsFrame, title='CTC Length', entryvar=self.ctcLengthVar, caption='', unitvar=self.unitLengthShort)
        ctcAng = EntryFrame(ctcParamsFrame, title='CTC Rest Angle', entryvar=self.ctcRestAngleVar, caption='', unit='(deg)')
        ctcTrav = EntryFrame(ctcParamsFrame, title='Max. Angle (+/-)', entryvar=self.maxAngleVar, caption='', unit='(deg)')
        currentAngle = EntryFrame(ctcParamsFrame, title='Current Angle', entryvar = self.currentAngleVar, unit='(Deg)', caption = EMPTY)
        mTorque = EntryFrame(ctcParamsFrame, title='Motor Torque', entryvar=self.motorTorqueVar, caption='', unitvar=self.unitTorque)
        mSpeed = EntryFrame(ctcParamsFrame, title='Motor Speed', entryvar=self.motorSpeedVar, caption='', unit='(rpm)')

        #grid input entries
        ctcLenGeo.grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ctcAng.grid(row=2, column=0, sticky='w', padx=5, pady=5)
        ctcTrav.grid(row=3, column=0, sticky='w', padx=5, pady=5)
        mTorque.grid(row=4, column=0, sticky='w', padx=5, pady=5)
        mSpeed.grid(row=5, column=0, sticky='w', padx=5, pady=5)
        currentAngle.grid(row=6, column=0, sticky='w', padx=5, pady=5)

        outputsSubContainer = ttk.Frame(geometryCalcFrame.outputsFrame)
        outputsSubContainer.grid(row=1, column=0, ipadx=10, ipady=10, padx=5)

        #calculation button
        geoCalcButt = ttk.Button(outputsSubContainer, text='Calculate', command = self.geometryPageCalculate)
        geoCalcButt.grid(row=0, column=2, padx=10, pady=15)

        #output question tooltip button
        geoQbutt = ttk.Button(outputsSubContainer, text='Help', command = self.geometryPageTooltip)
        geoQbutt.grid(row=1, column=2, padx=10, pady=15)

        #define output entries
        pitchEffEntry = OutputFrame(outputsSubContainer, title='Pitch Efficiency', variable=self.pitchEff, unit='%')
        pitchForceEntry = OutputFrame(outputsSubContainer, title='Pitch Force', variable=self.pitchForce, unitvar=self.unitForce)
        rollEffEntry = OutputFrame(outputsSubContainer, title='Roll Efficiency', variable=self.rollEff, unit='%')
        rollForceEntry = OutputFrame(outputsSubContainer, title='Roll Force', variable=self.rollForce, unitvar=self.unitForce)
        pitchRotEntry = OutputFrame(outputsSubContainer, title='Pitch Rotation', variable=self.pitchRot, unit='(Deg)')
        rollRotEntry = OutputFrame(outputsSubContainer, title='Roll Rotation', variable=self.rollRot, unit='(Deg)')
        pitchTravEntry = OutputFrame(outputsSubContainer, title='Pitch Travel', variable=self.pitchTrav, unitvar=self.unitLengthShort)
        rollTravEntry = OutputFrame(outputsSubContainer, title='Roll Travel', variable=self.rollTrav, unitvar=self.unitLengthShort)
        pitchTorqueEntry = OutputFrame(outputsSubContainer, title='Pitch Torque', variable=self.pitchTorque, unitvar=self.unitTorque)
        rollTorqueEntry = OutputFrame(outputsSubContainer, title='Roll Torque', variable=self.rollTorque, unitvar=self.unitTorque)

        #grid output entries
        pitchEffEntry.grid(row=0, column=0, padx=15, pady=3)
        rollEffEntry.grid(row=0, column=1, padx=15, pady=3)
        pitchForceEntry.grid(row=1, column=0, padx=15, pady=3)
        rollForceEntry.grid(row=1, column=1, padx=15, pady=3)
        pitchRotEntry.grid(row=2, column=0, padx=15, pady=3)
        rollRotEntry.grid(row=2, column=1, padx=15, pady=3)
        pitchTravEntry.grid(row=3, column=0, padx=15, pady=3)
        rollTravEntry.grid(row=3, column=1, padx=15, pady=3)
        pitchTorqueEntry.grid(row=4, column=0, padx=15, pady=3)
        rollTorqueEntry.grid(row=4, column=1, padx=15, pady=3)
        

        
#SPEED AND FORCE CALCULATION FRAME
#########################################################################        
        #speedForceCalFrame = ttk.Frame(calcFrame)

#MOTOR CALCULATION FRAME
#########################################################################               
        #motorCalFrame = ttk.Frame(calcFrame)



#ADD CALCULATION TABS TO THE TTK NOTEBOOK        
        calcTabs.add(ctcCalcFrameContainer, text="CTC Lever Calculations")
        calcTabs.add(geometryCalcFrameContainer, text="Sim Geometry Calculations")
        #calcTabs.add(geometryCalcFrameContainerb, text="Sim Geometry Calculations (Roll)")
        #calcTabs.add(speedForceCalFrame, text = "Speed and Force Calculations")
        #calcTabs.add(motorCalFrame, text="Motor Calculations")
        


#INITIALIZE MAINLOOP
##########################################################################               
if __name__ == "__main__":


    
    root = tk.Tk()
    root.geometry(defaultWindowSize)
    #root.minsize(minWindowSizex, minWindowSizeY)
    app = Application(root)
    root.iconbitmap(r'data\favicon.ico')
    root.mainloop()
    




