###############################################################################
#                             Python modules                                  #
###############################################################################

from math import * 
import numpy as np 
from scipy.interpolate import splrep, BSpline
import matplotlib.pyplot as plt 
from prettytable import PrettyTable
from stl import mesh 
import os 
import time 
from rich.console import Console 
from rich.table import Table
from rich.progress import Progress 
from rich.traceback import install 
import copy
import quaternion 
import shutil 
import zipfile
from pathlib import Path
import time

###############################################################################
#                         Braid class & functions                             #
###############################################################################


def director_vector(Point1,Point2):
    return (Point2 - Point1)/np.linalg.norm(Point2 - Point1)

class Braid :

    def __init__(self, Mandrel_Diameter, N_HG, N_slot,                                             # Braider Properties
     N_BY, N_AX, Pitch, N_Layers, N_Fiber_Per_Yarn, Yarn_Compaction, a_BY, b_BY, a_AX, b_AX,        # Braid Properties
     Reinf_Density, Matrix_Density, Ff_BY, Fm_BY, Ff_AX, Fm_AX,                                     # Material Properties
     Space_Mandrin, Space_BY_AX, Start_AX, N_Pitch_AX, N_Pitch_BY,                                  # Braid Adjustments
     NofInterpolPts, Smoothing):                                                                    # Spline Parameters

     # Mandrel Properties
     self.Mandrel_Diameter = Mandrel_Diameter

     # Braider Properties
     self.N_HG = N_HG
     self.N_slot = N_slot

     # Braid Properties
     self.N_BY = N_BY
     self.N_AX = N_AX
     self.Pitch = Pitch
     self.Type = Type # Structure of the braid, 1 for diamond, 2 for regular

     # Yarn Properties
     ## Shared
     self.N_Layers = N_Layers
     self.N_Fiber_Per_Yarn = N_Fiber_Per_Yarn
     self.Yarn_Compaction = Yarn_Compaction
     ## Section geometry axis
     self.a_BY = a_BY
     self.b_BY = b_BY
     self.a_AX = a_AX
     self.b_AX = b_AX

     # Material Properties
     ## Fiber Density
     self.Reinf_Density = Reinf_Density
     self.Matrix_Density = Matrix_Density

     self.rhof_BY = Reinf_Density
     self.rhof_AX = Reinf_Density 
     self.rhom_BY = Matrix_Density 
     self.rhom_AX = Matrix_Density
     ## Fiber Fineness
     self.Ff_BY = Ff_BY
     self.Fm_BY = Fm_BY
     self.Ff_AX = Ff_AX
     self.Fm_AX = Fm_AX

     ## Yarn matrix material areas (10-6 m2)
     # In the Brainding Yarns (BY)
     self.Sf_BY = Ff_BY/(self.rhof_BY)
     self.Sm_BY = Fm_BY/(self.rhom_BY)

     # In the Axial Yarns (AX)
     self.Sf_AX = Ff_AX/(self.rhom_AX)
     self.Sm_AX = Fm_AX/(self.rhom_AX)

     # Space Parameters
     self.Space_Mandrin = Space_Mandrin
     self.Space_BY_AX = Space_BY_AX
     self.Start_AX = Start_AX
     self.N_Pitch_AX = N_Pitch_AX
     self.N_Pitch_BY = N_Pitch_BY


     # Spline Parameters
     self.NofInterpolPts = NofInterpolPts
     self.Smoothing = Smoothing

    def Init_BY(self):
        "This method initialize the BY dictionnary values and fields"
        self.BY = {} # Initialize the BY dictionnary
        Length = ((self.N_HG)*2)*self.N_Pitch_BY
        for k in range(self.N_BY):  # For all braiding yarns 
            self.BY[k] = {'HG': np.zeros(Length),  # horngear index of the yarn 
                        'HGPos': [],    # where is the yarn on the horngear? (in, mid, out)
                        'YARNa' : self.a_BY[:,k], # set dimensions for yarn (ellipse shape)
                        'YARNb' : self.b_BY[:,k],
                        'R' : np.zeros(Length), # radius for the braid (a voir which one it is)
                        'Theta' : np.zeros(Length), # crown angle
                        'X' : np.zeros(Length), # coordonees x,y,z d'un p
                        'Y' : np.zeros(Length),
                        'Z' : np.zeros(Length),
                        'Braid_angle' : np.zeros(Length)} 

            for j in range(Length): # For all positions along pitch*Npitch
                self.BY[k]['HGPos'].append('n/d') 
    
    def Init_AX(self):
        "This method initialize the AX dictionnary values and fields"
        self.AX = {} # Initialize the AX dictionnary

        Start = (self.Start_AX)*(self.Pitch)
        End = (self.N_Pitch_AX)*(self.Pitch)
        Num =  2*(self.N_HG) + 1

        for k in range(self.N_AX): # For all axial yarns
            self.AX[k] = {'YARNa' : 1,
                        'YARNb' : self.b_AX[k],
                        'X' : 0,
                        'Y' : 0,
                        'Z' : np.linspace(Start, End, Num), 
                        'phi' : 2*np.pi*k/(self.N_HG), #crown angle
                        'R' : 0 }
    
    def Init_BY_Position(self):
        "This method creates braid yarns with respect to the selected type"
        
        # Create type 1 (1 over 1) BY
        if int(self.Type) == 1 :
            start = 0
            stop = self.N_HG
            step = 2
            for i in range(start, stop, step ): # For all even numbered horngears attribute starting horngear and position for counterclockwise yarns (even number : 0,2,4,6,...)
                k = i # The yarn number k is the same as the horngear number
                self.BY[k]['HG'][0] = i
                self.BY[k]['HGPos'][0] = 'OUT'

                k = i + 1 # The yarn number k is a unit higher than the horngear number
                self.BY[k]['HG'][0] = i
                self.BY[k]['HGPos'][0] = 'IN'

        # Create type 2 (2 over 2) BY
        if int(self.Type) == 2 :
            start = 0
            stop = self.N_HG
            step = 2
            for i in range(start,stop,step) : # for all even numbered horngears attribute starting horngear and position for counterclockwise yarns (even number : 0,2,4,6,..., NHG)
                k = i # The yarn number k is the same as the horngear number
                self.BY[k]['HG'][0] = i
                self.BY[k]['HGPos'][0] = 'OUT'

                k = i + 1 # The yarn number k is a unit higher than the horngear number
                self.BY[k]['HG'][0] = i
                self.BY[k]['HGPos'][0] = 'MID'

            for i in range(start, stop, step): # For all even numbered horngears attribute starting horngear and position for clockwise yarns(odd numbers : NHG+1, NHG+3, NHG+4, NHG+5,...,2*NHG-1)
                k = self.N_HG + i  # The yarn number is odd 
                self.BY[k]['HG'][0] = i + 1 
                self.BY[k]['HGPos'][0] = 'IN'
                
                # attribute starting horngear and position for counterclockwise yarns 
                # (even numbers: NHG+2, NHG+4, NHG+6, ..., 2*NHG)
                k = self.N_HG + 1 + i # The yarn number is even
                self.BY[k]['HG'][0] = i + 1
                self.BY[k]['HGPos'][0] = 'MID'
    
    def Propagation_Yarn_Position(self):
        "This method propagate the yarn position for the pitch"
        CW_start = 0
        CCW_start = 1
        End = self.N_BY
        Step = 2
        Length = ((self.N_HG)*2)*self.N_Pitch_BY
        for k in range(CW_start, End, Step): # For clockwise yarns, they have odd numbers
            for j in range(1, Length): # For all yarn positions along pitch after the first one

                if not(self.BY[k]['HGPos'][j-1] == 'MID'): # if the yarn at previous location is at 'IN' or 'OUT' position on the horngear
                    self.BY[k]['HG'][j] = self.BY[k]['HG'][j-1] # the yarn stays on that horngear
                    self.BY[k]['HGPos'][j] = 'MID' # the yarn moves to the 'MID' position

                else : # if it is not
                    if self.BY[k]['HG'][j-1] == self.N_HG - 1 : # if it is on the last horngear
                        self.BY[k]['HG'][j] = 0 # it returns on the horngear 0
                    else :
                        self.BY[k]['HG'][j] = self.BY[k]['HG'][j-1] + 1
                    
                    if self.BY[k]['HG'][j] % 2 == 1 : # the carrier is on a horngear numbered as odd
                        self.BY[k]['HGPos'][j] = 'IN'
                    else : # the carrier is on a horngear numbered as even
                        self.BY[k]['HGPos'][j] = 'OUT'

        for k in range(CCW_start, End, Step): # For CounterClockWise yarns, they have even numbers 
            for j in range(1, Length): # For all yarn positions along pitch after the first one 
                
                if not(self.BY[k]['HGPos'][j-1] == 'MID'): # if the yarn at the previous location is at in or out position on a horngear 
                    
                    if self.BY[k]['HG'][j-1] == 0: # if the yarn at previous position is on the first horngear 
                        self.BY[k]['HG'][j] = self.N_HG-1 # put it on the last one 
                    else:
                        self.BY[k]['HG'][j] = self.BY[k]['HG'][j-1]-1 # otherwise put it on the previous one 
                    self.BY[k]['HGPos'][j] = 'MID' # The position will be the mid one 

                else: # if the yarn is on the mid position 
                    self.BY[k]['HG'][j] = self.BY[k]['HG'][j-1] # it stays on the same horngear 
                    if self.BY[k]['HG'][j] % 2 == 1: # the carrier is on a horngear numbered as odd
                        self.BY[k]['HGPos'][j] = 'OUT'
                    else: # the carrier is placed on a horngear numbered as even
                        self.BY[k]['HGPos'][j] = 'IN'

    def Map_Yarn_Positions_Horngears(self):
        self.HGMap = {}
        "This method creates a map of the yarn positions on horngears"
        Length = ((self.N_HG)*2)*self.N_Pitch_BY
        # Initialize the horngear map
        for j in range(Length):
            self.HGMap[j] = {'HG' : {}}
            for i in range(self.N_HG):
                self.HGMap[j]['HG'][i] = {'HGNumPos' : -1*np.ones(3)}

        # Assign Yarn Position into the Horngear Map
        for j in range(Length):
            for i in range(self.N_HG):
                for k in range(self.N_BY):
                    if self.BY[k]['HG'][j] == i and self.BY[k]['HGPos'][j] == 'OUT' :
                        self.HGMap[j]['HG'][i]['HGNumPos'][0] = k

                    if self.BY[k]['HG'][j] == i and self.BY[k]['HGPos'][j] == 'MID' :
                        self.HGMap[j]['HG'][i]['HGNumPos'][1] = k

                    if self.BY[k]['HG'][j] == i and self.BY[k]['HGPos'][j] == 'IN' :
                        self.HGMap[j]['HG'][i]['HGNumPos'][2] = k
    
    def Assign_Coordinates(self):
        "This method assigns the coordinates to BY and AX yarns"
        Length = ((self.N_HG)*2)*self.N_Pitch_BY
        D_IN = (self.Mandrel_Diameter)/2 + self.Space_Mandrin 
        # AX Yarn
        for i in range(self.N_HG): # for all horngear
            self.AX[i]['R'] = D_IN
            self.AX[i]['X'] = self.AX[i]['R'] * cos(self.AX[i]['phi'])
            self.AX[i]['Y'] = self.AX[i]['R'] * sin(self.AX[i]['phi'])
        
        # BY Yarn
        ## IN 
        for k in range(self.N_BY): # for all the braided yarns from 0 to N_BY-1
            for j in range(Length): # for all pitch 
                if self.BY[k]['HGPos'][j] == 'IN' : # if the yarn is inside the braid
                    # what is the horngear before ?
                    PrevHG = self.BY[k]['HG'][j] # it is the same HG

                    # what id the horngear after ?
                    if self.BY[k]['HG'][j] == self.N_HG - 1 :
                        NextHG = 0
                    else :
                        NextHG = self.BY[k]['HG'][j] + 1

                    self.BY[k]['R'][j] = D_IN - self.Space_BY_AX - (self.BY[k]['YARNb'][j])/2 - self.AX[self.BY[k]['HG'][j]]['YARNb']/2
                    self.BY[k]['Braid_angle'][j] = atan(2*pi*(self.BY[k]['R'][j])/(self.Pitch))
                    self.BY[k]['Theta'][j] = 2*pi*(self.BY[k]['HG'][j])/(self.N_HG)

        ## OUT
        for k in range(self.N_BY): # for all yarns from 1 to N_BY
            for j in range(Length): # for all position along the pitch
                if self.BY[k]['HGPos'][j] == 'OUT': # It the yarn is outside the braid
                                        # what is the horngear before ?
                    PrevHG = self.BY[k]['HG'][j] # it is the same HG

                    # what id the horngear after ?
                    if self.BY[k]['HG'][j] == self.N_HG - 1 :
                        NextHG = 0
                    else :
                        NextHG = self.BY[k]['HG'][j] + 1
                    self.BY[k]['R'][j] = D_IN + self.Space_BY_AX  + (self.BY[k]['YARNb'][j])/2 + self.AX[self.BY[k]['HG'][j]]['YARNb']/2
                    self.BY[k]['Braid_angle'][j] = atan(2*pi*(self.BY[k]['R'][j])/(self.Pitch))
                    self.BY[k]['Theta'][j] =2*pi*(self.BY[k]['HG'][j])/(self.N_HG)
        
        ## MID
        for k in range(self.N_BY): # for all the braided yarns from 1 to N_BY
            for j in range(Length): # for all positions along the pitch
                if self.BY[k]['HGPos'][j] == 'MID' : # if the yarn is at the 'MID' positon of the Braid
                    # what is the horngear before ?
                    PrevHG = self.BY[k]['HG'][j] # it is the same HG

                    # what id the horngear after ?
                    if self.BY[k]['HG'][j] == self.N_HG - 1 :
                        NextHG = 0
                    else :
                        NextHG = self.BY[k]['HG'][j] + 1
                    
                    self.BY[k]['R'][j] = D_IN
                    self.BY[k]['Braid_angle'][j] = atan(2*pi*(self.BY[k]['R'][j])/(self.Pitch))
                    self.BY[k]['Theta'][j] = 2*pi/self.N_HG*(self.BY[k]['HG'][j] + self.AX[PrevHG]['YARNa']/(self.AX[PrevHG]['YARNa']+self.AX[NextHG]['YARNa']))
                
        for k in range(self.N_BY): # for all the braided yarns from 1 to N_BY 
            for j in range(Length): # for all positions along the pitch 
                self.BY[k]['X'][j] = self.BY[k]['R'][j]*cos(self.BY[k]['Theta'][j])
                self.BY[k]['Y'][j] = self.BY[k]['R'][j]*sin(self.BY[k]['Theta'][j])
                self.BY[k]['Z'][j] = self.Pitch/(self.N_HG*2)*(j) 

    def Data_Structure_Init(self):
        "Creation of a data structure for the spine points"
        self.SplineDataStruct = {'BY' : {}, 'AX' : {}}
        fields = ['YARNa', 'YARNb', 'HGPos', 'X', 'Y', 'Z', 'phi', 'theta', 'alpha']
        # Theta => Braiding angle
        # Alpha => Crimp angle
        # Phi => Rotation around Z axis

        ## Go fetch the Yarn numbers
        #BY

        for key in self.BY.keys():
            self.SplineDataStruct['BY'][key] = {}
            #Creating fields
            for field in fields:
                self.SplineDataStruct['BY'][key][field] = None

        #AX
        for key in self.AX.keys():
            self.SplineDataStruct['AX'][key] = {}
            for field in ['YARNa', 'YARNb', 'X', 'Y', 'Z']:
                self.SplineDataStruct['AX'][key][field] = None

    def Create_BY(self):
        "This method create the coordinates of the ellipse's point"
        # BY
        for i in range(self.N_BY):
            x = np.hstack((self.BY[i]['X'][:], self.BY[i]['X'][0])) # Add last X at j = 2*NHG+1
            x = np.hstack((x[-3:-1], x, x[1:3]))
            y = np.hstack((self.BY[i]['Y'][:], self.BY[i]['Y'][0])) # Add last Y at j = 2*NHG+1
            y = np.hstack((y[-3:-1], y, y[1:3]))
            z = np.hstack((self.BY[i]['Z'][:],self.Pitch))
            z = np.hstack((z[-3:-1] - self.Pitch, z, z[1:3] + self.Pitch))
            t = np.linspace(1 ,2*(self.N_HG)*self.N_Pitch_BY + 1, 2*(self.N_HG)*self.N_Pitch_BY + 5)
            tt = np.linspace(t[0],t[-1],(2*(self.N_HG)*self.N_Pitch_BY + 5) + (2*(self.N_HG)*self.N_Pitch_BY + 4)*(self.NofInterpolPts))

            tck_x = splrep(t,x,s=self.Smoothing)
            tck_y = splrep(t,y,s=self.Smoothing)
            tck_z = splrep(t,z,s=self.Smoothing)

            xx = BSpline(*tck_x)(tt)
            xx = xx[2*(self.NofInterpolPts+1):-2*(self.NofInterpolPts+1)]
            yy = BSpline(*tck_y)(tt)
            yy = yy[2*(self.NofInterpolPts+1):-2*(self.NofInterpolPts+1)]
            zz = BSpline(*tck_z)(tt)
            zz = zz[2*(self.NofInterpolPts+1):-2*(self.NofInterpolPts+1)]

            # Saving data into the SplineDataStruct
            self.SplineDataStruct['BY'][i]['X'] = xx 
            self.SplineDataStruct['BY'][i]['Y'] = yy 
            self.SplineDataStruct['BY'][i]['Z'] = zz 
            self.SplineDataStruct['BY'][i]['YARNa'] = np.zeros(len(zz))
            self.SplineDataStruct['BY'][i]['YARNa'] = np.interp(self.SplineDataStruct['BY'][i]['Z'],np.hstack((self.BY[i]['Z'],self.Pitch)),\
                                                        np.hstack((self.BY[i]['YARNa'],self.BY[i]['YARNa'][0])))
            self.SplineDataStruct['BY'][i]['YARNb'] = np.zeros(len(zz))
            self.SplineDataStruct['BY'][i]['YARNb'] = np.interp(self.SplineDataStruct['BY'][i]['Z'],np.hstack((self.BY[i]['Z'],self.Pitch)),\
                                                        np.hstack((self.BY[i]['YARNb'],self.BY[i]['YARNb'][0])))
            
            self.SplineDataStruct['BY'][i]['phi'] = np.zeros(len(zz)) # all rotation around Z axis are set to 0 rad 
            self.SplineDataStruct['BY'][i]['theta'] = np.zeros(len(zz)) # braiding angles are set to 0 rad
            self.SplineDataStruct['BY'][i]['alpha'] = np.zeros(len(zz)) # crimp angles 
            
        # AX   
        for i in range(self.N_AX): 
            x = self.AX[i]['X']*np.ones(2*self.N_HG+5)
            y = self.AX[i]['Y']*np.ones(2*self.N_HG+5)
            z = np.hstack((self.AX[i]['Z'][-3:-1]-6*self.Pitch,self.AX[i]['Z'],self.AX[i]['Z'][1:3]+6*self.Pitch))
            t = np.linspace(1,2*self.N_HG+1,2*self.N_HG+5)
            tt = np.linspace(t[0],t[-1],(2*self.N_HG + 5) + (2*self.N_HG + 4)*6*self.NofInterpolPts)
            tck_x = splrep(t,x,s=self.Smoothing)
            tck_y = splrep(t,y,s=self.Smoothing)
            tck_z = splrep(t,z,s=self.Smoothing)
            xx = BSpline(*tck_x)(tt)[6*2*(self.NofInterpolPts+1):-6*2*(self.NofInterpolPts+1)]
            yy = BSpline(*tck_y)(tt)[6*2*(self.NofInterpolPts+1):-6*2*(self.NofInterpolPts+1)]
            zz = BSpline(*tck_z)(tt)[6*2*(self.NofInterpolPts+1):-6*2*(self.NofInterpolPts+1)]
            self.SplineDataStruct['AX'][i]['X'] = xx 
            self.SplineDataStruct['AX'][i]['Y'] = yy
            self.SplineDataStruct['AX'][i]['Z'] = zz
            self.SplineDataStruct['AX'][i]['YARNa'] = self.AX[i]['YARNa']
            self.SplineDataStruct['AX'][i]['YARNb'] = self.AX[i]['YARNb']
                    
    def Ellipse(self, aby, bby, npts):
        "Method that creates a cloud of points to plot an ellipse"
        t = np.linspace(0, 2*np.pi, npts)
        points = np.zeros((npts,3))
        points[:,0] = bby/2 * np.cos(t)
        points[:,1] = aby/2 * np.sin(t)
        return points
        
    def SquareSpacking(self,S,YarnCompaction):
        #thickness = .1 mm, width = 13mm
        "Method that creates the points corresponding to each fiber in a Yarn AS A SQUARE"
        N_layers = self.N_Layers
        NFiberperYarn = self.N_Fiber_Per_Yarn
        points = np.zeros((NFiberperYarn,3))
        #Tolerance between layers to match with compaction
        Tolerance = (1/float(N_layers)) * ((1/np.sqrt(YarnCompaction)) - (2 * N_layers + 1)/(np.sqrt(3*N_layers *(N_layers +1 ) + 1))) * np.sqrt(S/np.pi)
        MaxCompaction = (NFiberperYarn)/(float((2 * N_layers + 1)**2)) # Maximun compaction possible

        if (Tolerance < 0): # if compaction not reachable
            print("Compaction too high, reduced to " + str(100 * MaxCompaction) + " %")
            Tolerance = 0 # reduce compaction to maximun possible
        
        # Generate the fibers of each layer 
        NFibersGenerated = 0 # number of fiber placed yet        
        for i in range(1, N_layers + 1):
            NFiberperLayer = 4 * i
            length = 0.1e-3 # thickness
            width = 6e-3 #S/length 
            perimeter = 2 * (length + width)
            points_per_side = i + 1
            # because we are in a rectangle: points aren't evenly distributed on the whole perimeter
            # decision: create 4 distinct linspaces (1 for each side of the rectangle)
            t_bottom = np.linspace(0, length, points_per_side)
            t_right = np.linspace(length, length + width, points_per_side)
            t_top = np.linspace(length+width, 2 * length + width, points_per_side)
            t_left = np.linspace(2*length + width, perimeter, points_per_side)

            t_common = np.concatenate([t_bottom, t_right, t_top, t_left])

            # remove common points (by removing the last point of every linspace)
            t = t_common[::2]
            for j in range(NFiberperLayer):
                # we start at the bottom left corner and go counter-clockwise
                if t[j] <= length:  # Bottom side
                    #points (0) = x, (1) = y
                    points[NFibersGenerated + j, 0] = -length/2 + t[j]
                    points[NFibersGenerated + j, 1] = -width/2
                elif t[j] <= (width + length):  # Right side
                    distanceOnSide = t[j] - length
                    points[NFibersGenerated + j, 0] = length/2
                    points[NFibersGenerated + j, 1] = -width/2 + distanceOnSide
                elif t[j] <= (width + 2 * length): # Top side
                    distanceOnSide = t[j] - width - length
                    points[NFibersGenerated + j, 0] = length/2 - distanceOnSide
                    points[NFibersGenerated + j, 1] = width/2
                else:  # Left side
                    distanceOnSide = t[j] - width - 2 * length
                    points[NFibersGenerated + j, 0] = -length/2
                    points[NFibersGenerated + j, 1] = width/2 - distanceOnSide
                
            NFibersGenerated += NFiberperLayer
        
        #First Fiber in the center
        return points
    def CircleSpacking(self,S,YarnCompaction):
        "Method that creates the points corresponding to each fiber in a Yarn"
        
        N_layers = self.N_Layers
        NFiberperYarn = self.N_Fiber_Per_Yarn
        points = np.zeros((NFiberperYarn,3))
        FiberRadius = np.sqrt(S/(NFiberperYarn*np.pi)) # radius of each fiber
        #Tolerance between layers to match with compaction
        Tolerance = (1/float(N_layers)) * ((1/np.sqrt(YarnCompaction)) - (2 * N_layers + 1)/(np.sqrt(3*N_layers *(N_layers +1 ) + 1))) * np.sqrt(S/np.pi)
        MaxCompaction = (NFiberperYarn)/(float((2 * N_layers + 1)**2)) # Maximun compaction possible

        if (Tolerance < 0): # if compaction not reachable
            print("Compaction too high, reduced to " + str(100 * MaxCompaction) + " %")
            Tolerance = 0 # reduce compaction to maximun possible
        
        # Generate the fibers of each layer 
        NFibersGenerated = 0 # number of fiber placed yet        
        for i in range(1, N_layers + 1):
            NFiberperLayer = 6 * i
            t = np.linspace (0, 2 * np.pi, NFiberperLayer + 1)
            LayerRadius = i * (2 * FiberRadius + Tolerance)
            for j in range(NFiberperLayer): # goes from point 0 to point 2pi - 1 (thus ignoring NFiberperLayer + 1)
                points[NFibersGenerated + j, 0] = LayerRadius * np.cos(t[j])
                points[NFibersGenerated + j, 1] = LayerRadius * np.sin(t[j])
            NFibersGenerated += NFiberperLayer
        
        #First Fiber in the center
        points[NFiberperYarn - 1, 0] = 0
        points[NFiberperYarn - 1, 1] = 0

        return points

    def Rotate(self, points, vdir, vini, axe):
        "Method that rotates points around an axis"
        #Normalisation of vectors
        vdir = vdir/np.linalg.norm(vdir)
        vini = vini/np.linalg.norm(vini)
        theta = np.arccos(np.dot(vdir,vini))

        if axe.lower() == 'x': 
            Rx = np.array([
                [1, 0, 0],
                [0, np.cos(theta), np.sin(theta)],
                [0, -np.sin(theta), np.cos(theta)]
            ])
            points = np.dot(points, Rx)
        
        elif axe.lower() == 'y': 
            Ry = np.array([
                [np.cos(theta), 0, -np.sin(theta)],
                [0, 1, 0],
                [np.sin(theta), 0, np.cos(theta)]
            ])
            points = np.dot(points, Ry)
        
        elif axe.lower() == 'z': 
            theta = np.arctan2(vdir[1],vdir[0])
            Rz = np.array([
                [np.cos(theta), np.sin(theta), 0],
                [-np.sin(theta), np.cos(theta), 0],
                [0, 0, 1]
            ])
            points = np.dot(points,Rz)
        
        return points

    def Rotate_Arbitrary(self, points, vdir): # This was vidr but there were only three instances of it here and lots of vdir variables that seemed to be related. (typo?)
        "Allow the rotation of points around an arbitrary axis"
        vnorm = np.array([0., 0., 1.])
        vdir = vdir/np.linalg.norm(vdir)
        axis = np.cross(vnorm, vdir)
        axis = axis/np.linalg.norm(axis)
        theta = -np.arccos(np.dot(vnorm, vdir))

        a_skew = np.array([
            [0, - axis[2], axis[1]],
            [axis[2], 0, -axis[0]], 
            [-axis[1], axis[0], 0]
        ])

        R = np.cos(theta) * np.eye(3) + (1 - np.cos(theta)) * axis.reshape(3, 1) * axis + np.sin(theta) * a_skew
        points = np.dot(points, R)

        return points

    def Translate(self, points, vtrans):
        "Method that allows the translations of points in a particular direction"
        u = np.zeros(points.shape)
        for i in range(points.shape[1]):
            u[:,i] = points[:,i] + vtrans[i]
        points = u

        return u

    def Create_Wireframe_Model(self) :
            "This method creates a wireframe model for the BY and AX yarns"
            self.filDict = {'BY':{}, 'AX': {}}

            # BY
            for i in self.SplineDataStruct['BY'].keys(): # looping through the different yarns
                xx = self.SplineDataStruct['BY'][i]['X']
                yy = self.SplineDataStruct['BY'][i]['Y']
                zz = self.SplineDataStruct['BY'][i]['Z']
                self.filDict['BY'][i] = {}

                # We compute phi angle wich is the rotation around Z axis for all points of the yarn i
                self.SplineDataStruct['BY'][i]['phi'] = np.arctan2(yy, xx)
                table = PrettyTable(['Point No', 'x', 'y', 'z'])
                
                for j in range(xx.shape[0]):
                    if j == 0: 
                        vdir= director_vector(
                            np.array([xx[-2], yy[-2], zz[-2]-self.Pitch]),
                            np.array([xx[j+1],yy[j+1],zz[j+1]])
                        )
                        splinePoints = np.array([
                            [xx[-2], yy[-2], zz[-2]-self.Pitch],
                            [xx[j+1],yy[j+1],zz[j+1]]
                        ])
                    elif j == xx.shape[0]-1: 
                        vdir = director_vector( 
                            np.array([xx[j-1], yy[j-1], zz[j-1]]),
                            np.array([xx[1],yy[1],self.Pitch + zz[1]])
                        )
                        splinePoints = np.array([
                            [xx[j-1], yy[j-1], zz[j-1]],
                            [xx[1],yy[1],self.Pitch + zz[1]]
                        ])
                    else:
                        vdir = director_vector(
                            np.array([xx[j-1],yy[j-1],zz[j-1]]),
                            np.array([xx[j+1],yy[j+1],zz[j+1]])
                        )
                        splinePoints = np.array([
                            [xx[j-1], yy[j-1], zz[j-1]],
                            [xx[j+1],yy[j+1],zz[j+1]]
                        ])
                    
                    aby = self.SplineDataStruct['BY'][i]['YARNa'][j]
                    bby = self.SplineDataStruct['BY'][i]['YARNb'][j]

                    # Compute the equation of the tangent plan
                    # The equation of the plan is ax + by + cz + d = 0 where :
                    a = xx[j]
                    b = yy[j]
                    c = 0 
                    d = -(a * xx[j] + b * yy[j] + c * zz[j])

                    # We compute the projection of the splinePoints j-1 (point 0), and j+1 (point 2) on the tangent plan 
                    # Through this process we obtain the projections of splinePoints on the tangent plan
                    lambda_0 = (a * splinePoints[0, 0] + b * splinePoints[0, 1] + c * splinePoints[0, 2] + d)/((a**2)+(b**2)+(c**2))
                    lambda_1 = (a * splinePoints[1, 0] + b * splinePoints[1, 1] + c * splinePoints[1, 2] + d)/((a**2)+(b**2)+(c**2))

                    splinePointsProjection = splinePoints - np.array([
                        [lambda_0 * a, lambda_0 * b, lambda_0 * c],
                        [lambda_1 * a, lambda_1 * b, lambda_1 * c]
                    ])

                    # We compute the projection of the vdir vecor on the tangent plan
                    vdirTangent = splinePointsProjection[1, :] - splinePointsProjection[0, :]

                    # We compute the angle between the projection of the director vector on the tangent plan
                    # and the vector k = [0, 0 1]. It is the braiding angle
                    arg = np.dot(vdirTangent, np.array([0,0,1]))/np.linalg.norm(vdirTangent)
                    arg = min([arg, 1]) if arg >= 0 else max([arg, 1]) # sometimes arg = 1.0000x which can lead to error so we round it to 1
                    if i % 2 == 0 : #CW
                        self.SplineDataStruct['BY'][i]['theta'][j] = np.arccos(arg)
                    else :
                        self.SplineDataStruct['BY'][i]['theta'][j] = -np.arccos(arg)
                    
                    # We create the ellipse points
                    self.EllipseCloud = self.SquareSpacking(S = self.Sf_BY + self.Sm_BY, YarnCompaction = self.Yarn_Compaction) # Draw the fibers

                    # we create a dictionnary containing the coordinates of x,y,z axis
                    # we want the rotation axis to change after each rotation we perform
                    axis_vectors = {'x' : np.array([1, 0, 0]),
                                'y' : np.array([0, 1, 0]),
                                'z' : np.array([0, 0, 1])}
                    
                    # we create the first quaternion for the rotation around z axis
                    qz = quaternion.from_rotation_vector(self.SplineDataStruct['BY'][i]['phi'][j]*axis_vectors['z'])
                    
                    # We use qz to rotate the EllipseCloud around z axis
                    self.EllipseCloud = quaternion.rotate_vectors(qz,self.EllipseCloud)
                    
                    # we update the local axis
                    axis_vectors['z'] = quaternion.rotate_vectors(qz,axis_vectors['z'])
                    axis_vectors['y'] = quaternion.rotate_vectors(qz,axis_vectors['y'])
                    axis_vectors['x'] = quaternion.rotate_vectors(qz,axis_vectors['x'])

                    #we create the second quaternion for the rotation around x axis
                    qx = quaternion.from_rotation_vector(-self.SplineDataStruct['BY'][i]['theta'][j]*axis_vectors['x'])

                    # we use qx to rotate the EllipseCloud around x axis
                    self.EllipseCloud = quaternion.rotate_vectors(qx, self.EllipseCloud)

                    # we update the local axis
                    axis_vectors['z'] = quaternion.rotate_vectors(qx,axis_vectors['z'])
                    axis_vectors['y'] = quaternion.rotate_vectors(qx,axis_vectors['y'])
                    axis_vectors['x'] = quaternion.rotate_vectors(qx,axis_vectors['x'])

                    # We compute the equation of the ZX plan using the coordinates in the axis_vectors
                    # The y axis acts as the normal vector
                    a,b,c = axis_vectors['y'][0],axis_vectors['y'][1],axis_vectors['y'][2]
                    d = -(a*axis_vectors['x'][0]+b*axis_vectors['x'][1]+c*axis_vectors['x'][2])
                    lambda_0 = (a*splinePoints[0,0]+b*splinePoints[0,1]+c*splinePoints[0,2]+d)/(a**2+b**2+c**2)
                    lambda_1 = (a*splinePoints[1,0]+b*splinePoints[1,1]+c*splinePoints[1,2]+d)/(a**2+b**2+c**2)
                    
                    splinePointsProjection =  splinePoints - np.array([
                        [lambda_0 * a, lambda_0*b, lambda_0*c],
                        [lambda_1 * a, lambda_1*b, lambda_1*c]
                    ])

                    vdirProj = splinePointsProjection[1,:] - splinePointsProjection[0,:]

                    # if the x component of vdirProj is in the x axis direction (angle inferior to pi/2), the alpha is positive
                    # otherwise, he is negative. 
                    arg = np.dot(vdirProj,axis_vectors['x'])/(np.linalg.norm(vdirProj)*np.linalg.norm(axis_vectors['x']))
                    arg = min([arg,1]) if arg >= 0 else max([arg,-1])
                    arg2 = np.dot(vdirProj,axis_vectors['z'])/(np.linalg.norm(vdirProj)*np.linalg.norm(axis_vectors['z']))
                    arg2 = min([arg2,1]) if arg >= 0 else max([arg2,-1])
                    if np.arccos(arg) <= np.pi/2: 
                        # Crimp angle positive
                        self.SplineDataStruct['BY'][i]['alpha'][j] = np.arccos(arg2)

                    else:
                        # Crimp angle negative
                        self.SplineDataStruct['BY'][i]['alpha'][j] = -np.arccos(arg2)
                    
                    # We create the third quaternion for the rotation around y axis
                    qy = quaternion.from_rotation_vector(self.SplineDataStruct['BY'][i]['alpha'][j]*axis_vectors['y'])

                    # We use qy to rotate the EllipseCloud around y axis
                    self.EllipseCloud = quaternion.rotate_vectors(qy,self.EllipseCloud)

                    # we update the local axis
                    axis_vectors['x'] = quaternion.rotate_vectors(qy,axis_vectors['x'])
                    axis_vectors['y'] = quaternion.rotate_vectors(qy,axis_vectors['y'])
                    axis_vectors['z'] = quaternion.rotate_vectors(qy,axis_vectors['z'])
                    
                    # We translate the rotated ellipseCloud
                    self.EllipseCloud = self.Translate(points=self.EllipseCloud, vtrans=np.array([xx[j], yy[j], zz[j]]))
                    
                    # We save the ellipseCloud in the filDict structure
                    self.filDict['BY'][i][j] = self.EllipseCloud

            # AX Yarns 
            for i in self.SplineDataStruct['AX'].keys(): 
                xx = self.SplineDataStruct['AX'][i]['X']
                yy = self.SplineDataStruct['AX'][i]['Y']
                zz = self.SplineDataStruct['AX'][i]['Z']
                self.filDict['AX'][i] = {}
                for j in range(xx.shape[0]):
                    if j < xx.shape[0]-1 :
                        vdir = np.array([xx[j+1] - xx[j], yy[j+1] - yy[j], zz[j+1] - zz[j]])
                    else:
                        vdir = np.array([xx[j] - xx[j-1], yy[j] - yy[j-1], zz[j] - zz[j-1]]) # for the last point
                    a = self.SplineDataStruct['AX'][i]['YARNa']
                    b = self.SplineDataStruct['AX'][i]['YARNb']
                    Sf_AX=self.Sf_AX[i]
                    Sm_AX=self.Sm_AX[i]
                    self.EllipseCloud = self.SquareSpacking(S = Sf_AX+Sm_AX, YarnCompaction=self.Yarn_Compaction) # Draw the fibers 
                    self.EllipseCloud = self.Rotate(points=self.EllipseCloud,
                                                    vdir = np.array([xx[j], yy[j], 0.]), 
                                                    vini = np.array([1.,0.,0.]), 
                                                    axe = 'z') # rotate the ellipse around z axis
                    self.EllipseCloud = self.Translate(points=self.EllipseCloud, vtrans=np.array([xx[j], yy[j], zz[j]]))
                    self.filDict['AX'][i][j] = self.EllipseCloud

    def Create_SplinePoints_Files(self,filename = 'it_0'):
        """Method that creates spline points files for each fiber"""
        yarnType_ = ['BY','AX']

        iterationFilePath = os.path.join('Splines_Points_Files/', filename)
        os.makedirs(iterationFilePath, exist_ok=True)

        # Now the code is storing .txt files from previous runs of Braid_Class so this line resets the iterationFilePath folder
        for f in Path(iterationFilePath).glob('*.txt'):
            f.unlink()

        N_Fiber_Per_Yarn = self.N_Fiber_Per_Yarn

        for yarnType in yarnType_:
            if yarnType == 'BY': #number of yarns differs on each type of yarn
                N = self.N_BY
            else :
                N = self.N_AX
            for i in range(N): #for each yarn
                for k in range(self.N_Fiber_Per_Yarn): #for each fiber                    
                    FiberNumber = i * N_Fiber_Per_Yarn + k
                    globals()[f'{yarnType}_fiber'] = open(f'{iterationFilePath}/{filename}_{yarnType}_fiber_{FiberNumber}.txt','w') # open .txt file

                    lignes = []
                    for j in range(len(self.SplineDataStruct[yarnType][i]['X'])): #coordinates at each zone of each fiber
                        x = self.filDict[yarnType][i][j][k][0]
                        y = self.filDict[yarnType][i][j][k][1]
                        z = self.filDict[yarnType][i][j][k][2]
                        lignes.append(f"{x:.10f} {y:.10f} {z:.10f}\n")

                # write in file and close 
                globals()[f'{yarnType}_fiber'].writelines(lignes)
                globals()[f'{yarnType}_fiber'].close()                
        
        print(f'Each yarn has {N_Fiber_per_Yarn} fibers') # For understanding how many fibers are in the yarns

        # Now to store the file in the desired directory (Temp cause Abaqus likes that)
        output_dir = Path('c:/Temp') # Desired location
        output_dir.mkdir(parents=True, exist_ok=True) # Make sure it exists
        archive_path = output_dir / 'Fiber_Points' # Define file name and path

        #Put all coordinate files in a single archive
        shutil.make_archive(str(archive_path),'zip', iterationFilePath)
        
        # Clean up old fiber folder and extract to the temp folder (same as the zip file)
        zip_path = Path('c:/Temp/Fiber_Points.zip') # Note the .zip
        extract_to = zip_path.parent / zip_path.stem # Folder we want extracted files in
        if extract_to.exists() and extract_to.is_dir(): # Check if it exists 
           shutil.rmtree(extract_to) # Delete it
           print('Deleted existing folder:', extract_to) # Confirm it was deleted
        
        with zipfile.ZipFile(zip_path,'r') as zip_ref: # Extract the contents
            zip_ref.extractall(extract_to)
            print('Extracted files into', archive_path)     

    def Create_STL_File(self, file_name = 'it_0', path=' '):
        # from stl import mesh
        """Method that creates a visualizable STL file from the points from self.FilDict dictionary, updated from Pierre-Richard's code"""
        os.makedirs(os.path.join('STL_Files/',file_name),exist_ok=True) 
        combinedFilesNames = []

        yarnTypes = ['BY','AX']
        for yarnType in yarnTypes : 
            for yarn in self.filDict[yarnType].keys(): # yarn number (for every yarn)
                nbr_rectangles = int(len(self.filDict[yarnType][yarn].keys())) # cross section number
                nbr_points = 4
                nbr_vertices = nbr_rectangles * nbr_points
                vertices = np.zeros((nbr_vertices,3))

                # Define the vertices 
                for i in range(nbr_rectangles): 
                    vertices[i*nbr_points:(i+1)*nbr_points,:] = self.filDict[yarnType][yarn][i]
                
                # Define faces 
                nbr_faces = 2*(nbr_rectangles-1) * nbr_points
                faces = np.zeros((nbr_faces,3))
                #k = 0 
                faces_list = []
                for j in range(nbr_rectangles-1): # from 0 to before last 
                    for i in range(nbr_points): # range 0,1,2,3
                        i_next = (i + 1) % nbr_points
                        v00 = j * nbr_points + i
                        v01 = j * nbr_points + i_next
                        v10 = (j+1) * nbr_points + i
                        v11 = (j+1) * nbr_points + i_next
                        
                        # triangle 1
                        faces_list.append([v00, v01, v10])
                        # triangle 2
                        faces_list.append([v01, v11, v10])
                
                # version ancienne: plus rapide, plus gros triangles
                        """faces[k, :] = [i + j*nbr_points, i_next + j*nbr_points, i + (j+1)*nbr_points]
                        k += 1
                        faces[k, :] = [i + (j+1)*nbr_points, i_next + (j+1)*nbr_points, i_next + j*nbr_points]
                        k += 1"""
                faces = np.array(faces_list, dtype=int)
                #faces = faces.astype(int)
                
                # Create the mesh 
                globals()[f'{yarnType}_yarn_{yarn}'] = mesh.Mesh(np.zeros(faces.shape[0],dtype=mesh.Mesh.dtype))
                globals()[f'{yarnType}_yarn_{yarn}'].vectors[:] = vertices[faces]
                """for i, f in enumerate(faces): 
                    for j in range(3): 
                        globals()[f'{yarnType}_yarn_{yarn}'].vectors[i][j] = vertices[f[j],:]"""

                # Write the mesh to file with stl extension 
                globals()[f'{yarnType}_yarn_{yarn}'].save(f'STL_Files/{file_name}/{file_name}_{yarnType}_yarn_{yarn}.stl')

        
            # Combine all BY yarns and all AX yarns  
            filenames = []
            for yarn in self.filDict[yarnType].keys():
                filenames.append(f'STL_Files/{file_name}/{file_name}_{yarnType}_yarn_{yarn}.stl')
            globals()[f'{file_name}_{yarnType}_yarns'] = mesh.Mesh.from_files(filenames=filenames)
            globals()[f'{file_name}_{yarnType}_yarns'].save(f'STL_Files/{file_name}/{file_name}_{yarnType}_yarns.stl')
            combinedFilesNames.append(f'STL_Files/{file_name}/{file_name}_{yarnType}_yarns.stl')
            
        # Combine both BY yarns and AX yarns 
        globals()[f'{file_name}_BY_AX_yarns'] = mesh.Mesh.from_files(filenames=combinedFilesNames)
        globals()[f'{file_name}_BY_AX_yarns'].save(f'STL_Files/{file_name}/{file_name}_BY_AX_yarns.stl')

    def Create_STEP_File(self, file_name = 'it_0', path = ' '):
        """Method that creates STEP files from the self.FilDict dictionary"""
        os.makedirs(os.path.join('STEP_Files/',file_name),exist_ok=True) 
        import cadquery as cq
        # goal: create rectangular wires for every cross-section, combining them using cadquery.loft
        #print(self.filDict)
        yarnTypes = ['BY','AX']
        for yarnType in yarnTypes : 
            for yarn in self.filDict[yarnType].keys(): #for every yarn number
                sorted_sections = [self.filDict[yarnType][yarn][k] for k in sorted(self.filDict[yarnType][yarn].keys())]

                # Create cross sections at the ends of the yarn (solves singularity issues)
                delta = 1e-8
                debut = sorted_sections[0] - delta
                fin = sorted_sections[-1] + delta
                sorted_sections = [debut] + sorted_sections + [fin]
                sorted_sections_tuple = [[tuple(coordinates.tolist()) for coordinates in section] for section in sorted_sections]
                wp = cq.Workplane("XY")

                #testing purposes                
                """if yarn == 11: 
                    print(f"debut:{debut}")
                    print(sorted_sections[0])"""

                # Add all wires sequentially on the stack
                for section in sorted_sections_tuple:
                    pts = [cq.Vector(*p) for p in section]
                    wp = wp.polyline(pts).close()

                # export individual stp 
                solid = wp.loft(combine=True)
                combined_list.append(solid)
                file_path = f'STEP_Files/{file_name}/{file_name}_{yarnType}_yarn_{yarn}.step' 
                cq.exporters.export(solid, file_path)
                combined_solid = cq.Compound.makeCompound([s.val() for s in combined_list])

            combined_solid.export(f'STEP_Files/{file_name}/{file_name}_{yarnType}_yarns.step')

    def __run__(self):
        """This method calls the other in order to realize all the required operations"""
        self.Init_BY()
        self.Init_AX()
        self.Init_BY_Position()
        self.Propagation_Yarn_Position()
        self.Map_Yarn_Positions_Horngears()
        self.Assign_Coordinates()
        self.Data_Structure_Init()
        self.Create_BY()
        self.Create_Wireframe_Model()
        self.Create_STL_File()
        self.Create_SplinePoints_Files()
        #self.Create_STEP_File()
              



###############################################################################
#                              Parameters                                     #
###############################################################################

# Material Constants
Reinf_Density = 2460
Matrix_Density = 1080

# Shared yarn properties
Yarn_Compaction = 0.7

# BraiderProperties
N_HG = 12 # Number of horngears
N_slot = 4  # Number of slots per horngear

# Braiding properties 
Mandrel_Perimeter = 79.9e-3 #, changed for Abaqus Mandrel
Mandrel_Diameter = Mandrel_Perimeter/pi
Pitch = 80e-3
N_Pitch = 1
N_BY = 12 # Amount of braiding yarns
N_AX = 12 # Amount of YAXial Yarns

Type = N_BY/N_HG # Structure of the braid, 1 for diamond, 2 for regular, ...

# Yarn Fineness in Tex
# Fineness of the reinforcing fibres (f) and the polymer (m) in the braiding yarns (BY)
Ff_BY = 2640*1e-6
Fm_BY = 1520*1e-6  
# Fineness of the reinforcing fibres (f) and the polymer (m) in the axial yarns (AX)
# The vector contains the same amount of elements than horngears since axial yarns passes though horngear centers
# Ff_AX = np.array([1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584]) *1e-6
Ff_AX = np.array([584,584,584,584,584,584,584,584,584,584,584,584]) *1e-6
Fm_AX = np.array([912,912,912,912,912,912,912,912,912,912,912,912]) *1e-6
#Ff_AX = np.array([2640,528,2640,528,2640,528,2640,528,2640,528,2640, 528])*1e-6 # 2640, 528
#Fm_AX = np.array([1520,304,1520,304,1520,304,1520,304,1520,304, 1520, 304])*1e-6 # 1520,304

# Assign material properties to each yarn components (f for fibers, m for matrix)
rhof_BY = Reinf_Density
rhof_AX = Reinf_Density 

rhom_BY = Matrix_Density 
rhom_AX = Matrix_Density 
      
# Assign yarn geometrical properties
RpYBY = Yarn_Compaction
RpYAX = Yarn_Compaction

# Yarn fiber material areas     
Sf_BY = Ff_BY/rhof_BY
Sf_AX = Ff_AX/rhof_AX

# Yarn matrix material area     
SmYBY = Fm_BY/rhom_BY
SmYAX = Fm_AX/rhom_AX

# Yarn total fineness
FYBY = Ff_BY + Fm_BY
FYAX = Ff_AX + Fm_AX

# Yarn total material area
SYBY = Sf_BY + SmYBY
SYAX = Sf_AX + SmYAX

# Yarn unconsolidated Area
AYBY = SYBY/RpYBY 
AYAX = SYAX/RpYAX 
InitialYarnAspectRatio = 4

N_Layers = 1 #Number of fiber layers in a Yarn, 0 = one fiber per Yarn
N_Fiber_per_Yarn = 0
for i in range(N_Layers):
    NFiberperLayer = 4*(i+1)
    N_Fiber_per_Yarn+=NFiberperLayer #Number of Fibers per yarn


Space_Mandrin = 6e-3
Space_BY_AX = 1e-3
StartAX = -2
NPitchAX = 4
N_Pitch_BY =1
Smoothing = 0
NofInterpolPts = 3


# Compute initial aby and bby values 
# BY
#edited to get circles
Dby = np.sqrt((4*AYBY/np.pi)) * np.ones((2*N_HG * N_Pitch_BY,N_BY))
a_BY = Dby*1e-3
b_BY = Dby*1e-3

# AX
Dax = np.sqrt((4*AYAX/np.pi)) * np.ones(N_AX)
a_AX = Dax*1e-3
b_AX = Dax*1e-3

#functioncallcounter = 0
#global functioncallcounter
#functioncallcounter += 1
#print('Function was called')
        
start_time = time.perf_counter()

tresse_initiale = Braid(Mandrel_Diameter, N_HG, N_slot,
                N_BY, N_AX, Pitch, N_Layers, N_Fiber_per_Yarn, Yarn_Compaction, a_BY, b_BY, a_AX, b_AX,
                Reinf_Density, Matrix_Density, Ff_BY, Fm_BY, Ff_AX, Fm_AX,
                Space_Mandrin, Space_BY_AX, StartAX, NPitchAX, N_Pitch_BY, NofInterpolPts, Smoothing)
tresse_initiale.__run__()

#fonction pour voir le temps
"""import cProfile
profiler = cProfile.Profile()
profiler.enable()

tresse_initiale.__run__()

profiler.disable()
profiler.print_stats(sort="tottime")"""

end_time = time.perf_counter()
print(f"Program time: {end_time - start_time:.3f} s")



