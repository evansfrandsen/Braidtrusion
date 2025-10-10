from abaqus import *
from abaqusConstants import *
import regionToolset
import numpy as np
import math
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
import displayGroupOdbToolset as dgo
import connectorBehavior

######################### Initialisation #######################################

DieInputFile = 'C:/temp/Dies_Assembly_Step2/PultusionSetupAssy_ForAbacus.CATProduct'

FiberCoordinatesFolder = 'C:/Temp/Fiber_Points_step2/'

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

########################Constants of the simulation############################

#Constants of the yarn
FfYBY = 2640e-6 #Kg/m  #for every Yarn
FfYAX = np.array([2640,528,2640,528,2640,528,2640,528,2640,528,2640,528])*1e-6 #Kg/m
#FfYAX = np.array([1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584,1584])*1e-6 #Kg/m
FmYBY=1520e-6 #Kg/m
FmYAX = np.array([1520,304,1520,304,1520,304,1520,304,1520,304,1520,304])*1e-6 #Kg/m
#FmYAX = np.array([912,912,912,912,912,912,912,912,912,912,912,912])*1e-6 #Kg/m
ReinfDensity = 2460 #Kg/m3
MatrixDensity = 1080 #Kg/m3
E_fiber = 72.3E9 #Pa #Young Modulus of fibers
nu_fiber=0.2 #Poisson ratio of fibers
YarnCompaction=0.7 #Yarn Compaction Ratio
N_layers=1 #0 for one fiber per Yarn
Vf=0.55 #Fiber Volume Fraction
Lfiber=0.5
NFiberperYarn = 1 #Number of fibers per Yarn
for i in range(0,N_layers): #Calculation of Fibers Number
    NFiberperLayer = ((i+1)*6)
    NFiberperYarn += NFiberperLayer

NReinf=int(round(NFiberperYarn*Vf)) #Number of Reinforcement Fibers
NPolym=int(NFiberperYarn-NReinf) #Number of polymer Fibers

#Material properties of each Yarn
rhofYBY = ReinfDensity
rhofYAX = ReinfDensity 
rhomYBY = MatrixDensity 
rhomYAX = MatrixDensity 

#Compaction of the Yarns
RpYBY = YarnCompaction
RpYAX = YarnCompaction

# Yarn fiber material areas     
SfYBY = FfYBY/rhofYBY #m2
SmYBY = FmYBY/rhomYBY #m2
SfYAX = FfYAX/rhofYAX #m2
SmYAX = FmYAX/rhomYAX #m2
# Yarn total material area
SYBY = SfYBY + SmYBY
SYAX = SfYAX + SmYAX

#Braid properties
Pitch1 = 50E-3
Nby = 24
Nax=12

#Theoretical Fiber Volume fraction
Vf_th = (Nby*SfYBY+sum(SfYAX))/(Nby*SYBY+sum(SYAX))

print('Theoretical Vf equals to '+str(100*Vf_th)+' % with the current finenesses for '+str(100*Vf)+' % expected')

#Constants of the fiber
BraidType = 1
N_Curved = 5
PitchN = 2

# Contact Properties
FricCoef_Braid=0.05
FricCoef_Pultrusion=0.3
ContactStiffness=1e11 #for linear Pressure Overclosure

###############################################################################
#                  Modelisation of Dies and Mandrel                           #
###############################################################################

#Generation of the Die and the Mandrel from Catia Part

mdb.openCatia(fileName=DieInputFile, topology=SOLID, useServer=True) #Import Geometry from CATIAV5
m=mdb.models['Model-1']

#Generation of Solid Parts
Mandrel=m.PartFromGeometryFile(bodyNum=1, combine=False, 
    dimensionality=THREE_D, geometryFile=mdb.acis, name='Mandrel', 
    scale=0.001, type=DISCRETE_RIGID_SURFACE) 
NewDie=m.PartFromGeometryFile(bodyNum=3, combine=False, 
    dimensionality=THREE_D, geometryFile=mdb.acis, name='NewStraightDie', 
    scale=0.001, type=DISCRETE_RIGID_SURFACE) 
NewMandrel=m.PartFromGeometryFile(bodyNum=2, combine=False, 
    dimensionality=THREE_D, geometryFile=mdb.acis, name='NewStraightMandrel', 
    scale=0.001, type=DISCRETE_RIGID_SURFACE) 
N_MultiDie=5 #Number of Multi-Dies
MultiDies = {} #Dictionnary containing all the Dies
for i in range(1,N_MultiDie):
    MultiDies["PartMultieDie-"+str(5-i)]=m.PartFromGeometryFile(bodyNum=9-i, combine=False, dimensionality=
        THREE_D, geometryFile=mdb.acis, name='MultiDie-'+str(5-i), scale=0.001, type=
        DISCRETE_RIGID_SURFACE)

#for the last longer Die
MultiDies["PartMultieDie-5"]=m.PartFromGeometryFile(bodyNum=4, combine=False, dimensionality=
    THREE_D, geometryFile=mdb.acis, name='MultiDie-5', scale=0.001, type=
    DISCRETE_RIGID_SURFACE)


#Mesh the Dies
for i in range(1,N_MultiDie+1): #For all Dies
    Part=MultiDies.get("PartMultieDie-"+str(i))
    faces=Part.faces[:]
    PickedRegions=(faces, )
    Part.setMeshControls(elemShape=QUAD, 
    regions=faces)
    Part.setElementType(regions=PickedRegions, elemTypes=(ElemType(elemCode=R3D4, elemLibrary=EXPLICIT), 
    ElemType(elemCode=R3D3, elemLibrary=EXPLICIT)), )
    Part.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.001)
    Part.generateMesh()

#Mesh of the mandrel with removed part
MandFaces=Mandrel.faces[:]
PickedRegions=(faces, )
Mandrel.setElementType(elemTypes=(ElemType(elemCode=R3D4, elemLibrary=EXPLICIT), 
    ElemType(elemCode=R3D3, elemLibrary=EXPLICIT)), regions=PickedRegions)
Mandrel.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.001)
Mandrel.generateMesh()

#Definition of Reference Points
Mandrel.ReferencePoint(point=(0.0, 0.0, 0.0))
for i in range(1,N_MultiDie+1):
    Part=MultiDies.get("PartMultieDie-"+str(i))
    Part.ReferencePoint(point=(0.0, 0.0, 0.12233+0.0254*(i-1)-0.435))

#Assembly of Dies and Mandrel
m.rootAssembly.DatumCsysByDefault(CARTESIAN)
MandrelAs=m.rootAssembly.Instance(dependent=ON, name='Mandrel-1', part=Mandrel)

#Set of Mandrel faces to apply Boundary conditions
MandrelFaces=MandrelAs.faces
SetMandBC = m.rootAssembly.Set(faces=MandrelFaces, name='Set-Mandrel_BC')

MultiDiesAs={}
SetMultiDiesBC={}
for i in range(1,N_MultiDie+1):
    Part=MultiDies.get("PartMultieDie-"+str(i))
    MultiDiesAs["MultiDieAs-"+str(i)]=m.rootAssembly.Instance(dependent=ON, name='MultiDie-'+str(i)+'-1', part=Part)
    #Set of Die faces to apply Boundary Conditions
    PartAs=MultiDiesAs.get("MultiDieAs-"+str(i))
    MultiDieFaces=PartAs.faces
    SetMultiDiesBC['SetMultiDieBC-'+str(i)]=m.rootAssembly.Set(faces=MultiDieFaces, name='Set-MultiDie-'+str(i)+'_BC')

#New reference point for Mandrel    
del Mandrel.features['RP']
Mandrel.ReferencePoint(point=(0.0, 0.0, 0.0))       
m.rootAssembly.regenerate()
Mandrel.regenerate()
del Mandrel.features['RP']
Mandrel.ReferencePoint(point=(0.0, 0.0, 0.09-0.435))
m.rootAssembly.regenerate() #to update assembly with the modifications of the parts

#Surfaces Definition
#Mandrel
Mandrel.Surface(name='Surf-MandrelMesh', side1Elements=Mandrel.elements[:])
#Dies
DiesSurf={}
for i in range(1,N_MultiDie+1):
    
    if i==N_MultiDie: #for the last bigger Die
        Part=MultiDies.get("PartMultieDie-"+str(i))
        DiesSurf["DieSurf-"+str(i)]=Part.Surface(name='Surf-'+str(i)+'_Mesh', 
            side2Elements=Part.elements[:])
    else : #For the 4 others Dies
        Part=MultiDies.get("PartMultieDie-"+str(i))
        DiesSurf["DieSurf-"+str(i)]=Part.Surface(name='Surf-'+str(i)+'_Mesh', 
            side2Elements=Part.elements[:])


########################### Mass Reference Points #############################

Mandrel_RP=Mandrel.Set(name='Set-Mandrel_RP', referencePoints=(Mandrel.referencePoints[6], ))
Mandrel.engineeringFeatures.PointMassInertia(
    alpha=0.0, composite=0.0, mass=0.001, name='Inertia-1', region=Mandrel_RP)
DiesRP={}
for i in range(1,N_MultiDie+1): #for the Dies
    Part=MultiDies.get("PartMultieDie-"+str(i))
    DiesRP["DieRP-"+str(i)]=Part.Set(name='Set-Die'+str(i)+'_RP', referencePoints=(
        Part.referencePoints[5], ))
    DieRP=DiesRP.get("DieRP-"+str(i))
    Part.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, 
        mass=0.001, name='Inertia-1', region=DieRP)

#############Creation of a long Die to simulate effect of Polymer###############

#Mesh of the new Die
#Select every surface
faces=NewDie.faces[:]
PickedRegions=(faces, )
NewDie.setMeshControls(elemShape=QUAD, regions=faces)
NewDie.setElementType(elemTypes=(ElemType(elemCode=R3D4, elemLibrary=STANDARD), 
    ElemType(elemCode=R3D3, elemLibrary=STANDARD)), regions=PickedRegions)
NewDie.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.001)
NewDie.generateMesh()

m.rootAssembly.regenerate()

#Assign a mass to the New Die
NewDie.ReferencePoint(point=(0.0, 0.0, 0.12233+0.0254*4-0.435))
NewDie.Set(name='SetRPMass_DieStraight', referencePoints=(NewDie.referencePoints[5], ))
NewDie.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, mass=0.001, 
    name='DieStraightMass', region=NewDie.sets['SetRPMass_DieStraight'])

########## Creation of a long Mandrel to simulate effect of Polymer ############

#Mesh the New Mandrel
faces=NewMandrel.faces[:]
PickedRegions=(faces, )
NewMandrel.setMeshControls(elemShape=QUAD, regions=faces)
NewMandrel.setElementType(elemTypes=(ElemType(elemCode=R3D4, elemLibrary=EXPLICIT), 
    ElemType(elemCode=R3D3, elemLibrary=EXPLICIT)), regions=PickedRegions)
NewMandrel.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.001)
NewMandrel.generateMesh()

#Assign a mass to the New Mandrel
NewMandrel.ReferencePoint(point=(0.0, 0.0, 0.09-0.435))
NewMandrel_RP=NewMandrel.Set(name='Set_RP_StraightMandrel', 
    referencePoints=(NewMandrel.referencePoints[5], ))
NewMandrel.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, 
    mass=0.001, name='Inertia-1', region=NewMandrel_RP) 
 
#################### Assembly New Straight Parts ###############################

NewDieAs=m.rootAssembly.Instance(dependent=ON, name='NewStraightDie-1', part=NewDie)
NewMandrelAs=m.rootAssembly.Instance(dependent=ON, name='NewStraightMandrel-1', part=NewMandrel)
m.rootAssembly.regenerate()

#To move Parts at the right place for the Braid
m.rootAssembly.translate(instanceList=('Mandrel-1', 'MultiDie-1-1', 'MultiDie-2-1', 
    'MultiDie-3-1', 'MultiDie-4-1', 'MultiDie-5-1', 'NewStraightDie-1', 'NewStraightMandrel-1'), 
    vector=(0.0, 0.0, 0.066))

#Generation of sets and surfaces from New Straight Parts
#For straight Mandrel
StraightMandSurf=m.rootAssembly.Surface(name='Surface_StraightMandrel', 
    side1Elements=NewMandrelAs.elements[:])
NewMandFaces=NewMandrelAs.faces
SetNewMandBC=m.rootAssembly.Set(faces=NewMandFaces, name = 'Set-NewMandrel_BC')
StraightMandRP=m.rootAssembly.Set(name='Set_RP_StraightMandrel_Assembly', 
    referencePoints=(NewMandrelAs.referencePoints[5], ))

#For Straight Die
NewDieFaces=NewDieAs.faces
SetNewDieBC=m.rootAssembly.Set(faces=NewDieFaces, name = 'Set-NewDie_BC')
StraightDieRP=m.rootAssembly.Set(name='Set_RP_StraightDie_Assembly', 
    referencePoints=(NewDieAs.referencePoints[5], ))
StraightDieSurf=m.rootAssembly.Surface(name='Surface_StraightDie', 
    side2Elements=NewDieAs.elements[:])    


###############################################################################
#                       Modelisation of the Braid                             #
###############################################################################

#################Choice of Polymer and Reinforcement Fibers####################

index_Reinf=[] #indexes of Reinforcement Fibers
index_Polym=[] #indexes of Polymer Fibers
if NReinf<NPolym:
    for i in range(1,NReinf+1):
        index_Reinf.append(2*i) #pairs indexes until reaching the target number of Reinforcement fibers
    for j in range(1,NFiberperYarn+1):
        if j not in index_Reinf:
            index_Polym.append(j) #the others indexes are Polymer fibers
else :
    for i in range(1,NPolym+1):
        index_Polym.append(2*i) #pairs indexes until reaching the target number of Polymer fibers
    for j in range(1,NFiberperYarn+1):
        if j not in index_Polym:
            index_Reinf.append(j) #the others indexes are Reinforcement fibers

index_Reinf_AX=[] #indexes of Reinforcement AX Fibers
index_Reinf_BY=[] #indexes of Reinforcement BY Fibers

for j in range(Nax): #for each Braiding Yarn
    for l in index_Reinf: #for Reinforcement Fibers
        k=j*NFiberperYarn+l #link with k above
        index_Reinf_AX.append(k)
for j in range(Nby): #for each Braiding Yarn
    for l in index_Reinf: #for Reinforcement Fibers
        k=j*NFiberperYarn+l #link with k above
        index_Reinf_BY.append(k)

#Get Node coordinates from the .txt files

#for Braiding Yarns
xyzBY=[] #coordinates of each node of each BY fiber
for i in  index_Reinf_BY :
    FiberCoords = open(FiberCoordinatesFolder+"it_0_BY_fiber_"+str(i-1)+".txt", "r") #Open the file
    Coords=[line.split(' ') for line in FiberCoords.read().splitlines()] #Get values inside the file
    Coords=np.array(Coords,dtype=float) #Convert strings into floats
    xyzBY.append(Coords)
    FiberCoords.close()

#for Axial Yarns
xyzAx=[] #coordinates of each node of each Ax fiber
for i in index_Reinf_AX :
    FiberCoords = open(FiberCoordinatesFolder+"it_0_AX_fiber_"+str(i-1)+".txt", "r") #Open the file
    Coords=[line.split(' ') for line in FiberCoords.read().splitlines()] #Get values inside the file
    Coords=np.array(Coords,dtype=float)
    xyzAx.append(Coords) #to convert from mm to m
    FiberCoords.close()


#Create the Parts for each BY fiber
FibersBY = {} #Dictionnary containing all the BY Fibers			
for k in range(1,len(xyzBY)+1):
	FibersBY["PartFiberBY-"+str(k)] = m.Part(name='fiberBY-'+str(k), dimensionality=THREE_D, type=DEFORMABLE_BODY)
	xxyyzzBY=xyzBY[k-1] #coordinates of each section of fiber
	fiberModel=FibersBY.get("PartFiberBY-"+str(k))
	for j in range(1,len(xxyyzzBY)-1): #for each node
		dtmBY = fiberModel.DatumPointByCoordinate(coords=tuple(xxyyzzBY[j]))
	d = fiberModel.datums
    #plot a line to create the fiber
	fiberModel.WirePolyLine(mergeType=IMPRINT, meshable=ON, points=d.values())	

#Create the Parts for each Ax fiber	
FibersAx = {} #Dictionnary containing all the Ax Fibers									
for k in range(1,len(xyzAx)+1):
	FibersAx["PartFiberAx-"+str(k)] = m.Part(name='fiberAx-'+str(k), dimensionality=THREE_D, type=DEFORMABLE_BODY)
	xxyyzzAx=xyzAx[k-1]
	fiberModel=FibersAx.get("PartFiberAx-"+str(k))	
	for j in range(1,len(xxyyzzAx)-1): #for each node
		dtmAx = fiberModel.DatumPointByCoordinate(coords=tuple(xxyyzzAx[j]))
	d = fiberModel.datums
    #plot a line to create the fiber
	fiberModel.WirePolyLine(mergeType=IMPRINT, meshable=ON, points=d.values())

###########################Creation of a Dummy Node############################
DummyNode=m.Part(name='DummyNode')
DummyNode.Node(coordinates=(0.0,0.0,-0.1), label=999)
DummyNode.SetFromNodeLabels(name='DummyNodeSetPart', nodeLabels=(999,))
m.rootAssembly.Instance(dependent=ON, name='DummyNode-1', part=DummyNode)
DummyNodeSet=m.rootAssembly.SetFromNodeLabels(name='DummyNodeSet', nodeLabels=(('DummyNode-1', (999, )), ))
#Create the dummy node for BC
DummynodeAX = {} #Dictionnary containing all the Ax dummy node
DummynodeBY = {} #Dictionnary containing all the BY dummy node										
for k in range(1,len(xyzAx)+1):
    xxyyzzAx=xyzAx[k-1]
    DummynodeAX["DummynodeAx1-"+str(k)]=m.Part(name='DummyNodeAX1-'+str(k))
    dummynode1=DummynodeAX.get("DummynodeAx1-"+str(k))
    dummynode1.Node(coordinates=tuple(xxyyzzAx[0]), label=1000+k)
    Dummy =dummynode1.SetFromNodeLabels(name='DummyNodeSetPartAX1-'+str(k), nodeLabels=(1000+k,))
    m.rootAssembly.Instance(dependent=ON, name='DummyNodeSetPartAX1-1-'+str(k), part=dummynode1)
    #m.rootAssembly.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, mass=1e-10, i11=1e-10, i22=1e-10, i33=1e-10, name='mass_dummynodeAX1'+str(k), region=Dummy)

    DummynodeAX["DummynodeAx2-"+str(k)]=m.Part(name='DummyNodeAX2-'+str(k))
    dummynode2=DummynodeAX.get("DummynodeAx2-"+str(k))
    dummynode2.Node(coordinates=tuple(xxyyzzAx[-1]), label=2000+k)
    Dummy =dummynode2.SetFromNodeLabels(name='DummyNodeSetPartAX2-'+str(k), nodeLabels=(2000+k,))
    m.rootAssembly.Instance(dependent=ON, name='DummyNodeSetPartAX2-1-'+str(k), part=dummynode2)
    #m.rootAssembly.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, mass=1e-10, i11=1e-10, i22=1e-10, i33=1e-10, name='mass_dummynodeAX2'+str(k), region=Dummy)
    
for k in range(1,len(xyzBY)+1):
    xxyyzzBY=xyzBY[k-1]
    DummynodeBY["DummynodeBY1-"+str(k)]=m.Part(name='DummyNodeBY1-'+str(k))
    dummynode1=DummynodeBY.get("DummynodeBY1-"+str(k))
    dummynode1.Node(coordinates=tuple(xxyyzzBY[0]), label=3000+k)
    Dummy =dummynode1.SetFromNodeLabels(name='DummyNodeSetPartBY1-'+str(k), nodeLabels=(3000+k,))
    m.rootAssembly.Instance(dependent=ON, name='DummyNodeSetPartBY1-1-'+str(k), part=dummynode1)
    #m.rootAssembly.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, mass=1e-10, i11=1e-10, i22=1e-10, i33=1e-10, name='mass_dummynodeBY1'+str(k), region=Dummy)

    DummynodeBY["DummynodeBY2-"+str(k)]=m.Part(name='DummyNodeBY2-'+str(k))
    dummynode2=DummynodeBY.get("DummynodeBY2-"+str(k))
    dummynode2.Node(coordinates=tuple(xxyyzzBY[-1]), label=4000+k)
    Dummy =dummynode2.SetFromNodeLabels(name='DummyNodeSetPartBY2-'+str(k), nodeLabels=(4000+k,))
    m.rootAssembly.Instance(dependent=ON, name='DummyNodeSetPartBY2-1-'+str(k), part=dummynode2)
    #m.rootAssembly.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, mass=1e-10, i11=1e-10, i22=1e-10, i33=1e-10, name='mass_dummynodeBY2'+str(k), region=Dummy)

DummyNode_BY_1 = []
DummyNode_Ax_1 = []
DummyNode_BY_2 = []
DummyNode_Ax_2 = []
for k in range(1,len(xyzBY)+1):
    DummyNode_BY_1.append( ('DummyNodeSetPartBY1-1-'+str(k), (3000+k,)) )
for k in range(1,len(xyzBY)+1):
    DummyNode_BY_2.append( ('DummyNodeSetPartBY2-1-'+str(k), (4000+k,)) )
for k in range(1,len(xyzAx)+1):
    DummyNode_Ax_1.append( ('DummyNodeSetPartAX1-1-'+str(k), (1000+k,)) )
for k in range(1,len(xyzAx)+1):
    DummyNode_Ax_2.append( ('DummyNodeSetPartAX2-1-'+str(k), (2000+k,)) )

DummyNodeBYSet=m.rootAssembly.SetFromNodeLabels(name='DummyNodeSetBY', nodeLabels=DummyNode_BY_1)
DummyNodeBYSet2=m.rootAssembly.SetFromNodeLabels(name='DummyNodeSetBY2', nodeLabels=DummyNode_BY_2)
DummyNodeAXSet=m.rootAssembly.SetFromNodeLabels(name='DummyNodeSetAX', nodeLabels=DummyNode_Ax_1)
DummyNodeAXSet2=m.rootAssembly.SetFromNodeLabels(name='DummyNodeSetAX2', nodeLabels=DummyNode_Ax_2)

DummyNodeforce_all=m.rootAssembly.SetByBoolean(name='DummyNodeforce_all', sets=(DummyNodeBYSet, DummyNodeAXSet))
DummyNodeBC_all=m.rootAssembly.SetByBoolean(name='DummyNodeBC_all', sets=(DummyNodeBYSet2, DummyNodeAXSet2))

##############################Material properties##############################

# Create Fiber material
FiberMaterial=m.Material(name='Material-1')
FiberMaterial.Density(table=((ReinfDensity, ), ))
FiberMaterial.Elastic(table=((E_fiber,nu_fiber), ))

#Braiding Yarns section
#for Reinforcement fibers
FiberRadiusBYReinf=np.sqrt(SfYBY/(NReinf*np.pi)) #Radius of each Fiber
ABYReinf=np.pi*FiberRadiusBYReinf**2
#Section of Fiber
m.TrussSection(area=ABYReinf, material='Material-1', name='fiberSectionBY_Reinf')


#Axial Yarns section
#for Reinforcement fibers
FiberRadiusAxReinf=np.sqrt(SfYAX/(NReinf*np.pi)) #Radius of each Fiber
AAXReinf=np.pi*FiberRadiusAxReinf**2
#Section of Fiber
for i in range(len(AAXReinf)):
    m.TrussSection(area=AAXReinf[i], material='Material-1', name='fiberSectionAx_Reinf_'+str(i))

###############################Mesh of BY fibers###############################

for k in range(1,len(xyzBY)+1):
    FiberBY=FibersBY.get("PartFiberBY-"+str(k))
    FiberBY.setElementType(elemTypes=(ElemType(elemCode=T3D2, elemLibrary=EXPLICIT), ), 
    regions=(FiberBY.edges[:], ))
    FiberBY.seedEdgeByNumber(constraint=FINER, edges=FiberBY.edges[:], number=1)
    FiberBY.generateMesh()

#Mesh of Ax fibers
for k in range(1,len(xyzAx)+1):
    FiberAx=FibersAx.get("PartFiberAx-"+str(k))	
    FiberAx.setElementType(elemTypes=(ElemType(elemCode=T3D2, elemLibrary=EXPLICIT), ), 
        regions=(FiberAx.edges[:], ))
    FiberAx.seedEdgeByNumber(constraint=FINER, edges=FiberAx.edges[:], number=1)
    FiberAx.generateMesh()

#############Creation of an Orphan Mesh, copies of the Fibers###################
#Dictionaries with Copies
M_FibersBY={}
M_FibersAx={}

#Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    FiberBY=FibersBY.get("PartFiberBY-"+str(k))
    M_FibersBY["M_PartFiberBY-"+str(k)]=FiberBY.PartFromMesh(copySets=True, name='M_fiberBY-'+str(k))
    
#Axial Yarn Fibers
for k in range(1,len(xyzAx)+1):
    FiberAx=FibersAx.get("PartFiberAx-"+str(k))
    M_FibersAx["M_PartFiberAx-"+str(k)]=FiberAx.PartFromMesh(copySets=True, name='M_fiberAx-'+str(k))

#Renumber nodes
mdb.meshEditOptions.setValues(enableUndo=True, maxUndoCacheElements=0.5)
#for Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    M_FiberBY=M_FibersBY.get("M_PartFiberBY-"+str(k))
    M_FiberBY.renumberNode(increment=1, nodes=(M_FiberBY.nodes[1], 
		M_FiberBY.nodes[0]), startLabel=1)
#for Axial Yarn fibers
for k in range(1,len(xyzAx)+1):
    M_FiberAx=M_FibersAx.get("M_PartFiberAx-"+str(k))
    M_FiberAx.renumberNode(increment=1, nodes=(M_FiberAx.nodes[1], 
		M_FiberAx.nodes[0]), startLabel=1)
	
# Assign Section to the Orphan Mesh
#for Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    M_FiberBY=M_FibersBY.get("M_PartFiberBY-"+str(k))
    M_FiberBY.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(elements=M_FiberBY.elements[:]), sectionName='fiberSectionBY_Reinf', thicknessAssignment=FROM_SECTION)
    
#for Axial Yarn fibers
for k in range(1,len(xyzAx)+1):
    M_FiberAx=M_FibersAx.get("M_PartFiberAx-"+str(k))
    j = (k-1)/(len(index_Reinf))
    M_FiberAx.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(elements=M_FiberAx.elements[:]), sectionName='fiberSectionAx_Reinf_'+str(j), thicknessAssignment=FROM_SECTION)
		
# Assign Orientation to the Fibers
#for Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    M_FiberBY=M_FibersBY.get("M_PartFiberBY-"+str(k))
    M_FiberBY.assignBeamSectionOrientation(method=N1_COSINES, n1=(0.0, 0.0, -1.0), 
        region=Region(elements=M_FiberBY.elements[:]))
#for Axial Yarn fibers
for k in range(1,len(xyzAx)+1):
    M_FiberAx=M_FibersAx.get("M_PartFiberAx-"+str(k))
    M_FiberAx.assignBeamSectionOrientation(method=N1_COSINES, n1=(0.0, 0.0, -1.0), 
        region=Region(elements=M_FiberAx.elements[:]))

##############################Assembly of fibers################################

As_M_FibersBY={}
As_M_FibersAx={}

#for Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    m.rootAssembly.DatumCsysByDefault(CARTESIAN)
    M_FiberBY=M_FibersBY.get("M_PartFiberBY-"+str(k))
    As_M_FibersBY["As_M_FiberBY-"+str(k)]=m.rootAssembly.Instance(dependent=ON, 
        name='M_fiberBY-'+str(k), part=M_FiberBY)
    
#for Axial Yarn fibers
for k in range(1,len(xyzAx)+1):
    m.rootAssembly.DatumCsysByDefault(CARTESIAN)
    M_FiberAx=M_FibersAx.get("M_PartFiberAx-"+str(k))
    As_M_FibersAx["As_M_FiberAx-"+str(k)]=m.rootAssembly.Instance(dependent=ON, 
        name='M_fiberAx-'+str(k), part=M_FiberAx)

SurfacesBY={}
SurfacesAx={}

#Get the Surfaces of every Fiber
#for Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    M_FiberBY=M_FibersBY.get("M_PartFiberBY-"+str(k))
    #Surface of the whole fiber
    SurfacesBY['SurfaceBY'+str(k)]=M_FiberBY.Surface(circumElements=M_FiberBY.elements[:], 
        name='SurfBY-'+str(k))
    #Surface at the end of the fiber
    M_FiberBY.Surface(end2Elements=M_FiberBY.elements[1:3], name='SurfEndBYEl1-'+str(k))
#for Axial Yarn fibers
for k in range(1,len(xyzAx)+1):
    M_FiberAx=M_FibersAx.get("M_PartFiberAx-"+str(k))
    #Surface of the whole fiber
    SurfacesAx['SurfaceAx'+str(k)]=M_FiberAx.Surface(circumElements=M_FiberAx.elements[:], 
        name='SurfAx-'+str(k))
    #Surface at the end of the fiber
    M_FiberAx.Surface(end2Elements=M_FiberAx.elements[1:3], name='SurfEndAxEl1-'+str(k))

###############################################################################
#                         Build Abaqus simulations                            #
###############################################################################


#############################Creation of steps##################################
#Second Step for pultrusion beginning
Step1=m.ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=0.2)
m.SmoothStepAmplitude(data=((0.0, 0.0), (0.2,1.0)), name='Amp-Die', timeSpan=STEP)


############################Create contact properties###########################

#In Step 1
IntProp1=m.ContactProperty('IntProp-1')
IntProp1.TangentialBehavior(dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, table=((FricCoef_Pultrusion, ), ), temperatureDependency=OFF)
IntProp1.NormalBehavior(pressureOverclosure=LINEAR, contactStiffness=ContactStiffness, constraintEnforcementMethod=DEFAULT)
	
Int1=m.ContactExp(createStepName='Step-1', name='Int-1')
Int1.includedPairs.setValuesInStep(stepName='Step-1', useAllstar=ON)
Int1.contactPropertyAssignments.appendInStep(assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Step-1')	

# Smooth step
m.SmoothStepAmplitude(data=((0.0, 0.0), (0.2, 1.0)), name='Amp-1', timeSpan=STEP)

####################Create Sets to assign BCs and loads#########################

#For Boundary Conditions
nodesBC_BY = []
nodesBC_Ax = []
#for Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    Node=len(xyzBY[k-1])-2 #take last node of each fiber
    nodesBC_BY.append( ('M_fiberBY-'+str(k), (Node, )) )
#for Axial Yarn fibers
for k in range(1,len(xyzAx)+1):
    Node=len(xyzAx[k-1])-2 #take last node of each fiber
    nodesBC_Ax.append( ('M_fiberAx-'+str(k), (Node, )) )

BC_BY=m.rootAssembly.SetFromNodeLabels(nodeLabels=nodesBC_BY, name='BoundaryConditionsBY')
BC_Ax=m.rootAssembly.SetFromNodeLabels(nodeLabels=nodesBC_Ax, name='BoundaryConditionsAx')

#For loads
nodesLoad_BY = []
nodesLoad_Ax = []
#for Braiding Yarn fibers
for k in range(1,len(xyzBY)+1):
    #take first Node of each fiber
    nodesLoad_BY.append( ('M_fiberBY-'+str(k), (1, )) )
#for Axial Yarn fibers
for k in range(1,len(xyzAx)+1):
    nodesLoad_Ax.append( ('M_fiberAx-'+str(k), (1, )) )
Load_BY=m.rootAssembly.SetFromNodeLabels(nodeLabels=nodesLoad_BY, name='ForceBY')
Load_Ax=m.rootAssembly.SetFromNodeLabels(nodeLabels=nodesLoad_Ax, name='ForceAx')

BCSet=m.rootAssembly.SetByBoolean(name='BoundaryConditionsSET', sets=(BC_BY, BC_Ax))
ForceSet=m.rootAssembly.SetByBoolean(name='ForceSET', sets=(Load_BY, Load_Ax))


############################# History Outputs #################################

#Global outputs
m.historyOutputRequests['H-Output-1'].setValues(variables=
    ('ALLDC', 'ALLFD', 'ALLIE', 'ALLKE', 'ALLSE', 'ALLVD', 'ALLWK', 'ALLCW', 'ALLPW', 'ETOTAL'))
m.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'LE', 'UT', 'CSTRESS', 'CFORCE', 'COORD'))

#Forces and Area of Contact Outputs
m.rootAssembly.regenerate() 

#For mandrel
m.HistoryOutputRequest(createStepName='Step-1', filter=ANTIALIASING, 
    name='H-Output-ContactMandrel', numIntervals=200, rebar=EXCLUDE, 
    region=MandrelAs.surfaces['Surf-MandrelMesh'], sectionPoints=DEFAULT, 
    variables=('CFNM', 'CFSM', 'CFTM','CAREA'))

#For Dies
for i in range(1,N_MultiDie+1):
    DieSurf=m.rootAssembly.allInstances['MultiDie-'+str(i)+'-1'].surfaces['Surf-'+str(i)+'_Mesh']
    m.HistoryOutputRequest(createStepName='Step-1', filter=
        ANTIALIASING, name='H-Output-ContactDie'+str(i), rebar=EXCLUDE, region=DieSurf, 
        sectionPoints=DEFAULT, variables=('CFNM', 'CFSM', 'CFTM', 'CAREA'))

#For Straight Die and Mandrel
m.HistoryOutputRequest(createStepName='Step-1', filter=ANTIALIASING, 
    name='H-Output-StraightDie', rebar=EXCLUDE, region=StraightDieSurf, 
    sectionPoints=DEFAULT, variables=('CFNM', 'CFSM', 'CFTM', 'CAREA')) 
m.HistoryOutputRequest(createStepName='Step-1', filter=ANTIALIASING, 
    name='H-Output-StraightMandrel', rebar=EXCLUDE, region=StraightMandSurf, 
    sectionPoints=DEFAULT, variables=('CFNM', 'CFSM', 'CFTM', 'CAREA'))

#Get Reaction Forces in direction 3
#For mandrel
m.HistoryOutputRequest(createStepName='Step-1', filter=ANTIALIASING, 
    name='H-Output-Mandrel_RF', rebar=EXCLUDE, region=m.rootAssembly.sets['Mandrel-1.Set-Mandrel_RP'], 
    sectionPoints=DEFAULT, variables=('RF3', ))

#For Dies
for i in range(1,N_MultiDie+1):
    DieRP=m.rootAssembly.sets['MultiDie-'+str(i)+'-1.Set-Die'+str(i)+'_RP']
    m.HistoryOutputRequest(createStepName='Step-1', filter=
        ANTIALIASING, name=
        'H-Output-Die'+str(i)+'_RF', rebar=EXCLUDE, region=DieRP, 
        sectionPoints=DEFAULT, variables=('RF3', ))

#For Straight Mandrel and Die
m.HistoryOutputRequest(createStepName='Step-1', name='H-Output-StraightMandrel_RF', 
    rebar=EXCLUDE, region=StraightMandRP, filter=ANTIALIASING, 
    sectionPoints=DEFAULT, variables=('RF3', ))
m.HistoryOutputRequest(createStepName='Step-1', name='H-Output-StraightDie_RF', 
    rebar=EXCLUDE, region=StraightDieRP, filter=ANTIALIASING, 
    sectionPoints=DEFAULT, variables=('RF3', ))
################### Multi-point Constraints btw set of DOF ####################


Name_connector = 'Spring'
Connector = m.ConnectorSection(name= Name_connector, translationalType=AXIAL)
elastic_0 = connectorBehavior.ConnectorElasticity(components=(1, ), table=((100.0, ), ))
m.sections[Name_connector].setValues(behaviorOptions =(elastic_0, ) )
m.sections[Name_connector].behaviorOptions[0].ConnectorOptions()


datum1 = mdb.models['Model-1'].rootAssembly.datums[1]
for i in range(1,len(xyzAx)+1):
    n1 = m.rootAssembly.instances['M_fiberAx-'+str(i)].nodes
    r11 = m.rootAssembly.instances["DummyNodeSetPartAX2-1-"+str(i)].nodes
    Wire = m.rootAssembly.WirePolyLine(points=((n1[1], r11[0]), ), mergeType=IMPRINT, meshable=OFF)
    edges1 = mdb.models['Model-1'].rootAssembly.edges.findAt((n1[1].coordinates, ))
    m.rootAssembly.Set(edges=edges1, name='WireAXR-'+str(i)+'-Set-1')
    region=m.rootAssembly.sets['WireAXR-'+str(i)+'-Set-1']
    csa = m.rootAssembly.SectionAssignment(sectionName=Name_connector, region=region)
    #m.rootAssembly.ConnectorOrientation(orient2sameAs1=OFF, region=csa.getSet(), localCsys1=datum1)

for i in range(1,len(xyzBY)+1):
    n1 = m.rootAssembly.instances['M_fiberBY-'+str(i)].nodes
    r11 = m.rootAssembly.instances["DummyNodeSetPartBY2-1-"+str(i)].nodes
    Wire = m.rootAssembly.WirePolyLine(points=((n1[1], r11[0]), ), mergeType=IMPRINT, meshable=OFF)
    edges1 = mdb.models['Model-1'].rootAssembly.edges.findAt((n1[1].coordinates, ))
    m.rootAssembly.Set(edges=edges1, name='WireBYR-'+str(i)+'-Set-1')
    region=m.rootAssembly.sets['WireBYR-'+str(i)+'-Set-1']
    csa = m.rootAssembly.SectionAssignment(sectionName=Name_connector, region=region)
    #m.rootAssembly.ConnectorOrientation(orient2sameAs1=OFF, region=csa.getSet(), localCsys1=datum1)

for i in range(1,len(xyzAx)+1):
    n1 = m.rootAssembly.instances['M_fiberAx-'+str(i)].nodes
    r11 = m.rootAssembly.instances["DummyNodeSetPartAX1-1-"+str(i)].nodes
    Wire = m.rootAssembly.WirePolyLine(points=((n1[-1], r11[0]), ), mergeType=IMPRINT, meshable=OFF)
    edges1 = mdb.models['Model-1'].rootAssembly.edges.findAt((n1[-1].coordinates, ))
    m.rootAssembly.Set(edges=edges1, name='WireAXL-'+str(i)+'-Set-1')
    region=m.rootAssembly.sets['WireAXL-'+str(i)+'-Set-1']
    csa = m.rootAssembly.SectionAssignment(sectionName=Name_connector, region=region)
    #m.rootAssembly.ConnectorOrientation(orient2sameAs1=OFF, region=csa.getSet(), localCsys1=datum1)

for i in range(1,len(xyzBY)+1):
    n1 = m.rootAssembly.instances['M_fiberBY-'+str(i)].nodes
    r11 = m.rootAssembly.instances["DummyNodeSetPartBY1-1-"+str(i)].nodes
    Wire = m.rootAssembly.WirePolyLine(points=((n1[-1], r11[0]), ), mergeType=IMPRINT, meshable=OFF)
    edges1 = mdb.models['Model-1'].rootAssembly.edges.findAt((n1[-1].coordinates, ))
    m.rootAssembly.Set(edges=edges1, name='WireBYL-'+str(i)+'-Set-1')
    region=m.rootAssembly.sets['WireBYL-'+str(i)+'-Set-1']
    csa = m.rootAssembly.SectionAssignment(sectionName=Name_connector, region=region)
    #m.rootAssembly.ConnectorOrientation(orient2sameAs1=OFF, region=csa.getSet(), localCsys1=datum1)

BC_all=m.rootAssembly.SetByBoolean(name='BC_all', sets=(DummyNodeBC_all, ForceSet))
m.Equation(name='Equation2L', terms=( (-1.0, 'BC_all', 3), (1.0, 'DummyNodeSet', 3), ))
############################# Loads on Fibers #################################
#Force on Fibers on Step 1
LoadStp1_ReinfR=m.ConcentratedForce(amplitude='Amp-Die', cf3=0.4, createStepName='Step-1', distributionType=UNIFORM, field='', localCsys=None, name='Load_stp1_ReinfR', region=BCSet)

#Move of Reinforcement fibers to the Die
BraidBCs_Stp1=m.DisplacementBC(amplitude='Amp-Die', createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BraidBCsL_stp1', region=DummyNodeforce_all, u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
BraidBCs_Stp1=m.DisplacementBC(amplitude='Amp-Die', createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BraidBCsL2_stp1', region=ForceSet, u1=UNSET, u2=UNSET, u3=-580E-3, ur1=UNSET, ur2=UNSET, ur3=UNSET)
BraidBCs_Stp1=m.DisplacementBC(amplitude='Amp-Die', createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BraidBCsL3_stp1', region=DummyNodeBC_all, u1=0.0, u2=0.0, u3=-580E-3, ur1=UNSET, ur2=UNSET, ur3=UNSET)

BraidBCs_Stp1=m.DisplacementBC(amplitude='Amp-Die', createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BCsPBC', region=DummyNodeSet, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

################# Boundary Conditions on Dies and Mandrels ####################

#Reference Points of the Dies fixed
for i in range(1,N_MultiDie+1): 
    DieRP=m.rootAssembly.sets['MultiDie-'+str(i)+'-1.Set-Die'+str(i)+'_RP']
    m.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, 
        fieldName='', fixed=OFF, localCsys=None, name='Die'+str(i)+'_RP', region=DieRP, 
        u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)

#Reference point of the Mandrel fixed
Mandrel_RP=m.rootAssembly.sets['Mandrel-1.Set-Mandrel_RP']
m.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Mandrel_RP', 
    region=Mandrel_RP, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)

#Mandrel fixed in both Steps
m.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, 
    fieldName='', fixed=OFF, localCsys=None, name='Mandrel', region=SetMandBC, 
    u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)


#No displacements admitted for reference points of Straight Parts
m.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, 
    fieldName='', fixed=OFF, localCsys=None, name='StraightMandrel_RP', 
    region=StraightMandRP, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)
m.DisplacementBC(amplitude=UNSET, createStepName='Step-1', distributionType=UNIFORM, 
    fieldName='', fixed=OFF, localCsys=None, name='StraightDie_RP', 
    region=StraightDieRP, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)

