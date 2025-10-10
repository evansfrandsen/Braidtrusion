from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior


session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

N_layers=1 #0 for one fiber per Yarn
NFiberperYarn = 1 #Number of fibers per Yarn
for i in range(0,N_layers): #Calculation of Fibers Number
    NFiberperLayer = ((i+1)*6)
    NFiberperYarn += NFiberperLayer
Nby = 24
Nax = 12

Test = 'ActualTestEvans' #Ligne a modifier
odb = session.odbs['C:/Temp/ActualTestEvans.odb'] #Ligne a modifier


input_filename = 'C:/temp/FieldOutputReport/FieldOutput_'+ Test +'.txt'

def Extract_Lines(input_file, output_file, start_marker, end_marker):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        copy = False
        skip_lines = 4
        for line in infile:
            if start_marker in line:
                copy = True
                skip_lines = 4
                continue
            if end_marker in line:
                break

            if copy :
                if skip_lines > 0 :
                    skip_lines -= 1
                    continue
                line2 = line.split()
                if line2 == []:
                    break
                outfile.write(str(line2[1])+ ' ' + str(line2[2]) + ' ' + str(line2[3])+ "\n")

def Extract_Yarn_coord(Type_Yarn, N_Yarn, N_Fiber_per_Yarn):

    for i in range (N_Yarn*N_Fiber_per_Yarn):
        output_filename = "C:/temp/Fiber_Points_step2/it_0_"+ Type_Yarn +"_fiber_"+ str(i) +".txt" 
        start_text = "Field Output reported at nodes for part: M_FIBER"+ Type_Yarn +"-" +str(i+1)
        end_text = "Field Output reported at nodes for part: M_FIBER"+ Type_Yarn +"-" +str(i+2)

        Extract_Lines(input_filename, output_filename, start_text, end_text)

odbName=session.viewports[session.currentViewportName].odbDisplay.name
session.odbData[odbName].setValues(activeFrames=(('Step-1', (20, )), ))


session.fieldReportOptions.setValues(reportFormat=NORMAL_ANNOTATED, separateTables= OFF)
session.writeFieldReport(fileName= input_filename, append=OFF, sortItem='Node Label', odb=odb, step=0, frame=20, outputPosition=NODAL,  variable=(('COORD', NODAL, ((COMPONENT, 'COOR1'), (COMPONENT, 'COOR2'), (COMPONENT, 'COOR3'), )), ), stepFrame=SPECIFY)

Extract_Yarn_coord(Type_Yarn='AX', N_Yarn=Nax, N_Fiber_per_Yarn=NFiberperYarn)
Extract_Yarn_coord(Type_Yarn='BY', N_Yarn=Nby, N_Fiber_per_Yarn=NFiberperYarn)
