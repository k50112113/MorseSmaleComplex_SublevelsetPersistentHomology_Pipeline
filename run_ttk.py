import inspect
import sys
import os
import vtk 
from vtk import (
    vtkDataObject,
    vtkTableWriter,
    vtkThreshold,
    #vtkXMLPolyDataWriter,
    #vtkXMLUnstructuredGridReader,
    #vtkXMLStructuredGridReader,
    #vtkXMLUnstructuredGridWriter,
    vtkXMLImageDataReader,
    vtkDelimitedTextWriter,
    vtkDataObjectToTable,
    vtkPolyData,
    vtkTable
)
from topologytoolkit import (
    ttkMorseSmaleComplex,
    ttkPersistenceCurve,
    ttkPersistenceDiagram,
    ttkTopologicalSimplification,
    ttkScalarFieldSmoother,
)

# python run_ttk.py butane_2d_log_prob_example.vti 1 1 0 results

if len(sys.argv) == 6:
    input_file_path = sys.argv[1]
    output_folder_path = sys.argv[2]
    persistent_theshold = float(sys.argv[3])
    smooth_itr = int(sys.argv[4])
    connector = int(sys.argv[5])
    os.system(f"mkdir -p {output_folder_path}")
else:
    print("Invalid arguments")
    exit()

# 1. loading the input data
#reader = vtkXMLUnstructuredGridReader()
reader = vtkXMLImageDataReader()
reader.SetFileName(input_file_path)
#reader.Update()

smooth = ttkScalarFieldSmoother()
smooth.SetInputConnection(reader.GetOutputPort())
smooth.SetInputArrayToProcess(0, 0, 0, 0, "scalars")
smooth.SetNumberOfIterations(smooth_itr)
smooth.SetDebugLevel(3)
smooth.Update()

# 2. computing the persistence curve
curve = ttkPersistenceCurve()
curve.SetInputConnection(smooth.GetOutputPort())
curve.SetInputArrayToProcess(0, 0, 0, 0, "scalars")
curve.SetDebugLevel(3)
curve.Update()

outputpc = curve.GetOutputDataObject(0)
table0 = vtkDataObjectToTable()
table0.SetInputConnection(curve.GetOutputPort(0))
table0.Update()
table0.GetOutput().AddColumn(outputpc.GetColumn(0))
table0.GetOutput().AddColumn(outputpc.GetColumn(1))
table0.Update()
writer1 = vtkDelimitedTextWriter()
writer1.SetFileName(f"{output_folder_path}/pc.txt")
writer1.SetInputConnection(table0.GetOutputPort())
writer1.Update()
writer1.Write()

# 3. computing the persitence diagram
diagram = ttkPersistenceDiagram()
diagram.SetInputConnection(smooth.GetOutputPort())
diagram.SetInputArrayToProcess(0, 0, 0, 0, "scalars")
diagram.SetDebugLevel(3)
#diagram.Update()

# 4. selecting the critical point pairs
criticalPairs = vtkThreshold()
criticalPairs.SetInputConnection(diagram.GetOutputPort())
criticalPairs.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "PairIdentifier")
criticalPairs.ThresholdBetween(0.0, 999999)
#criticalPairs.Update()

outputcp = criticalPairs.GetOutputDataObject(0)
table0 = vtkDataObjectToTable()
table0.SetInputConnection(criticalPairs.GetOutputPort(0))
table0.Update()
table0.GetOutput().AddColumn(outputcp.GetPoints().GetData())
table0.Update()
writer1 = vtkDelimitedTextWriter()
writer1.SetFileName(f"{output_folder_path}/cp.txt")
writer1.SetInputConnection(table0.GetOutputPort())
writer1.Update()
writer1.Write()

# 5. selecting the most persistent pairs
persistentPairs = vtkThreshold()
persistentPairs.SetInputConnection(criticalPairs.GetOutputPort())
persistentPairs.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "Persistence")
persistentPairs.ThresholdBetween(persistent_theshold, 999999)
# persistentPairs.Update()

outputpp = persistentPairs.GetOutputDataObject(0)
table0 = vtkDataObjectToTable()
table0.SetInputConnection(persistentPairs.GetOutputPort(0))
table0.Update()
table0.GetOutput().AddColumn(outputpp.GetPoints().GetData())
table0.Update()
writer1 = vtkDelimitedTextWriter()
writer1.SetFileName(f"{output_folder_path}/pp.txt")
writer1.SetInputConnection(table0.GetOutputPort())
writer1.Update()
writer1.Write()

# compute MS complex without filtration
morseSmaleComplex = ttkMorseSmaleComplex()
morseSmaleComplex.SetInputConnection(smooth.GetOutputPort())
morseSmaleComplex.SetInputArrayToProcess(0, 0, 0, 0, "scalars")
morseSmaleComplex.SetDebugLevel(3)
morseSmaleComplex.Update()

output0 = morseSmaleComplex.GetOutputDataObject(0)
table0 = vtkDataObjectToTable()
table0.SetInputConnection(morseSmaleComplex.GetOutputPort(0))
table0.Update()
table0.GetOutput().AddColumn(output0.GetPoints().GetData())
table0.Update()
writer1 = vtkDelimitedTextWriter()
writer1.SetFileName(f"{output_folder_path}/nofilter-critical-points.txt")
writer1.SetInputConnection(table0.GetOutputPort())
writer1.Update()
writer1.Write()

output1 = morseSmaleComplex.GetOutputDataObject(1)
table1 = vtkDataObjectToTable()
table1.SetInputConnection(morseSmaleComplex.GetOutputPort(1))
table1.Update()
table1.GetOutput().AddColumn(output1.GetPoints().GetData())
table1.Update()
writer1 = vtkDelimitedTextWriter()
writer1.SetFileName(f"{output_folder_path}/nofilter-1-sep-points.txt")
writer1.SetInputConnection(table1.GetOutputPort())
writer1.Update()
writer1.Write()

table11 = vtkTable()
for i in range(output1.GetCellData().GetNumberOfArrays()):
    table11.AddColumn(output1.GetCellData().GetArray(i))

writer11 = vtkDelimitedTextWriter()
writer11.SetFileName(f"{output_folder_path}/nofilter-1-sep-cells.txt")
#writer11.SetInputConnection(table11.GetOutputPort())
writer11.SetInputData(table11)
writer11.Update()
writer11.Write()

# simplifying the input data to remove non-persistent pairs
topologicalSimplification = ttkTopologicalSimplification()
topologicalSimplification.SetInputConnection(0, smooth.GetOutputPort())
topologicalSimplification.SetInputArrayToProcess(0, 0, 0, 0, "scalars")
topologicalSimplification.SetInputConnection(1, persistentPairs.GetOutputPort())
topologicalSimplification.SetDebugLevel(3)
#topologicalSimplification.Update()

# compute MS complex with filtration
morseSmaleComplex = ttkMorseSmaleComplex()
morseSmaleComplex.SetInputConnection(topologicalSimplification.GetOutputPort())
morseSmaleComplex.SetInputArrayToProcess(0, 0, 0, 0, "scalars")
if connector > 0:
    morseSmaleComplex.SetSaddleConnectorsPersistenceThreshold(persistent_theshold)
    morseSmaleComplex.SetReturnSaddleConnectors(connector)
morseSmaleComplex.SetDebugLevel(3)
morseSmaleComplex.Update()

output0 = morseSmaleComplex.GetOutputDataObject(0)
table0 = vtkDataObjectToTable()
table0.SetInputConnection(morseSmaleComplex.GetOutputPort(0))
table0.Update()
table0.GetOutput().AddColumn(output0.GetPoints().GetData())
table0.Update()
writer1 = vtkDelimitedTextWriter()
writer1.SetFileName(f"{output_folder_path}/critical-points.txt")
writer1.SetInputConnection(table0.GetOutputPort())
writer1.Update()
writer1.Write()

output1 = morseSmaleComplex.GetOutputDataObject(1)
table1 = vtkDataObjectToTable()
table1.SetInputConnection(morseSmaleComplex.GetOutputPort(1))
table1.Update()
table1.GetOutput().AddColumn(output1.GetPoints().GetData())
table1.Update()
writer1 = vtkDelimitedTextWriter()
writer1.SetFileName(f"{output_folder_path}/1-sep-points.txt")
writer1.SetInputConnection(table1.GetOutputPort())
writer1.Update()
writer1.Write()

table11 = vtkTable()
for i in range(output1.GetCellData().GetNumberOfArrays()):
    table11.AddColumn(output1.GetCellData().GetArray(i))

writer11 = vtkDelimitedTextWriter()
writer11.SetFileName(f"{output_folder_path}/1-sep-cells.txt")
#writer11.SetInputConnection(table11.GetOutputPort())
writer11.SetInputData(table11)
writer11.Update()
writer11.Write()
