#--Modules
import numpy as np
import vtk

#--Data
section = np.loadtxt("section.txt")
nsec    = len(np.unique(section[:,0]))
nz      = int(len(section[:,0])/nsec)

#--Create Polygon

#Create Polygon
points = vtk.vtkPoints()
scalar = vtk.vtkDoubleArray()
pointScalar = vtk.vtkDoubleArray()
cells  = vtk.vtkCellArray()

gridMesh = vtk.vtkUnstructuredGrid()
polyMesh = vtk.vtkPolyData()
      
#Filling Points
for i in range(nsec*nz):
    points.InsertNextPoint(section[i,0], section[i,1], section[i,2])
    pointScalar.InsertNextTuple([section[i,3]])
            
#Insert Cells
quad = vtk.vtkQuad()
        
for i in range(nsec-1):
    for h in range(nz-1):
        quad.GetPointIds().SetId(0, h + nz*i + 0)
        quad.GetPointIds().SetId(1, h + nz*i + 1)
        quad.GetPointIds().SetId(2, h + nz*(i+1) + 1)
        quad.GetPointIds().SetId(3, h + nz*(i+1) + 0)
                
        gridMesh.InsertNextCell(quad.GetCellType(), quad.GetPointIds())
                
        #Cell Scalar Value
        value = (section[ h + nz*i + 0, 3] + section[ h + nz*i + 1, 3] + section[ h + nz*(i+1) + 1, 3] + section[ h + nz*(i+1) + 0, 3])/4
        scalar.InsertNextTuple([(value)])
        
        #Polydata Cell
        cells.InsertNextCell(quad)

#Grid Mesh
gridMesh.SetPoints(points)
gridMesh.GetCellData().SetScalars(scalar)

gridRange = gridMesh.GetScalarRange()

Range = polyMesh.GetScalarRange()
gridMapper = vtk.vtkDataSetMapper()
gridMapper.SetInputData(gridMesh)
gridMapper.SetScalarRange(gridRange)

gridActor = vtk.vtkActor()
gridActor.SetMapper(gridMapper)
        
#Poly Mesh
polyMesh.SetPoints(points)
polyMesh.SetPolys(cells)
polyMesh.GetPointData().SetScalars(pointScalar)

polyMapper = vtk.vtkPolyDataMapper()
polyMapper.SetInputData(polyMesh)
#polyMapper.SetScalarRange(Range)

polyActor = vtk.vtkActor()
polyActor.SetMapper(polyMapper)

#--Contour Making
nCont = 15

#Contour Bands
Range = polyMesh.GetScalarRange()
dR    = (Range[1] - Range[0])/nCont

bands = []

for i in range(nCont):
    bands += [Range[0] + dR*i]

#Contour Filter
contour = vtk.vtkBandedPolyDataContourFilter()
contour.SetInputData(polyMesh)

for i in range(nCont):
    contour.SetValue(i, bands[i])
    
contour.SetScalarModeToIndex()
contour.GenerateContourEdgesOn()

contourMapper = vtk.vtkPolyDataMapper()
contourMapper.SetInputConnection(contour.GetOutputPort())
contourMapper.SetScalarRange(Range)
#contourMapper.SetScalarModeToUseCellData() #Untuk categorical Data

contourActor = vtk.vtkActor()
contourActor.SetMapper(contourMapper)
contourActor.GetProperty().SetOpacity(0.0)

#Contour Line
edgeMapper = vtk.vtkPolyDataMapper()
edgeMapper.SetInputData(contour.GetContourEdgesOutput()) #Extrack Edge/Line from contour Filter

edgeMapper.SetResolveCoincidentTopologyToPolygonOffset()

edgeActor = vtk.vtkActor()
edgeActor.SetMapper(edgeMapper)
edgeActor.GetProperty().SetColor(0, 0, 0)

#--Renderer and Window
renderer = vtk.vtkRenderer()
renderer.AddActor(contourActor) #smooth section
renderer.AddActor(edgeActor) #contour line
renderer.AddActor(gridActor) #grid section

#Renderer Properties
renderer.SetBackground(255, 255, 255)

renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderWindow.Render()
renderWindowInteractor.Start()