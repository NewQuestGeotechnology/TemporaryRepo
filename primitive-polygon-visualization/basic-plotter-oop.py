#Object Oriented Program

#Modules
import vtk
import numpy as np

#Class Initiation

#Surface/Topography
class dem_surface:
    def __init__(self, name, data, main):
        self.name = name
        self.data = data
        self.xyz = np.loadtxt(self.data)
        
        #Data Selection
        self.xyz = self.xyz[(self.xyz[:,0] > main.xmin) & (self.xyz[:,0] < main.xmax) & (self.xyz[:,1] > main.ymin) & (self.xyz[:,1] < main.ymax)]

        nx = len(np.unique(self.xyz[:,0]))
        ny = len(np.unique(self.xyz[:,0]))
        
        #Data Sorting        
        sort_data = np.zeros((1,3))
        y = np.unique(self.xyz[:,1])
        
        for i in range(len(y)):
            temp_data = self.xyz[self.xyz[:,1] == y[i]]
            sort_data = np.concatenate((sort_data,temp_data), axis=0)

        self.xyz = sort_data[1:-1,:]
            
        #Polygon Components
        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        #Fill Points
        for i in range(ny):
            for j in range(nx):
                coord = self.xyz[nx*i + j, 0], self.xyz[nx*i + j, 1], self.xyz[nx*i + j, 2]
                points.InsertNextPoint(coord)

        #Fill Rectangle
        quad = vtk.vtkQuad()

        for h in range(ny):
            if h < ny - 1:
                for i in range(nx):

                    if i < nx - 1:
                        quad.GetPointIds().SetId(0, nx*h + i + nx + 0)
                        quad.GetPointIds().SetId(1, nx*h + i + nx + 1)
                        quad.GetPointIds().SetId(2, nx*h + i + 1)
                        quad.GetPointIds().SetId(3, nx*h + i + 0)

                        cells.InsertNextCell(quad)

                    else:
                        break

            else:
                break
        
        #Mesh Data
        self.mesh = vtk.vtkPolyData()
        self.mesh.SetPoints(points)
        self.mesh.SetPolys(cells)
        main.surface_renderer(self.mesh)
        
#Fault
class fault:
    def __init__(self, name, data, main):
        self.name = name
        self.data = data
        self.xyz = np.loadtxt(data)

        #Parameter
        ndata = 51

        #Polygon
        point = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        #Filling Points
        for i in range(ndata*2):
            coord = self.xyz[i,0], self.xyz[i,1], self.xyz[i,2]
            point.InsertNextPoint(coord) 

        #Create Polygon
        quad = vtk.vtkQuad()

        for i in range(ndata-2):
            quad.GetPointIds().SetId(0, ndata+i+0)
            quad.GetPointIds().SetId(1, ndata+i+1)
            quad.GetPointIds().SetId(2, i+1)
            quad.GetPointIds().SetId(3, i+0)

            cells.InsertNextCell(quad)

        #Create Mesh
        self.mesh = vtk.vtkPolyData()
        self.mesh.SetPoints(point)
        self.mesh.SetPolys(cells)
        main.fault_renderer(self.mesh)

#Well
class well:
    def __init__(self, name, data, main):
        self.name = name
        self.xyz = np.loadtxt(data)

        #parameter
        ndata = len(self.xyz[:,0])

        #Polygon
        point = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        #Filling Points
        for i in range(ndata):
            coord = self.xyz[i,0], self.xyz[i,1], self.xyz[i,2]
            point.InsertNextPoint(coord) 

        #Create Polyline
        line = vtk.vtkPolyLine()
        
        line.GetPointIds().SetNumberOfIds(ndata)
        for i in range(ndata):
            line.GetPointIds().SetId(i,i)

        cells.InsertNextCell(line)

        #Create Mesh
        self.mesh = vtk.vtkPolyData()
        self.mesh.SetPoints(point)
        self.mesh.SetLines(cells)
        
        #Point Data
        point_data = vtk.vtkDoubleArray()
        
        for i in range(ndata):
            point_data.InsertNextTuple([self.xyz[i,2]*0.0005])
            
        self.mesh.GetPointData().SetScalars(point_data)
        main.well_renderer(self.mesh)
        
#Section
class resistivity_section: #Input must already a section file
    def __init__(self, name, data, main):
        self.data = np.loadtxt(data)
        
        #Parameter
        ndata = len(self.data[:,0])
        
        #Create Polygon
        points = vtk.vtkPoints()
        cells  = vtk.vtkCellArray()
        
        #Filling Points
        for i in range(ndata):
            coord = self.data[i,0], self.data[i,1], self.data[i,2]
            points.InsertNextPoint(coord)
            
        #Creating Cells
        triangle = vtk.vtkTriangle()
        quad = vtk.vtkQuad()
        
        #Data Identification
        block_number = np.unique(self.data[:,1])
        block_sum = []
        
        for i in range(len(block_number)):
            block_sum = block_sum + [len(self.data[self.data[:,1] == block_number[i]])]    
        
        #--Filling CellArray
        
        #Triangle Cells
        triangle.GetPointIds().SetId(0,0) 
        triangle.GetPointIds().SetId(1,1)
        triangle.GetPointIds().SetId(2,block_sum[0])
        
        cells.InsertNextCell(triangle)
        
        #Quad Cell
        for i in range(block_sum[0]-2):
            quad.GetPointIds().SetId(0,i+2)
            quad.GetPointIds().SetId(1,block_sum[0]+i+1)
            quad.GetPointIds().SetId(2,block_sum[0]+i+0)
            quad.GetPointIds().SetId(3,i+1)
            
            cells.InsertNextCell(quad)
        
        #PointValue
        point_data = vtk.vtkDoubleArray()
        
        for i in range(ndata):
            point_data.InsertNextTuple(np.log10([self.data[i,3]]))
        
        #Create Mesh & Set Value
        self.mesh = vtk.vtkPolyData()
        self.mesh.SetPoints(points)
        self.mesh.SetPolys(cells)
        self.mesh.GetPointData().SetScalars(point_data)
        
        #Renderer
        main.surface_renderer(self.mesh)
        
#Visualization
class visualize:

    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        self.renderer = vtk.vtkRenderer()
        self.colors = vtk.vtkNamedColors()
        self.xmin = int(xmin)
        self.ymin = int(ymin)
        self.xmax = int(xmax)
        self.ymax = int(ymax)
        self.zmin = int(zmin)
        self.zmax = int(zmax)

    def surface_renderer(self, point):    
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(point)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetPointSize(10)

        self.renderer.AddActor(actor)
        
    def fault_renderer(self, point):
        color = self.colors.GetColor3d("red")

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(point)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetOpacity(0.5)
        actor.GetProperty().SetPointSize(10)

        self.renderer.AddActor(actor)

    def well_renderer(self, point):    
        #Testing Tube
        tube = vtk.vtkTubeFilter()
        tube.SetInputData(point)
        tube.SetRadius(100)
        tube.SetNumberOfSides(15)
        tube.Update()

        tubeMapper = vtk.vtkPolyDataMapper()
        tubeMapper.SetInputConnection(tube.GetOutputPort())

        tubeActor = vtk.vtkActor()
        tubeActor.SetMapper(tubeMapper)

        color = self.colors.GetColor3d("blue")
        tubeActor.GetProperty().SetColor(color)

        self.renderer.AddActor(tubeActor)        

    def show_window(self):
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(self.renderer)
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)

        renderWindow.Render()
        renderWindowInteractor.Start()

#Main Window
main = visualize(797000,806000,9194000,9206000,-3000,3000)

#Add Fault
#fault1 = fault("fault_1", "fault.txt", main)
#fault1 = fault("fault_2", "fault_2.txt", main)

#Add DEM
#dem = dem_surface("surface", "DEM.xyz", main)

#Add Well
#well1 = well("well_1", "d14trj.txt", main)

#Add Section
section = resistivity_section("section_1", "section.txt", main)

#Show Window
main.show_window()