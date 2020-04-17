#--Modules
import sys
import vtk
import numpy as np
from PyQt5 import Qt
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QPushButton, QAction, QLineEdit, QMessageBox, QLabel
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot

#--Main Window Start
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

class MainWindow(Qt.QMainWindow):

    def __init__(self, parent = None):
        Qt.QMainWindow.__init__(self, parent)
        
        #--Window and Renderer Setting
        
        #Main Window Initialization
        self.frame = Qt.QFrame()
        self.vl = Qt.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)

        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(255, 255, 255)
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        
        #Renderer Setting
        self.ren.ResetCamera()
        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)
        
        #Start Renderer
        self.iren.Initialize()
        self.show()      
        #----------------------------------
        
        #--Menu & Functions 
        
        #Menu Bar to Wrap All Menu
        mainMenu = self.menuBar()
        
        #Menu on Menu Bar
        fileMenu = mainMenu.addMenu('Meshgrid')
        
        #Submenu Show Meshgrid
        addmeshgrid = Qt.QAction('Add Meshgrid', self)
        addmeshgrid.setShortcut('Ctrl+M')
        addmeshgrid.triggered.connect(self.add_mesh)
        fileMenu.addAction(addmeshgrid)
        
        #Submenu Add Section
        add_section = Qt.QAction("Add Section", self)
        add_section.triggered.connect(self.add_section)
        fileMenu.addAction(add_section)           
        #-----------------------------------
        
        #--Data        
        self.meshgrid = np.loadtxt("meshgrid.txt")
        
    
    #--Functions    
    def add_mesh(self):  
        data = self.meshgrid
        
        #Parameter
        ndata = len(data[:,0])
        nx    = len(np.unique(data[:,0]))
        ny    = len(np.unique(data[:,1]))
        nz    = int(ndata/nx/ny)
        
        #Polygon
        points = vtk.vtkPoints()
        cells  = vtk.vtkCellArray()       
        point_data = vtk.vtkDoubleArray()
        scalar = vtk.vtkDoubleArray()
        
        #Filling Points
        for i in range(ndata):
            coord = data[i,0], data[i,1], data[i,2]
            points.InsertNextPoint(coord)
            point_data.InsertNextTuple([np.log10(data[i,3])*0.5])
            
        #Create Polygon
        pixel = vtk.vtkHexahedron()
        mesh = vtk.vtkUnstructuredGrid()
        
        for g in range(ny-1):
            for h in range(nx-1):
                for i in range(nz-1):
                    pixel.GetPointIds().SetId(0, nx*nz*g + nz*h + i)
                    pixel.GetPointIds().SetId(1, nx*nz*g + nz*(h+1) + i)
                    pixel.GetPointIds().SetId(2, nx*nz*(g+1) + nz*(h+1) + i)
                    pixel.GetPointIds().SetId(3, nx*nz*(g+1) + nz*h + i)
                    
                    pixel.GetPointIds().SetId(4, nx*nz*g + nz*h + i+1)
                    pixel.GetPointIds().SetId(5, nx*nz*g + nz*(h+1) + i+1)
                    pixel.GetPointIds().SetId(6, nx*nz*(g+1) + nz*(h+1) + i+1)
                    pixel.GetPointIds().SetId(7, nx*nz*(g+1) + nz*h + i+1)
                
                    cells.InsertNextCell(pixel)
                    mesh.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())
                    
                    scalar_value = (data[nx*nz*g + nz*h + i,3] + data[nx*nz*g + nz*(h+1) + i,3] + data[nx*nz*(g+1) + nz*(h+1) + i,3] + data[nx*nz*(g+1) + nz*h + i,3] + data[nx*nz*g + nz*h + i+1,3] + data[nx*nz*g + nz*(h+1) + i+1,3] + data[nx*nz*(g+1) + nz*(h+1) + i+1,3] + data[nx*nz*(g+1) + nz*h + i+1,3])/8
                    scalar.InsertNextTuple([np.log10(scalar_value)*0.5])
                    
        #Mesh
        mesh.SetPoints(points)
        mesh.GetCellData().SetScalars(scalar)
        
        #Create Mapper, Actor
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(mesh)
        
        self.actor  = vtk.vtkActor()
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().EdgeVisibilityOn()
        self.ren.AddActor(self.actor) 
        
        #Renderer Setting
        self.ren.ResetCamera()
        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)

        self.iren.Initialize()
        self.show()
    
    def remove_mesh(self):
        #Renderer Setting
        self.ren.RemoveActor(self.actor)
        
    def add_section(self):

        #Data
        meshgrid = self.meshgrid
        
        #Test Input
        print("section coordinate (x0, x1, y0, y1)")
        section_coor = input()
        
        section_coor = section_coor.split(",")
        
        x0 = float(section_coor[0])
        x1 = float(section_coor[1])
        y0 = float(section_coor[2])
        y1 = float(section_coor[3])
        
        #Data Parameter
        ndata = len(meshgrid[:,0])
        nx = len(np.unique(meshgrid[:,0]))
        ny = len(np.unique(meshgrid[:,1]))
        nz = int(ndata/nx/ny)
        
        minx = min(meshgrid[:,0])
        miny = min(meshgrid[:,1])
        
        delx = meshgrid[nz,0] - meshgrid[nz-1,0]
        dely = meshgrid[nz*nx,1] - meshgrid[nz*nx-1,1]
        
        #Mesh Creation
        d = (delx+dely)/2
        r = ((x1-x0)**2 + (y1-y0)**2)**0.5
        dr = r/(int(r/d))
        
        mesh = np.zeros((4,4))
        nsec = int(r/d)+1
        section = np.zeros((nsec*nz, 4))
        
        for i in range(nsec):
            x = x0 + dr*i*(x1-x0)/r
            y = y0 + dr*i*(y1-y0)/r
            
            nmx = int((x0 + d*i*(x1-x0)/r - minx)/delx)
            nmy = int((y0 + d*i*(y1-y0)/r - miny)/dely)
            
            #--Interpolasi Liniar
            
            for h in range(nz):
                section[i*nz+h, 0] = x
                section[i*nz+h, 1] = y
                
                #Neighbors
                mesh[0,:] = meshgrid[(nx*(nmy+0)+nmx+0)*nz + h, :]
                mesh[1,:] = meshgrid[(nx*(nmy+0)+nmx+1)*nz + h, :]        
                mesh[2,:] = meshgrid[(nx*(nmy+1)+nmx+0)*nz + h, :]
                mesh[3,:] = meshgrid[(nx*(nmy+1)+nmx+1)*nz + h, :]
                
                #Elevation Interpolation
                a0 = ((mesh[1,2] - mesh[0,2])/(mesh[1,0] - mesh[0,0])*(section[i*nz+h,0] - mesh[0,0])) + mesh[0,2]
                b0 = ((mesh[3,2] - mesh[2,2])/(mesh[3,0] - mesh[2,0])*(section[i*nz+h,0] - mesh[2,0])) + mesh[2,2]
                section[i*nz+h,2]  = ((b0 - a0) / (mesh[2,1]-mesh[0,1])*(section[i*nz+h,1] - mesh[0,1])) + a0
                
                a0 = ((mesh[1,3] - mesh[0,3])/(mesh[1,0] - mesh[0,0])*(section[i*nz+h,0] - mesh[0,0])) + mesh[0,3]
                b0 = ((mesh[3,3] - mesh[2,3])/(mesh[3,0] - mesh[2,0])*(section[i*nz+h,0] - mesh[2,0])) + mesh[2,3]
                section[i*nz+h,3]  = np.log10(((b0 - a0) / (mesh[2,1]-mesh[0,1])*(section[i*nz+h,1] - mesh[0,1])) + a0)
        
        self.test = section
        np.savetxt("section.txt", section, "%.2f")
        
        #--Creating Vtk Object
        #Create Polygon
        points = vtk.vtkPoints()
        mesh   = vtk.vtkUnstructuredGrid()
        scalar = vtk.vtkDoubleArray()
        
        #Filling Points
        for i in range(nsec*nz):
            points.InsertNextPoint(section[i,0], section[i,1], section[i,2])
            
        #Insert Cells
        quad = vtk.vtkQuad()
        
        for i in range(nsec-1):
            for h in range(nz-1):
                quad.GetPointIds().SetId(0, h + nz*i + 0)
                quad.GetPointIds().SetId(1, h + nz*i + 1)
                quad.GetPointIds().SetId(2, h + nz*(i+1) + 1)
                quad.GetPointIds().SetId(3, h + nz*(i+1) + 0)
                
                mesh.InsertNextCell(quad.GetCellType(), quad.GetPointIds())
                
                #Cell Scalar Value
                value = (section[ h + nz*i + 0, 3] + section[ h + nz*i + 1, 3] + section[ h + nz*(i+1) + 1, 3] + section[ h + nz*(i+1) + 0, 3])/4
                scalar.InsertNextTuple([(value)*0.5])
                
        mesh.SetPoints(points)
        mesh.GetCellData().SetScalars(scalar)
        
        #Create Mapper, Actor
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(mesh)
        
        self.section_actor  = vtk.vtkActor()
        self.section_actor.SetMapper(mapper)
        self.section_actor.GetProperty().EdgeVisibilityOn()
        
        self.ren.AddActor(self.section_actor) 
        
        #Renderer Setting
        self.ren.ResetCamera()
        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)

        self.iren.Initialize()
        self.show()

#--Widget

#Pop Up Window        
class examplePopup(QWidget):
    def __init__(self, name):
        super().__init__()

        self.name = name

        self.initUI()

    def initUI(self):
        lblName = QLabel(self.name, self)
       
if __name__ == "__main__":
    app = Qt.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())