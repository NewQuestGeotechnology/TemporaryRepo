#Making Meshgrid Class
import numpy as np
import utilites
from datetime import datetime

print("start: " + str(datetime.now().time()))
class meshgridc:
    
    def __init__(self, data, station, topo, xmin, xmax, ymin, ymax, zmin, zmax, d):
        self.data = np.loadtxt(data)
        self.station = np.loadtxt(station, usecols=(1,2,3))
        self.topo = np.loadtxt(topo)
        
        #Data & Topo Filtering
        self.topo = self.topo[(self.topo[:,0] > xmin-100) & (self.topo[:,0] < xmax+100) & (self.topo[:,1] > ymin-100) & (self.topo[:,1] < ymax+100)]
        self.data = self.data[(self.data[:,0] > xmin) & (self.data[:,0] < xmax) & (self.data[:,1] > ymin) & (self.data[:,1] < ymax) & (self.data[:,2] > zmin) & (self.data[:,2] < zmax)]
        
        #Meshgrid
        self.meshgrid = utilites.meshgrid(self.topo, self.data, xmin, xmax, ymin, ymax, zmin, d)        
     
#Data
mesh    = "GeoSystem2020-topo.txt"
topo    = "DEM.xyz"
station = "Koordinat_Stasiun.txt"

d    = 300
xmin = 795727.0
xmax = 805323.0
ymin = 9197285.0
ymax = 9206100.0
zmin = -3000
zmax = 3000

#Object 
meshgrid = meshgridc(mesh, station, topo, xmin, xmax, ymin, ymax, zmin, zmax, d)
print("start: " + str(datetime.now().time()))