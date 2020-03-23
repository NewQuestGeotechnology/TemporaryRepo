#Tes DEM
import numpy as np

#Data Loading
data = np.loadtxt("GeoSystemTopo01.txt")
stasiun = np.loadtxt("Koordinat_Stasiun.txt", usecols=(1,2,3))
#dem = np.loadtxt("DEM_SORT.txt")

#Parameter
ndata = len(data[:,0])
xmin = min(data[:,0])
ymin = min(data[:,1])
xmax = max(data[:,0])
ymax = max(data[:,1])

#Center of Station
x_center = (max(stasiun[:,0]) + min(stasiun[:,0]))/2 
y_center = (max(stasiun[:,1]) + min(stasiun[:,1]))/2 

#Center of Data (Dummy)
data_xcenter = (xmin + xmax)/2
data_ycenter = (ymin + ymax)/2

#Center of Data (Real)
x0 = x_center - (data_xcenter - xmin)
y0 = y_center - (data_ycenter - ymin)

#Coordinate Conversion
for i in range(ndata):
    data[i,0] = data[i,0] - xmin + x0
    data[i,1] = data[i,1] - ymin + y0

np.savetxt("DATA_CENTER.txt", data, "%.0f")