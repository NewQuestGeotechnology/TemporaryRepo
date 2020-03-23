#Plot Resistivity
import numpy as np
import matplotlib.pyplot as plt

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
    
#View area (xmin, xmax, ymin, ymax)
view_area = np.zeros((4))
view_area[0] = min(stasiun[:,0] - 2000)
view_area[1] = max(stasiun[:,0] + 2000)
view_area[2] = min(stasiun[:,1] - 2000)
view_area[3] = max(stasiun[:,1] + 2000)

#Data selection
def data_selection(data, xmin, xmax, ymin, ymax, zmin, zmax):
    selected_data = data[(data[:,0]>xmin) & (data[:,0]<xmax) & (data[:,1]>ymin) & (data[:,1]<ymax) & (data[:,2]>zmin) & (data[:,2]<zmax)]
    return selected_data

data = data_selection(data, view_area[0], view_area[1], view_area[2], view_area[3], -3000, 3000)

#Tes section
section = data[data[:,1] == data[0,1]]

def data_sorting(data, axis):
    index = np.unique(data[:,axis])
    
    sorted_data = np.array(([[0,0,0,0]]))
    
    for i in range(len(index)):
        temp_data = data[data[:,axis] == index[i]]
        sorted_data = np.concatenate((sorted_data,temp_data), axis=0)
        
    sorted_data = sorted_data[1:-1,:]
    
    return sorted_data

section = data_sorting(section,0)

section0 = section[(section[:,0] == section[0,0]) & (section[:,1] == section[0,1])]
section1 = section[(section[:,0] == section[47,0]) & (section[:,1] == section[47,1])]
section2 = section[(section[:,0] == section[98,0]) & (section[:,1] == section[98,1])]

#Data Plotting
plt.scatter(section0[:,0], section0[:,2])
plt.scatter(section1[:,0], section1[:,2], c="red")
plt.scatter(section2[:,0], section2[:,2])
#np.savetxt("section_test.txt", section, "%.0f")