#Resistivity Section

#Modules
import numpy as np

#Function
def resistivity_section(data, xmin, xmax, ymin, ymax, zmin, zmax):
    
    #Data import & Parameter
    data = np.loadtxt(data)
    
    #Data Selection (Section)
    x = np.unique(data[:,0])
    temp_data = data[(data[:,0] == x[int(len(x)/2)])]
    
    #Data Selection (Wide)
    temp_data = temp_data[temp_data[:,1] > 28000]
    temp_data = temp_data[temp_data[:,1] < 38000]
    
    #Data Selection (Depth)
    temp_data = temp_data[temp_data[:,2] > zmin]
    temp_data = temp_data[temp_data[:,2] < zmax]
    
    #Sorting Data
    y = np.unique(temp_data[:,1])    
    
    temporary = np.array([[0,0,0,0]])
    for i in range(len(y)):
        temporary_data = temp_data[temp_data[:,1] == y[i]]
        temporary = np.concatenate((temporary,temporary_data), axis=0)
        
    temporary = temporary[1:-1,:]
    
    section = np.zeros((len(temp_data),2))
    section[:,0] = temp_data[:,1]
    section[:,1] = temp_data[:,2]
    
    return temporary
     
#Execution
section = resistivity_section("GeoSystemTopo01.txt", 797000,806000,9194000,9206000,-3000,3000)
print(len(section))
np.savetxt("section.txt", section, "%.0f")