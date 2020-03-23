#Tes DEM
import numpy as np

data = np.loadtxt("DEM_SORT.xyz")
ndata = len(data[:,0])

#Sorting Data
y = np.unique(data[:,1])

#Sorting Data
sort_data = np.array(([[0,0,0]]))

for i in range(len(y)):
    temp_data = data[data[:,1] == y[i]]
    sort_data = np.concatenate((sort_data,temp_data), axis=0)

sort_data = sort_data[1:-1,:]

np.savetxt("DEM_SORT.txt", sort_data, "%.1f")