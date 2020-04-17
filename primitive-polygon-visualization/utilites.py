#Utilities

import numpy as np
import vtk

#Data Sorter
def data_sorter(data):
    
    #--Parameter
    y = np.unique(data[:,1])
    
    #--Data Sorting
    
    #Data Template
    sorted_data = np.array([[0]*len(data[0,:])])
    
    #Sorting Process
    for i in range(len(y)):
        temp_data = data[data[:,1] == y[i]]
                   
        sorted_data = np.concatenate((sorted_data, temp_data))
     
    sorted_data = sorted_data[1:-1,:]
    
    return sorted_data

#xyz interpolator
def xyz_grid_interpolator(data, point):
    
    #--Parameters
    x = np.unique(data[:,0])
    y = np.unique(data[:,1])
    
    nx = len(x)
    ny = len(y)
    
    xmin_data = min(data[:,0])
    ymin_data = min(data[:,1])
    
    #Average Grid Spacing
    dx = (np.sum(x - x[0]))/((nx-1)/2*(2 + (nx-2)))
    dy = (np.sum(y - y[0]))/((ny-1)/2*(2 + (ny-2)))
    
    #--Point Location
    nxp = int((point[0] - xmin_data)/dx)
    nyp = int((point[1] - ymin_data)/dy)
        
    #neighbor points
    ng = np.zeros((4,3))
        
    ng[0,:] = data[nx*nyp + nxp + 0, :]
    ng[1,:] = data[nx*nyp + nxp + 1, :]
    ng[2,:] = data[nx*(nyp+1) + nxp + 0, :]
    ng[3,:] = data[nx*(nyp+1) + nxp + 1, :]
        
    #Linear Interpolation
    x1 = (ng[1,2]-ng[0,2]) / (ng[1,0]-ng[0,0])*(point[0]-ng[0,0]) + ng[0,2]
    x2 = (ng[3,2]-ng[2,2]) / (ng[3,0]-ng[2,0])*(point[0]-ng[2,0]) + ng[2,2]
    x3 = (x2-x1) / (ng[2,1]-ng[0,1])*(point[1]-ng[0,1]) + x1
        
    return x3
    
#Pythagoras
def pythagoras(data1,data2):
    ndata = len(data1)
    
    if ndata == 4:
        return ((data1[0]-data2[0])**2 + (data1[1]-data2[1])**2 + (data1[2]-data2[2])**2)**0.5
    
    elif ndata == 3:
        return ((data1[0]-data2[0])**2 + (data1[1]-data2[1])**2)**0.5
    
#Krigging For Meshgrid (Ordinary)
def xyzv_krigging(data, grid):
    
    #Data Range
    data_xrange = max(data[:,0]) - min(data[:,0])
    data_yrange = max(data[:,1]) - min(data[:,1])
    data_zrange = max(data[:,2]) - min(data[:,2])
    
    #--Neighboring Data Identification
    def krigcal(point, data, data_xrange, data_yrange, data_zrange):
        alpha = 0
        neighbors = []
        
        while len(neighbors) < 4:
            alpha = alpha + 0.05
            neighbors = data[(data[:,0] > (point[0]-data_xrange*alpha)) & (data[:,0] < (point[0]+data_xrange*alpha)) & (data[:,1] > (point[1]-data_yrange*alpha)) & (data[:,1] < (point[1]+data_yrange*alpha)) & (data[:,2] > (point[2]-data_zrange*alpha)) & (data[:,2] < (point[2]+data_zrange*alpha))]
        
        sort_index = np.argsort(np.apply_along_axis(pythagoras, 1, neighbors, point))
        
        neighborsb = np.zeros((4,4))
        neighborsb[0,:] = neighbors[sort_index[0],:]
        neighborsb[1,:] = neighbors[sort_index[1],:]
        neighborsb[2,:] = neighbors[sort_index[2],:]
        neighborsb[3,:] = neighbors[sort_index[3],:]

        neighbors = neighborsb
            
        nn = len(neighbors[:,0])
        
        if np.average(neighbors[:,3]) == neighbors[0,3]:
            output = neighbors[0,3]
        
        else:  
            
            #--Creating Semivariogram
            #Lag Tolerance
            nlag = 10
            
            rmax = ((min(neighbors[:,0])-max(neighbors[:,0]))**2+(min(neighbors[:,1])-max(neighbors[:,1]))**2+(min(neighbors[:,2])-max(neighbors[:,2]))**2)**0.5
            dlag = rmax/(nlag)
            lag  = np.arange(0, rmax+dlag, dlag)
            
            #Calculating Semivariogram
            svariogram = np.zeros((len(lag),3))
            svariogram[:,0] = lag
            
            #Semivariance
            start = 1
            
            for h in range(nn):
                for i in range(start,nn):
                    index = int(((neighbors[i,0]-neighbors[h,0])**2+(neighbors[i,1]-neighbors[h,1])**2+(neighbors[i,2]-neighbors[h,2])**2)**0.5/dlag)
                    
                    svariogram[index,1] += 1
                    svariogram[index,2] += (neighbors[i,3] - neighbors[h,3])**2
                    
                start += 1 
            
            #Semivariance Calculation
            svariogram = svariogram[svariogram[:,1] != 0]
            
            def semvar(data):
                return data / [1, 1, (2*data[1])]
            
            svariogram = np.apply_along_axis(semvar, 1, svariogram)
            
            #Linear Semivariogram Coef
            sill = max(svariogram[:,2])*0.9
            rng  = ((data_xrange*alpha)**2 + (data_yrange*alpha)**2 + (data_zrange*alpha)**2)**2
            
            #K Matrix
            K = np.zeros((nn,nn))
            
            for h in range(nn):
                for i in range(nn):
                    K[h,i] = ((neighbors[i,0]-neighbors[h,0])**2 + (neighbors[i,1]-neighbors[h,1])**2 + (neighbors[i,2]-neighbors[h,2])**2)**0.5
                    K[h,i] = sill - sill*(1.5*(K[h,i]/rng) - 0.5*(K[h,i]/rng)**2)   
                    
            #lam Matrix
            def cov(data, b, rng, sill):
                r = ((data[0]-b[0])**2 + (data[1]-b[1])**2 + (data[2]-b[2])**2)**0.5
                return sill - sill*(1.5*(r/rng) - 0.5*(r/rng)**2)
            
            k = np.apply_along_axis(cov, 1, neighbors, point, rng, sill)
                
            lam = np.linalg.solve(K,k)
            
            #Regional & Residual
            average = (neighbors[0,3] + neighbors[1,3] + neighbors[2,3] + neighbors[3,3])/4
            res = neighbors[:,3] - average
            
            #Result
            output = average + np.dot(res.T, lam)
        
        return output
    
    return np.apply_along_axis(krigcal, 1, grid, data, data_xrange, data_yrange, data_zrange)

#Krigging For Meshgrid (Ordinary)
def xyzv_linear_gridding(data, grid):
    
    #Data Range
    data_xrange = max(data[:,0]) - min(data[:,0])
    data_yrange = max(data[:,1]) - min(data[:,1])
    data_zrange = max(data[:,2]) - min(data[:,2])
    
    #--Neighboring Data Identification
    def krigcal(point, data, data_xrange, data_yrange, data_zrange):
        alpha = 0
        neighbors = []
        
        while len(neighbors) < 4:
            alpha = alpha + 0.05
            neighbors = data[(data[:,0] > (point[0]-data_xrange*alpha)) & (data[:,0] < (point[0]+data_xrange*alpha)) & (data[:,1] > (point[1]-data_yrange*alpha)) & (data[:,1] < (point[1]+data_yrange*alpha)) & (data[:,2] > (point[2]-data_zrange*alpha)) & (data[:,2] < (point[2]+data_zrange*alpha))]
        
        sort_index = np.argsort(np.apply_along_axis(pythagoras, 1, neighbors, point))
        
        neighborsb = np.zeros((4,4))
        neighborsb[0,:] = neighbors[sort_index[0],:]
        neighborsb[1,:] = neighbors[sort_index[1],:]
        neighborsb[2,:] = neighbors[sort_index[2],:]
        neighborsb[3,:] = neighbors[sort_index[3],:]

        neighbors = neighborsb
        
        #Linear Approach
        A = np.zeros((4,4))
        b = np.zeros((4))
            
        A[:,0] = neighbors[:,0]
        A[:,1] = neighbors[:,1]
        A[:,2] = neighbors[:,2]
        A[:,3] = 1
            
        b[:] = neighbors[:,3]
        
        if np.linalg.det(np.dot(A.T,A)) == 0:
            output =neighbors[0,3]
            
        else:  
            #Solving
            coef = np.linalg.solve(np.dot(A.T,A), np.dot(A.T,b))
            output = point[0]*coef[0] + point[1]*coef[1] + point[2]*coef[2] + coef[3]
            
        return output
    
    return np.apply_along_axis(krigcal, 1, grid, data, data_xrange, data_yrange, data_zrange)

#Meshgrid Creation
def meshgrid(topo, data, xmin, xmax, ymin, ymax, zmin, d):
    
    #topo filtering
    topo = topo[(topo[:,0] > xmin-300) & (topo[:,0] < xmax+300) & (topo[:,1] > ymin-300) & (topo[:,1] < ymax+300)]
    topo = data_sorter(topo)
    
    #meshgrid X and Y
    dx = (xmax - xmin)/int((xmax - xmin)/d)
    x  = np.concatenate((np.arange(xmin, xmax, dx), np.array([xmax])))
    nx = len(x)
    
    dy = (ymax - ymin)/int((ymax - ymin)/d)
    y  = np.concatenate((np.arange(ymin, ymax, dy), np.array([ymax])))
    ny = len(y)
    
    #Meshgrid z
    zmax = min(topo[:,2])
    zblock = 10
    multp  = 1.5
    z_bm = [0]
    
    #Z Benchmark
    i = 0
    
    while z_bm[i] < zmax - (zmin):
        z_bm += [z_bm[i] + zblock*(multp**(i+1))]
        i += 1
    
    z_bm[-1] = z_bm[-1] - (z_bm[-1] - (zmax - (zmin)))
    
    z = -np.array(z_bm)
    nz = len(z)
    
    #Meshgrid Iniziation
    meshgrid = np.array([[0,0,0,0]])
    temp     = np.zeros((nz,4))
    
    #Meshgrid Creation
    for i in range(ny):
        for h in range(nx):
            elev = xyz_grid_interpolator(topo, np.array([x[h], y[i]]))
            ztemp = z + elev
            ztemp[-1] = zmin
            
            temp[:,0] = x[h]
            temp[:,1] = y[i]
            temp[:,2] = ztemp
            
            meshgrid = np.concatenate((meshgrid, temp), 0)
            
    meshgrid = meshgrid[1::,:]
    
    #Krigging Process
    def xyzv_krigging(data, grid):
    
        #Data Range
        data_xrange = max(data[:,0]) - min(data[:,0])
        data_yrange = max(data[:,1]) - min(data[:,1])
        data_zrange = max(data[:,2]) - min(data[:,2])
    
        #--Neighboring Data Identification
        def krigcal(point, data, data_xrange, data_yrange, data_zrange):
            alpha = 0
            neighbors = []
        
            while len(neighbors) < 4:
                alpha = alpha + 0.05
                neighbors = data[(data[:,0] > (point[0]-data_xrange*alpha)) & (data[:,0] < (point[0]+data_xrange*alpha)) & (data[:,1] > (point[1]-data_yrange*alpha)) & (data[:,1] < (point[1]+data_yrange*alpha)) & (data[:,2] > (point[2]-data_zrange*alpha)) & (data[:,2] < (point[2]+data_zrange*alpha))]
            
            sort_index = np.argsort(((neighbors[:,0]-point[0])**2 + (neighbors[:,1]-point[1])**2 + (neighbors[:,2]-point[2])**2)**0.5)
            
            neighborsb = np.zeros((4,4))
            neighborsb[0,:] = neighbors[sort_index[0],:]
            neighborsb[1,:] = neighbors[sort_index[1],:]
            neighborsb[2,:] = neighbors[sort_index[2],:]
            neighborsb[3,:] = neighbors[sort_index[3],:]

            neighbors = neighborsb
            
            nn = len(neighbors[:,0])
            
            if np.average(neighbors[:,3]) == neighbors[0,3]:
                output = neighbors[0,3]
        
            else:  
            
            #--Creating Semivariogram
            #Lag Tolerance
                nlag = 10
            
                rmax = ((min(neighbors[:,0])-max(neighbors[:,0]))**2+(min(neighbors[:,1])-max(neighbors[:,1]))**2+(min(neighbors[:,2])-max(neighbors[:,2]))**2)**0.5
                dlag = rmax/(nlag)
                lag  = np.arange(0, rmax+dlag, dlag)
            
                #Calculating Semivariogram
                svariogram = np.zeros((len(lag),3))
                svariogram[:,0] = lag
            
                #Semivariance
                start = 1
            
                for h in range(nn):
                    for i in range(start,nn):
                        index = int(((neighbors[i,0]-neighbors[h,0])**2+(neighbors[i,1]-neighbors[h,1])**2+(neighbors[i,2]-neighbors[h,2])**2)**0.5/dlag)
                    
                        svariogram[index,1] += 1
                        svariogram[index,2] += (neighbors[i,3] - neighbors[h,3])**2
                    
                    start += 1 
            
                #Semivariance Calculation
                svariogram = svariogram[svariogram[:,1] != 0]
            
                def semvar(data):
                    return data / [1, 1, (2*data[1])]
            
                svariogram = np.apply_along_axis(semvar, 1, svariogram)
            
                #Linear Semivariogram Coef
                sill = max(svariogram[:,2])*0.9
                rng  = ((data_xrange*alpha)**2 + (data_yrange*alpha)**2 + (data_zrange*alpha)**2)**2
            
                #K Matrix
                K = np.zeros((nn,nn))
            
                for h in range(nn):
                    for i in range(nn):
                        K[h,i] = ((neighbors[i,0]-neighbors[h,0])**2 + (neighbors[i,1]-neighbors[h,1])**2 + (neighbors[i,2]-neighbors[h,2])**2)**0.5
                        K[h,i] = sill - sill*(1.5*(K[h,i]/rng) - 0.5*(K[h,i]/rng)**2)   
                    
                    #lam Matrix
                def cov(data, b, rng, sill):
                    r = ((data[0]-b[0])**2 + (data[1]-b[1])**2 + (data[2]-b[2])**2)**0.5
                    return sill - sill*(1.5*(r/rng) - 0.5*(r/rng)**2)
            
                k = np.apply_along_axis(cov, 1, neighbors, point, rng, sill)
                
                lam = np.linalg.solve(K,k)
            
                #Regional & Residual
                average = (neighbors[0,3] + neighbors[1,3] + neighbors[2,3] + neighbors[3,3])/4
                res = neighbors[:,3] - average
            
                #Result
                output = average + np.dot(res.T, lam)
        
            return output
    
        return np.apply_along_axis(krigcal, 1, grid, data, data_xrange, data_yrange, data_zrange)
    
    c = xyzv_krigging(data, meshgrid)
    meshgrid[:,3] = c
    
    return meshgrid
    

#Krigging For Meshgrid (Ordinary)
def xyzv_linear_gridding(data, grid):
    
    #Data Range
    data_xrange = max(data[:,0]) - min(data[:,0])
    data_yrange = max(data[:,1]) - min(data[:,1])
    data_zrange = max(data[:,2]) - min(data[:,2])
    
    #--Neighboring Data Identification
    def krigcal(point, data, data_xrange, data_yrange, data_zrange):
        alpha = 0
        neighbors = []
        
        while len(neighbors) < 4:
            alpha = alpha + 0.05
            neighbors = data[(data[:,0] > (point[0]-data_xrange*alpha)) & (data[:,0] < (point[0]+data_xrange*alpha)) & (data[:,1] > (point[1]-data_yrange*alpha)) & (data[:,1] < (point[1]+data_yrange*alpha)) & (data[:,2] > (point[2]-data_zrange*alpha)) & (data[:,2] < (point[2]+data_zrange*alpha))]
        
        sort_index = np.argsort(np.apply_along_axis(pythagoras, 1, neighbors, point))
        
        neighborsb = np.zeros((4,4))
        neighborsb[0,:] = neighbors[sort_index[0],:]
        neighborsb[1,:] = neighbors[sort_index[1],:]
        neighborsb[2,:] = neighbors[sort_index[2],:]
        neighborsb[3,:] = neighbors[sort_index[3],:]

        neighbors = neighborsb
        
        #Linear Approach
        A = np.zeros((4,4))
        b = np.zeros((4))
            
        A[:,0] = neighbors[:,0]
        A[:,1] = neighbors[:,1]
        A[:,2] = neighbors[:,2]
        A[:,3] = 1
            
        b[:] = neighbors[:,3]
        
        if np.linalg.det(np.dot(A.T,A)) == 0:
            output =neighbors[0,3]
            
        else:  
            #Solving
            coef = np.linalg.solve(np.dot(A.T,A), np.dot(A.T,b))
            output = point[0]*coef[0] + point[1]*coef[1] + point[2]*coef[2] + coef[3]
            
        return output
    
    return np.apply_along_axis(krigcal, 1, grid, data, data_xrange, data_yrange, data_zrange)
    
#Meshgrid Slicing
def meshgrid_slice(meshgrid, x0, x1, y0, y1, main):
    
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
            quad.GetPointIds().SetId(0, nz*i + 0)
            quad.GetPointIds().SetId(1, nz*i + 1)
            quad.GetPointIds().SetId(2, nz*(i+1) + 1)
            quad.GetPointIds().SetId(3, nz*(i+1) + 0)
            
            mesh.InsertNextCell(quad.GetCellType(), quad.GetPointIds())
            
            #Cell Scalar Value
            value = (section[ nz*i + 0, 3] + section[ nz*i + 1, 3] + section[ nz*(i+1) + 1, 3] + section[ nz*(i+1) + 0, 3])/4
            scalar.InsertNextTuple([np.log10(value)*0.5])
            
    mesh.SetPoints(points)
    mesh.GetCellData().SetScalars(scalar)
    main.mesh_renderer(mesh)
            
    return section
#point = np.array([[795727, 9197285, 0, 0], [795727, 9197285, 0, 0]])
#data  = np.loadtxt("D:\WORK\MODEL_WNI\data\mesh_adj.txt")

#a = xyzv_krigging(data, point)