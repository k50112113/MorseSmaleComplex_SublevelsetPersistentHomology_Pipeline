import sys
from _utils import read_grid_data, write_vtk

grid_data_filename = sys.argv[1]
vtk_filename = sys.argv[2]
data, dim, nbin = read_grid_data(grid_data_filename)
if dim == 2:
    x = data[0].reshape((nbin,nbin))
    y = data[1].reshape((nbin,nbin))
    data = data[2].reshape((nbin,nbin,1))
    spacing = [x[1][0]-x[0][0],y[0][1]-y[0][0],1]
    origin = [-1,-1,0]
elif dim == 3:
    x = data[0].reshape((nbin,nbin,nbin))
    y = data[1].reshape((nbin,nbin,nbin))
    z = data[2].reshape((nbin,nbin,nbin))
    data = data[3].reshape((nbin,nbin,nbin))
    spacing = [x[1][0][0]-x[0][0][0],y[0][1][0]-y[0][0][0],z[0][0][1]-z[0][0][0]]
    origin = [-1,-1,-1]
else:
    print("TTK only supports 2 or 3 dimensions currently, no outputs are created")
    exit()

write_vtk(data, spacing, origin, filename = vtk_filename)