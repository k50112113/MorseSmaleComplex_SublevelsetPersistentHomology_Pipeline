import numpy as np

def read_grid_data(filename):
    # file format: the file contains N^D lines (D = 2 or 3), each line represents a sample and the number of columns should be D+1. The first D columns represents the grid coordinates, the D+1 column represents the value in that grid cell. For instance, there are 3,600 lines with 3 columns in the file, therefore, the data represents a 2D grid with 60 x 60 bins.
    #
    # 2D grid data example:
    # X1 Y1 v1
    # X2 Y2 v2
    # X3 Y3 v3
    # ...
    # 
    # 3D grid data example:
    # X1 Y1 Z1 v1
    # X2 Y2 Z2 v2
    # X3 Y3 Z3 v3
    # ...
    #
    #
    # Returns:
    #   data (numpy.ndarray): float, shape = (D x N^D)
    #   dim: int
    #   nbin: int
    #
    # note: use data[i].reshape((nbin, nbin, nbin, ...)) to recover the N-dimensional data for the i^th feature

    data = []
    with open(filename,"r") as fin:
        for aline in fin:
            linelist = aline.strip().split()
            data.append([float(i) for i in linelist])
    data = np.array(data).transpose()
    dim = data.shape[0] - 1
    nbin = round(data.shape[1]**(1/dim))
    return data, dim, nbin

def write_vtk(data, spacing, origin, filename):
    # write VTK file for TTK
    from tvtk.common import configure_input
    from tvtk.api import tvtk
    data = data.astype('f')
    dims = data.shape
    spoints = tvtk.StructuredPoints(spacing=spacing, origin=origin, dimensions=dims)
    s = data.transpose().copy()
    spoints.point_data.scalars = np.ravel(s)
    #spoints.point_data.scalars = np.ravel(data,order='F')
    spoints.point_data.scalars.name = 'scalars'

    # Writes legacy ".vtk" format if filename ends with "vtk", otherwise
    # this will write data using the newer xml-based format.
    # write_data(spoints, filename+'.vti')
    # Uncomment the next two lines to save the dataset to a VTK XML file.
    writer = tvtk.XMLImageDataWriter(file_name = filename)
    # writer = tvtk.XMLStructuredGridWriter(file_name=filename+'.vti')
    configure_input(writer, spoints) # <== will work
    writer.write()

def read_ttk_morse_smale_raw_output(filename):
    output = []
    with open(filename, "rb") as fin:
        fin.readline()
        for aline in fin:
            tmp = ""
            output.append([])
            for achar in aline:
                if achar == ord('\n'):
                    output[-1].append(tmp)
                    tmp = ""
                    break
                elif achar==ord(','):
                    output[-1].append(tmp)
                    tmp = ""
                else:
                    if achar <= 3:
                        tmp = str(achar)
                    else:
                        tmp += chr(achar)
    return output

def read_ttk_critical_pairs(filename):
    dim = 3
    critical_pair_ = []
    pair_type_ = []
    ttk_vertexscalarfield = []
    with open(filename) as fin:
        fin.readline()
        for aline in fin:
            linelist = aline.strip().split(',')
            critical_pair_.append([float(i) for i in linelist[5:5+dim]])
            pair_type_.append(int(linelist[1]))
            ttk_vertexscalarfield.append(int(linelist[0]))
    critical_pair_ = np.array(critical_pair_)
    ttk_vertexscalarfield = np.array(ttk_vertexscalarfield)
    pair_type_ = np.array(pair_type_)
    if critical_pair_[2].all() == 0:
        dim = 2
        critical_pair_ = critical_pair_[:,:2]
    return critical_pair_, pair_type_, dim, ttk_vertexscalarfield

def read_ttk_morse_smale(results_folder, sep_interval=2):
    dim = 3
    ms_critical_point_  = [] # (N x 2) coordinate of each critical point
    ms_critical_type_   = [] # (N x 1) 2: max, 1: saddle, 0: min
    ms_critical_cellid_ = [] # (N x 1) id of each critical point
    ttk_vertexscalarfield = []
    '''
    ms_max_point_       = [] # (N x 1) coordinate of each maximum point
    ms_sad_point_       = [] # (N x 1) coordinate of each saddle point
    ms_max_cellid_      = [] # (N x 1) id of each maximum point
    ms_sad_cellid_      = [] # (N x 1) id of each saddle point
    '''

    ms_1sep_point_         = [] # (M x 2) coordinate of each 1-separation point
    ms_1sep_cellid_        = [] # (M x 1) id of each 1-separation point
    ms_1sep_celldim        = [] # (M x 1) 3: body, 2: face, 1: edge, 0: point
    ms_1sep_sep_type_      = [] # (M x 1) dim-1: accending separation, 0: decending separation
    ms_1sep_startid_       = [] # (K x 1) source id of each 1-separation point
    ms_1sep_id_            = [] # (K x 1) id of each 1-separation point
    ms_1sep_endid_         = [] # (K x 1) destination id of each 1-separation point

    output = read_ttk_morse_smale_raw_output("%s/critical-points.txt"%(results_folder))
    for linelist in output:
        if len(linelist) > 0:
            ms_critical_point_.append([float(i) for i in linelist[6:6+dim]])
            ms_critical_type_.append(int(linelist[0]))
            ms_critical_cellid_.append(int(linelist[1]))
            ttk_vertexscalarfield.append(int(linelist[4]))
       
    ms_critical_point_ = np.array(ms_critical_point_)
    ms_critical_type_ = np.array(ms_critical_type_)
    ms_critical_cellid_ = np.array(ms_critical_cellid_)
    ttk_vertexscalarfield = np.array(ttk_vertexscalarfield)
    if ms_critical_point_[2].all() == 0:
        dim = 2
        ms_critical_point_ = ms_critical_point_[:,:2]

    output = read_ttk_morse_smale_raw_output("%s/1-sep-points.txt"%(results_folder))
    for linelist in output:
        if len(linelist) > 0:
            ms_1sep_point_.append([float(i) for i in linelist[3:3+dim]])
            ms_1sep_cellid_.append(int(linelist[2]))
            ms_1sep_celldim.append(int(linelist[1]))
                
    output = read_ttk_morse_smale_raw_output("%s/1-sep-cells.txt"%(results_folder))
    for linelist in output:
        if len(linelist) > 0:
            ms_1sep_startid_.append(int(linelist[0]))
            ms_1sep_endid_.append(int(linelist[1]))
            ms_1sep_id_.append(int(linelist[2]))
            ms_1sep_sep_type_.append(int(linelist[3]))

    sep_point_ = []                   # (W x m x 2) W = the total number of separation lines (accending + decending)
    sep_celldim = []
    sep_point_start_end_type_ = []    # (W x 2) (source id, destination id, separation type)
    startid = -1
    endid = -1
    sepid = -1
    offset = 0
    for index in range(len(ms_1sep_sep_type_)):
        if sepid != ms_1sep_id_[index]:
            sepid = ms_1sep_id_[index]
            startid = ms_1sep_startid_[index]
            endid = ms_1sep_endid_[index]
            atype = ms_1sep_sep_type_[index]
            sep_point_start_end_type_.append([startid,endid,atype])
            sep_point_.append([])
            sep_celldim.append([])
            offset += 1
        if ms_1sep_celldim[index+offset] >= 0:
            sep_point_[len(sep_point_)-1].append(ms_1sep_point_[index+offset])
            sep_celldim[len(sep_celldim)-1].append(ms_1sep_celldim[index+offset])    

    for index in range(len(sep_point_)-1,-1,-1):
        _,_,atype = sep_point_start_end_type_[index]
        if atype == -1:
            del sep_point_[index]
            del sep_point_start_end_type_[index]

    '''
    ms_max_point_ = np.array(ms_max_point_)
    ms_sad_point_ = np.array(ms_sad_point_)
    ms_max_cellid_ = np.array(ms_max_cellid_)
    ms_sad_cellid_ = np.array(ms_sad_cellid_)
    '''
    ms_1sep_point_ = np.array(ms_1sep_point_)
    ms_1sep_cellid_ = np.array(ms_1sep_cellid_)
    ms_1sep_startid_ = np.array(ms_1sep_startid_)
    ms_1sep_endid_ = np.array(ms_1sep_endid_)

    for sep_index in range(len(sep_point_)):
        sep_point_[sep_index] = np.array(sep_point_[sep_index])
        if sep_point_[sep_index].shape[0] > 2:
            sep_point_[sep_index] = sep_point_[sep_index][1::sep_interval]
    ms_critical_cellid_tmp_ = list(ms_critical_cellid_)
    for sep_index in range(len(sep_point_)):
        sep_point_start_end_type_[sep_index][0] = ms_critical_cellid_tmp_.index(sep_point_start_end_type_[sep_index][0])
        sep_point_start_end_type_[sep_index][1] = ms_critical_cellid_tmp_.index(sep_point_start_end_type_[sep_index][1])

    return ms_critical_point_, ms_critical_type_, sep_point_, sep_point_start_end_type_, sep_celldim, ttk_vertexscalarfield

def get_significant_value(value, significant_digits=1, round_method = "round"):
    import math
    k = 0
    while abs(value*(10**k)) < 10**significant_digits and value != 0:
        k += 1
    if round_method == "round":
        return round(value*(10**k))/(10**k)
    elif round_method == "floor":
        return math.floor(value*(10**k))/(10**k)
    elif round_method == "ceil":
        return math.ceil(value*(10**k))/(10**k)
    
def plot(CV, colorcode=None, cmap=None, elev=13., azim=0):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    dim = CV.shape[1]
    if dim == 3:
        axsize = 0.8
    else:
        axsize = 0.7

    fig = plt.figure(figsize = (10, 8))
    if dim == 3:
        ax = fig.add_axes([.5*(1-axsize), .5*(1-axsize), axsize, axsize], projection='3d')
        ax.clear()
        ax.view_init(elev=elev, azim=azim)
        ax.tick_params(axis="x",which ="major",length=9,width=2,labelsize=20, pad=10)
        ax.tick_params(axis="y",which ="major",length=9,width=2,labelsize=20, pad=10)
        ax.tick_params(axis="x",which ="minor",length=6,width=2,labelsize=20, pad=10)
        ax.tick_params(axis="y",which ="minor",length=6,width=2,labelsize=20, pad=10)
        ax.tick_params(axis="z",which ="major",length=9,width=2,labelsize=20, pad=10)
        ax.tick_params(axis="z",which ="minor",length=6,width=2,labelsize=20, pad=10)
        ax.set_xlabel("$CV_1$",fontsize=26,labelpad=24)
        ax.set_ylabel("$CV_2$",fontsize=26,labelpad=24)
        ax.set_zlabel("$CV_3$",fontsize=26,labelpad=24)
    elif dim == 2:
        ax = fig.add_axes([.7*(1-axsize), .5*(1-axsize), axsize, axsize])
        ax.tick_params(axis="x",which ="major",length=9,width=2,labelsize=24, pad=10)
        ax.tick_params(axis="y",which ="major",length=9,width=2,labelsize=24, pad=10)
        ax.tick_params(axis="x",which ="minor",length=6,width=2,labelsize=24, pad=10)
        ax.tick_params(axis="y",which ="minor",length=6,width=2,labelsize=24, pad=10)
        ax.set_xlabel("$CV_1$",fontsize=26)
        ax.set_ylabel("$CV_2$",fontsize=26)
    
    ax.spines["top"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["right"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.tick_params(axis="both", direction="in", which ="both", top=True, right=True)

    vmin = get_significant_value(colorcode.min(),significant_digits = 1,round_method="ceil")
    vmax = get_significant_value(colorcode.max(),significant_digits = 1,round_method="floor")
    yticklabels = [vmin + i/4*(vmax-vmin) for i in range(5)]
    
    nbin = round(CV.shape[0]**(1/dim))
    if dim == 2:
        
        cv1 = CV[:,0].reshape(nbin,nbin)
        cv2 = CV[:,1].reshape(nbin,nbin)
        colorcode = colorcode.reshape(nbin,nbin)
        
        surf = ax.contourf(cv1,cv2,colorcode, vmin=vmin, vmax=vmax, cmap=cmap, levels=500, extend='both', norm = None)
        for c in surf.collections:
            c.set_edgecolor("face")
        cb = fig.colorbar(surf,ax=ax,pad=0.08,ticks=yticklabels,extend='both')
        cb.ax.tick_params(labelsize=20)
        
        return fig, ax, plt
    
    elif dim == 3:
        cv1 = CV[:,0].reshape(nbin,nbin,nbin)
        cv2 = CV[:,1].reshape(nbin,nbin,nbin)
        cv3 = CV[:,2].reshape(nbin,nbin,nbin)
        colorcode_flat = np.exp(colorcode)
        colorcode = colorcode.reshape(nbin,nbin,nbin)
        min_tmp = colorcode_flat.min()
        for index in range(len(colorcode_flat)):
            if colorcode_flat[index] == min_tmp:
                colorcode_flat[index] = 0

        iso_prob = 0.98
        spacing = CV[:,2][1]-CV[:,2][0]
        left = colorcode_flat.min()
        right = colorcode_flat.max()
        
        mid = (left+right)/2
        sum_tmp = np.sum(colorcode_flat[np.where(colorcode_flat>mid)[0]])*spacing**3
        itrmax = 500
        itr = 0
        while abs(sum_tmp - iso_prob)/iso_prob > 1E-4 and itr < itrmax:
            itr += 1
            if sum_tmp > iso_prob:
                left = mid
            else:
                right = mid
            mid = (left+right)/2
            sum_tmp = np.sum(colorcode_flat[np.where(colorcode_flat>mid)[0]])*spacing**3
        iso_val=np.log(mid)
        
        from skimage import measure
        # verts, faces, norm, val = measure.marching_cubes_lewiner(volume = colorcode, level= iso_val, spacing=(spacing, spacing, spacing)) # deprecated
        verts, faces, norm, val = measure.marching_cubes(volume = colorcode, level= iso_val, spacing=(spacing, spacing, spacing))
        surf = ax.plot_trisurf(verts[:, 0]-1, verts[:,1]-1, faces, verts[:, 2]-1,color='k', lw=0,alpha=0.2)

        return fig, ax, plt
                
        