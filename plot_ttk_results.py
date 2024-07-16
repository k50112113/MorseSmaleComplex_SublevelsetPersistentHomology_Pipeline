import sys
from _utils import read_grid_data, read_ttk_morse_smale, plot

grid_data_filename = sys.argv[1]
output_folder_path = sys.argv[2]

# plot density contour
data, dim, nbin = read_grid_data(grid_data_filename)
CV = data[:dim].transpose()
colorcode = data[dim]
fig, ax, plt = plot(CV, colorcode=colorcode,cmap='jet',elev=14,azim=37)
# ax.clear()

# plot ttk path
if dim == 2:
    sep_linewidth = [2,5]
    sep_linestyle = ['-','-']
    sep_color = ['red','black']
    cmap = 'jet'
elif dim == 3:
    sep_linewidth = [0,0,2]
    sep_linestyle = ['-','-',':']
    sep_color = ['red','blue','black']
    cmap = 'coolwarm'
ms_critical_point_, \
ms_critical_type_, \
sep_point_, \
sep_point_start_end_type_, \
sep_celldim, \
ttk_vertexscalarfield = read_ttk_morse_smale(results_folder = output_folder_path)
for i_sep in range(len(sep_point_)):
    sep_start, sep_end, sep_type = sep_point_start_end_type_[i_sep]
    sep_point = sep_point_[i_sep].transpose()
    ax.plot(*sep_point,color=sep_color[sep_type],lw=sep_linewidth[sep_type],ls=sep_linestyle[sep_type],zorder=10)
ms_critical_point = ms_critical_point_.transpose()
ax.scatter(*ms_critical_point,s=300,linewidths=2,edgecolors='white',c=ms_critical_type_,cmap=cmap,zorder=20)

plt.savefig("%s/plot.png"%(output_folder_path), transparent = False)

# plt.show()


