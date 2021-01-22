
import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import interpolate
from scipy.interpolate import splev, splrep, griddata



icolour = ["b", "r", "k", "g", "c", "m", "y", "grey", "olive", "brown", "pink"] ## n=11
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
iline = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)),  (0, (3, 1, 1, 1, 1, 1))]
#iline = ["solid", 'dashed', 'dashdot', "dotted", 'dashdotdotted', 'densely dashed',  'densely dashdotted', 'densely dashdotdotted']

def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
					bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)


def Plot_2D(x, data):
	for i in range(2, len(data[0])):
# ORIGINAL
#		plt.plot(x, [j for j in data[:, i]], "o")
# INTERPOLATION
#	    y = interpolate.interp1d(x, [j for j in data[:, i]], kind="linear")
#   	 plt.plot(np.linspace(min(x), max(x), 150), y(np.linspace(min(x), max(x), 150)))
# SPLINES
		spl = splrep(x, [j for j in data[:, i]], k=3)
		y = splev(np.linspace(min(x), max(x), 150), spl)
		plt.plot(np.linspace(min(x), max(x), 150), y, linestyle=iline[i-2], color=icolour[i-2], label="$Au_{" + str(i) + "}$")
# ANNOTATE
	plt.annotate("$\\theta_{Au_{1}}^{t=0}=0.5ML$", xy=(40, 0.001), xycoords="data", color=icolour[2], size=14,
				 xytext=(40, 0.001), textcoords="data", horizontalalignment="right", verticalalignment="top")

	plt.ylabel("$\\theta_{Au_{i}}$ (ML)", fontsize=14)
#	plt.xticks(np.arange(0,Xmax, 0.25))
	plt.xlim([min(x), max(x)*1.1])
#	plt.ylim([min(y), 0.1])
#	plt.tick_params(axis='both', labelrotation=0, labelsize=16)               # custimise tick labels
#	plt.grid(True)
	plt.legend(loc='best')
#	plt.title("$TM_{" + str(int(data[0,0])) + "}$ $dist_{" + str(sys.argv[2]) + "}$")
	plt.tight_layout()
	plt.ion()
	plt.show()
	SaveFig()


def Plot_3D(x, y, data):
	z = data[:, 2]

	figure = plt.figure(figsize=(12, 10), clear=True)
	ax = figure.add_subplot(111, projection='3d')

	z = list(map(float, z))
	grid_x, grid_y = np.mgrid[min(x):max(x):200j, min(y):max(y):200j]
	grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')	# better "cubic" but it generates Nan values
	surface = ax.plot_surface(grid_x, grid_y, grid_z, cmap="plasma", alpha=0.7, edgecolor='none', linewidth=0, antialiased=False)
	surface.set_clim(min(z), max(z))

	figure.colorbar(surface, shrink=0.5, aspect=8)
	ax.set_xlabel(r'Temperature (K)', rotation=0, fontsize=16)
	ax.set_ylabel(r'time (s)', rotation=0, fontsize=16)
	ax.set_zlabel("$\\theta_{Au_{10}}$ (ML)", fontsize=16) #  % labels_species[int(plot_species[0])], fontsize=16)
#	ax.set_zlim([0, 0.11])
#	ax.zaxis.set_ticklabels([])
#	ax.zaxis.set_major_formatter(FormatStrFormatter("%.2g"))
	ax.view_init(azim=55, elev=10)
	plt.ion()
	plt.show()
	SaveFig()




###############################################################################################################

name = sys.argv[1][:-4]
data = np.loadtxt(sys.argv[1], comments="#")

temp = data[:, 0]
time = data[:, 1]
if len(set(temp)) > 1 and len(set(time)) > 1:
	plt.xlabel("Temperature (K)", fontsize=14)
	plt.xlabel("time (s)", fontsize=14)
	Plot_3D(temp, time, data)
elif len(set(temp)) > 1 and len(set(time)) == 0:
	plt.xlabel("Temperature (K)", fontsize=14)
	Plot_2D(temp, data)
elif len(set(temp)) == 0 and len(set(time)) > 1:
	plt.xlabel("time (s)", fontsize=14)
	Plot_2D(time, data)

