
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



def i_Axis():
	answer = int(input("Which axis would you like for this data, 1 or 2? One number.\n"))
	while 1 > int(answer) > 2:
		answer = int(input("Select one number, axis 1 or axis 2.\n"))
	return answer

def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
					bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)

def Plot_2D(x, ylabel, data, axis, i_label):
	if axis == 1:
		line = ax1.plot(x, [j for j in data[:, 2]], "k-", label="$Au_{" + str(i_label) + "}$")
		ax1.set_ylabel("$\\theta_{Au_{" + str(sys.argv[i_sys][-5]) + str("}}$ (ML)"), fontsize=14)
		ax1.tick_params(axis='y')
	else:
		ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
		line = ax2.plot(x, [j for j in data[:, 2]], linestyle=iline[1], color=icolour[0], label="$Au_{" + str(i_label) + "}$")
		ax2.set_ylabel("$\\theta_{Au_{" + str(sys.argv[i_sys][-5]) + str("}}$ (ML)"), fontsize=14, color=icolour[0])
		ax2.tick_params(axis='y', color=icolour[0], labelcolor=icolour[0])
		ax2.spines['right'].set_color(icolour[0])
	return line

###############################################################################################################

name = sys.argv[1][:-4]
data = np.loadtxt(sys.argv[1], comments="#")
temp = data[:, 0]
time = data[:, 1]


if len(set(temp)) > 1 and len(set(time)) == 1:
	fig, ax1 = plt.subplots()
	ax1.set_xlabel("Temperature (K)", fontsize=14)
	line_legend = []
	for i_sys in range(1, len(sys.argv)):
		data = np.loadtxt(sys.argv[i_sys], comments="#")
		temp = data[:, 0]
		time = data[:, 1]
		xlabel = "Temperature (K)"
		ylabel = str(sys.argv[i_sys][-5])
		line = Plot_2D(temp, ylabel, data, i_Axis(), sys.argv[i_sys][-5])
		line_legend += line
	ax1.legend(line_legend, [l.get_label() for l in line_legend], loc=0)
	fig.tight_layout()  # otherwise the right y-label is slightly clipped
	plt.ion()
	plt.show()
	SaveFig()



