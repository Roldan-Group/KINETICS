



import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import splev, splrep


def SaveFig():
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg", figsize=(12, 10), clear=True,
					bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)


name = sys.argv[1][:-4]
data = np.loadtxt(sys.argv[1], comments="#")

temp = data[:, 0]
time = data[:, 1]
if len(set(temp)) > 1:
    x = temp
    plt.xlabel("Temperature (K)", fontsize=14)
else:
    x = time
    plt.xlabel("time (s)", fontsize=14)
print(len(data[:,2]))
for i in range(2, len(data[0])-2):
    plt.plot(x, [np.log(j) for j in data[:, i]], "o")
    y = interpolate.interp1d(x, [np.log(j) for j in data[:, i]], kind="linear")
    x2 = np.linspace(min(x), max(x), 15)
    y2 = y(x2)
#    plt.plot(np.linspace(min(x), max(x), 5), y(np.linspace(min(x), max(x), 5)))
# splines
    spl = splrep(x2, y2, k=3)
    y = splev(np.linspace(min(x), 510, 500), spl)
    plt.plot(np.linspace(min(x), 510, 500), y)



plt.ylabel("$ln(\\theta_{i})$", fontsize=14)
#plt.xticks(np.arange(0,Xmax, 0.25))
plt.xlim([0, max(x)*1.1])
plt.ylim([-30, -10])
#plt.tick_params(axis='both', labelrotation=0, labelsize=16)               # custimise tick labels
#plt.grid(True)
#plt.legend(loc='best')
#plt.title("$TM_{" + str(int(data[0,0])) + "}$ $dist_{" + str(sys.argv[2]) + "}$")
plt.ion()
plt.show()
SaveFig()


