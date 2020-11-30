



import sys
import numpy as np
import matplotlib.pyplot as plt


file_name = [str(i) for i in sys.argv[1:]]
output = file_name[0][:-4]
for i in file_name[1:]:
    output = str(output + "_" + i[:-4])


xmin = 500
xmax = 5000
xstep = 500
ymin = 0
ymax = 5000
for name in file_name:
    data = np.loadtxt(name, comments="#")                            # import and reads data
    x = [data[i, 0] for i in range(len(data))]
    if min(x) > xmin:
        xmin = min(x)
    if max(x) < xmax:
        xmax = max(x)


for n, name in enumerate(file_name):
    data = np.loadtxt(name, comments="#")                    	# import and reads data
    x = [data[i, 0] for i in range(len(data)) if data[i, 0] > xmin and data[i, 0] < xmax]

    if name[-3:] == "dat":
        y0 = [data[i, 1] for i in range(len(data)) if data[i, 0] > xmin and data[i, 0] < xmax]
        y = [i/abs(max(y0)) for i in y0]
    elif name[-3:] == "txt":
        y0 = [(data[i, 1] - max(data[:,1]))*(-1) for i in range(len(data)) if data[i,0] > xmin and data[i, 0] < xmax]
        y = [i/abs(max(y0)) for i in y0]

    if len(file_name) < 4:
        if n == 0:
            plt.plot(x, y, "k-", lw=1.5, label=str(name[:-4]))
        elif n == 1:
            plt.plot(x, y, "b--", lw=1.5, label=str(name[:-4]))
        elif n == 2:
            plt.plot(x, y, "r:", lw=1.5, label=str(name[:-4]))
        elif n == 3:
            plt.plot(x, y, "g-.", lw=1.5, label=str(name[:-4]))
    else:
        plt.plot(x, y,  lw=1.5, label=str(name[:-4]))


plt.xlabel(r'ω /$cm^{-1}$', fontsize=14)
plt.ylabel("Intensity /a.u.", fontsize=14)
plt.xticks(np.arange(xmin, xmax, xstep))   # Xmin,Xmax,Xstep
plt.xlim([xmin, xmax])
plt.yticks([])
plt.tick_params(axis='both', labelrotation=0, labelsize=12)    		# custimise tick labels
plt.legend(loc='best')

#plt.annotate("Bulk", xy=(12,0.99), xycoords="data", size=14,
#                     xytext=(11.75,0.85), textcoords="data",
#                     arrowprops=dict(arrowstyle="->", fc="0.6"),
#                     horizontalalignment="right", verticalalignment="top")
plt.tight_layout()
plt.savefig("Frequencies.png", dpi=300, orientation='landscape', transparent=True)
plt.show()


