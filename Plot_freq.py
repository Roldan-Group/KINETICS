



import sys
import numpy as np
import matplotlib.pyplot as plt


file_name = [str(i) for i in sys.argv[1:]]
output = file_name[0][:-4]
for i in file_name[1:]:
    output = str(output + "_" + i[:-4])


Xmin = 0
Xmax = 5000
Xstep = 500
Ymin = 0
Ymax = 5000
for name in file_name:
    DataSet = np.loadtxt(name, comments="#")                            # import and reads data
    X = [ DataSet[i,0] for i in range(len(DataSet))]
    if min(X) > Xmin:
        Xmin = min(X)
    if max(X) < Xmax:
        Xmax = max(X)


for n, name in enumerate(file_name):
    DataSet = np.loadtxt(name, comments="#")                    	# import and reads data
    X = [ DataSet[i,0] for i in range(len(DataSet)) if DataSet[i,0] > Xmin and DataSet[i,0] < Xmax]

    if name[-3:] == "dat":
        y0 = [ DataSet[i,1] for i in range(len(DataSet)) if DataSet[i,0] > Xmin and DataSet[i,0] < Xmax]
        Y = [ (i)/abs(max(y0)) for i in y0 ]
    elif name[-3:] == "txt":
        y0 = [ (DataSet[i,1] - max(DataSet[:,1]))*(-1)  for i in range(len(DataSet)) if DataSet[i,0] > Xmin and DataSet[i,0] < Xmax]
        Y = [ (i)/abs(max(y0)) for i in y0 ]

    if len(file_name) < 4:
        if n == 0:
            plt.plot(X, Y, "k-", lw=1.5, label=str(name[:-4]))
        elif n == 1:
            plt.plot(X, Y, "b--", lw=1.5, label=str(name[:-4]))
        elif n == 2:
            plt.plot(X, Y, "r:", lw=1.5, label=str(name[:-4]))    
        elif n == 3:
            plt.plot(X, Y, "g-.", lw=1.5, label=str(name[:-4]))
    else:
        plt.plot(X, Y,  lw=1.5, label=str(name[:-4]))


plt.xlabel(r'ω /$cm^{-1}$',fontsize=14)
plt.ylabel("Intensity /a.u.", fontsize=14)
plt.xticks(np.arange(Xmin,Xmax,Xstep))   # Xmin,Xmax,Xstep
plt.xlim([Xmin,Xmax])
plt.yticks([])
plt.tick_params(axis='both',labelrotation=0,labelsize=12)    		# custimise tick labels
plt.legend(loc='best')

#plt.annotate("Bulk", xy=(12,0.99), xycoords="data", size=14,
#                     xytext=(11.75,0.85), textcoords="data",
#                     arrowprops=dict(arrowstyle="->", fc="0.6"),
#                     horizontalalignment="right", verticalalignment="top")


plt.tight_layout()
plt.savefig("G_comparison.png", dpi=300, orientation='portrait',transparent=True)
plt.show()


