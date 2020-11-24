'''
    Alberto Roldan

    Plots results from Data.dat

'''

import sys
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d                                   # Spyder said it is not used

from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter    # Spyder said it is not used

try:
    data = sys.argv[1]
    Xcolumn = int(sys.argv[2])
    Ycolumn = int(sys.argv[3])
    Zcolumn = int(sys.argv[4])
except IOError:
    raise Exception("   Provide 4 arguments: DataSet and x y z columns number (starting from zero)")

DataSet = np.loadtxt(data, comments="#")                    # import and reads data
maxColumns=len(DataSet[1])                                  # determines the number of columns in data
columnLabels=open(data,"r").readlines()[0]                  # labels are the firts raw of data
label=columnLabels.split()                                  # makes the line a list
label.pop(0)                                                # remove "#" from the list

xlabel = label[Xcolumn]
ylabel = label[Ycolumn]
zlabel = label[Zcolumn]

X=[]; Y=[]; Z=[]
for i in range(len(DataSet[:,Xcolumn])):
    if DataSet[i,Xcolumn] > 0 and DataSet[i,Ycolumn] > 0 and DataSet[i,Zcolumn] < 0:
        X.append(DataSet[i,Xcolumn]/DataSet[i,0])
        Y.append(DataSet[i,Ycolumn])
        Z.append(DataSet[i,Zcolumn])


figure = plt.figure(figsize=(11.69,16.53),clear=True)       # prepares a figure
ax = figure.gca(projection='3d')
ax.scatter3D(X, Y, Z, cmap=cm.coolwarm)
ax.set_xlabel(xlabel, rotation=0, fontsize=10)
ax.set_ylabel(ylabel, rotation=0, fontsize=10)
ax.set_zlabel(zlabel, rotation=0, fontsize=10)

plt.show()

ax.plot_trisurf(X, Y, Z, cmap="viridis", linewidth=0, antialiased=False)
ax.view_init(azim=-135, elev=10)
plt.savefig(xlabel[:4]+"_"+ylabel[:4]+"_Norm"+zlabel[:4]+".png", dpi=300, orientation='landscape',transparent=True)

#plt.show()                                                  # opens the plot in python








