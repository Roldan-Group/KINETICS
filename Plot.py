'''
    Alberto Roldan

    Plots results from Data.dat

'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter    # Spyder said it is not used


arguments = sys.argv
data = arguments[1]

dataset = np.loadtxt(data, comments="#")                    # import and reads data
max_columns = len(dataset[1])                                # determines the number of columns in data
column_labels = open(data,"r").readlines()[max_columns + 7]     # labels are the firsts raw of data
#column_labels = open(data,"r").readlines()[0]
label = column_labels.split()                                # makes the line a list
label.pop(0)                                                # remove "#" from the list

plt.figure(figsize=(16*2, 16*2), clear=True)        # prepares a A4 figure
k=1
for n in range(max_columns):                             # scan through columns; avoids the N of cluster atoms
    x = dataset[:, n]
    for m in range(max_columns):
        y = dataset[:, m] # for i in range(len(dataset[:, m])) if dataset[i, m] < 0.1]
        plt.subplot(max_columns, max_columns, k)              # each row is a column, +1 for the combined plots below
        plt.plot(x, y, "bo", markersize=2)                     # plots with blue dots (bo)
        plt.xlabel(label[n], fontsize=8)
        plt.ylabel(label[m], fontsize=8)
        plt.tick_params(axis='both', labelrotation=0, labelsize=6)    # customise tick labels
#        plt.axis("tight")
        plt.grid(True)                                      # adds lighter lines across the plots
        k+=1

plt.subplots_adjust(top=0.99, bottom=0.05, left=0.05, right=0.99, hspace=0.7, wspace=0.7)
plt.savefig("RelationPlots.png", dpi=300, orientation='landscape', transparent=True)
#plt.show()                                                 # opens the plot in python







i=0; x_columns=[]; y_columns=[]; z_columns=[];
for parameter in label:
    if parameter.startswith('dist_Mg') is True:                 # as a function of cs_i
#        x_columns.append(i)
        x_columns = i
        i+=1
    elif parameter.startswith('dist_O') is True:            # as a function of the average distance to different sites // perpendicular distance :: Z  --> Zdist
#        y_columns.append(i)
        y_columns = i
        i+=1
    elif parameter.startswith('Eadh') is True:              # as a function of cluster coordination (cc)
        z_columns = i
        i+=1
    else:
        i+=1

x = np.array([dataset[i,x_columns] for i in range(len(dataset[:,x_columns])) if dataset[i,z_columns] < 0.1])
y = np.array([dataset[i,y_columns] for i in range(len(dataset[:,y_columns])) if dataset[i,z_columns] < 0.1])
xlabel = label[x_columns]
ylabel = label[y_columns]

z = np.array([z for z in dataset[:,z_columns] if z < 0.1])
zlabel = label[z_columns]

names = [open(data,"r").readlines()[max_columns + 8 + i].split()[-1] for i in range(len(dataset[:, x_columns]))
          if dataset[i,z_columns] < 1.0]
for i in range(len(z)):
    if z[i] == min(z):
        x_text = x[i]
        y_text = y[i]

figure = plt.figure(figsize=(11.69, 16.53), clear=True)       # prepares a figure
ax = figure.add_subplot(111, projection='3d')

surface = ax.plot_trisurf(x,y,z, cmap="viridis", edgecolor='none', linewidth=0, antialiased=False)
figure.colorbar(surface, shrink=0.5, aspect=8)
ax.scatter3D(x, y, z, cmap="viridis")#, edgecolor='none', linewidth=0, antialiased=False)

ax.text(x_text, y_text, min(z), str("(" + str(x_text) + "," + str(y_text) + ")"))
ax.set_xlabel(xlabel+" /$\\AA$", rotation=0, fontsize=10)
ax.set_ylabel(ylabel+" /$\\AA$", rotation=0, fontsize=10)
ax.set_zlabel(zlabel+" /eV", rotation=0, fontsize=10)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter("%.02f"))
ax.view_init(azim=-135, elev=10)
#plt.savefig("Eadh3D.png", dpi=300, orientation='landscape', transparent=True)
#plt.show()                                                  # opens the plot in python


