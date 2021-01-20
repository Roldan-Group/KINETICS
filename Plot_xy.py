



import sys
import numpy as np
import matplotlib.pyplot as plt


file_name = sys.argv[1]
output = str(file_name[:-4])
data = np.loadtxt(file_name, comments="#")                            # import and reads data
x = data[:, 0]
y = data[:, 1]

plt.plot(x, y, "ko", ms=1, label=output)
plt.plot(np.linspace(min(x), max(x)), np.linspace(min(x), max(x)), "b:")
#plt.xlabel(r'ω /$cm^{-1}$',fontsize=14)

#plt.ylabel("Intensity /a.u.", fontsize=14)
#plt.xticks(np.arange(Xmin,Xmax,Xstep))   # Xmin,Xmax,Xstep
#plt.xlim([Xmin,Xmax])
#plt.yticks([])
#plt.tick_params(axis='both', labelrotation=0, labelsize=12)    		# custimise tick labels
plt.legend(loc='best')
#plt.annotate("Bulk", xy=(12,0.99), xycoords="data", size=14,
#                     xytext=(11.75,0.85), textcoords="data",
#                     arrowprops=dict(arrowstyle="->", fc="0.6"),
#                     horizontalalignment="right", verticalalignment="top")
plt.title(output)
plt.tight_layout()
#plt.savefig(output + ".png", dpi=300, orientation='landscape', transparent=True)
plt.show()


