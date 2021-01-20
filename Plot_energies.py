



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


file_name = sys.argv[1]
output = str(file_name[:-4])
data = np.loadtxt(file_name, comments="#")                            # import and reads data
n = data[:, 0]
x = data[:, 2]
y = data[:, 1]

data2 = np.loadtxt(sys.argv[2], comments="#")                            # import and reads data
x2 = data2[:, 1]
y2 = data2[:, 0]


def trendline(x, a, b):
    return a*x + b  # lineal trend

popt, pcov = curve_fit(trendline, x, y)#, bounds=([d0_eq*0.8, 0, r0_eq*0.8], [d0_eq*1.2, 50, r0_eq*1.2]))
a, b = popt
r2 = 1 - np.sqrt(sum([(y[i] - trendline(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
print(popt, r2)
#label_trend = "$E_{Adh} = %.3f \cdot e^{\minus %.3f (r - %.3f)} \minus 2 \cdot " \
#              "e^{\minus %.3f (r \minus %.3f)}$ ; $R^{2}$= %.2f" %(d_eq, 2*a, r_eq, a, r_eq, r2)
plt.plot(np.linspace(min(x)-np.abs(min(x))*0.1, max(x)+np.abs(max(x))*0.1, 101),
         trendline(np.linspace(min(x)-np.abs(min(x))*0.1, max(x)+np.abs(max(x))*0.1, 101), *popt), "b:")# , label=label_trend)

plt.plot(x, y, "ko")
plt.plot(x2, y2, "ro", fillstyle="none")
#plt.plot(np.linspace(min(x), max(x)), np.linspace(min(x), max(x)), "b:")

for i, m in enumerate(n):
    plt.annotate("$Au_{" + str(int(m)) + "} + 1$", xy=(x[i], y[i]), xycoords="data",# size=14,
                        xytext=(x[i]-0.07, y[i]+0.01), textcoords="data",
#                        arrowprops=dict(arrowstyle="->", fc="0.6"),
                        horizontalalignment="right", verticalalignment="top")
plt.annotate("$Au_{2} + 1$", xy=(x2[0], y2[0]), xycoords="data",# size=14,
                    xytext=(x2[0]-0.03, y2[0]+0.01), textcoords="data",
                    horizontalalignment="right", verticalalignment="top")
plt.annotate("$Au_{4} + 1$", xy=(x2[1], y2[1]), xycoords="data",# size=14,
                    xytext=(x2[1]-0.03, y2[1]+0.01), textcoords="data",
                    horizontalalignment="right", verticalalignment="top")
plt.annotate("$Au_{6} + 1$", xy=(x2[2], y2[2]), xycoords="data",# size=14,
                    xytext=(x2[2]-0.03, y2[2]+0.01), textcoords="data",
                    horizontalalignment="right", verticalalignment="top")


plt.xlabel("$E_{r}$ (eV)", size=14)
#plt.xlabel("$Au_{n}$", size=14)
plt.ylabel("$E_{a}$ (eV)", size=14)
#plt.xticks(np.arange(min(x), max(x)+1, 1))   # Xmin,Xmax,Xstep
#plt.xlim([min(x)-1.5, max(x)+0.5])
#plt.yticks([])
#plt.ylim([min(y)-np.abs(min(y))*0.1, max(y)+np.abs(max(y))*0.1])
#plt.tick_params(axis='both', labelrotation=0, labelsize=12)    		# custimise tick labels
#plt.legend(loc='best')
#plt.title(output)
plt.tight_layout()
plt.savefig(output + ".svg", figsize=(12, 10), dpi=300, orientation='landscape', transparent=True)
plt.show()


