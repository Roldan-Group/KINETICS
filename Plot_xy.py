



import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


file_name = sys.argv[1]
output = str(file_name[:-4])
data = np.loadtxt(file_name, comments="#")                            # import and reads data
#x = data[:, 1]
y = data[:, 0]
x = np.arange(0, len(y), 1)+2

x_range = np.abs(max(x) - min(x))
y_range = np.abs(max(y) - min(y))

def function(x, a, b, c, d):
    return a + b * np.exp(-((np.log(x)-np.log(c))**2)/(2*d))

def fitting(x, y):
    popt, pcov = curve_fit(function, x, y)  #, bounds=([0, 0, 0, 0], [1, 1, 1, 1]))
    return popt

popt = fitting(x, y)
root_squared = 1-np.sqrt(sum([(y[i] - function(x[i], *popt))**2 for i in range(len(x))])/sum(i*i for i in y))
print(popt, root_squared)
popt = np.array([0.19404975, 0.81379323, 1.13040037, 0.25108869])
a, b, c, d = popt
root_squared = 0.916335584466851
print(popt, root_squared)
label_trend = "$E_{a} = %.3f \plus %.3f \cdot e^{\minus \\frac{(log(n) \minus log(%.3f))^{2}}{2 \cdot %.3f}}$" \
                  " ; $R^{2}$= %.2f" % (a, b, c, d, root_squared)

plt.plot(x, y, "ko")    #, ms=1, label=output)
plt.plot(x, y, "k:")
plt.plot(np.linspace(min(x)-x_range*0.01, max(x)+x_range*0.05),
         function(np.linspace(min(x)-x_range*0.01, max(x)+x_range*0.05), *popt), "b--", lw=0.5, label=label_trend)

plt.xlabel("$Au_{n}$", fontsize=14)    #$E_{r}$ $(eV)$', fontsize=14)
plt.ylabel("$E_{a}$ $(eV)$", fontsize=14)
plt.xticks(np.arange(min(x), max(x)+1))   # Xmin,Xmax,Xstep
#plt.xlim([min(x)-x_range*0.15, max(x)+x_range*0.15])
#plt.yticks([])
plt.ylim([min(y)-y_range*0.15, max(y)+y_range*0.15])
plt.tick_params(axis='both', labelrotation=0, labelsize=12)    		# custimise tick labels

for i in range(len(x)):
    plt.annotate("$Au_{"+str(i+1)+"\plus 1}$", xy=(x[i], y[i]), xycoords="data",  # size=14,
                     xytext=(x[i]-x_range*0.05, y[i]-y_range*0.08), textcoords="data" #,
#                     arrowprops=dict(arrowstyle="->", fc="0.6"),
#                     horizontalalignment="right", verticalalignment="top"
                 )

plt.legend(loc='best')
#plt.title(output)
plt.tight_layout()
#plt.savefig(output + ".png", dpi=300, orientation='landscape', transparent=True)
plt.show()


