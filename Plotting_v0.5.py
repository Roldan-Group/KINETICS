''' 
	Plots XY columns from KINETICS ouputs
	read input.mk.in to know the species
		by Alberto Roldan
			08/2020


	USAGE: Plot_v.py file.mk.in
		 : Plot_v.py file.mk.in  n n 1 0 0 all 1.0 0


	** molecues and surface_species are only those systems with NO negative frequencies

'''


import sys, os
import subprocess
import numpy as np
from scipy.interpolate import griddata
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from termcolor import colored

Joules = 1.6021766208e-19
Av = 6.0221409e23

# Plot styles
imarker = ['o',"s","v","H","X","*","^","<",">","p","P","h","1","2","3","4","d","+"]
icolour = ["k", "b", "r", "g", "c", "m"]


def Line_colour():
	icolour = ["black", "blue", "red", "green", "cyan", "magenta", "yellow"]
	print("Which line colour would you like? One number.\n")
	n = 0
	for i, c in enumerate(icolour):
		if i/3 < 1*(1+n):
			print(" [%d] %-15s" % (i, c), end='', flush=True)
		else:
			n += 1
			print(" [%d] %-15s" % (i, c))
	answer = input("\n")
	while int(answer) > len(icolour):
		print("Select a single line colour.")
		answer = input("Which line colour would you like? One number.\n")

	return str(icolour[int(answer)])

# Library containing the parameters fitting the Shomate equations from NIST
def Library (system):
	shomate_equation = {
		"O2":    [100, 700,  31.32234,  -20.23531,  57.86644,  -36.50624, -0.007374,   -8.903471, 246.7945,    0.0   ],
		"H2":    [298, 1000, 33.066178, -11.363417, 11.432816, -2.772874, -0.158558,   -9.980797, 172.707974,  0.0   ],
		"H2O":   [500, 1700, 30.09200,    6.832514,  6.793435, -2.534480,  0.082139, -250.8810,   223.3967, -241.8264],
		"NH3":   [298, 1400, 19.99563,   49.77119, -15.37599,   1.921168,  0.189174,  -53.30667,  203.8591,  -45.89806],
		"N2H4":  [800, 2000, 35.18240,   96.05260, -40.50130,   6.668070, -0.874233,   77.99150,  249.4250,   95.35340],
		"NH2NH2":[800, 2000, 35.18240,   96.05260, -40.50130,   6.668070, -0.874233,   77.99150,  249.4250,   95.35340],
		"CO":    [298, 1300, 25.56759,    6.096130,  4.054656, -2.671301,  0.131021,  -118.0089,  227.3665, -110.5271],
		"CO2":   [298, 1200, 24.99735,   55.18696,  -33.69137,  7.948387, -0.136638,  -403.6075,  228.2431, -393.5224],
						 }
	return shomate_equation[system]

def Thermodynamics(syst):
	# define the Thermodynamic functions
	def function_Cp(T, a, b, c, d, e, f, g, h):
		T = T / 1000
		return a + b * T + c * T ** 2 + d * T ** 3 + e / T ** 2

	def function_S(T, a, b, c, d, e, f, g, h):
		T = T / 1000
		return a * np.log(T) + b * T + 1 / 2 * c * T ** 2 + 1 / 3 * d * T ** 3 - e / (2 * T ** 2) + g

	def function_H(T, a, b, c, d, e, f, g, h):
		T = T / 1000
		return (a * T + 1 / 2 * b * T ** 2 + 1 / 3 * c * T ** 3 + 1 / 4 * d * T ** 4 - e / T + f - h) * 1000  # from KJ/mol to J/mol

	def function_G(T, a, b, c, d, e, f, g, h):
		T = T / 1000
		return (a * T + 1 / 2 * b * T ** 2 + 1 / 3 * c * T ** 3 + 1 / 4 * d * T ** 4 - e / T + f - h) * 1000 - \
			   (T * 1000) * (a * np.log(T) + b * T + 1 / 2 * c * T ** 2 + 1 / 3 * d * T ** 3 - e / (2 * T ** 2) + g)


# get parameters from NIST in Library to compare
	try:
		X2min, X2max, a,b,c,d,e,f,g,h = Library(syst)
	except:
		X2min = 0
		pass

	H0 = 0                                   # reference enthalpy at 298 K
	variables_thermo = [ "Cp", "S", "H", "G" ]
	for var in variables_thermo:
		datainput = "./THERMODYNAMICS/DATA/" + str(syst) + "/" + str(var) + str(syst) + ".dat"
		data = np.loadtxt(datainput, comments="#")
		X = data[:,0]
		Y = data[:,1]*Joules*Av
		if var is "H":
			for i in range(len(X)):
				if X[i] == 298:
					H0 = Y[i]
				elif X[i] > 298 and H0 == 0:
					H0 = Y[i]
		if var is "Cp" and 0 < X2min < max(X):
			plt.plot(np.linspace(X2min, max(X), 101),
					 function_Cp(np.linspace(X2min, max(X), 101), a,b,c,d,e,f,g,h), "k", label= "Shomate")
		if var is "S" and 0 < X2min < max(X):
			plt.plot(np.linspace(X2min, max(X), 101),
					 function_S(np.linspace(X2min, max(X), 101), a,b,c,d,e,f,g,h), "k", label= "Shomate")
		if var is "H" and 0 < X2min < max(X):
			plt.plot(np.linspace(X2min, max(X), 101),
					 function_H(np.linspace(X2min, max(X), 101), a,b,c,d,e,f,g,h)/1000, "k", label= "Shomate")
		if var is "G" and 0 < X2min < max(X):
			plt.plot(np.linspace(X2min, max(X), 101),
					 function_G(np.linspace(X2min, max(X), 101), a,b,c,d,e,f,g,h)/1000, "k", label= "Shomate")

		if var is "H" or var is "G":
			plt.ylabel(var + " ($KJ \cdot mol^{-1}$)")
			plt.plot(X, (Y-H0)/1000, "bo-", label="Computed")
		else:
			plt.ylabel(var + " ($J \cdot mol^{-1}$)")
			plt.plot(X, Y, "bo-", label="Computed")

		plt.xlabel(r'Temperature /K')
		plt.title(syst)
		plt.subplots_adjust(left=0.15, right=0.8, top=0.9, bottom=0.1)
		plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
#		plt.savefig("./THERMODYNAMICS/PLOTS/" + syst + "/" + var + syst + ".svg",
#				dpi=300, orientation='portrait', transparent=True)
		plt.show()


def Reaction_Constants(processes):
	sticky = []
	arrhenius = []
	r_constant = []
	temperature = np.loadtxt("./KINETICS/PROCESS/ReactionParameters1.dat", comments="#")[:,0]
	for n, process in enumerate(processes):
		datainput = "./KINETICS/PROCESS/ReactionParameters" + str(n+1) + ".dat"
		data = np.loadtxt(datainput, comments="#")
		if process == "A":
			sticky.append([n, data[:,1]])
			arrhenius.append([n, data[:,2]])
			r_constant.append([n, data[:,3]])
		elif process == "R" or process == "D":
			arrhenius.append([n, data[:,1]])
			r_constant.append([n, data[:,2]])

	plt.xlabel(r'Temperature (K)')
	for n, Y in sticky:
		plt.plot(temperature, np.log(Y), label=str("process" + str(n+1)))
		plt.ylabel("$log(\\sigma)$")
	plt.subplots_adjust(left=0.15, right=0.8, top=0.9, bottom=0.1)
	plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
#	plt.savefig("./KINETICS/PROCESS/Sticky.svg", dpi=300, orientation='landscape', transparent=True)
	plt.show()
	for n, Y in arrhenius:
		plt.plot(temperature, np.log(Y), label=str("process" + str(n+1)))
		plt.ylabel("$log(A)$")
	plt.subplots_adjust(left=0.15, right=0.8, top=0.9, bottom=0.1)
	plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
#	plt.savefig("./KINETICS/PROCESS/Arrhenius.svg", dpi=300, orientation='landscape', transparent=True)
	plt.show()
	for n, Y in r_constant:
		plt.plot(temperature, np.log(Y), label=str("process" + str(n+1)))
		plt.ylabel("$log(K)$")
	plt.subplots_adjust(left=0.15, right=0.8, top=0.9, bottom=0.1)
	plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
#	plt.savefig("./KINETICS/PROCESS/RConstants.svg", dpi=300, orientation='landscape', transparent=True)
	plt.show()


def	Species_2D(experiment, labels_species, conditions, species, time_range,
					   plot_v, plot_ph, plot_temp, plot_time, plot_species, plot_ini_species, ir):
# plot_species is plot_systems
	y = {}
	if type(plot_temp) is list:
		x = plot_temp
		if type(plot_time) is list:
			plot_time = float(plot_time[-1])
		name = str("V=%.2f_pH=%.2f_t=%.2f_%s" % (plot_v, plot_ph, plot_time, plot_ini_species))
		plt.xlabel(r'Temperature (K)', size="16")
	elif type(plot_time) is list:
		x = plot_time
		if type(plot_temp) is list:
			plot_temp = float(plot_temp[-1])
		name = str("V=%.2f_pH=%.2f_T=%.2f_%s" % (plot_v, plot_ph, plot_temp, plot_ini_species))
		plt.xlabel(r'time (s)', size="16")
	for i in range(len(conditions)):
		v, ph, temp, time = conditions[i]
		if v == plot_v and ph == plot_ph:
			if experiment != "TPR":
				if time == min(time_range) and experiment == "const_TEMP" or\
					time == min(time_range) and temp == min(temp_range) and experiment == "variable_TEMP":
					label_comment = ''
					for j in range(len(labels_species[:-1])):
						if species[i, j] != 0:
							label_comment += labels_species[j] + "=" + str(round(species[i, j], 2))
			else:
				if temp == min(temp_range[:-1]):
					label_comment = ''
					if type(plot_species) is int:
						label_comment = labels_species[plot_species] + "=" + str(round(species[i, plot_species], 2))
					else:
						for j in range(len(plot_species)):
							if species[i, j] != 0:
								label_comment += labels_species[j] + "=" + str(round(species[i, j], 2))
			if label_comment in plot_ini_species:
# multiple temperatures, one time
				if type(plot_temp) is list and temp in plot_temp:
					if time == plot_time:
						if type(plot_species) is list:
							for spec in plot_species:
								spec = int(spec)
								if labels_species[spec] not in y:
									y[labels_species[spec]] = [species[i, spec]]
								else:
									y[labels_species[spec]].append(species[i, spec])
						else:
							spec = int(plot_species)
							if labels_species[spec] not in y:
								y[labels_species[spec]] = [species[i, spec]]
							else:
								y[labels_species[spec]].append(species[i, spec])
# one temperature, multiple times
				elif type(plot_temp) is float and temp == plot_temp:
					if time in plot_time:
						if type(plot_species) is list:
							for spec in plot_species:
								spec = int(spec)
								if labels_species[spec] not in y:
									y[labels_species[spec]] = [species[i, spec]]
								else:
									y[labels_species[spec]].append(species[i, spec])
						else:
							spec = int(plot_species)
							if labels_species[spec] not in y:
								y[labels_species[spec]] = [species[i, spec]]
							else:
								y[labels_species[spec]].append(species[i, spec])
	if experiment != "TPR":
		plt.ylabel("$\\theta_{i}$ (ML)", size="16")
		plt.yticks(np.arange(0, 1, step=0.1))
		plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.15)
		for spec in y:
			if len(y) > 1:
				plt.plot(x, y[spec], ls="-", label=spec + "$_{" + plot_ini_species + "}$")
			else:
				line_color = Line_colour()
				plt.plot(x, y[spec], ls="-", color=line_color, label=spec + "$_{" + plot_ini_species + "}$")
	else:
		plt.ylabel("$\\frac{\\delta P_{i}}{\\delta T}$ (a.u.)", size="16")
		plt.yticks([])
		plt.subplots_adjust(left=0.2, right=0.75, top=0.9, bottom=0.15)
		for spec in y:
			y_tpr = []
			comment = str(spec) + "= " + str(round(y[spec][0], 2))
			n = 1
			for j in range(len(y[spec])+1):
				if j < len(x) * n:
					y_tpr.append((y[spec][j]-sum(y_tpr)))
				else:
					plt.plot(x, y_tpr, ls="-", label=comment)
					n += 1
					if j < len(y[spec]):
						y_tpr = []
						y_tpr.append(0)
						comment = str(spec) + "= " + str(round(y[spec][j], 2))
#	plt.ylim(-0.01, 1.01)
#	plt.title(name)
	plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
#	plt.savefig("./KINETICS/PLOTS/"+experiment+"/"+name+".svg",
#				bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)
	plt.ion()
	plt.show()
	SaveFig(experiment)
	Extract_numeric_data(experiment, labels_species, conditions, species, time_range,
						 plot_v, plot_ph, plot_temp, plot_time, plot_species, plot_ini_species)
	if experiment != "TPR":
		ir_z_spectra = []
		frequencies = ir[spec][:, 0]
		i_intensity = {}
		for spec in y:
			i_intensity[spec] = []
			if len(y[spec]) > 50:
# 50 lines in the IR plot
				new_x = [x[i] for i in range(0, 50, int(len(y[spec])/50))]
				y[spec] = [y[spec][i] for i in range(0, 50, int(len(y[spec])/50))]
			else:
				new_x = x
			for j in range(len(new_x)):
					i_intensity[spec].append([y[spec][j] * ir[spec][w, 1] for w in range(len(frequencies))])
#		print(i_intensity["H2"][-1][242], i_intensity["N2"][-1][2352], i_intensity["NH3"][-1][3262])    # species ; variable (T or t) ; frequency
		for j in range(len(new_x)):
			for w in range(len(frequencies)):
				ir_z_spectra.append(sum([i_intensity[spec][j][w] for spec in y]))
		ir_z_spectra = np.reshape(ir_z_spectra, (len(frequencies), len(new_x)), order="F")

		IR_3D(experiment, plot_temp, plot_time, frequencies, new_x, ir_z_spectra, [i for i in y])
		Extract_numeric_data(experiment, labels_species, conditions, species, time_range,
						 plot_v, plot_ph, plot_temp, plot_time, plot_species, plot_ini_species)



def	Species_3D(experiment, labels_species, conditions, species, time_range,
					   plot_v, plot_ph, plot_temp, plot_time, plot_species, plot_ini_species):
# plot_species is plot_systems
	x = []
	y = []
	z = []
	name = str("V=%.2f_pH=%.2f_%s" % (plot_v, plot_ph, plot_ini_species))
	for i in range(len(conditions)):
		v, ph, temp, time = conditions[i]
		if v == plot_v and ph == plot_ph:
			if time == min(time_range) and experiment == "const_TEMP" or\
					time == min(time_range) and temp == min(temp_range) and experiment == "variable_TEMP":
				label_comment = ''
				for j in range(len(labels_species)-1):
					if species[i, j] != 0:
						label_comment += labels_species[j] + "=" + str(round(species[i, j], 2))
			if label_comment == plot_ini_species:
				if temp in plot_temp and time in plot_time:
						spec = int(plot_species[0])
						x.append(temp)
						y.append(time)
						z.append(species[i, spec])

	figure = plt.figure(figsize=(12, 10), clear=True)
	ax = figure.add_subplot(111, projection='3d')

	z = list(map(float, z))
	grid_x, grid_y = np.mgrid[min(x):max(x):200j, min(y):max(y):200j]
	grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')	# better "cubic" but it generates Nan values
	surface = ax.plot_surface(grid_x, grid_y, grid_z, cmap="plasma", alpha=0.7, edgecolor='none', linewidth=0, antialiased=False)
	surface.set_clim(min(z), max(z))

#	surface = ax.plot_trisurf(x, y, z, cmap="viridis", edgecolor='none', linewidth=0, antialiased=False)
	figure.colorbar(surface, shrink=0.5, aspect=8)
#	ax.scatter3D(x, y, z, cmap="viridis")
	ax.set_xlabel(r'Temperature (K)', rotation=0, fontsize=16)
	ax.set_ylabel(r'time (s)', rotation=0, fontsize=16)
	ax.set_zlabel("$\\theta_{i}$ (ML)", fontsize=16) #  % labels_species[int(plot_species[0])], fontsize=16)
#	ax.set_zlim([0, 0.11])
#	ax.zaxis.set_ticklabels([])
#	ax.zaxis.labelpad = 5
	ax.zaxis._axinfo['label']['space_factor'] = 1.0
#	ax.zaxis.set_major_locator(LinearLocator(5))
	ax.zaxis.set_major_formatter(FormatStrFormatter("%.16g"))
	ax.view_init(azim=-155, elev=15)
#	ax.set_title(name)
#	plt.savefig("./KINETICS/PLOTS/"+experiment+"/"+name+".svg",
#				bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)
	plt.ion()
	plt.show()
	SaveFig(experiment)
	Extract_numeric_data(experiment, labels_species, conditions, species, time_range,
						 plot_v, plot_ph, plot_temp, plot_time, plot_species, plot_ini_species)

def SaveFig(experiment):
	answer = str(input("Would you like to save the figure (y/n)?\n"))
	if answer == "y":
		figure_out_name = "KINETICS/PLOTS/" + experiment + "/" + str(input("What would it be the figure name (a word & no format)?\n"))
		plt.savefig(figure_out_name + ".svg",
					bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)


def IR_3D(experiment, plot_temp, plot_time, x, y, z, spec_names):
	z_max = 0
	for i in range(len(y)):
		if max(z[:, i]) > z_max:
			z_max = max(z[:, i])
	name = "IR"
	for i in spec_names:
		name += "_" + i
	figure = plt.figure(figsize=(11.69, 16.53), clear=True)
	ax = figure.add_subplot(111, projection='3d')
	for i in range(len(y)):
		ax.plot(x, np.full((len(x), 1), y[i]), z[:, i]/z_max)
	ax.set_xlabel("$\\omega$ $(cm^{-1})$", rotation=0, fontsize=16)
	if type(plot_temp) is list:
		ax.set_ylabel(r'Temperature (K)', rotation=0, fontsize=16)
	elif type(plot_time) is list:
		ax.set_ylabel(r'time (s)', rotation=0, fontsize=16)
	ax.zaxis.set_ticklabels([])
	ax.set_zlabel("Intensity (a.u.)", fontsize=16)
#	ax.zaxis.set_major_locator(LinearLocator(5))
#	ax.zaxis.set_major_formatter(FormatStrFormatter("%.2g"))
	ax.view_init(azim=-45, elev=20)
	ax.set_title(name)
#	plt.savefig("./KINETICS/PLOTS/"+experiment+"/"+name+".svg",
#				bbox_inches='tight', dpi=300, orientation='landscape', transparent=True)
	SaveFig(experiment)
	plt.show()

def Species(molecules, surface_species):
	labels_spe = []
	for i in molecules:
		labels_spe.append(i)
	for i in surface_species:
		labels_spe.append(i)
	print("\nThe species in the system are:")
	n = 1
	for i in range(1, len(labels_spe)):
		if i/4 < n:
			print(" [%d] %-15s\t" % (i-1, labels_spe[i-1]), end='', flush=True)
		else:
			n += 1
			print(" [%d] %-15s" % (i-1, labels_spe[i-1]))
	print(" [%d] %-15s" % (len(labels_spe)-1, labels_spe[-1]))
	print(colored("\n [%d] %s" % (len(labels_spe), "ALL Molecules"), "green"))
	print(colored(" [%d] %s" % (len(labels_spe)+1, "ALL Surface Species"), "red"))
	print(colored(" [%d] %s" % (len(labels_spe)+2, "ALL"), "blue"))
	answer = input("Which species do you want to include in the plot? Numbers space separated. Note that for multiple T"
				   " and t, only ONE species (value) is accepted.\n")
	answer = [int(i) for i in answer.split(" ")]
	if answer[0] == len(labels_spe):
		plot_spe = molecules
	elif answer[0] == len(labels_spe)+1:
		plot_spe = [i for i in labels_spe if i not in molecules]
	elif answer[0] == len(labels_spe)+2:
		plot_spe = [i for i in labels_spe]
	else:
		for i in answer:
			if i > len(labels_spe)+2 or i < 0:
				print("Select the right species")
				exit()
		plot_spe = [labels_spe[i] for i in answer]
	return plot_spe


def Species_spe(molecules, surface_species, spe):
	labels_spe = []
	for i in molecules:
		labels_spe.append(i)
	for i in surface_species:
		labels_spe.append(i)
	if spe == len(labels_spe):
		plot_spe = molecules
	elif spe == len(labels_spe)+1:
		plot_spe = [i for i in labels_spe if i not in molecules]
	elif spe == len(labels_spe):
		plot_spe = [i for i in labels_spe]
	else:
		plot_spe = [labels_spe[spe]]
		while plot_spe[0] not in labels_spe:
			print("Select the right species")
			answer = int(input("Which species do you want to include in the plot?\n"))
			plot_spe = [labels_spe[answer]]
	return plot_spe


def IRs(freq_path, plot_species):
	ir = {}
	if not os.path.exists("./IRs"):
		os.makedirs("./IRs")
	if not os.path.exists("./IRs/originals/"):
		os.makedirs("./IRs/originals")
	for i in plot_species:
		if not os.path.exists("./IRs/originals/" + i + ".dat"):
			p = subprocess.Popen(["/home/alberto/Software/KINETICS/vasp2ir_DFT.pl", freq_path[i]])
			p.wait()
			if os.path.exists("./IRSPECTRA"):
				ir[i] = np.loadtxt("./IRSPECTRA", comments="#")
				os.rename("./IRSPECTRA", "./IRs/originals/" + i + ".dat")
				os.remove("./IRCAR")
		else:
			ir[i] = np.loadtxt("./IRs/originals/" + i + ".dat")
	return ir


def	Extract_numeric_data(experiment, labels_species, conditions, species, time_range,
					   plot_v, plot_ph, plot_temp, plot_time, plot_species, plot_ini_species):
	answer = str(input("Would you like to extract numeric data from the previous plot (y/n)?\n"))
	if answer == "y":
		data_out_name = "KINETICS/DATA/" + experiment + "/" + str(input("What would it be the output file name (a word)?\n"))
		data_out = open(data_out_name + ".dat", "w+")
		data_out.write("# v=" + str(plot_v) + " ph=" + str(plot_ph) + " " + experiment + "\n")
		data_out.write("#\t temp\ttime")
		for spec in plot_species:
			data_out.write("\t{:s}" .format(labels_species[spec]))
		data_out.write("\n")

		for i in range(len(conditions)):
			v, ph, temp, time = conditions[i]
			if v == plot_v and ph == plot_ph:
				if time == min(time_range) and experiment == "const_TEMP" or \
					time == min(time_range) and temp == min(temp_range) and experiment == "variable_TEMP":
					label_comment = ''
					for j in range(len(labels_species)-1):
						if species[i, j] != 0:
							label_comment += labels_species[j] + "=" + str(round(species[i, j], 2))
				if label_comment == plot_ini_species:
					if type(plot_temp) is list and temp in plot_temp or type(plot_temp) is float and temp == plot_temp:
						if type(plot_time) is list and time in plot_time or type(plot_time) is float and time == plot_time:
							data_out.write("{:> 9.2f} {:> 3.6f}" .format(temp, time))
							for spec in plot_species:
								data_out.write(" {: 2.16g}" .format(species[i, spec]))
							data_out.write("\n")
		data_out.close()


#####################################################################################################################
######################################################################################################################

# reads the initial file mk.in
input_file = sys.argv[1]
try:
	f = open(input_file)
	lines = f.readlines()
except IOError:
	raise Exception("   "+input_file+" is not a valid!")
if len(sys.argv) > 2:
	thermo = str(sys.argv[2])
	print("################\n# thermo=", thermo)
else:
	thermo = None
if len(sys.argv) > 3:
	krate = str(sys.argv[3])
	print("# krate=", krate)
else:
	krate = None
if len(sys.argv) > 4:
	exp = str(sys.argv[4])
	print("# exp=", exp)
else:
	exp = None
if len(sys.argv) > 5:
	v0 = str(sys.argv[5])
	print("# v=", v0)
else:
	v0 = None
if len(sys.argv) > 6:
	ph0 = str(sys.argv[6])
	print("# ph=", ph0)
else:
	ph0 = None
if exp != "3":
	if len(sys.argv) > 7:
		tem = str(sys.argv[7])
		print("# tem=", tem)
	else:
		tem = None
	if tem == "ALL" or tem == "All" or tem == "all":
		if len(sys.argv) > 8:
			tim = str(sys.argv[8])
			print("# tim=", tim)
		else:
			tim = None
		if len(sys.argv) > 9:
			spe = str(sys.argv[9])
			print("# spe=", spe)
		else:
			spe = None
		if len(sys.argv) > 10:
			comp = str(sys.argv[10])
			print("# comp=", comp, "\n################")
		else:
			comp = None
	else:
		if tem is not None:
			tim = "ALL"
			print("# tim= all")
		else:
			tim = None
		if len(sys.argv) > 8:
			spe = str(sys.argv[8])
			print("# spe=", spe)
		else:
			spe = None
		if len(sys.argv) > 9:
			comp = str(sys.argv[9])
			print("# comp=", comp, "\n################")
		else:
			comp = None
else:
	tem = "ALL"
	print("# tem= all")
	tim = None
	if len(sys.argv) > 7:
		spe = str(sys.argv[7])
		print("# spe=", spe)
	else:
		spe = None
	if len(sys.argv) > 8:
		comp = str(sys.argv[8])
		print("# comp=", comp, "\n################")
	else:
		comp = None

# Loops the initial file looking for system's names
systems = []
molecules = []
surface_species = []
processes = []
freq_path = {}
for line in lines:
	words = line.split(" ")
	words = [i.strip() for i in words if i]
# Looking for parameters in the InputFile
	if words[0] == "SYSTEM":
		systems.append(words[2])
	if words[0] == "FREQPATH":
		freq_path[systems[-1]] = words[2]
	if words[0] == "FREQ":
		freq = ["no" for i in words[2:] if float(i) < 0]
	if words[0] == "IPRESSURE" or words[0] == "RPRESSURE":
		molecules.append(systems[-1])
	if words[0] == "ICOVERAGE" or words[0] == "RCOVERAGE":
		if "no" not in freq:
			surface_species.append(systems[-1])
	if words[0] == "PROCESS":
		processes.append(words[2])

# Loops the initial systems to plot THERMODYNAMIC variables
if thermo is None:
	answer = input("Do you want to plot Cp, S, H and G?\n")
else:
	answer = thermo
if answer.startswith('y') is True or answer.startswith('Y') is True:
	for syst in systems:
		Thermodynamics(syst)

# Plots all the PROCESS log(parameters) as a function of T
if krate is None:
	answer = input("Do you want to plot the log of sticky coefficients,"
				   " Arrhenius pre-exponentials and reaction constants?\n")
else:
	answer = krate
if answer.startswith('y') is True or answer.startswith('Y') is True:
	Reaction_Constants(processes)

if spe is None:
	plot_species = Species(molecules, surface_species)
else:
	plot_species = Species_spe(molecules, surface_species, int(spe))

# Prepares to plot the experiments
if exp is None:
	answer = input("Chose the experiment from which do you want to plot results:\n"
				   " [1] Constant Temperature\n"
				   " [2] Variable Temperature\n"
				   " [3] Temperature Programmed Reaction (TPR / TPD)\n")
else:
	answer = exp
if answer == "1":
	experiment = "const_TEMP"
	experiment_path = "./KINETICS/DATA/" + experiment + "/solution" + experiment + ".dat"
elif answer == "2":
	experiment = "variable_TEMP"
	experiment_path = "./KINETICS/DATA/" + experiment + "/solution" + experiment + ".dat"
elif answer == "3":
	experiment = "TPR"
	experiment_path = "./KINETICS/DATA/" + experiment + "/TPR_" + plot_species[0] + ".dat"

try:
	datainput = experiment_path
	labels = open(datainput, "r").readline()				    # labels are the first raw of data
	labels = labels.split()                                 # makes the line a list
	labels.pop(0)
	labels_species = labels[4:]
	conditions = np.loadtxt(datainput, comments="#")[:, :4]
	species = np.loadtxt(datainput, comments="#")[:, 4:]
except IOError:
	raise Exception("   " + experiment_path + " results not found!")
	pass
# conditions to plot
v_range = [float(i) for i in sorted(set(conditions[:, 0]))]
if v0 is None:
	if len(v_range) > 1:
		print("\nThe applied external potentials are: ", v_range)
		plot_v = float(input("Which external applied voltage do you want to set in the plot? One value."))
		while plot_v not in v_range:
			print("Select a right potential")
			plot_v = float(input("Which external applied voltage do you want to set in the plot? One value."))
	else:
		plot_v = float(v_range[0])
else:
	plot_v = float(v0)

ph_range = [float(i) for i in sorted(set(conditions[:, 1]))]
if ph0 is None:
	if len(ph_range) > 1:
		print("\nThe simulated pH are: ", ph_range)
		plot_ph = float(input("Which pH do you want to set in the plot? One value."))
		while plot_ph not in ph_range:
			print("Select a right pH")
			plot_ph = float(input("Which pH do you want to set in the plot? One value."))
	else:
		plot_ph = float(ph_range[0])
else:
	plot_ph = float(ph0)

temp_range = [float(i) for i in sorted(set(conditions[:, 2]))]
if tem is None:
	if len(temp_range) > 1 and experiment != "TPR":
		print("\nThe simulated temperatures are:")
		temp_range.append("ALL")
		n = 1
		for i, t in enumerate(temp_range):
			if i/6 < n:
				if type(t) is float:
					print(" %.2f\t" % (t), end='', flush=True)
				else:
					print(colored(" %s" % (t), "green"))
			else:
				n += 1
				if type(t) is float:
					print(" %.2f" % (t))
				else:
					print(colored(" %s" % (t), "green"))
		answer = input("Which temperature(s) do you want to set in the plot? Note that a single T will produce a "
						  "2D species vs. t plot and multiple T will ask you for an specific time\n")
		if answer == "ALL" or answer == "all" or answer == "All":
			plot_temp = temp_range[:-1]
		else:
			plot_temp = []
			answer = [float(i) for i in answer.split()]
			if len(answer) > 1:
				for t in answer:
					if t not in temp_range:
						print("Select the right temperature")
						exit()
					else:
						plot_temp.append(t)
			else:
				if answer[0] not in temp_range:
					print("Select the right temperature")
					exit()
				else:
					plot_temp = float(answer[0])
	elif len(temp_range) == 1 and experiment != "TPR":
		plot_temp = float(temp_range[0])
	elif experiment == "TPR":
		plot_temp = temp_range[:-1]
else:
	if tem == "ALL" or tem == "all" or tem == "All":
		plot_temp = temp_range[:]
	else:
		plot_temp = float(tem)

time_range = [float(i) for i in sorted(set(conditions[:, 3]))]
if type(plot_temp) is list and experiment != "TPR":
	if tim is None:
		print("\nThe simulated times are:")
		time_range.append("ALL")
		n = 1
		for i, t in enumerate(time_range):
			if i/6 < n:
				if type(t) is float:
					print(" %.6f\t" % (t), end='', flush=True)
				else:
					print(colored(" %s" % (t), "green"))
			else:
				n += 1
				if type(t) is float:
					print(" %.6f" % (t))
				else:
					print(colored(" %s" % (t), "green"))
		answer = input("Which time(s) do you want to set in the plot? Note that a single t value will produce a "
						  "2D species vs. T plot and multiple t will generate a 3D plot for a single species.\n")
		if answer == "ALL" or answer == "all" or answer == "All":
			plot_time = time_range[:-1]
		else:
			plot_time = []
			answer = [float(i) for i in answer.split()]
			if len(answer) > 1:
				for t in answer:
					if t not in time_range:
						print("Select the right time")
						exit()
					else:
						plot_time.append(t)
			else:
				if answer[0] not in time_range:
					print("Select the right time")
					exit()
				else:
					plot_time = float(answer[0])
	else:
		if tim == "ALL" or tim == "all" or tim == "All":
			plot_time = time_range[:-1]
		else:
			plot_time = float(tim)
elif type(plot_temp) is float and experiment != "TPR":
	plot_time = time_range[:-1]
elif experiment == "TPR":
	plot_time = float(time_range[0])

# multiple initial species' compositions
initial_compositions = []
for i in range(len(conditions)):
	v, ph, temp, time = conditions[i]
	if experiment != "TPR":
		if v == min(v_range) and ph == min(ph_range) and temp == min(temp_range[:-1]) and time == min(time_range[:-1]):
			label_comment = ''
			for j in range(len(labels_species[:-1])):
				if species[i, j] != 0:
					label_comment += labels_species[j] + "=" + str(round(species[i, j], 2))
			if label_comment != '':
				initial_compositions.append(label_comment)
	else:
		if v == min(v_range) and ph == min(ph_range) and temp == min(temp_range[:-1]):
			label_comment = ''
			if type(plot_species) == int:
				label_comment = labels_species[plot_species] + "=" + str(round(species[i, plot_species], 2))
			elif type(plot_species) == list:
				for j in range(len(plot_species)):
					if species[i, j] != 0:
						label_comment += labels_species[j] + "=" + str(round(species[i, j], 2))
			if label_comment != '':
				initial_compositions.append(label_comment)
if experiment != "TPR":
	if comp is None:
		if len(initial_compositions) > 1:
			print("\nThere are different initial compositions of the system:")
			for i, c in enumerate(initial_compositions):
				print(" [%d] %-15s" % (i, c))
			answer = input("Which initial composition do you want to include in the plot? One number.\n")
			while int(answer) > len(initial_compositions):
				print("Select a single initial composition")
				answer = input("Which conditions do you want to include in the plot? One number.\n")
			plot_ini_composition = initial_compositions[int(answer)]
		else:
			plot_ini_composition = initial_compositions[0]
	else:
		plot_ini_composition = initial_compositions[int(comp)]
else:
	plot_ini_composition = initial_compositions

plot_systems = [i for i in range(len(labels_species)) if labels_species[i] in plot_species]

if experiment == "const_TEMP" or experiment == "variable_TEMP":
	if type(plot_temp) is float and type(plot_time) is list:
		ir = IRs(freq_path, plot_species)
		Species_2D(experiment, labels_species, conditions, species, time_range,
					plot_v, plot_ph, plot_temp, plot_time, plot_systems, plot_ini_composition, ir)

	elif type(plot_temp) is list and type(plot_time) is float:
		ir = IRs(freq_path, plot_species)
		Species_2D(experiment, labels_species, conditions, species, time_range,
					plot_v, plot_ph, plot_temp, plot_time, plot_systems, plot_ini_composition, ir)
	elif type(plot_temp) is list and type(plot_time) is list:
		if len(plot_systems) > 1:
			print("Select the right conditions T, t, and number of species")
			exit()
		Species_3D(experiment, labels_species, conditions, species, time_range,
					plot_v, plot_ph, plot_temp, plot_time, plot_systems, plot_ini_composition)
elif experiment == "TPR":
	for spe in plot_systems:
		ir = IRs(freq_path, [labels_species[spe]])
		Species_2D(experiment, labels_species, conditions, species, time_range,
					plot_v, plot_ph, plot_temp, plot_time, spe, plot_ini_composition, ir)



