"""
	This script builds on the perl version by A.Roldan.

"""

import os, pathlib
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp


def printdata(rconditions, species, experiment, arguments, values):
	data = [[*rconditions.keys(), *species]]
	for row in values:
		data.append([*arguments, *row])
	maxlen = [max([len(f"{data[r][c]}") + 2 for r in range(len(data))]) for c in
			  range(len(data[0]))]  # max length per column

	folder = './KINETICS/DATA'
	outputfile = folder + "/" + str(experiment) + ".dat"
	if not pathlib.Path(folder).exists():
		pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
	output = open(outputfile, "a")
	output.write("#")
	for i in range(len(data[0])):
		output.write(" {val:>{wid}s}".format(wid=maxlen[i], val=data[0][i]))  # headings
	output.write("\n")
	for row in data[1:]:
		for i in range(len(row)):
			output.write(" {val:>{wid}.3{c}}".format(wid=maxlen[i], val=row[i],
													 c='e' if row[i] > 1e3 or np.abs(row[i]) < 1e-2 else 'f'))
	output.write("\n")
	output.close()

def ode_solver(time, species, rhs, ics, arguments):
	t, temp = sp.symbols('time temperature', positive=True, real=True)
	t_span = (time[:1])  # time grid
	t_eval = np.arange(*time)
	conditions = [t, temp]  # conditions required to lambdify
	f_ode = sp.lambdify((*conditions, *species), rhs, "numpy")  # Convert symbolic into numerical
	def ode_system(t, species, temp):    # define ode function compatible with solve_ivp
		return f_ode(t, temp, *species)

	print(*arguments)

	sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval, args=(*arguments,), method='BDF', rtol=1e-8, atol=1e-10)
	return list(zip(sol.t, *sol.y))


class ConsTemperature:
	def __init__(self, rconditions, systems, equations):
		t , temp = sp.symbols('time temperature', positive=True, real=True)
		ics, species = self.initial_species(systems)
		rhs = []
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				rhs.append(equations[name])     # rhs for SciPy
		for temp in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
			data = ode_solver(rconditions["time"], species, rhs, ics, [temp])
			printdata(dict(rconditions), species, str("Cons_Temperature"), [temp], data)


        '''    A, B, C = sol.y
        rate_overall = k2 * B[-1]  # e.g. formation rate of C at steady state
        return rate_overall
        degree_of_rate_control(k_list, solver_func, delta=0.05)
        '''



	@staticmethod
	def initial_species(systems):  # process is processes[process]
		ics = []    # list of initial concentrations in the order of systems[names]
		species = []    # list of species in the order on systems
		surf0 = {}
		sites_name = []
		for name in systems.keys():
			if systems[name]['kind'] == 'surface':
				sites_name.extend(systems[name]['sites'])
		for s in list(set(sites_name)):
			surf0[s] = 1
		for name in systems.keys():
			if systems[name]['kind'] == "molecule":
				ics.append(systems[name]["pressure0"])
				species.append(name)
			elif systems[name]['kind'] == "adsorbate":
				ics.append(systems[name]["coverage0"])
				''' adsorbates have only one kind of adsorption site per system '''
				surf0[systems[name]['sites'][0]] -= systems[name]["coverage0"] * systems[name]["nsites"]
				species.append(name)
		for s in surf0.keys():
			ics.append(surf0[s])
			species.append(s)
		return ics, species
	'''
	@staticmethod
	def rate_control(species, equations): 
	def degree_of_rate_control(k_list, solver_func, delta):
		\''' The Degree of Rate Control (DRC), introduced by C. T. Campbell (J. Catal. 204, 520, 2001), quantifies how
		   sensitive the overall reaction rate is to each elementary step’s rate constant. \'''
		k_list = np.array(k_list, dtype=float)
		#Provided?  r0 = solver_func(k_list)  # steady-state overall rate with nominal k
	
		X = np.zeros_like(k_list)
	
		for i in range(len(k_list)):
			k_up = k_list.copy()
			k_down = k_list.copy()

			k_up[i] *= (1 + delta)
			k_down[i] *= (1 - delta)

			r_up = solver_func(k_up)
			r_down = solver_func(k_down)

			X[i] = (np.log(r_up) - np.log(r_down)) / (np.log(k_up[i]) - np.log(k_down[i]))

		# Normalize to sum = 1 (optional)
		if np.sum(X) != 0:
			X /= np.sum(X)

		#   return X        ploted instead <<< make it general also for SELECTIVITY CONTROL
		import matplotlib.pyplot as plt
        labels = [f"k{i+1}" for i in range(len(X))]
        plt.bar(labels, X)
        plt.ylabel("Degree of Rate Control (X_i)")
        plt.title("Rate Control Analysis (Campbell)")
        plt.show()
		    x_limits = [min(xmin)-0.5, max(xmax)+1]
        ax1.plot(x_limits, [0, 0], linestyle=":", lw=0.5, color=icolour[0])
        ax1.tick_params(labelsize=14)
        ax1.set_xticks([])
        ax1.set_xticks(xtick_location) 
        ax1.set_xticklabels(xtick_label, rotation=25, ha="right")    # rotation=0, ha="center")
        ax1.set_xlim(x_limits)
        ax1.set_ylabel(y_label, fontsize=18)
        ax1.set_ylim(y_limit)
        #nyticks = 5 
        #ax1.set_yticks(np.round(np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[1], nyticks), 1))
        legend = ax1.legend(legend_lines, legend_labels, loc='best', fontsize=14) #upper left   #best
        fig.tight_layout()
        plt.ion()
        plt.show()
        SaveFig()
        def SaveFig():
        answer = str(input("Would you like to save the figure (y/n)?\n"))
        if answer == "y":
            figure_out_name = str(input("What would it be the figure's name (a word & no format)?\n"))
            plt.savefig(figure_out_name + ".svg", dpi=300, orientation='landscape', transparent=True)




	# Selectivity control
	'''

class TPR:
	def __init__(self, rconditions, systems, processes, equations):
		t , temp = sp.symbols('time temperature', positive=True, real=True)
		species = self.initial_species(systems)     # list of species in the order on systems
		rhs = []
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				rhs.append(equations[name])     # rhs for SciPy
		arguments = [temp]
		''' In Temperature Desorption experiments, the initial conditions consider only adsorbed molecules.
		Therefore for Reaction type A, the coverage of the product will be 1.
		There will be a TPD/TPR for each of the molecules in process[kind]=A'''
		tpd_done = []
		for process in processes.keys():
			if processes[process]['kind'] == 'A':
				ics = []  # list of initial concentrations in the order of species
				out_file = ''
				for name in species:
					for i, product in enumerate(processes[process]['products']):
						if name == product:
							ics.append(np.floor([1/processes[process]['pstoichio'][i]])[0])
							out_file = 'TPD_' + str(name)
						else:
							ics.append(0.)
				if out_file not in tpd_done:
					# T rate =   0.01, 0.1,  1, 5, 10,    15, 20 K/s
					for t in [1000, 100, 10, 2, 1, 0.6666, 0.5]:
						time = [0, 10/t, 10/t]       # temperature rate in k/s
						data = []
						for temp in np.arange(50, 1000, 10):     # integrate at different temperatures
							data.append(ode_solver(time, species, rhs, ics, [temp]))
						printdata(dict(rconditions), species, out_file + str(time[-1]), [temp], data)
					tpd_done.append(out_file)

	@staticmethod
	def initial_species(systems):  # process is processes[process]
		species = []    # list of species in the order on systems
		surf0 = {}
		sites_name = []
		for name in systems.keys():
			if systems[name]['kind'] == 'surface':
				sites_name.extend(systems[name]['sites'])
		for s in list(set(sites_name)):
			surf0[s] = 1
		for name in systems.keys():
			if systems[name]['kind'] == "molecule":
				species.append(name)
			elif systems[name]['kind'] == "adsorbate":
				species.append(name)
		for s in surf0.keys():
			species.append(s)
		return species





