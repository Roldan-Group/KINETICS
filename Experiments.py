"""
	This script builds on the perl version by A.Roldan.

"""

import os, pathlib
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from Symbols_def import t, temp, constants



def printdata(experiment, data):
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
	''' time: tuple or list, e.g. (0, 10, 0.1)
		t_num is the numeric valu of time of integration (solve_ivp)
		species: list of sympy symbols [theta_A, theta_B, ...]
		rhs: list of sympy expressions for d(species)/dt
		ics: list of initial conditions
		arguments: tuple or list of non-time parameters (e.g. temp, pressures, etc.)
		constants: dict of physical constants {h: ..., kB: ..., ...} '''
	''' substitute the symbolic constants, e.g. h, kb, ..., by its values '''
	rhs = [i.subs(constants) for i in rhs]
	conditions = [t, temp]  # conditions required to lambdify
	f_ode = sp.lambdify((*conditions, *species), rhs, ["numpy", 'sympy'])  # Convert symbolic into numerical

	def ode_system(t_num, y, *args):    # define ode function compatible with solve_ivp
		temp_num = args[0] if len(args) > 0 else None
		dydt = f_ode(t_num, temp_num, *y)
		return np.array(dydt, dtype=float).flatten()

	''' substitute rconditions'''
	t_span = (time[:2])  # time grid
	t_eval = np.arange(*time)
	sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval, args=(*arguments,), method='BDF', rtol=1e-8, atol=1e-10)
	data = []
	for t in sol.t:
		data.append([temp, t, *sol.y.T])

	print(data)

	return sol, data


class ConsTemperature:
	def __init__(self, rconditions, systems, processes, equations):
		''' in arguments, the temperature numeric value should be the first entry (arg[0] '''
		ics, species = self.initial_species(systems)
		rhs = []
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				rhs.append(equations[name])     # rhs for SciPy
		''' evaluate the Reaction Rates '''
		data = [[*rconditions.keys(), *species]]  # basic: Temp, time, species
		rates_ss = [[*rconditions.keys(), *species]]
		rates_avg = [[*rconditions.keys(), *species]]

		if isinstance(rconditions['temperature'], (int, float)): # single temperature
			sol, sol_T = ode_solver(rconditions["time"], species, rhs, ics, [rconditions['temperature']])
			data.append(sol_T)
			ss, avg = self.rates(processes, species, temp, sol)
			rates_ss.append(ss)
			rates_avg.append(avg)
		else:      # temperature ramp
			for temp in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
				sol, sol_T = ode_solver(rconditions["time"], species, rhs, ics, [temp])
				data.append(sol_T)  # transposed solution: (time x species)
				ss, avg = self.rates(processes, species, temp, sol)
				rates_ss.append(ss)
				rates_avg.append(avg)
		start = time.time()
		print("\t ... Generating Concentrations and Rates ...")
		printdata(str("Cons_Temperature"), data)
		printdata("SteadyState_Rates", rates_ss)    # temperature x processes
		self.barplot(processes, "Steady-state", rates_ss, 0.5)
		printdata("Average_Rates", rates_avg)       # temperature x processes
		self.barplot(processes, "Average", rates_avg, 0.5)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")
		start = time.time()
		print("\t ... Generating Degree of Rate and Selectivity Control ...")
		self.degree_of_rate_control(rconditions, processes, systems, species, ics, rates_ss)

		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")



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

	@staticmethod
	def rates(processes, species, temp, sol):
		''' - equations: list of sympy expressions r1(..), r2(..), ... using species symbols
		 - species: list of sympy symbols [theta_A, theta_B, ...]
		substitute the symbolic constants, e.g. h, kb, ..., by its values '''
		''' evaluating the rates as a function of time '''
		rates_ss = [temp]   # rates at the steady-state, i.e. over the last 10% of the time points
		rates_avg = [temp]    # rate averages along t_span using the trapezoidal rule to integrate the rate curve.
		for process in processes:
			ss, avg = rki_value(species, processes[process]['krate0'], processes[process]['reactants'], temp, sol)
			rates_ss += ss
			rates_avg += avg
		return rates_ss, rates_avg

	@staticmethod
	def rki_value(species, krate, reactants, temp, sol):
		''' - equations: list of sympy expressions r1(..), r2(..), ... using species symbols
		 - species: list of sympy symbols [theta_A, theta_B, ...]
		substitute the symbolic constants, e.g. h, kb, ..., by its values '''
		''' evaluating the rates as a function of time '''
		args = [sol.y[i] for i in range(len(species)) if species[i] in reactants]
		args.append(temp)
		rate_time = sp.lambdify((*args), krate.subs(constants), ['numpy', 'sympy'])
		n_tail = int(0.1 * len(rate_time))  # 10% of the last points
		rate_ss = rate_time[-n_tail:, :].mean(axis=0)   # rates at the steady-state, i.e. over the last 10% of the time points
		rate_avg = np.trapezoid(rate_time, sol.t) / (sol.t[-1] - sol.t[0])   # rate averages along t_span using the trapezoidal rule to integrate the rate curve.
		return rate_ss, rate_avg

	@staticmethod
	def barplot(processes, experiment, data, bar_width):
		icolour = ["b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki",
		   "teal"]
		ipatterns = ["///", "...", "xx", "**", "\\", "|", "--", "++", "oo", "OO"]

		y_label = experiment + "TOF ($molecule \cdot site^{-1} \cdot s^{-1}$)"
		labels = [process for process in processes]   # processes
		data_label = [data[0][0], data[-1][0]]  # initial and final temperature
		x = [len(data_label) * i for i in range(len(labels))]   # 2 temp per specie
		y = []
		for i in data[0]:   # for every column
			y.append(data[0][i])
			y.append(data[-1][i])
		# GroupedBars(labels, "DACs", np.array(y), "$E_{Ads}$ (eV)", data_label)

		fig, ax1 = plt.subplots(figsize=(8, 6), clear=True)
		y = np.reshape(y, (len(x), len(data_label)))
		for n in range(len(x)):
			for i in range(len(data_label)):
				offset = bar_width * i
				ax1.bar(x[n] + offset, y[n][i], width=bar_width, hatch=ipatterns[i],
					color=icolour[i], label=data_label[i])#, alpha=0.7)
		x_limit = [ax1.get_xlim()[0], ax1.get_xlim()[1]]
		ax1.plot(x_limit, [0, 0], "k:", lw=1.5)
		# ax1.set_xlabel(x_label, fontsize=18)
		ax1.set_xticks([i + ((len(data_label) -1) * bar_width)/2 for i in x])
		ax1.set_xticklabels([labels[i] for i in range(1, len(labels), len(data_label))], rotation=0, ha="right")
		ax1.set_xlim(x_limit)
		ax1.tick_params(axis='x', rotation=45, labelsize=16)
		ax1.set_ylabel(y_label, fontsize=20)
		#ax1.set_ylim([-0.18, 0.4])
		ax1.tick_params(axis='y', rotation=0, labelsize=16)
		# ax1.set_ylim(0, max(list(np.array(y[-2]) + y_top) + y[-1])*1.05)
		# legend = ax1.legend(loc="best", fontsize=14)
		leg_lines, leg_labels = ax1.get_legend_handles_labels()
		legend = ax1.legend(leg_lines[0:len(data_label)], leg_labels[0:len(data_label)], loc='best', fontsize=14)
		fig.tight_layout()
		plt.ion()
		plt.savefig('./KINETICS/DATA/' + experiment + "_Rates.svg", dpi=300, orientation='landscape', transparent=True)

	@staticmethod
	def degree_of_rate_control(rconditions, processes, systems, species, ics, rates_ss):
		''' The Degree of Rate Control (DRC), introduced by C. T. Campbell (J. Catal. 204, 520, 2001), quantifies how
		   sensitive the overall reaction rate is to each elementary step’s rate constant. '''
		time = rconditions['time']
		k_list = [processes[process]['krate0'] for process in processes]
		r0 = rates_ss[:][1:]  # steady-state rate for each k (temperature x krates)
		eps = 1e-3  # constant perturbation factor, small enough to retain linearity

		drc = [] # the degree of rate control will have krates for initial temperature
		''' because the forward and backward reactions are distinguishable BUT both have to be perturbed to 
		maintain the total Krate, the loop goes in 2 by 2. The result (drc) will have 1/2 len(Krate) as it is for 
		processes not for forward and backward constants'''
		for i in range(0, len(processes), 2):
			processes[str(i)]['krate0'] *= (1 + eps)    # the forward constant (1) + the perturbation
			processes[str(i+1)]['krate0'] *= (1 - eps)  # the backward constant (1) - the perturbation
			''' because rhs derives from the ConsTemperature REquations (in Kintics.py), processes[process][krate0] 
			has to be modified and REquations recalled to generate the perturbed rhs '''
			equations = REquations(dict(processes), dict(systems)).constemperature
			''' get the perturbed rhs '''
			rhs = []
			for name in species:  # lists of names with the order of ics
				if name in equations.keys():
					rhs.append(equations[name])  # rhs for SciPy
			for temp in [rates_ss[1][0], rates_ss[-1][0]]:  # initial (labels!) and final temperautes
				sol, _ = ode_solver(time, species, rhs, ics, [temp])
				rates, _ = rates(processes, species, temp, sol)   # r_f contains temp x krates, no labels
				r_f = rates[i+1]  # because the first column is temperature
				r_b = rates[i+1+1]    # temperature + backward constant
				k_f, _ = rki_value(species, processes[str(i)]['krate0'], processes[str(i)]['reactants'], temp, sol)
				k_b, _ = rki_value(species, processes[str(i+1)]['krate0'], processes[str(i+1)]['reactants'], temp, sol)
				drc.append( (np.log(r_f) - np.log(r_b)) / (np.log(k_f) - np.log(k_b)) )
		''' reset the reaction constants to the original form'''
		for i, process in enumerate(processes):
			processes[process]['krate0'] = k_list[i]
		# Normalize to sum = 1 (optional)
		if np.sum(drc) != 0:
			drc /= np.sum(drc)
		barplot(processes, "Degree_of_Rate_Control", drc, 0.5)
		data = [[*rconditions.keys(), *processes]]  # basic: Temp, time, processes
		for i in np.array(drc).reshape(-1, 2):
			data.append(i)
		printdata("Degree_of_Rate_Control", data)



	# Selectivity control


class TPR:
	def __init__(self, rconditions, systems, processes, equations):
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





