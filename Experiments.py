"""
	This script builds on the perl version by A.Roldan.

"""

import pathlib
import time
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from Symbols_def import t, temp, constants, sym_equation
from Kinetics import REquations
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt




def printdata(experiment, data):
	maxlen = 12 #[max([len(f"{data[r][c]}") + 2 for r in range(len(data))]) for c in range(len(data[0]))]  # max
	# length per column
	folder = './KINETICS/DATA'
	outputfile = folder + "/" + str(experiment) + ".dat"
	if not pathlib.Path(folder).exists():
		pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
	output = open(outputfile, "a")
	output.write("#")
	for i in range(len(data[0])):
		output.write(" {val:>{wid}s}".format(wid=maxlen, val=data[0][i]))  # headings
	output.write("\n")
	for row in data[1:]:
		for i in range(len(row)):
			output.write(" {val:>{wid}.3{c}}".format(wid=maxlen+1, val= float(row[i]),
												c='e' if float(row[i]) > 1e3 or np.abs(float(row[i])) < 1e-2 else 'f'))
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
	''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species'''
	f_ode = sp.lambdify((*conditions, *species), rhs, ["numpy", 'sympy'])  #
	''' Old versions of scipy do not accept args, so temp_num and other could be unwrapped here.'''
	def ode_system(t_num, y, *args):    # define ode function compatible with solve_ivp
		if len(args) < 1:
			raise ValueError("ERROR! Temperature not passed into ODE solver.")
		temp_num = args[0]
		dydt = f_ode(t_num, *args, *y)   # There is ONLY t and Temp, but more can be added
		return np.array(dydt, dtype=float).flatten()
	''' substitute rconditions'''
	t_span = (time[:2])  # time grid
	t_eval = np.arange(*time)
	sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval,
					args=arguments if isinstance(arguments, tuple) else (arguments,),
					method='BDF', rtol=1e-8, atol=1e-10)
	# Combine outputs (arguments + time + species)
	n_points = len(sol.t)
	data = np.column_stack([
		np.tile(arguments, (n_points, 1)),  # repeat all args per row
		sol.t,
		sol.y.T
	])
	return sol, data.flatten().tolist()


class ConsTemperature:
	def __init__(self, rconditions, systems, processes, equations):
		start = time.time()
		print("\t ... Generating Concentrations and Rates ...")
		''' in arguments, the temperature numeric value should be the first entry (arg[0] '''
		ics, species = self.initial_species(systems)
		rhs = []
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				rhs.append(equations[name])     # rhs for SciPy
		''' evaluate the Reaction Rates '''
		data = [*rconditions.keys(), *species]  # basic: Temp, time, species
		data_elements = len(data)
		rates_ss = [[key for key in rconditions.keys() if key != "time"] + [*processes.keys()]]  # time not required as it is at steady-state
		rates_avg = rates_ss.copy()   # time not required as it is an average along time
		rates_elements = len(rates_ss[0])
		if isinstance(rconditions['temperature'], (int, float)): # single temperature
			sol, sol_T = ode_solver(rconditions["time"], species, rhs, ics, (rconditions['temperature'],))
			data += sol_T
			ss, avg = ConsTemperature.rates(processes, species, rconditions['temperature'], sol)
			rates_ss.append(ss)
			rates_avg.append(avg)
		else:      # temperature ramp
			for temp_num in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
				sol, sol_T = ode_solver(rconditions["time"], species, rhs, ics, (temp_num,))
				data += sol_T  # transposed solution: (time x species)
				ss, avg = ConsTemperature.rates(processes, species, temp_num, sol)
				rates_ss.append(ss)
				rates_avg.append(avg)
		data = np.array(data).reshape(-1, data_elements)  # n columns, inferred rows
		printdata(str("Cons_Temperature"), data.tolist())

		labels = [process for process in processes]
		printdata("SteadyState_Rates", rates_ss)    # temperature x processes
		ConsTemperature.barplot(processes, "SteadyState Rates", rates_ss, labels, 0.5)
		printdata("Average_Rates", rates_avg)  # temperature x processes
		ConsTemperature.barplot(processes, "Average Rates", rates_avg, labels, 0.5)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")

		start = time.time()
		print("\t ... Generating Degree of Rate and Selectivity Control ...")
		ConsTemperature.degree_of_rate_control(rconditions, processes, systems, species, ics)
		ConsTemperature.degree_of_selectivity_control(rconditions, systems, processes)
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
	def rates(processes, species, temp_num, sol):
		''' - equations: list of sympy expressions r1(..), r2(..), ... using species symbols
		 - species: list of sympy symbols [theta_A, theta_B, ...]
		substitute the symbolic constants, e.g. h, kb, ..., by its values '''
		''' evaluating the rates as a function of time '''
		rates_ss = [temp_num]   # rates at the steady-state, i.e. over the last 10% of the time points
		rates_avg = [temp_num]    # rate averages along t_span using the trapezoidal rule to integrate the rate curve.
		for process in processes:
			ss, avg = ConsTemperature.rki_value(species, processes[process]['krate0'], temp_num, sol)
			rates_ss.append(ss)
			rates_avg.append(avg)
		return rates_ss, rates_avg

	@staticmethod
	def rki_value(species, krate, temp_num, sol):
		''' - equations: list of sympy expressions r1(..), r2(..), ... using species symbols
		 - species: list of sympy symbols [theta_A, theta_B, ...]
		substitute the symbolic constants, e.g. h, kb, ..., by its values '''
		''' evaluating the rates as a function of time '''
		rate_fn = sp.lambdify((temp, *species), krate.subs(constants), ['numpy'])
		rate_time = []
		for tidx in range(sol.y.shape[1]):    # species values at time tidx
			rate_time.append(rate_fn(temp_num, *sol.y[:, tidx]))
		n_tail = int(0.1 * len(rate_time))  # 10% of the last points
		rate_ss = np.array(rate_time[-n_tail:]).mean(axis=0)   # rates at the steady-state, i.e. over the last 10% of
		# the time points
		''' Numpy versions > 2.0 uses "trapezoid" instead of "trapz" '''
		rate_avg = np.trapz(rate_time, sol.t) / (sol.t[-1] - sol.t[0])   # rate averages along t_span using the trapezoidal rule to integrate the rate curve.
		return rate_ss, rate_avg

	@staticmethod
	def barplot(processes, experiment, data, labels, bar_width):
		icolour = ["b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki",
		   "teal"]
		ipatterns = ["///", "...", "xx", "**", "\\", "|", "--", "++", "oo", "OO"]

		y_label = "TOF ($molecule \cdot site^{-1} \cdot s^{-1}$)"
		x_label = experiment
		temps = [row[0] for row in data[1:]]  # temperatures; first row is for labels
		gap = 0.7   # gap between group of columns, e.g. processes
		x = np.arange(len(labels)) * (1 + gap)
		y = []
		for row in range(len(data)-1):   # for every process | column 0 is temperature
			y.append([float(data[row+1][column+1]) for column in range(len(labels))])
		# Convert to array and apply positive floor for log-scale
		y = np.array(y, dtype=float)
		y[y <= 0] = 1e-20    # to prevent log crash

		fig, ax1 = plt.subplots(figsize=(8, 6), clear=True)
		ax1.set_ylim(bottom=np.min(y) * 0.5, top=np.max(y) * 5)  # Prevents from hitting zero during autoscaling
		for n in range(len(y)):    # processes | the first column is temperature
			offset = (n - len(temps)/2) * bar_width + bar_width/2
			ax1.bar(x + offset, y[n], log=True,
					width=bar_width, hatch=ipatterns[n], color=icolour[n], label=f"{temps[n]} K", alpha=0.7)

		x_limit = [ax1.get_xlim()[0], ax1.get_xlim()[1]]
		ax1.plot(x_limit, [1e-20, 1e-20], "k-", lw=1.5)
		ax1.set_xlim(x_limit)
		ax1.set_xlabel(x_label, fontsize=16)
		ax1.set_xticks(x)
		ax1.tick_params(axis='x', rotation=0, labelsize=16)
		ax1.set_xticklabels(labels, rotation=0, ha="center")

		ax1.set_ylabel(y_label, fontsize=18)
		ax1.tick_params(axis='y', rotation=0, labelsize=16)
		leg_lines, leg_labels = ax1.get_legend_handles_labels()
		legend = ax1.legend(leg_lines[0:len(temps)], leg_labels[0:len(temps)], loc='best', fontsize=16)
		fig.tight_layout()
		plt.ion()
		plt.savefig('./KINETICS/DATA/' + "_".join(experiment.split()) + ".svg", dpi=300, orientation='landscape',
					transparent=True)

	@staticmethod
	def degree_of_rate_control(rconditions, processes, systems, species, ics):
		''' The Degree of Rate Control (DRC), introduced by C. T. Campbell (J. Catal. 204, 520, 2001), quantifies how
		   sensitive the overall reaction rate is to each elementary step’s rate constant. '''
		time = rconditions['time']
		k_list = [processes[process]['krate0'] for process in processes]
		eps = 1e-3  # constant perturbation factor, small enough to retain linearity
		''' because the forward and backward reactions are distinguishable BUT both have to be perturbed to 
		maintain the total Krate, the loop goes in 2 by 2. The result (drc) will have 1/2 len(Krate) as it is for 
		processes not for forward and backward constants'''
		labels = [f"{i}/{i+1}" for i in range(1, len(processes), 2)]
		data = [[*list(rconditions.keys())[:-1], *labels]]    # no time because it is at steady-state.
		for temp_num in np.arange(*rconditions['temperature']): # temperatures
			data_row = [temp_num]
			drc = []
			for i in range(1, len(processes), 2):
				processes[str(i)]['krate0'] *= (1 + eps)    # the forward constant (1) + the perturbation
				processes[str(i+1)]['krate0'] *= (1 - eps)  # the backward constant (1) - the perturbation
				''' because rhs derives from the ConsTemperature REquations (in Kinetics.py), processes[process][
				krate0] has to be modified and REquations recalled to generate the perturbed rhs '''
				equations = REquations(dict(processes), dict(systems)).constemperature
				''' get the perturbed rhs '''
				rhs = []
				for name in species:  # lists of names with the order of ics
					if name in equations.keys():
						rhs.append(equations[name])  # rhs for SciPy

				sol, _ = ode_solver(time, species, rhs, ics, (temp_num,))
				rates, _ = ConsTemperature.rates(processes, species, temp_num, sol)   # r_f contains temp x krates,
				# no labels
				r_f = rates[i]  #  starts from 1 because the first column is temperature
				r_b = rates[i+1]    # backward constant right after the forward process (i)
				k_f, _ = ConsTemperature.rki_value(species, processes[str(i)]['krate0'], temp_num, sol)
				k_b, _ = ConsTemperature.rki_value(species, processes[str(i+1)]['krate0'], temp_num, sol)
				drc.append( (np.log(r_f) - np.log(r_b)) / (np.log(k_f) - np.log(k_b)) )
				''' reset the reaction constants to the original form'''
				processes[str(i)]['krate0'] = k_list[i-1]   # k_list starts from 0 but processes from 1
				processes[str(i+1)]['krate0'] = k_list[i]   # k_list starts from 0 but processes from 1
			# Normalize drc to sum = 1 (optional)
			if np.sum(drc) != 0:
				drc /= np.sum(drc)
			data_row.extend(drc)
			data.append(data_row)
		printdata("Degree_of_Rate_Control", data)
		ConsTemperature.barplot(labels, "Degree of Rate Control", data, labels, 0.5)

	@staticmethod
	def degree_of_selectivity_control(rconditions, systems, processes):
		'''	Compute Campbell's Degree of Selectivity Control (DSC) for a product. '''
		# symbolic rate constants, one per reaction
		k_symbols = {process: sp.symbols(f'k_{process}') for process in processes}
		# Substitute actual expressions
		k_exprs = {process: processes[process]['krate0'].subs(constants) for process in processes}

		# products: molecules in systems without initial pressure
		products = []
		for name in systems:
			if systems[name]['kind'] == 'molecule':
				if systems[name]['pressure0'] == 0.0:
					products.append(name)

		# symbolic rates
		r_symbolic = {}
		r_sum = 0
		for name in products:
			r_symbolic[name] = sym_equation(processes, name)
			r_sum += r_symbolic[name]
		# Symbolic selectivity expression
		s_expr = {}
		for name in products:
			s_expr[name] = r_symbolic[name] / r_sum
		
		ds_dk = {name: {} for name in s_expr}  # nested dictionary of products and process
		dsc = {name: {} for name in s_expr}    # nested dictionary of products and process
		dsc_num = {name: {} for name in s_expr}    # nested dictionary with numeric dsc
		dsc_func = {name: {} for name in s_expr}    # nested lambdify function
		for name in products:
			# Symbolic derivatives dS/dk for each product and each k
			for process in processes:
				ds_dk[name][process] = sp.diff(s_expr[name], k_symbols[process])
				dsc[name][process] = ds_dk[name][process] * k_symbols[process] / s_expr[name]
				dsc_num[name][process] = dsc[name][process].subs(k_exprs)
				dsc_func[name][process] = sp.lambdify(temp, dsc_num[name][process], "numpy")

		labels = [f"{i}" for i in range(1, len(processes))]
		for name in products:
			data = [[*list(rconditions.keys())[:-1], *labels]]    # no time because it is at steady-state
			for temp_num in np.arange(*rconditions['temperature']): # temperatures
				row = [temp_num]
				row.extend([dsc_func[name][process](temp_num) for process in processes])
				data.append(row)

			print("name", name, "data", data)

			printdata(name + "_Degree_of_Selectivity_Control", data)
			ConsTemperature.barplot(labels, name +" Degree of Rate Control", data, labels, 0.5)


class TPR:
	def __init__(self, rconditions, systems, processes, equations):
		species = self.initial_species(systems)     # list of species in the order on systems
		rhs = []
		start = time.time()
		print("\t ... Generating TPR ...")
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				rhs.append(equations[name])     # rhs for SciPy
		''' In Temperature Desorption experiments, the initial conditions consider only adsorbed molecules.
		Therefore for Reaction type A, the coverage of the product will be 1.
		There will be a TPD/TPR for each of the molecules in process[kind]=A'''
		tpd_done = []
		for process in processes.keys():
			if processes[process]['kind'] == 'D':
				ics = []  # list of initial concentrations in the order of species
				out_file = ''
				for name in species:
					if name in processes[process]['reactants']:
						for i, reactant in enumerate(processes[process]['reactants']):    # starting point of adsorbed species
							if name == reactant:
								ics.append(np.floor([1/processes[process]['rstoichio'][i]]))
								out_file = 'TPD_' + str(name)
					else:
						ics.append(0.)
				n_elements = len([*rconditions.keys(), *species])
				if out_file not in tpd_done:
					# t_rate =   0.01, 0.1,  1, 5, 10  K/min
					for t_num in [60]:  #################[60000, 6000, 600, 120, 60]:
						t_eval = [0, 2*(10/t_num), 10/t_num]       # temperature rate in k/s --- from 0 to 10/t_num
						data = [*rconditions.keys(), *species]  # basic: Temp, time, species
						for temp_num in np.arange(50, 1010, 10): # integrate at different temperatures
							sol, sol_T = ode_solver(t_eval, species, rhs, ics, (temp_num,))
							data += sol_T[-n_elements:]  # transposed solution: (time x species)
							''' the initial concentration for the next temperature is the same as the last time 
							on the current temperature '''
							ics = sol.y.T[-1]
						data = np.array(data).reshape(-1, n_elements)  # n columns,
						# inferred rows
						printdata(out_file + "_" + str(round(t_eval[-1], 2)), data.tolist())
					tpd_done.append(out_file)
		print("\t\t\t\t", round((time.time() - start)/60, 3), " minutes")

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





