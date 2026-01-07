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
from types import SimpleNamespace




def printdata(experiment, data):
	maxlen = 12 #[max([len(f"{data[r][c]}") + 2 for r in range(len(data))]) for c in range(len(data[0]))]  # 12/ length
	# per column
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
												c='e' if float(row[i]) > 1e3 or
														1e-20 < np.abs(float(row[i])) < 1e-2 else 'f'))
		output.write("\n")
	output.close()

def ode_solver(systems, time, species, surfaces_s, adsorbates_by_surface, rhs, ics, arguments):
	''' time: tuple or list, e.g. (0, 10, 0.1)
		species: list of sympy symbols [theta_A, theta_B, ...]
		rhs: list of sympy expressions for d(species)/dt
		ics: list of initial conditions
		arguments: tuple or list of non-time parameters (e.g. temp, pressures, etc.) '''
	''' substitute the symbolic constants, e.g. h, kb, ..., by its values '''
	ics = np.asarray(ics, dtype=float)
	rhs = [i.subs(constants) for i in rhs]
	conditions = [t, temp]  # conditions required to lambdify

	''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species'''
	f_ode = sp.lambdify((*conditions, *species, *surfaces_s), rhs, ["numpy", 'sympy'])  #
	''' substitute rconditions'''
	t_span = (time[:2])  # time grid
	t_eval = np.arange(*time)

	''' Old versions of scipy do not accept args, so temp_num and other could be unwrapped here.'''
	def ode_system(t_num, y, *args):    # define ode function compatible with solve_ivp
		if len(args) < 1:
			raise ValueError("ERROR! Temperature not passed into ODE solver.")
		temp_num = args[0]
		y = np.clip(y, 0.0, None)     # Enforce physical bounds on adsorbates
		''' reconstruct algebraic expression to define Surface sites in the rhs'''
		surface_values = []
		for s_sym in surfaces_s:
			s_name = str(s_sym)
			coverage = 1.0    # Start with one full site
			for idx in adsorbates_by_surface[s_name]:	# Subtract adsorbates occupying this surface
				name = species[idx]     # idx is the position of the adsorbate in species
				coverage -= y[idx] * systems[name]["nsites"]
			if coverage > 1:
				raise ValueError(f"ERROR! Coverage {s_sym} is bigger than 1.")
			coverage = np.clip(coverage, 0.0, 1.0) # Prevent zero / negative free coverage
			surface_values.append(coverage)
		dydt = f_ode(t_num,  *[temp_num], *y, *surface_values)	# Evaluate RHS | There is ONLY t and Temp, for now
		if not np.all(np.isfinite(dydt)):	# Safety check
			raise RuntimeError(f"Non-finite RHS at t={t_num}\n", f"y={y}\n",
							   f"surface={surface_values}\n", f"dydt={dydt}")
		return np.array(dydt, dtype=float)  # .flatten()

	''' The jacobian provides stability and accelerate the ODE convergence '''
	def jac_numeric(t, y, temp_num):
		jacobian = np.zeros((len(y), len(y)))
		f0 = ode_system(t, y, temp_num)
		for i in range(len(y)):
			h = 1e-8 * max(1.0, abs(y[i]))  # “Perturb species i by a very small amount (10⁻⁸), recompute the RHS.
			y2 = y.copy()
			y2[i] += h
			jacobian[:, i] = (ode_system(t, y2, temp_num) - f0) / h
		return jacobian

	sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval,
					args=arguments if isinstance(arguments, tuple) else (arguments,),
					method='BDF', jac=jac_numeric, rtol=1e-6, atol=1e-10)	#, events=steady_state_event)  # Alternative: "BDF" or "Radau"
	if not sol.success:
		raise RuntimeError(f"ODE solver failed: {solution.message}")
	''' Build the augmented solution (gasses + adsorbates + surfaces) '''
	surface_values = ConsTemperature.surface_coverages(sol, species, surfaces_s, adsorbates_by_surface, systems)
	solution = SimpleNamespace()
	solution.t = sol.t
	solution.y = np.vstack([sol.y, *surface_values])
	solution.species = species + [str(s) for s in surfaces_s]
	solution.success = sol.success
	solution.message = sol.message
	data = np.column_stack([np.tile(arguments, (len(sol.t), 1)), solution.t,	solution.y.T])

	return sol, solution, data.flatten().tolist()   # gas + ads, gas + ads + surfaces


class ConsTemperature:
	def __init__(self, rconditions, systems, processes, equations):
		start = time.time()
		print("\t... Generating Concentrations and Rates ...")
		''' in arguments, the temperature numeric value should be the first entry (arg[0] '''
		ics, species, surfaces, adsorbates_by_surface = self.initial_species(processes, systems)
		rhs = []
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				rhs.append(equations[name])     # rhs for SciPy
		''' evaluate the Reaction Rates '''
		data = [*rconditions.keys(), *species, *surfaces]  # basic: Temp, time, species
		data_elements = len(data)
		rates_ss = [[key for key in rconditions.keys() if key != "time"] + [*processes.keys()]]  # time not required as it is at steady-state
		rates_avg = rates_ss.copy()   # time not required as it is an average along time
		dsc_data = {}   # nested dictionary of temperatures and concentrations at time[-1]
		if isinstance(rconditions['temperature'], (int, float)): # single temperature
			sol, solution, sol_T = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface,  rhs, ics,
									(rconditions['temperature'],))
			data += sol_T
			dsc_data[rconditions['temperature']] = solution.y[:, -1]
			ss, avg = ConsTemperature.rates(processes, species, rconditions['temperature'], sol)
			rates_ss.append(ss)
			rates_avg.append(avg)
		else:      # temperature ramp
			for temp_num in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
				sol, solution, sol_T = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface, rhs,
										ics, (temp_num,))
				data += sol_T  # transposed solution: (time x species)
				dsc_data[temp_num] = solution.y[:, -1]
				ss, avg = ConsTemperature.rates(processes, species, temp_num, sol)
				rates_ss.append(ss)
				rates_avg.append(avg)
		data = np.array(data).reshape(-1, data_elements)  # n columns, inferred rows
		printdata(str("Cons_Temperature"), data.tolist())

		labels = [process for process in processes]
		printdata("SteadyState_Rates", rates_ss)    # temperature x processes
		printdata("Average_Rates", rates_avg)  # temperature x processes
		ylabel = "TOF ($molecule \cdot site^{-1} \cdot s^{-1}$)"
		ConsTemperature.barplot(processes, "SteadyState Rates", ylabel, rates_ss, labels, 0.5)
		ConsTemperature.barplot(processes, "Average Rates", ylabel, rates_avg, labels, 0.5)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")

		start = time.time()
		print("\t... Generating Degree of Rate Control ...")
		ConsTemperature.degree_of_rate_control(rconditions, processes, systems, species, surfaces, adsorbates_by_surface, rhs, ics)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")
		start = time.time()
		print("\t... Generating Degree of Selectivity Control ...")

		print("DSC_DATA\n", dsc_data)

		ConsTemperature.degree_of_selectivity_control(rconditions, systems, processes, species, surfaces, dsc_data)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")

	@staticmethod
	def initial_species(processes, systems):  # process is processes[process]
		no_ts = []
		for process in processes:   # to avoid including Transition States
			no_ts.extend(processes[process]['reactants'])
			no_ts.extend(processes[process]['products'])
		ics_gases = []    # list of initial concentrations in the order of systems[names]
		ics_ads = []
		gases = []    # list of species in the order on systems
		ads = []
		surfaces = []
		for name in list(set(no_ts)):
			if systems[name]['kind'] == "molecule":
				gases.append(name)
			elif systems[name]['kind'] == "adsorbate":
				ads.append(name)
			elif systems[name]['kind'] == 'surface':  # Identify surface symbols and their adsorbates
				surfaces.append(name)
		for name in sorted(gases):
			ics_gases.append(systems[name]["pressure0"])
		for name in sorted(ads):
			ics_ads.append(systems[name]["coverage0"])
		species = sorted(gases) + sorted(ads)
		adsorbates_by_surface = {}
		for s in sorted(surfaces):
			adsorbates_by_surface[str(s)] = [i for i, name in enumerate(species) if
											 systems[name]["kind"] == "adsorbate"
											 and systems[name]["sites"] == systems[str(s)]["sites"]]

		return ics_gases + ics_ads, species, sorted(surfaces), adsorbates_by_surface

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
			- substitute the symbolic constants, e.g. h, kb, ..., by its values '''
		''' evaluating the rates as a function of time '''
		rate_fn = sp.lambdify((temp, *species), krate.subs(constants), ['numpy', 'sympy'])
		rate_time = []
		for tidx in range(sol.y.shape[1]):    # species values at time tidx
			rate_time.append(rate_fn(temp_num, *sol.y[:, tidx]))
		n_tail = max(1, int(0.1 * len(rate_time)))  # if rate_time !=0, get 10% of the last points
		rate_ss = np.array(rate_time[-n_tail:]).mean(axis=0)   # rates at the steady-state, i.e. over the last 10% of the time points
		# Defensive conditions
		if len(sol.t) < 2:
			raise RuntimeError("ODE solver returned too few time points in Experiments.rki_value.")
		dt = sol.t[-1] - sol.t[0]
		if dt == 0 or not np.isfinite(dt):
			raise RuntimeError("Invalid time interval in ODE solution in Experiments.rki_value.")
		if np.any(~np.isfinite(rate_time)):
			raise RuntimeError("rate_time contains NaN or inf -- Experiments.rki_value.")
		''' Numpy versions > 2.0 uses "trapezoid" instead of "trapz" '''
		rate_time = np.asarray(rate_time)
		rate_avg = np.trapz(rate_time, sol.t, axis=0) / dt  # rate averages along t_span using the trapezoidal rule to integrate the rate curve.
		return rate_ss, rate_avg

	@staticmethod
	def barplot(processes, experiment, y_label, data, labels, bar_width):
		icolour = ["b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki",
		   "teal"]
		ipatterns = ["///", "...", "xx", "**", "\\", "|", "--", "++", "oo", "OO"]

		x_label = experiment
		temps = [row[0] for row in [data[1], data[-1]]]  # temperatures; first row is for labels
		gap = 0.7   # gap between group of columns, e.g. processes
		x = np.arange(len(labels)) * (1 + gap)
		temp0 = [float(data[1][column+1]) for column in range(len(labels))] # for every process | column 0 is
		temp1 = [float(data[-1][column+1]) for column in range(len(labels))]   # lables is the list of processes
		y = [temp0, temp1]
		y = np.array(y, dtype=float) # Convert to array and apply positive floor for log-scale
		y[y <= 0] = 1e-20    # to prevent log crash

		fig, ax1 = plt.subplots(figsize=(10, 6), clear=True)
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
		ax1.tick_params(axis='x', rotation=0, labelsize=14)
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
	def surface_coverages(sol, species, surfaces, adsorbates_by_surface, systems):
		n_points = len(sol.t)
		surface_values = []
		for s_sym in surfaces:
			s_name = str(s_sym)
			coverage = np.ones(n_points)
			for idx in adsorbates_by_surface[s_name]:  # idx is the position of the adsorbate in species
				name = species[idx]  # name is the system, e.g, H2O
				coverage -= np.abs(sol.y[idx]) * systems[name]["nsites"]
			if np.any(coverage < -1e-20) or np.any(coverage > 1 + 1e-20):  # the 1e-20 gives some wiggle
				raise ValueError(f"ERROR! Coverage {s_sym} is beyon the [0,1] limit.")
			surface_values.append(np.maximum(coverage, 0.0))  # clipping the lower limit
		return surface_values

	@staticmethod
	def degree_of_rate_control(rconditions, processes, systems, species, surfaces, adsorbates_by_surface, rhs, ics):
		''' The Degree of Rate Control (DRC), introduced by C. T. Campbell (J. Catal. 204, 520, 2001), quantifies how
		   sensitive the overall reaction rate is to each elementary step’s rate constant. '''
		time = rconditions['time']
		k_list = [processes[process]['krate0'] for process in processes]    # the first process is "1"
		eps = 1e-3  # constant perturbation factor, small enough to retain linearity
		''' Find the forming product to define the net rate '''
		products = []
		for name in systems:
			if systems[name]['kind'] == 'molecule':
				if systems[name]['pressure0'] == 0.0:
					products.append(name)
		''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species'''
		ics = np.asarray(ics, dtype=float)
		rhs = [i.subs(constants) for i in rhs]
		conditions = [t, temp]  # conditions required to lambdify
		f_ode = sp.lambdify((*conditions, *species, *surfaces), rhs, ["numpy", 'sympy'])  #

		labels = [f"{i}" for i in range(1, len(processes))]
		data = [[*list(rconditions.keys())[:-1], *labels]]    # no time because it is at steady-state.
		for temp_num in [rconditions['temperature'][0], rconditions['temperature'][1]]: # initial and final temperatures
			sol_base, _, _ = ode_solver(systems, time, species, surfaces, adsorbates_by_surface, rhs, ics, (temp_num,))
			ics_local = sol_base.y[:, -1]
			''' Find the net_rate and surface converages at temp_num using sol_base'''
			surface_values = ConsTemperature.surface_coverages(sol_base, species, surfaces, adsorbates_by_surface, systems)
			dydt_ss = f_ode(sol_base.t[-1], temp_num, *sol_base.y[:, -1],  *surface_values)
			net_rate = dydt_ss[species.index(products[0])]

			data_row = [temp_num]
			drc = []
			for i in range(1, len(processes)):
				processes[str(i)]['krate0'] *= (1 + eps)    # the forward constant (1) + the perturbation
				''' because rhs derives from the ConsTemperature REquations (in Kinetics.py), processes[process][
				krate0] has to be modified and REquations recalled to generate the perturbed rhs '''
				equations = REquations(dict(processes), dict(systems)).constemperature
				''' get the perturbed rhs '''
				rhs_local = []
				for name in species:  # lists of names with the order of ics
					if name in equations.keys():
						rhs_local.append(equations[name])  # rhs for SciPy
				sol, _, _ = ode_solver(systems, time, species, surfaces, adsorbates_by_surface, rhs_local, ics_local, (temp_num,))
				ics_local = sol.y[:, -1]
				''' Find the net_rate and surface converage at temp_num using sol_base'''
				surface_values = ConsTemperature.surface_coverages(sol, species, surfaces, adsorbates_by_surface,
				                                                   systems)
				dydt_ss = f_ode(sol.t[-1], temp_num, *sol.y[:, -1], *surface_values)
				net_rate1 = dydt_ss[species.index(products[0])]
				''' compute the DRC '''
				drc.append((np.log(np.abs(net_rate1)) - np.log(np.abs(net_rate))) / np.log(1 + eps))
				''' reset the reaction constants to the original form'''
				processes[str(i)]['krate0'] = k_list[i-1]   # k_list starts from 0 but processes from 1
			data_row.extend(drc)
			data.append(data_row)
		printdata("Degree_of_Rate_Control", data)
		ylabel = "$DRC_{i}$"
		ConsTemperature.barplot(labels, "Degree of Rate Control", ylabel, data, labels, 0.5)

	@staticmethod
	def degree_of_selectivity_control(rconditions, systems, processes, species, surfaces, dsc_data):
		'''	Compute Campbell's Degree of Selectivity Control (DSC) for a product. '''
		# symbolic rate constants, one per reaction
		k_symbols = {process: sp.symbols(f'k_{process}') for process in processes}
		# Substitute actual expressions
		k_exprs = {f'k_{process}': processes[process]['krate0'].subs(constants) for process in processes}

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

		ds_dk = {name: {} for name in s_expr}  # nested dictionary of differentials of products and process
		dsc = {name: {} for name in s_expr}    # nested dictionary of products and process
		dsc_num = {name: {} for name in s_expr}    # nested dictionary with numeric dsc
		dsc_func = {name: {} for name in s_expr}    # nested lambdify function
		for name in products:
			# Symbolic derivatives dS/dk for each product and each k
			for process in processes:
				ds_dk[name][process] = sp.diff(s_expr[name], k_symbols[process])
				dsc[name][process] = ds_dk[name][process] * k_symbols[process] / s_expr[name]
				dsc_num[name][process] = dsc[name][process].subs(k_exprs)
				dsc_func[name][process] = sp.lambdify((temp, *species, *surfaces), dsc_num[name][process], "numpy")




				#print("DSC_TYPE\n", type(dsc_func[name][process](450.0)))
				free_syms = dsc_num[name][process].free_symbols
				print("DSC SYMBOLS - ", f"NAME: {name} ||| PROCESS: {process}:", free_syms)



		labels = [f"{i}" for i in range(1, len(processes))]
		for name in products:
			data = [[*list(rconditions.keys())[:-1], *labels]]    # no time because it is at steady-state
			for temp_num in [rconditions['temperature'][0], rconditions['temperature'][1]]: # nitial and final temperatu
				row = [temp_num]
				row.extend([dsc_func[name][process](temp_num, *dsc_data[temp_num]) for process in processes])
				data.append(row)

			# print("SELECTIVITY: name", name, "data", data)

			printdata(name + "_Degree_of_Selectivity_Control", data)
			ylabel = "$DSC_{i}$"
			ConsTemperature.barplot(labels, name +" Degree of Rate Control", ylabel, data, labels, 0.5)


class TPR:
	def __init__(self, rconditions, systems, processes, equations):
		species, surfaces, adsorbates_by_surface = self.initial_species(processes, systems)     # list of species in the order on systems
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
								ics.append(float(np.floor(1/processes[process]['rstoichio'][i])))
								out_file = 'TPD_' + str(name)
					else:
						ics.append(0.)
				ics = np.asarray(ics, dtype=float)
				n_elements = len([*rconditions.keys(), *species])
				if out_file not in tpd_done:
					# t_rate =      1, 10  K/min
					for t_num in [600, 60]:
						t_step = (10/t_num) / 100	# it is recomended to have at least 100 stepts
						t_eval = [0, 10/t_num+t_step, t_step]       # temperature rate in k/s --- from 0 to 10/t_num
						data = [*rconditions.keys(), *species]  # basic: Temp, time, species
						for temp_num in np.arange(50, 1010, 10): # integrate at different temperatures
							sol, _, sol_T = ode_solver(systems, t_eval, species, surfaces, adsorbates_by_surface, rhs, ics, (temp_num,))
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
	def initial_species(processes, systems):
		no_ts = []
		for process in processes:   # to avoid including Transition States
			no_ts.extend(processes[process]['reactants'])
			no_ts.extend(processes[process]['products'])
		gases = []    # list of species in the order on systems
		ads = []
		surfaces = []
		for name in list(set(no_ts)):
			if systems[name]['kind'] == "molecule":
				gases.append(name)
			elif systems[name]['kind'] == "adsorbate":
				ads.append(name)
			elif systems[name]['kind'] == 'surface':	# Identify surface symbols and their adsorbates
				surfaces.append(name)
		species = sorted(gases) + sorted(ads)
		adsorbates_by_surface = {}
		for s in sorted(surfaces):
			adsorbates_by_surface[str(s)] = [i for i, name in enumerate(species) if
												 systems[name]["kind"] == "adsorbate"
												 and systems[name]["sites"] == systems[str(s)]["sites"]]
		return species, sorted(surfaces), adsorbates_by_surface
