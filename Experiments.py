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

	print(f"Experiment: {experiment}, Data shape: {np.array(data).shape}")

	maxlen = [max([len(f"{data[r][c]}") + 2 for r in range(len(data))]) for c in range(len(data[0]))]  # 12/ length
	# per column
	folder = './KINETICS/DATA'
	outputfile = folder + "/" + str(experiment) + ".dat"
	if not pathlib.Path(folder).exists():
		pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
	output = open(outputfile, "a")
	output.write("#")
	for i in range(len(data[0])):
		output.write(" {val:>{wid}s}".format(wid=maxlen[i], val=str(data[0][i])))  # headings
	output.write("\n")
	for row in data[1:]:
		for i in range(len(row)):
			output.write(" {val:>{wid}.3{c}}".format(wid=maxlen[i], val= float(row[i]),
			                                         c='f' if 1e-5 < np.abs(float(row[i])) < 1e3 or 0. else 'e'))
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
		surface_values = [] ## no worthy to pass it through the subroutine, too many patches.
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
					method='BDF', jac=jac_numeric, rtol=1e-6, atol=1e-8)	#, events=steady_state_event)  # Alternative: "BDF" or "Radau"
	if not sol.success:
		raise RuntimeError(f"ODE solver failed: {sol.message}")
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
	def __init__(self, rconditions, systems, processes, equations, equation_factors):
		start = time.time()
		print("\t... Generating Concentrations and Rates ...")
		''' in arguments, the temperature numeric value should be the first entry (arg[0] '''
		ics, species, surfaces, adsorbates_by_surface = self.initial_species(processes, systems)
		rhs = []
		for name in species:        # lists of names with the order of ics
			dydt = 0
			for i in range(len(equations[name])):
				dydt += equation_factors[name][i]*equations[name][i]
			rhs.append(dydt)     # rhs for SciPy
		''' evaluate the Reaction Rates '''
		data = [*rconditions.keys(), *species, *surfaces]  # basic: Temp, time, species
		data_elements = len(data)
		rates_ss = [key for key in rconditions.keys() if key != "time"] + [*species]  # time not
		rates_elements = len(rates_ss)
		# required as it is at steady-state
		rates_avg = rates_ss.copy()   # time not required as it is an average along time
		drc_data = {}   # nested dictionary of temperatures and concentrations at time[-1]
		dsc_data = {}   # nested dictionary of temperatures and concentrations at time[-1]
		if isinstance(rconditions['temperature'], (int, float)): # single temperature
			sol, solution, sol_T = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface,  rhs, ics,
									(rconditions['temperature'],))
			data += sol_T
			drc_data[rconditions['temperature']] = sol.y
			dsc_data[rconditions['temperature']] = solution.y[:, -1]
			ss, avg = ConsTemperature.rates(systems, species, surfaces, adsorbates_by_surface, rhs, sol,
											rconditions['temperature'])
			rates_ss += ss
			rates_avg += avg
		else:      # temperature ramp
			for temp_num in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
				sol, solution, sol_T = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface,
												  rhs, ics, (temp_num,))


				print("TEMP:", temp_num, "SOL:", sol.y)


				data += sol_T  # transposed solution: (time x species)
				drc_data[str(temp_num)] = sol
				dsc_data[str(temp_num)] = solution.y[:, -1]
				ss, avg = ConsTemperature.rates(systems, species, surfaces, adsorbates_by_surface, equations, sol,
												temp_num)
				rates_ss += ss
				rates_avg += avg
		data = np.array(data).reshape(-1, data_elements)  # n columns, inferred rows
		printdata("Cons_Temperature", data)
		rates_ss = np.array(rates_ss).reshape(-1, rates_elements)
		rates_avg = np.array(rates_avg).reshape(-1, rates_elements)
		printdata("SteadyState_Rates", rates_ss)    # temperature x processes
		printdata("Average_Rates", rates_avg)  # temperature x processes
		ylabel = "TOF ($molecule \cdot site^{-1} \cdot s^{-1}$)"
		ConsTemperature.barplot("SteadyState Rates", "Species", ylabel, rates_ss, species, 0.5)
		ConsTemperature.barplot("Average Rates", "Species", ylabel, rates_avg, species, 0.5)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")

		start = time.time()
		print("\t... Generating Degree of Rate Control ...")
		ConsTemperature.degree_of_rate_control(rconditions, processes, systems, species, surfaces,
											   adsorbates_by_surface, equations, equation_factors, drc_data)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")
		#start = time.time()
		#print("\t... Generating Degree of Selectivity Control ...")

		#print("DSC_DATA\n", dsc_data)
		#ConsTemperature.degree_of_selectivity_control(rconditions, systems, processes, species, surfaces, dsc_data)
		#print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")

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
	def rates(systems, species, surfaces, adsorbates_by_surface, equations, sol, temp_num):
		''' rhs: list of sympy expressions r1(..), r2(..), ... using species symbols
		species: list of sympy symbols [theta_A, theta_B, ...] '''
		''' Find the surface converages at temp_num using sol'''
		surface_values = ConsTemperature.surface_coverages(sol, species, surfaces, adsorbates_by_surface, systems)
		''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species'''
		conditions = [t, temp]  # conditions required to lambdify
		rate_equations = []
		for name in species:
			r_i = 0
			for i in range(len(equations[name])):
				r_i += equations[name][i]
			rate_equations.append(r_i.subs(constants))
		rate_fn = sp.lambdify((*conditions, *species, *surfaces), rate_equations,["numpy", 'sympy'])

		if len(sol.t) < 2:		# Defensive conditions
			raise RuntimeError("ODE solver returned too few time points in Experiments.rates.")
		dt = sol.t[-1] - sol.t[0]
		if dt == 0 or not np.isfinite(dt):
			raise RuntimeError("Invalid time interval in ODE solution in Experiments.rates.")

		''' evaluating the rates as a function of time '''
		rate_time = np.array(rate_fn(sol.t, temp_num, *sol.y,  *surface_values), dtype=float).T # (times x species)
		if rate_time.ndim != 2:
			print("rate_time", rate_time)
			raise ValueError(f"rate_time should be 2D (time, species), got {rate_time.shape}")
		n_tail = max(1, int(0.1 * len(rate_time)))  # if rate_time !=0, get 10% of the last points
		rate_ss = [float(temp_num)] + [i for i in np.array(rate_time[-n_tail:]).mean(axis=0).tolist()]   # rates at the
		# steady-state, i.e. over the last 10% of the time points
		if np.any(~np.isfinite(rate_time)):
			raise RuntimeError("rate_time contains NaN or inf -- Experiments.rates.")
		''' Numpy versions > 2.0 uses "trapezoid" instead of "trapz" '''
		avg = np.trapz(rate_time, sol.t, axis=0) / dt  # rate averages along t_span using the trapezoidal rule
		rate_avg = [float(temp_num)] + avg.tolist()
		return rate_ss, rate_avg

	@staticmethod
	def barplot(experiment, x_label, y_label, data, labels, bar_width):
		icolour = ["b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki",
		   "teal"]
		ipatterns = ["///", "...", "xx", "**", "\\", "|", "--", "++", "oo", "OO"]

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
		#ax1.plot(x_limit, [1e-20, 1e-20], "k-", lw=1.5)
		ax1.set_xlim(x_limit)
		ax1.set_xlabel(x_label, fontsize=16)
		ax1.set_xticks(x)
		ax1.tick_params(axis='x', rotation=0, labelsize=14)
		ax1.set_xticklabels(labels, rotation=45, ha="center")

		ax1.set_ylim([1e-10, ax1.get_ylim()[1]])
		ax1.set_ylabel(y_label, fontsize=18)
		ax1.tick_params(axis='y', rotation=0, labelsize=16)

		leg_lines, leg_labels = ax1.get_legend_handles_labels()
		legend = ax1.legend(leg_lines[0:len(temps)], leg_labels[0:len(temps)], loc='best', fontsize=16)
		plt.title(experiment)
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
				coverage -= np.abs(sol.y[idx, -1]) * systems[name]["nsites"]
			if np.any(coverage < -1e-20) or np.any(coverage > 1 + 1e-20):  # the 1e-20 gives some wiggle
				raise ValueError(f"ERROR! Coverage {s_sym} is beyon the [0,1] limit.")
			surface_values.append(np.clip(coverage, 0.0, 1.))     # Enforce physical bounds on adsorbates
			# np.maximum(coverage, 0.0))  # clipping the lower limit
		return surface_values

	@staticmethod
	def degree_of_rate_control(rconditions, processes, systems, species, surfaces, adsorbates_by_surface, equations,
							   equation_factors, sol_base):
		''' The Degree of Rate Control (DRC), introduced by C. T. Campbell (J. Catal. 204, 520, 2001), quantifies how
		   sensitive the overall reaction rate is to each elementary step’s rate constant. '''
		k_list = [processes[process]['krate0'] for process in processes]    # the first process is "1"
		temp_list = [float(i) for i in sol_base.keys()]
		eps = 1e-6  # constant perturbation factor, small enough to retain linearity
		''' Find the reactants and forming products to define the net rate '''
		reactants = []
		products = []
		for name in systems:
			if systems[name]['kind'] == 'molecule':
				if systems[name]['pressure0'] > 0.0:
					reactants.append(name)
				elif systems[name]['pressure0'] == 0.0:
					products.append(name)
		rhs = []
		for name in species:  # lists of names with the order of ics
			dydt = 0
			for i in range(len(equations[name])):
				dydt += equation_factors[name][i] * equations[name][i]
			rhs.append(dydt)  # rhs for SciPy
		''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species'''
		conditions = [t, temp]  # conditions required to lambdify
		#f_ode = sp.lambdify((*conditions, *species, *surfaces), [i.subs(constants) for i in rhs], ["numpy",
		# 'sympy'])  #

		''' the selection of +i and -i could be improve but cheching reactants and products.'''
		labels = [f"{i}/{i+1}" for i in range(1, len(processes), 2)]
		data_i = {f'{name}': [[*list(rconditions.keys())[:-1], *labels]] for name in reactants+products} # at steady-state
		for temp_num in [temp_list[0] ]: #, temp_list[-1]]: # initial and final temperatures
			data_row = {f'{name}': [temp_num] for name in reactants+products}
			sol = sol_base[str(temp_num)]
			ics_local = sol.y[:, -1]
			''' Find the species concentrations and surface converages at temp_num using sol_base'''
			coverages = {f"{species[i]}": np.abs(sol.y[i, -1]) for i in range(len(species))}
			surface_values = ConsTemperature.surface_coverages(sol, species, surfaces, adsorbates_by_surface, systems)
			surfaces_num = {f"{surfaces[i]}": np.abs(surface_values[i][-1]) for i in range(len(surfaces))}
			''' DRC must be computed using a reaction rate, not a species balance.
			Correct choices:	rate of product formation
								rate of reactant consumption
								sum of elementary fluxes forming product  '''
			net_rate0 = {}
			for name in products:
				r_i = 0
				for i in range(len(equations[name])):
					r_i += equations[name][i]
				rate_equation = r_i.subs(constants)
				net_rate0[name] = np.float64(sp.lambdify((*conditions[1:], *species, *surfaces),
			                                    rate_equation, "numpy")(temp_num, **coverages, **surfaces_num))

				print(net_rate0[name])

				if net_rate0[name] <= 0.:
					print(name, rate_equation)
					print("temp", temp_num)
					print("coverages", coverages)
					print("surfaces", surfaces_num)
					exit()


				print(name, "rate0", net_rate0[name])

			for i in range(1, len(processes), 2):
				''' The partial derivative is taken holding constant the rate constants, k_j, for all 
				other steps j ≠ i and the equilibrium constant, K_i, for step i (and all other steps too, 
				since their forward and reverse rate constants are held fixed). Note that keeping K_i constant
				means that the forward and reverse rate constants for step i, k_i and k_–i, both must be 
				varied by equal factors so that their ratio remains constant. 
				Campbel, ACS Catal. 2017, 7, 4, 2770–2779 https://doi.org/10.1021/acscatal.7b00115 '''

				processes[str(i)]['krate0'] *= (1 + eps)    # the forward constant (1) + the perturbation
				processes[str(i+1)]['krate0'] *= (1 + eps)    # the forward constant (1) + the perturbation
				''' get the perturbed rhs: 
				because rhs derives from the ConsTemperature REquations (in Kinetics.py), processes[process][krate0]
				has to be modified and REquations recalled to generate the perturbed rhs '''
				sol, _, _ = ode_solver(systems, rconditions['time'], species, surfaces, adsorbates_by_surface, rhs,
				                       ics_local, (temp_num,))
				ics_local = sol.y[:, -1]
				''' Find the coverages, surface coverage and net rates at temp_num using previous sol '''
				local_coverages = {f"{species[s]}": np.abs(sol.y[s, -1]) for s in range(len(species))}
				surface_values = ConsTemperature.surface_coverages(sol, species, surfaces, adsorbates_by_surface,
				                                                   systems)
				local_surfaces_num = {f"{surfaces[s]}": np.abs(surface_values[s][-1]) for s in range(len(surfaces))}
				for name in products:
					r_i = 0
					for e in range(len(equations[name])):
						r_i += equations[name][e]
					rate_equation = r_i.subs(constants)

					net_rate1 = (sp.lambdify((*conditions[1:], *species, *surfaces), rate_equation, "numpy")
								 (temp_num, **local_coverages, **local_surfaces_num))




					if net_rate0[name] <= 0 or net_rate1 <= 0:
						a = np.nan
					else:
						a = float(np.log(net_rate1/net_rate0[name]) / np.log(1. + eps))

					print("RATES for ", name,  net_rate0[name], net_rate1, f"process {i}/{i + 1} at {temp_num} K =", a)

					data_row[name].append(a)

				''' reset the reaction constants to the original form'''
				processes[str(i)]['krate0'] = k_list[i-1]   # k_list starts from 0 but processes from 1
				processes[str(i+1)]['krate0'] = k_list[i]   # k_list starts from 0 but processes from 1
			for name in products:
				data_i[name].append(data_row[name])

		for p in products:
			printdata(p + "_Degree_of_Rate_Control", data_i[p])
			ylabel = f"$DRC_{{i,\ {p} }}$"

			print(f"DRC to barplot for {p}")

			ConsTemperature.barplot("", "Process Reactions", ylabel, data_i[p], labels, 0.5)
			for r in reactants:
				if isinstance(data_i[p], float) and isinstance(data_i[r], float):
					dsc = (np.array(data_i[p]) - np.array(data_i[r])).tolist()
					printdata(f'{p}/{r}' + "_Degree_of_Selectivity_Control", dsc)
					ylabel = f"$DSC_{{i,\ {p}/{r} }}$"
					ConsTemperature.barplot("", f"Degree of Selectivity Control for {p}/{r}", ylabel, dsc, labels, 0.5)

	@staticmethod
	def degree_of_selectivity_control(rconditions, systems, processes, species, surfaces, dsc_data):
		'''	From the previous definition of the degree of rate control, the degree of selectivity control for step i,
		 whereby the net rate is replaced with the selectivity to the desired product P from the most
		 valuable reactant R, S = rP/rR, where rP is the rate of production of P and rR is the rate of
		 consumption of R (scaled by the stoichiometric ratio of P:R) --> DSC_i = DRC_i,P - DRC_i,R
		Campbel, ACS Catal. 2017, 7, 4, 2770–2779 https://doi.org/10.1021/acscatal.7b00115 '''

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



		labels = [f"{i}" for i in range(1, len(processes)+1)]
		for name in products:
			data = [[*list(rconditions.keys())[:-1], *labels]]    # no time because it is at steady-state
			temp_list = [float(i) for i in dsc_data.keys()]
			for temp_num in [temp_list[0], temp_list[-1]]:  # initial and final temperatures
				print(f"\t\tDSC of {name} at {temp_num} K")
				row = [temp_num]
				row.extend([dsc_func[name][process](temp_num, *dsc_data[temp_num]) for process in processes])
				data.append(row)

			# print("SELECTIVITY: name", name, "data", data)

			printdata(name + "_Degree_of_Selectivity_Control", data)
			ylabel = "$DSC_{i}$"
			ConsTemperature.barplot(labels, name +" Degree of Rate Control", ylabel, data, labels, 0.5)


class TPR:
	def __init__(self, rconditions, systems, processes, equations, equation_factors):
		species, surfaces, adsorbates_by_surface = self.initial_species(processes, systems)     # list of species in the order on systems
		start = time.time()
		print("\t ... Generating TPR ...")
		rhs = []
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				dydt = 0
				for i in range(len(equations)):
					dydt += equation_factors[name][i]*equations[name][i]
				rhs.append(dydt)     # rhs for SciPy
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
