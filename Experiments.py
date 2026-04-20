"""
	This script builds on the perl version by A.Roldan.

"""

import pathlib
import time
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from Symbols_def import t, temp, kb, JtoeV, constants
#from Kinetics import REquations
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
			a = np.abs(float(row[i]))
			output.write(" {val:>{wid}.3{c}}".format(wid=maxlen[i], val= a, c='f' if 1e-3 < a < 1e3 or a == 0. else
			'e'))
		output.write("\n")
	output.close()

def ode_solver(systems, time, species, surfaces_s, adsorbates_by_surface, gas_number, rhs, ics, arguments):
	''' time: tuple or list, e.g. (0, 10, 0.1)
		species: list of sympy symbols [theta_A, theta_B, ...]
		rhs: list of sympy expressions for d(species)/dt
		ics: list of initial conditions
		arguments: tuple or list of non-time parameters (e.g. temp, pressures, etc.) '''
	''' substitute the symbolic constants, e.g. h, kb, ..., by its values '''
	ics = np.asarray(ics, dtype=float)
	#rhs = [i.subs(constants) for i in rhs]	# already done in ConsTEMP
	conditions = [t, temp]  # conditions required to lambdify

	''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species
	surfaces have been substituted in the rhs '''
	f_ode = sp.lambdify((*conditions, *species), rhs, ["numpy", 'sympy'])  #
	jacobian = sp.Matrix(rhs).jacobian(species)
	f_jac = sp.lambdify((*conditions, *species), jacobian, ["numpy", "sympy"])
	''' substitute rconditions'''
	t_span = (time[:2])  # time grid
	t_eval = np.arange(*time)

	''' Old versions of scipy do not accept args, so temp_num and other could be unwrapped here.'''
	def ode_system(t_num, y, *args):    # define ode function compatible with solve_ivp
		if len(args) < 1:
			raise ValueError("ERROR! Temperature not passed into ODE solver.")
		dydt = f_ode(t_num, *args, *y)	# Evaluate RHS | There is ONLY t and Temp, for now
		if not np.all(np.isfinite(dydt)):	# Safety check
			raise RuntimeError(f"Non-finite RHS at t={t_num}\n", f"y={y}\n", f"dydt={dydt}")
		return np.array(dydt, dtype=float)

	''' The jacobian provides stability and accelerate the ODE convergence '''
	def jac_analytic(t_num, y, *args):
		jac_eval = f_jac(t_num, *args, *y)
		return np.array(jac_eval, dtype=float)

	sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval,
					args=arguments if isinstance(arguments, tuple) else (arguments,),
					method='LSODA', jac=jac_analytic, rtol=1e-6, atol=1e-8)	#, rtol <=1e-4 | BDF or Radau
	if not sol.success:
		raise RuntimeError(f"ODE solver failed: {sol.message}")

	''' Build the augmented solution (gasses + adsorbates + surfaces) '''
	sol.y = clip_species(surfaces_s, adsorbates_by_surface, gas_number, sol.y.copy())
	surface_values = surface_coverages(sol.y, species, surfaces_s,	adsorbates_by_surface, systems)
	solution = SimpleNamespace()
	solution.t = sol.t
	solution.y = np.vstack([sol.y, *surface_values])
	solution.species = species + [str(s) for s in surfaces_s]
	solution.success = sol.success
	solution.message = sol.message
	data = np.column_stack([np.tile(arguments, (len(sol.t), 1)), solution.t,	solution.y.T])

	return sol, solution, data.flatten().tolist()   # gas + ads, gas + ads + surfaces

def check_adsorbates(systems, species, surfaces_s, adsorbates_by_surface, y, t_check):
	''' Print adsorbates' coverage '''
	for s_sym in surfaces_s:
		s_name = str(s_sym)
		message = []
		a_sum = 0
		for idx in adsorbates_by_surface[s_name]:  # idx is the position of the adsorbate in species
			name = species[idx]  # name is the system, e.g, H2O
			message.append(f"{species[idx]}:{np.abs(y[idx, t_check]) * systems[name]["nsites"]}")
			a_sum += np.abs(y[idx, t_check]) * systems[name]["nsites"]
	print(f"ADSORBATES(t={t_check}): {message} --> SUM: {a_sum}")

def clip_species(surfaces_s, adsorbates_by_surface, gas_number, y):
	for s_sym in surfaces_s:
		s_name = str(s_sym)
		for idx in adsorbates_by_surface[s_name]:  # idx is the position of the adsorbate in species
			y[idx] = np.clip(y[idx], 0.0, 1.0)  # adsorbates between 0 and 1
	for idx in gas_number:
		y[idx] = np.clip(y[idx], 0.0, None)  # ensures gases > 0
	return y

def surface_coverages(y, species, surfaces, adsorbates_by_surface, systems):
	surface_values = []
	for s_sym in surfaces:
		s_name = str(s_sym)
		coverage = np.ones(y.shape[1])
		message = []
		for idx in adsorbates_by_surface[s_name]:  # idx is the position of the adsorbate in species
			name = species[idx]  # name is the system NOT symbol, e.g, H2O
			coverage -= y[idx, :] * systems[name]["nsites"]
			message.append(f"{species[idx]}:{y[idx] * systems[name]["nsites"]}")
		#if np.any(coverage < -1e-8) or np.any(coverage > 1. + 1e-8):  # the 1e-8 some wiggle for optimisation
		#	check_adsorbates(systems, species, surfaces, adsorbates_by_surface, y, -1)
		#	raise ValueError(f"ERROR! Coverage {s_name}: min={coverage.min()}, max={coverage.max()}")
		surface_values.append(np.clip(coverage, 0.0, 1.0))  # Enforce physical meaning
	return surface_values

class ConsTemperature:
	def __init__(self, rconditions, systems, processes, equations, equation_factors):
		start = time.time()
		print("\t... Generating Concentrations and Rates ...")
		''' in arguments, the temperature numeric value should be the first entry (arg[0] '''
		ics, species, surfaces, gas_number, adsorbates_by_surface = self.initial_species(processes, systems)

		''' reconstruct algebraic expression to define Surface sites '''
		surface_expr = {}	 ## no worthy to pass it through the subroutine, too many patches.
		for s_sym in surfaces:
			s_name = str(s_sym)
			coverage_expr = 1.0    # Start with one full site
			for idx in adsorbates_by_surface[s_name]:	# Subtract adsorbates occupying this surface
				name = species[idx]     # idx is the position of the adsorbate in species
				coverage_expr -=  sp.symbols(name) * systems[str(name)]["nsites"]	#	symbolic
			surface_expr[s_sym] = coverage_expr

		''' Defines rhs and substitutes the contants and Surface sites'''
		rhs = []
		for name in species:        # lists of names with the order of ics
			dydt = 0
			for i in range(len(equations[name])):
				dydt += equation_factors[name][i] * equations[name][i]
			rhs.append(dydt.subs(constants))     # rhs for SciPy
		rhs = [i.subs(surface_expr) for i in rhs]		# This replaces ALL surfaces simultaneously.

		''' evaluate the Reaction Rates '''
		data = [*rconditions.keys(), *species, *surfaces]  # basic: Temp, time, species
		data_elements = len(data)
		rates_ss = [key for key in rconditions.keys() if key != "time"] + [*species]  # time not
		rates_elements = len(rates_ss)		# required as it is at steady-state
		rates_avg = rates_ss.copy()   # time not required as it is an average along time
		drc_data = {}   # nested dictionary of temperatures and concentrations at time[-1]
		dsc_data = {}   # nested dictionary of temperatures and concentrations at time[-1]
		if isinstance(rconditions['temperature'], (int, float)): # single temperature
			sol, solution, sol_T = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface,
											  gas_number, rhs, ics,	(rconditions['temperature'],))
			data += sol_T
			drc_data[rconditions['temperature']] = sol.y
			dsc_data[rconditions['temperature']] = solution.y[:, -1]
			ss, avg = ConsTemperature.rates(species, rhs, sol, rconditions['temperature'])
			rates_ss += ss
			rates_avg += avg
		else:      # temperature ramp
			for temp_num in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
				sol, solution, sol_T = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface,
												  gas_number, rhs, ics, (temp_num,))
				data += sol_T  # transposed solution: (time x species)
				drc_data[str(temp_num)] = sol
				dsc_data[str(temp_num)] = solution.y[:, -1]
				ss, avg = ConsTemperature.rates(species, rhs, sol, temp_num)
				rates_ss += ss
				rates_avg += avg
		''' Cons_Temp Experiments '''
		data = np.array(data).reshape(-1, data_elements)  # n columns, inferred rows
		printdata("Cons_Temperature", data)
		''' Rates Experiments '''
		rates_ss = np.array(rates_ss).reshape(-1, rates_elements)
		rates_avg = np.array(rates_avg).reshape(-1, rates_elements)
		printdata("SteadyState_Rates", rates_ss)    # temperature x processes
		printdata("Average_Rates", rates_avg)  # temperature x processes
		ylabel = "TOF ($molecule \cdot site^{-1} \cdot s^{-1}$)"
		ConsTemperature.barplot("SteadyState Rates", "Species", ylabel, rates_ss, species, 0.5)
		ConsTemperature.barplot("Average Rates", "Species", ylabel, rates_avg, species, 0.5)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")

		print("DATA RATES", rates_avg)


		start = time.time()
		print("\t... Generating Degree of Rate Control ...")
		ConsTemperature.analytical_degree_of_rate_control(systems, processes, species, equation_factors, drc_data,
														  surfaces, adsorbates_by_surface, gas_number)
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
		gas_number = [i for i, name in enumerate(gases)]

		return ics_gases + ics_ads, species, sorted(surfaces), gas_number, adsorbates_by_surface

	@staticmethod
	def rates(species, rhs, sol, temp_num):
		''' rhs: list of sympy expressions r1(..), r2(..), ... using species symbols
		species: list of sympy symbols [theta_A, theta_B, ...] '''
		if len(sol.t) < 2:		# Defensive conditions
			raise RuntimeError("ODE solver returned too few time points in Experiments.rates.")
		dt = sol.t[-1] - sol.t[0]
		if dt == 0 or not np.isfinite(dt):
			raise RuntimeError("Invalid time interval in ODE solution in Experiments.rates.")

		''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species 
		 No need to include the surfaces because rhs has them substituted with adsorbates. '''
		conditions = [t, temp]  # conditions required to lambdify
		rate_fn = sp.lambdify((*conditions, *species), rhs,["numpy", 'sympy'])

		''' evaluating the rates as a function of time '''
		rate_time = np.array(rate_fn(sol.t, temp_num, *sol.y), dtype=float).T # (times x species)
		if rate_time.ndim != 2:
			print("rate_time", rate_time)
			raise ValueError(f"rate_time should be 2D (time, species), got {rate_time.shape}")
		n_tail = max(1, int(0.1 * len(rate_time)))  # if rate_time !=0, get 10% of the last points
		rate_ss = [float(temp_num)] + [i for i in np.array(rate_time[-n_tail:]).mean(axis=0).tolist()]
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
		ax1.set_xlabel(x_label, fontsize=18)
		ax1.set_xticks(x)
		ax1.tick_params(axis='x', rotation=0, labelsize=14)
		ax1.set_xticklabels(labels, rotation=45, ha="right")

		#ax1.set_ylim([1e-10, ax1.get_ylim()[1]])
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
	def analytical_degree_of_rate_control(systems, processes, species, equation_factors, sol_base, surfaces,
										  adsorbates_by_surface, gas_number):
		''' The Degree of Rate Control (DRC), introduced by C. T. Campbell (J. Catal. 204, 520, 2001), quantifies how
		   sensitive the overall reaction rate is to each elementary step’s rate constant. '''
		''' DRC must be computed using a reaction rate, not a species balance.
		Correct choices:	rate of product formation
							rate of reactant consumption
							sum of elementary fluxes forming product '''
		# --- Symbols ---
		species_syms = {s: sp.Symbol(s) for s in species}
		k_list = list(processes.keys())
		k_syms = {pr: sp.Symbol(f"k_{pr}") for pr in processes}
		k_sym_list = [k_syms[k] for k in k_list]
		concentrations_vec = sp.Matrix([species_syms[s] for s in species])
		all_variables = (temp, *concentrations_vec, *k_sym_list)

		# --- Reaction Step pairing ---
		step_pairs = [(k_list[i], k_list[i + 1]) for i in range(0, len(k_list), 2)]
		if not step_pairs:
			raise ValueError("No forward/backward step pairs found")

		data_heading = ["Temperature"] + [f"Step_{i}" for i in step_pairs]

		# --- Rate constants ---
		k_func = {k: sp.lambdify(temp, processes[k]['krate0'].subs(constants), 'numpy') for k in k_list}
		def keq_func(temp_num):	# Equilibriuma constants : DRC and DSC per forward/backward
			kf = np.array([k_func[j_f](temp_num) for j_f, _ in step_pairs])
			kb = np.array([k_func[j_b](temp_num) for _, j_b in step_pairs])
			return kf / kb

		# --- Surface expressions ---
		surface_expr = {}
		for s_sym in surfaces:
			s_name = str(s_sym)
			coverage_expr = 1.0
			for idx in adsorbates_by_surface[s_name]:
				name = species[idx]
				coverage_expr -= species_syms[name] * systems[str(name)]["nsites"]
			surface_expr[s_sym] = coverage_expr

		# --- Build RHS ---
		f_list = []
		for name in species:
			dydt = 0
			for j, pr in enumerate(k_list):
				expr = k_syms[pr]
				for i, r in enumerate(processes[pr]['reactants']):
					if r in species_syms:
						expr *= species_syms[r] ** processes[pr]['rstoichio'][i]
					else:
						expr *= surface_expr[r] ** processes[pr]['rstoichio'][i]
				dydt += equation_factors[name][j] * expr
			f_list.append(dydt)
		f_vec = sp.Matrix(f_list)
		f_func = sp.lambdify(all_variables, f_vec, "numpy")

		# --- Jacobian ---
		j_sym = f_vec.jacobian(concentrations_vec)
		j_func = sp.lambdify(all_variables, j_sym, "numpy")

		# --- Products ---
		products = [name for name in systems if systems[name]['kind'] == 'molecule' and systems[name]['pressure0'] == 0]
		rate_exprs = {p: f_vec[species.index(p)] for p in products}
		r_funcs = {p: sp.lambdify(all_variables, rate_exprs[p], "numpy") for p in products}

		# --- ∂r/∂θ : Matrix (Np × Ns) ---
		dr_dtheta_sym = sp.Matrix([[sp.diff(rate_exprs[p], s) for s in concentrations_vec] for p in products])
		dr_dtheta_func = sp.lambdify(all_variables, dr_dtheta_sym, "numpy")

		# --- ∂f/∂k : Matrix (Ns × Nk) with thermodynamic constraint ---
		df_dkf_cols = []
		df_dkb_cols = []
		for pr_f, pr_b in step_pairs:
			kf = k_syms[pr_f]
			kb = k_syms[pr_b]	# it is dependent of kf and link through keq
			term_f = kf * f_vec.diff(kf)
			term_b = kb * f_vec.diff(kb)
			df_dkf_cols.append(term_f)
			df_dkb_cols.append(term_b)
		df_dkf_sym = sp.Matrix.hstack(*df_dkf_cols)
		df_dkb_sym = sp.Matrix.hstack(*df_dkb_cols)
		df_dkf_func = sp.lambdify(all_variables, df_dkf_sym, "numpy")
		df_dkb_func = sp.lambdify(all_variables, df_dkb_sym, "numpy")

		# --- ∂r/∂k direct : Matrix (Np × Nk) ---
		dr_dkf_cols = []
		dr_dkb_cols = []
		for pr_f, pr_b in step_pairs:
			kf = k_syms[pr_f]
			kb = k_syms[pr_b]
			term_f = sp.Matrix([kf * sp.diff(rate_exprs[p], kf) for p in products])
			term_b = sp.Matrix([kb * sp.diff(rate_exprs[p], kb) for p in products])
			dr_dkf_cols.append(term_f)
			dr_dkb_cols.append(term_b)
		dr_dkf_direct_sym = sp.Matrix.hstack(*dr_dkf_cols)
		dr_dkb_direct_sym = sp.Matrix.hstack(*dr_dkb_cols)
		dr_dkf_direct_func = sp.lambdify(all_variables, dr_dkf_direct_sym, "numpy")
		dr_dkb_direct_func = sp.lambdify(all_variables, dr_dkb_direct_sym, "numpy")

		# --- Results ---
		drc_results = {p: {} for p in products}
		drc_data = {p: [] for p in products}
		dsc_results = {p: {} for p in products}
		dsc_data = {p: [] for p in products}

		# --- Loop for Temperatures ---
		temp_list = [float(t) for t in sol_base.keys()]
		for temp_num in temp_list:
			conc_guess = sol_base[str(temp_num)].y[:, -1]
			# --- thermodynamic constraint: forward rates are independent & backward rates follow thermodynamics
			keq_vals = keq_func(temp_num)
			k_values = np.array([k_func[k](temp_num) for k in k_func])
			for idx, (pr_f, pr_b) in enumerate(step_pairs):
				j_f = k_list.index(pr_f)
				j_b = k_list.index(pr_b)
				k_values[j_b] = k_values[j_f] / keq_vals[idx]

			# --- steady state ---
			def steady(theta):
				vals = (temp_num, *theta, *k_values)
				return np.array(f_func(*vals), dtype=float).flatten()

			sol = root(steady, conc_guess, method='lm')
			concentrations = clip_species(surfaces, adsorbates_by_surface, gas_number, sol.x)
			all_vals = (temp_num, *concentrations, *k_values)

			# --- Evaluate matrices ---
			jacobian = np.array(j_func(*all_vals), dtype=float)		# J
			df_dkf = np.array(df_dkf_func(*all_vals), dtype=float)	# B forwards
			df_dkb = np.array(df_dkb_func(*all_vals), dtype=float) 	# B backwards
			df_dk = df_dkf + (1.0 / keq_vals) * df_dkb
			dr_dkf = np.array(dr_dkf_direct_func(*all_vals), dtype=float)	# H forwards
			dr_dkb = np.array(dr_dkb_direct_func(*all_vals), dtype=float)	# H backwards
			dr_dk_direct = dr_dkf + (1.0 / keq_vals) * dr_dkb
			dr_dconc = np.array(dr_dtheta_func(*all_vals), dtype=float)	# G

			# --- Solve sensitivities ---
			dtheta_dk = -np.linalg.solve(jacobian, df_dk)

			# --- Rate vector ---
			r_vec = np.array([r_funcs[p](*all_vals) for p in products], dtype=float)

			# --- Total derivative ---
			dr_dk_total = dr_dk_direct + dr_dconc @ dtheta_dk  # Matrix (Np × Nk)

			# --- DRC ---
			x = dr_dk_total / r_vec[:, None]

			# --- diagnostic ---
			print("sum DRC:", np.sum(X, axis=1))
			print("r_vec:", r_vec)
			print("cond(J):", np.linalg.cond(jacobian))
			print("||df_dk||:", np.linalg.norm(df_dk))
			print("||dtheta_dk||:", np.linalg.norm(dtheta_dk))

			# --- Store DRC ---
			for i, p in enumerate(products):
				row = [temp_num]
				for j, (pr_f, _) in enumerate(step_pairs):
					drc_results[p][pr_f] = x[i, j]
					row.append(float(x[i, j]))
				drc_data[p].append(row)
				print(drc_data[p])


			# --- DSC ---
			r_tot = np.sum(r_vec)
			if abs(r_tot) < 1e-20:
				continue

			s_vec = r_vec / r_tot
			kdr = dr_dk_total
			kdr_sum = np.sum(kdr, axis=0)
			dsc = (kdr - s_vec[:, None] * kdr_sum) / (s_vec[:, None] * r_tot)
			for i, p in enumerate(products):
				for j, (pr_f, _) in enumerate(step_pairs):
					dsc_results[p][pr_f] = dsc[i, j]

		for p in products:
			printdata(f"{p}_Degree_of_Rate_Control", data_heading + drc_data[p])
			ylabel = f"$DRC_{p}$)"
			ConsTemperature.barplot(f"{p}_DRC", "Reaction Step", ylabel, drc_results[p], species, 0.5)
			printdata(f"{p}_Degree_of_Selectivity_Control", dsc_results[p])
			ylabel = f"$DSC_{p}$)"
			ConsTemperature.barplot(f"{p}_DSC", "Reaction Step", ylabel, drc_results[p], species, 0.5)


class TPR:
	def __init__(self, rconditions, systems, processes, equations, equation_factors):
		species, surfaces, adsorbates_by_surface, gas_number = self.initial_species(processes, systems)     # list of
		# species in the order on systems



		start = time.time()
		print("\t ... Generating TPR ...")
		rhs = []
		for name in species:        # lists of names with the order of ics
			dydt = 0
			for i in range(len(equations[name])):
				dydt += equation_factors[name][i] * equations[name][i]
			rhs.append(dydt)     # rhs for SciPy
		rhs = [i.subs(constants) for i in rhs]
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
							sol, _, sol_T = ode_solver(systems, t_eval, species, surfaces, adsorbates_by_surface,
							                                  gas_number, rhs, ics, (temp_num,))

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
		gas_number = [i for i, name in enumerate(gases)]
		return species, sorted(surfaces), adsorbates_by_surface, gas_number
