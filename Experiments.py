"""
	This script builds on the perl version by A.Roldan.

"""

import pathlib
import time
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

from Symbols_def import t, temp, constants
from types import SimpleNamespace
from Diagnostics import Diagnostics



def printdata(experiment, data):
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
	data = np.column_stack([np.tile(arguments, (len(sol.t), 1)), solution.t,	solution.y.T])

	return sol, data   # gas + ads, gas + ads + surfaces

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
		print("\t... Solving ODE System ...")
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
		data = np.array([*rconditions.keys(), *species, *surfaces])  # basic: Temp, time, species
		if isinstance(rconditions['temperature'], (int, float)): # single temperature
			sol, solution = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface,
											  gas_number, rhs, ics,	(rconditions['temperature'],))

			data = np.vstack([data, solution])
		else:      # temperature ramp
			for temp_num in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
				sol, solution = ode_solver(systems, rconditions["time"], species, surfaces, adsorbates_by_surface,
												  gas_number, rhs, ics, (temp_num,))
				data = np.vstack([data, solution])  # transposed solution: (time x species)
		''' Cons_Temp Experiments '''
		printdata("Cons_Temperature", data)
		print("\t\t\t\t", round((time.time() - start) / 60, 3), " minutes")

		print("\t... Generating Diagnostics ...")
		start = time.time()
		Diagnostics(systems, processes, species, equation_factors, surfaces, adsorbates_by_surface,
							gas_number, surface_expr, data)
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
