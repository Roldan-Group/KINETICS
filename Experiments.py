"""
	This script builds on the perl version by A.Roldan.

"""

import pathlib
import time
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from Symbols_def import t, temp, constants, chem_label
from types import SimpleNamespace
from Diagnostics import Diagnostics



def printdata(experiment, data):
	maxlen = [max([len(f"{data[r][c]}") + 2 for r in range(len(data))]) for c in range(len(data[0]))]  # 12/ length
	# per column
	folder = experiment[0]
	outputfile = folder + "/" + str(experiment[1]) + ".dat"
	if not pathlib.Path(folder).exists():
		pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
	output = open(outputfile, "a")
	output.write("#")
	for i in range(len(data[0])):
		output.write(" {val:>{wid}s}".format(wid=maxlen[i], val=str(data[0][i])))  # headings
	output.write("\n")
	for row in data[1:]:	# first column
		for i in range(len(row)):
			a = float(row[i])
			output.write(" {val:>{wid}.3{c}}".format(wid=maxlen[i], val= a, c='f' if 1e-3 < a < 1e3 or a == 0. else
			'e'))
		output.write("\n")
	output.close()

def safe_exp(x):
	return np.exp(np.clip(x, -700, 700))

def ode_solver(systems, t_span, t_eval, species, clipping, ics, rhs_func, jac_func, arguments):
	''' time: tuple or list, e.g. (0, 10, 0.1)
		species: list of sympy symbols [theta_A, theta_B, ...]
		rhs: list of sympy expressions for d(species)/dt
		ics: list of initial conditions
		arguments: tuple or list of non-time parameters (e.g. temp, pressures, etc.) '''
	''' substitute the symbolic constants, e.g. h, kb, ..., by its values '''
	ics = np.asarray(ics, dtype=float)
	surfaces_s, adsorbates_by_surface, gas_number = clipping

	''' Old versions of scipy do not accept args, so temp_num and other could be unwrapped here.'''
	def ode_system(t_num, y, *args):    # define ode function compatible with solve_ivp
		dydt = rhs_func(t_num, *args, *y)	# --- Evaluate RHS
		if not np.all(np.isfinite(dydt)):		# --- Safety check ---
			raise RuntimeError(f"Non-finite RHS at t={t_num}\n y={y}\n dydt={dydt}")
		return np.array(dydt, dtype=float)

	''' The jacobian provides stability and accelerate the ODE convergence '''
	def jac_analytic(t_num, y, *args):
		jac_eval = jac_func(t_num, *args, *y)
		return np.array(jac_eval, dtype=float)

	sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval,
					args=arguments if isinstance(arguments, tuple) else (arguments,),
					method='LSODA', jac=jac_analytic, rtol=1e-6, atol=1e-8) 	# BDF or LSODA
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
	theta_max = 1.0
	for s_sym in surfaces_s:
		s_name = str(s_sym)
		adsorbate_indices = adsorbates_by_surface[s_name]  # idx is the position of the adsorbate in species
		theta_total = np.sum(y[adsorbate_indices])
		if theta_total > theta_max:
			y[adsorbate_indices] *= theta_max / theta_total
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

class Isothermal:
	def __init__(self, rconditions, systems, processes, equations, equation_factors):
		start = time.time()
		print("\t... Solving ODE System ...")
		''' in arguments, the temperature numeric value should be the first entry (arg[0] '''
		conditions = [t, temp]  # conditions required to lambdify
		ics, species, surfaces, gas_number, adsorbates_by_surface = self.initial_species(processes, systems)
		clipping = (surfaces, adsorbates_by_surface, gas_number)

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
		''' Convert symbolic into numerical --- ics is the initial concentrations in the same order than species
		surfaces have been substituted in the rhs '''
		rhs_func = sp.lambdify((*conditions, *species), rhs,  [{"exp": safe_exp}, "numpy"] )
		jacobian = sp.Matrix(rhs).jacobian(species)
		jac_func = sp.lambdify((*conditions, *species), jacobian,  [{"exp": safe_exp}, "numpy"])
		''' substitute rconditions'''
		t_span = (rconditions["time"][:2])  # time grid
		t_eval = np.arange(*rconditions["time"])

		''' evaluate the Reaction Rates '''
		data = np.array([*rconditions.keys(), *species, *surfaces])  # basic: Temp, time, species
		if isinstance(rconditions['temperature'], (int, float)): # single temperature
			sol, solution = ode_solver(systems, t_span, t_eval, species, clipping, ics, rhs_func, jac_func,
									   (rconditions['temperature'],))
			data = np.vstack([data, solution])
		else:      # temperature ramp
			for temp_num in np.arange(*rconditions["temperature"]):     # integrate at different temperatures
				sol, solution = ode_solver(systems, t_span, t_eval, species, clipping, ics, rhs_func, jac_func,
										   (temp_num,))
				data = np.vstack([data, solution])  # transposed solution: (time x species)
		''' Isothermal Experiments '''
		printdata(['./KINETICS/ISOTHERMAL/', "Isothermal"], data)
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
	def __init__(self, systems, processes, equations, equation_factors):
		''' In Temperature Desorption experiments, the initial conditions consider only adsorbed molecules.
		Therefore for Reaction type A, the coverage of the product will be 1.
		There will be a TPD/TPR for each of the molecules in process[kind]=A
		> Integrate in TIME, not temperature
		 T(t) = T0 + beta*t
		 dtheta/dt = R(theta, T(t)) '''
		start = time.time()
		print("\t ... Generating TPR ...")

		conditions = [temp]  # conditions required to lambdify: no time
		species, surfaces, adsorbates_by_surface, gas_number = self.initial_species(processes, systems)
		adsorbates = [name for idx, name in enumerate(species) if idx not in gas_number]

		''' reconstruct algebraic expression to define Surface sites and RHS'''
		surface_expr = {}	 ## no worthy to pass it through the subroutine, too many patches.
		for s_sym in surfaces:
			s_name = str(s_sym)
			coverage_expr = 1.0    # Start with one full site
			for idx in adsorbates_by_surface[s_name]:	# Subtract adsorbates occupying this surface
				name = species[idx]     # idx is the position of the adsorbate in species
				coverage_expr -=  sp.symbols(name) * systems[str(name)]["nsites"]	#	symbolic
			surface_expr[s_sym] = coverage_expr
		# --- RHS generation
		rhs = []
		for name in adsorbates:     # ONLY adsorbates, no gaseous species
			expr = sp.Integer(0)
			for f, e in zip(equation_factors[name], equations[name]):
				expr += f * e
			rhs.append(expr.subs(surface_expr).subs(constants))
		rhs_func = sp.lambdify((*conditions, *adsorbates), rhs,  [{"exp": safe_exp}, "numpy"])
		gas_rate_expr = {}  # ONLY gases
		for idx in gas_number:
			name = species[idx]
			expr = sp.Integer(0)
			for f, e in zip(equation_factors[name], equations[name]):
				if f > 0:
					expr += f * e
			gas_rate_expr[name] = expr.subs(surface_expr).subs(constants)
		gas_rate_func = {gas: sp.lambdify((*conditions, *adsorbates), gas_rate_expr[gas], [{"exp": safe_exp}, "numpy"])
		                 for gas in gas_rate_expr}

		temp_i = 100.0
		temp_f = 1000.0
		dtemp = 1  # K resolution
		heating_rates = [1.0, 10.0]  # Already in K/min

		for gas_idx in gas_number:
			gas_name = species[gas_idx]
			ads_name, ads_idx = TPR.desorbing_gas(processes, adsorbates, gas_name)
			ics = np.zeros(len(adsorbates), dtype=float)
			ics[ads_idx] = 0.99 / float(systems[gas_name]['nmolsite'])	# initial for the adsorbate at 1 ML
			for beta_min in heating_rates:		# Heating rates 1 & 10 (K/min)
				beta = beta_min / 60       # heating rate in K/s

				temps, concentrations, rate = TPR.semi_implicit_temp(rhs_func, gas_rate_func, adsorbates, ics, temp_i,
				                                                     temp_f, dtemp, beta)
				data_conc = np.column_stack([temps, concentrations])
				conc_header = (["Temperature(K)"] + list(adsorbates))
				data_conc = np.vstack([conc_header, data_conc])
				printdata(['./KINETICS/TPR/', f"TPD_{gas_name}_concentrations_{int(beta_min)}K_min"], data_conc)

				tpd_data = np.column_stack([temps] + [rate[g] for g in gas_rate_expr.keys()])
				tpd_header = (["Temperature(K)"] + [g for g in gas_rate_expr.keys()])
				tpd_data = np.vstack([tpd_header, tpd_data])
				printdata(['./KINETICS/TPR/', f"TPD_{gas_name}_{int(beta_min)}K_min"], tpd_data)

				TPR.tpd_plot(gas_name, tpd_data, beta_min)

		elapsed = (time.time() - start) / 60
		print("\t\t\t\t", round(elapsed, 3), "minutes")

	@staticmethod
	def initial_species(processes, systems):
		"""
		Collect species preserving order and excluding transition states.
		"""
		seen = []
		for process in processes.values():
			seen.extend(process["reactants"])
			seen.extend(process["products"])
		ordered_species = list(dict.fromkeys(seen))		# preserve order

		gases = []
		ads = []
		surfaces = []
		for name in ordered_species:
			kind = systems[name]["kind"]
			if kind == "molecule":
				gases.append(name)
			elif kind == "adsorbate":
				ads.append(name)
			elif kind == "surface":
				surfaces.append(name)
		species = sorted(gases) + sorted(ads)
		adsorbates_by_surface = {s: [i for i, name in enumerate(species)
									 if systems[name]["kind"] == "adsorbate" and
									 systems[name]["sites"] == systems[s]["sites"]]	for s in sorted(surfaces)}
		# Correct indexing in species array
		gas_number = [i for i, name in enumerate(species) if systems[name]["kind"] == "molecule"]

		return species, sorted(surfaces), adsorbates_by_surface, gas_number

	@staticmethod
	def desorbing_gas(processes, adsorbates, name):
		''' Find the adsorbate species that desorbs into gas species `name`.
		Returns: 	ads_name : Adsorbate species name.
					ads_idx : Index of adsorbate in species list.'''
		for kind, proc in processes.items():
			if proc.get("kind") != "D":			# reaction type must be desorption
				continue
			reacts = proc.get("reactants", [])
			prods = proc.get("products", [])
			if name not in prods:			# check whether target gas is produced
					continue
			for r in reacts:		# find adsorbate in reactants
				ads_name = r
				ads_idx = adsorbates.index(r)
				return ads_name, ads_idx
		raise ValueError(f"No desorption reaction found for gas '{name}'")

	@staticmethod
	def semi_implicit_temp(rhs_func, gas_rate_func, adsorbates, ics, temp_i, temp_f, dtemp, beta, n_corrector=5,  nsub=1):
		''' For TPD/TPR, a semi-implicit (backward Euler / predictor-corrector) scheme is often much more stable
		than explicit Euler while remaining far simpler than a stiff global ODE solver.
		This approach avoids:
			- BDF instability
			- Jacobian singularities
			- LSODA internal failures
			- symbolic stiffness catastrophes
		gives:
			- deterministic behavior
			- fixed temperature spacing
			- easier debugging
			- physical clipping
			- robust TPD spectra '''
		temps = np.arange(temp_i, temp_f + dtemp, dtemp)
		rate = {g: np.zeros(len(temps), dtype=float) for g in gas_rate_func}
		nadsorbates = len(adsorbates)
		concentrations = np.zeros((len(temps), nadsorbates), dtype=float)
		concentrations[0] = np.asarray(ics, dtype=float)

		np.seterr(over='raise', invalid='raise')
		dtemp_sub = dtemp / nsub    # Internal temperature substep
		for i in range(len(temps) - 1):
			temp_old = temps[i]
			y_old = concentrations[i].copy()
			for m in range(nsub):	# Internal substepping
				temp_sub = temp_old + m * dtemp_sub
				temp_new = temp_sub + dtemp_sub
				y_old = np.clip(y_old, 0.0, 1.0)
				y_old = np.nan_to_num(y_old, nan=0.0, posinf=1.0, neginf=0.0)
				try:
					dydtemp_old = np.asarray(rhs_func(temp_sub, *y_old), dtype=float) / beta # ---- RHS in temperature space
				except FloatingPointError:
					print("Overflow at:")
					print("T =", temp_sub)
					for s, v in zip(adsorbates, y_old):
						print(f"{s:20s} {v: .5e}")
					raise
				y_trial = y_old + dtemp_sub * dydtemp_old
				y_trial = np.clip(y_trial, 0.0, 1.0)
				y_trial = np.nan_to_num(y_trial, nan=0.0, posinf=1.0, neginf=0.0)
				for k in range(n_corrector):	# Corrector iterations
					try:
						dydtemp_new = np.asarray(rhs_func(temp_new, *y_trial), dtype=float) / beta
					except FloatingPointError:
						print("Overflow at NEW:")
						print("T =", temp_sub)
						for s, v in zip(adsorbates, y_old):
							print(f"{s:20s} {v: .5e}")
						raise
					y_next = y_old + 0.5 * dtemp_sub * (dydtemp_old + dydtemp_new)	# Trapezoidal semi-implicit update
					y_next = np.clip(y_next, 0.0, 1.0)
					err = np.max(np.abs(y_next - y_trial))
					y_trial = 0.5 * y_trial + 0.5 * y_next		# damping
					if err < 1e-8:
						break
				y_old = y_trial.copy()
			concentrations[i + 1] = y_old

			for gas in gas_rate_func:# ---- TPD spectra
				rate[gas][i] = gas_rate_func[gas](temp_old, *y_old)

			if not np.all(np.isfinite(y_old)):
				raise RuntimeError(f"Non-finite solution at T={temps[i+1]}")

		return temps, concentrations, rate

	@staticmethod
	def tpd_plot(gas_name, data, heating_rate):
		icolour = ["k", "b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki", "teal"]
		iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1))]
		names = data[0, 1:] # first row are labels
		x = data[1:, 0].astype(float) - 273 # making it in centigrades
		ydata = data[1:, 1:].astype(float)
		y_max = np.max(ydata)

		fig, ax1 = plt.subplots(figsize=(10, 6), clear=True)
		for i in range(ydata.shape[1]):
			y = ydata[:, i] / y_max
			if np.max(y) > 0.01:	# arbitrary argument
				ax1.plot(x, y, color=icolour[i % len(icolour)], linestyle=iliner[i % len(iliner)], lw=1.5,
						 label=chem_label(names[i]))
				# --- find the TPD peaks
				peaks, properties = find_peaks(y, height=0.05,  # minimum peak height
												  prominence=0.02,  # minimum prominence
												  width=5,  # for broad desorption bands
												  distance=10)  # minimum separation (points)
				for peak_idx in peaks:
					ax1.annotate(f"{int(x[peak_idx])}", xy=(x[peak_idx], y[peak_idx]),
								 xytext=(0, 8), textcoords="offset points", color=icolour[i % len(icolour)], fontsize=14,
								 ha="center",  annotation_clip=False)

		ax1.set_xlabel("Temperature ($ \\degree C $)", fontsize=18)
		ax1.set_xlim([min(x), max(x)])
		ax1.tick_params(axis='x', rotation=0, labelsize=16)
		ax1.set_ylim([-0.01, 1.1])   # Normalised
		ax1.set_ylabel(f"$\\frac{{\\delta P_{{{chem_label(gas_name)}}}}}{{\\delta T}}$ (a.u.)", fontsize=18)
		ax1.tick_params(axis='y', rotation=0, labelsize=16)
		legend = ax1.legend(loc="best", fontsize=14)
		plt.title(f"TPD ({int(heating_rate)} $K \\cdot min^{{{-1}}}$)")
		fig.tight_layout()
		plt.ion()
		plt.savefig(f"KINETICS/TPR/TPD_{gas_name}_{int(heating_rate)}K_min.svg", dpi=300, orientation='landscape',
					transparent=True)




