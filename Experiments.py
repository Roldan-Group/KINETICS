"""
	This script builds on the perl version by A.Roldan.

"""

import pathlib
import time
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks, savgol_filter
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
	return np.exp(np.clip(x, -200, 200))

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
	"""
	Temperature Programmed Desorption / Reaction simulator.

	IMPORTANT DESIGN CHOICES
	------------------------
	1) Integrate ONLY adsorbates dynamically
	2) Surface sites reconstructed algebraically
	3) Gas production reconstructed from desorption rates
	4) Temperature ramp:
		:contentReference[oaicite:0]{index=0}
	5) Uses existing ODE solver infrastructure
	6) Preserves numerical stability and site conservation
	"""
	def __init__(self, systems, processes, equations, equation_factors):
		start = time.time()
		print("\t ... Generating TPR ...")

		# --- Species bookkeeping
		conditions = [t, temp]
		species, surfaces, adsorbates_by_surface, gas_number = self.initial_species(processes, systems)
		gas_species = [species[i] for i in gas_number]
		adsorbates = [s for s in species if systems[s]["kind"] == "adsorbate"]
		ads_symbols = [sp.symbols(a) for a in adsorbates]

		# --- Algebraic surface balances
		surface_expr = {}
		eps = 1e-30
		for s_sym in surfaces:
			s_name = str(s_sym)
			coverage_raw = 1.0
			for idx in adsorbates_by_surface[s_name]:
				name = species[idx]
				coverage_raw -= (sp.symbols(name) * systems[name]["nsites"])
			surface_expr[s_sym] = (coverage_raw + sp.sqrt(coverage_raw**2 + eps**2)) / 2 # like: max(A, 0) but smoothly

		# --- Build RHS ONLY for adsorbates
		rhs_ads = []
		for name in adsorbates:
			expr = sp.Integer(0)
			for f, e in zip(equation_factors[name],	equations[name]):
				#expr += f * e
				max_rate = 1e7   # Rate Regularisation
				expr += f * max_rate * sp.tanh(e / max_rate)    # no clipping, smoothening
			expr = expr.subs(surface_expr).subs(constants)
			rhs_ads.append(expr)
		rhs_func = sp.lambdify((*conditions, *ads_symbols), rhs_ads, [{"exp": safe_exp}, "numpy"])

		# --- Jacobian
		#jac_expr = sp.Matrix(rhs_ads).jacobian(ads_symbols)
		#jac_func = sp.lambdify((*conditions, *ads_symbols), jac_expr, [{"exp": safe_exp}, "numpy"])

		# --- Gas desorption rates
		gas_rate_expr = {}
		for gas in gas_species:
			expr = sp.Integer(0)
			for f, e in zip(equation_factors[gas], equations[gas]):
				if f > 0:
					expr += f * e
			expr = expr.subs(surface_expr).subs(constants)
			gas_rate_expr[gas] = expr
		gas_rate_func = {g: sp.lambdify((*conditions, *ads_symbols), gas_rate_expr[g], [{"exp": safe_exp},
																						"numpy"]) for g in gas_rate_expr}
		# --- TPR Conditions
		temp_i = 10.0
		temp_f = 1273.0
		dtemp = 1.0 	# K temperature accuracy

		heating_rates = [1.0] #[1.0, 10.0]  # K/min
		# --- Generate one TPD per desorbing gas
		for gas_name in gas_species:
			ads_name, ads_idx_global = self.desorbing_gas(processes, species, gas_name)
			ads_idx = adsorbates.index(ads_name)		# local adsorbate index

			# --- Initial conditions (ADSORBATES ONLY)
			ics = np.zeros(len(adsorbates), dtype=float)
			ics[ads_idx] = (0.95 / float(systems[gas_name]["nmolsite"]))    # <1 ML for stability

			for beta_min in heating_rates:

				print(f"\t TPD {gas_name} @ ",f"{beta_min} K/min")

				beta = beta_min / 60.0  # K/s

				# --- Time grid from temperature ramp
				temps = np.arange(temp_i, temp_f + dtemp, dtemp)
				t_final = (temp_f - temp_i) / beta
				t_eval = (temps - temp_i) / beta
				t_span = (0.0, t_final)


				def rhs_time(t_num, temp_num, *y):		# --- ODE wrappers
					return np.asarray(rhs_func(t_num, temp_num,	*y), dtype=float)
				#def jac_time(t_num, temp_num, *y):
				#	return np.asarray(jac_func(t_num, temp_num,	*y), dtype=float)
				def temperature_profile(t):		# --- Temperature argument
					return temp_i + beta * t
				def ode_rhs(t_num, y):		# --- solve_ivp wrappers
					temp_num = temperature_profile(t_num)
					dydt = rhs_time(t_num, temp_num, *y)
					if not np.all(np.isfinite(dydt)):
						print("NON-FINITE RHS")
						print("T =", temp_num)
						print("y =", y)
						print("dydt =", dydt)
						raise FloatingPointError
					return dydt
				#def ode_jac(t_num, y):
				#	temp_num = temperature_profile(t_num)
				#	return jac_time(t_num, temp_num, *y)

				try:
					sol = solve_ivp(ode_rhs, t_span, ics, t_eval=t_eval, method="BDF", rtol=1e-3, atol=1e-8) # dense_output=True, jac=ode_jac,
					if not sol.success:
						print(f"WARNING: solver stopped early for {gas_name}")
						print(sol.message)

					temps = np.asarray(temps) #temp_i + beta * sol.t)		# Recover temperatures grid
					concentrations = np.zeros((len(temps), len(adsorbates)))		# allocate full concentration matrix
					if len(sol.t) > 0:
						concentrations[:len(sol.t), :] = sol.y.T  # --- Build concentration matrix
					concentrations[concentrations < 0.0] = 0.0		# --- Tiny negative cleanup ONLY

					for row in concentrations:		# --- Surface conservation correction
						total = 0.0
						for i, ads in enumerate(adsorbates):
							total += (row[i]* systems[ads]["nsites"])
						if total > 1.0:
							row[:] /= total

					tpd = {}		# --- Reconstruct gas production
					for gas in gas_species:
						rates = np.zeros(len(temps))
						for itemp, (t_num, temp_num) in enumerate(zip(sol.t, temps)):
							rates[itemp] = gas_rate_func[gas](t_num, temp_num, *concentrations[itemp])
						rates = np.maximum(rates, 0.0)
						tpd[gas] = rates

					# --- Save concentrations
					data_conc = np.column_stack([temps, concentrations])
					conc_header = (["Temperature(K)"] + adsorbates)
					data_conc = np.vstack([conc_header,	data_conc])
					printdata(['./KINETICS/TPR/', f"TPD_{gas_name}_concentrations_", f"{int(beta_min)}K_min"],	data_conc)

					# --- Save TPD spectra
					tpd_matrix = np.column_stack([temps] + [tpd[g] for g in gas_species])
					tpd_header = (["Temperature(K)"] + gas_species)
					tpd_data = np.vstack([tpd_header, tpd_matrix])
					printdata(['./KINETICS/TPR/', f"TPD_{gas_name}_{int(beta_min)}K_min"], tpd_data)
					self.tpd_plot(gas_name,	tpd_data, beta_min)

				except Exception as err:
					print(f"TPR solver failed for {gas_name}: e{err}\n")

		elapsed = (time.time() - start) / 60
		print("\t\t\t\t", round(elapsed, 3), "minutes")

	@staticmethod
	def initial_species(processes, systems):
		seen = []
		for process in processes.values():
			seen.extend(process["reactants"])
			seen.extend(process["products"])
		ordered_species = list( dict.fromkeys(seen))
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
		adsorbates_by_surface = {s: [i for i, name in enumerate(species) if systems[name]["kind"] == "adsorbate"
									 and systems[name]["sites"]	== systems[s]["sites"]]	for s in sorted(surfaces)}
		gas_number = [i for i, name in enumerate(species) if systems[name]["kind"] == "molecule"]
		return (species, sorted(surfaces), adsorbates_by_surface, gas_number)

	@staticmethod
	def desorbing_gas(processes, species, gas_name):
		for _, proc in processes.items():
			if proc.get("kind") != "D":
				continue
			reacts = proc.get("reactants", [])
			prods = proc.get("products", [])
			if gas_name not in prods:
				continue
			for r in reacts:
				return r, species.index(r)
		raise ValueError(f"No desorption process ", f"for gas '{gas_name}'")

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
		ax1.set_ylabel("$\\frac{{\\delta P_{i}}}{{\\delta T}}$ (a.u.)", fontsize=18)
		ax1.tick_params(axis='y', rotation=0, labelsize=16)
		legend = ax1.legend(loc="best", fontsize=14)
		plt.title(f"TPD ({int(heating_rate)} $K \\cdot min^{{{-1}}}$)")
		fig.tight_layout()
		plt.ion()
		plt.savefig(f"KINETICS/TPR/TPD_{gas_name}_{int(heating_rate)}K_min.svg", dpi=300, orientation='landscape',
					transparent=True)


