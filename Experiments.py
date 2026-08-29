"""
Designed for Main.py, where stages are called as Routine(model).
Results are stored in memory under model.experiments.
Expected previous stages:
	RConstants(model)
	REquations(model)
Expected equation convention:
	rate equations are already signed contributions, i.e. dY/dt is obtained by
	summing equations[name]. equation_factors are retained only as metadata for
	production/consumption filtering.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace
from typing import Any

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from Symbols_def import t, temp, constants, pressure_std


@dataclass(slots=True)
class TimeSeriesResult:
	headers: list[str]
	values: np.ndarray
	species: list[str] = field(default_factory=list)
	surfaces: list[str] = field(default_factory=list)
	metadata: dict[str, Any] = field(default_factory=dict)

@dataclass(slots=True)
class TPRResult:
	gas: str
	heating_rate: float
	temperatures: np.ndarray
	adsorbates: list[str]
	gas_species: list[str]
	concentrations: np.ndarray
	spectra: dict[str, np.ndarray]
	metadata: dict[str, Any] = field(default_factory=dict)

@dataclass(slots=True)
class ExperimentResults:
	isothermal: TimeSeriesResult | None = None
	tpr: list[TPRResult] = field(default_factory=list)

def safe_exp(x):
	return np.exp(np.clip(x, -200, 200))

def kind(species) -> str:
	return str(species.kind).lower()

def is_molecule(species) -> bool:
	return kind(species) == "molecule"

def is_surface(species) -> bool:
	return kind(species) == "surface"

def is_adsorbate(species) -> bool:
	return kind(species) in {"adsorbate", "transitionstate", "ts"}

def species_symbol(name: str) -> sp.Symbol:
	return sp.Symbol(name)

def ensure_experiment_results(model) -> ExperimentResults:
	if not hasattr(model, "experiments") or model.experiments is None:
		model.experiments = ExperimentResults()
	return model.experiments

def scan_values(scan, *, default: list[float] | None = None) -> list[float]:
	if scan is None:
		if default is None:
			raise ValueError("Required scan/condition is missing.")
		return default
	if hasattr(scan, "values"):
		return list(scan.values())
	if isinstance(scan, (int, float)):
		return [float(scan)]
	if isinstance(scan, (list, tuple)):
		if len(scan) == 1:
			return [float(scan[0])]
		if len(scan) == 3:
			return list(np.arange(float(scan[0]), float(scan[1]) + 1e-12, float(scan[2])))
	raise TypeError(f"Cannot convert condition to scan values: {scan!r}")

def get_equation_tuple(model):
	"""
	Return (equations, equation_factors, tpd, tpd_factors).
	Preferred source is model.equations if a future refactor stores equations
	explicitly on the model. For compatibility with the previous refactored
	Kinetics.py, this also supports Kinetics.get_rate_equations(model).
	"""
	if hasattr(model, "equations") and model.equations is not None:
		equations = model.equations
		if isinstance(equations, tuple) and len(equations) == 4:
			return equations
		if all(hasattr(equations, attr) for attr in ("isothermal", "equation_factors", "tpd", "tpd_factors")):
			return (equations.isothermal, equations.equation_factors, equations.tpd, equations.tpd_factors,)
	try:
		from Kinetics import get_rate_equations  # type: ignore
	except Exception:
		get_rate_equations = None
	if get_rate_equations is not None:
		equations = get_rate_equations(model)
		if equations is not None:
			return equations
	raise ValueError("Rate equations are missing. Run Kinetics(model) before Experiments.")

def participating_species(model, *, include_ts: bool = False) -> list[str]:
	seen: list[str] = []
	for process in model.processes.values():
		for item in process.reactants:
			if item.species not in seen:
				seen.append(item.species)
		if include_ts:
			for item in process.ts:
				if item.species not in seen:
					seen.append(item.species)
		for item in process.products:
			if item.species not in seen:
				seen.append(item.species)
	return seen

def initial_species(model):
	"""
	Return initial conditions and species bookkeeping.
	Dynamic species are molecules plus adsorbates. Surfaces are reconstructed
	algebraically from site balances rather than integrated dynamically.
	"""
	names = participating_species(model, include_ts=False)
	gases = sorted(name for name in names if is_molecule(model.species[name]))
	adsorbates = sorted(name for name in names if is_adsorbate(model.species[name]))
	surfaces = sorted(name for name in names if is_surface(model.species[name]))
	dynamic = gases + adsorbates
	initial_conditions = [float(getattr(model.species[name], "pressure0", 0.0)) for name in gases]
	initial_conditions += [float(getattr(model.species[name], "coverage0", 0.0)) for name in adsorbates]
	gas_indices = [dynamic.index(name) for name in gases]
	adsorbates_by_surface: dict[str, list[int]] = {}

	for surface_name in surfaces:
		surface = model.species[surface_name]
		surface_site = getattr(surface, "sites", None)
		adsorbates_by_surface[surface_name] = [dynamic.index(name) for name in adsorbates
											   if getattr(model.species[name], "sites", None) == surface_site]
	return initial_conditions, dynamic, surfaces, gas_indices, adsorbates_by_surface

def nsites(model, name: str) -> float:
	value = getattr(model.species[name], "nsites", None)
	return float(value) if value is not None else 1.0

def build_surface_expressions(model, dynamic_species: list[str], surfaces: list[str],
							  adsorbates_by_surface: dict[str, list[int]], *, smooth: bool = True,
							  smoothing_epsilon: float = 1.0e-12,):
	"""Build algebraic surface-coverage expressions."""
	expressions: dict[sp.Symbol, sp.Expr] = {}
	eps = sp.Float(smoothing_epsilon) # smoothening value
	for surface_name in surfaces:
		coverage = sp.Float(1.0)
		for idx in adsorbates_by_surface[surface_name]:
			ads_name = dynamic_species[idx]
			coverage -= species_symbol(ads_name) * nsites(model, ads_name)
		if smooth:
			coverage = (coverage + sp.sqrt(coverage**2 + eps**2)) / 2
		expressions[species_symbol(surface_name)] = coverage
	return expressions

def clip_dynamic_species(model, dynamic_species: list[str], surfaces: list[str],
						 adsorbates_by_surface: dict[str, list[int]], gas_indices: list[int], y: np.ndarray) -> np.ndarray:
	"""Clip gases to non-negative values and enforce surface site conservation."""
	y = np.asarray(y, dtype=float).copy()
	for surface_name in surfaces:
		indices = adsorbates_by_surface[surface_name]
		if not indices:
			continue
		total = np.zeros(y.shape[1]) if y.ndim == 2 else 0.0
		for idx in indices:
			total += y[idx] * nsites(model, dynamic_species[idx])
		if y.ndim == 2:
			mask = total > 1.0
			if np.any(mask):
				for idx in indices:
					y[idx, mask] /= total[mask]
		elif total > 1.0:
			for idx in indices:
				y[idx] /= total
	if y.ndim == 2:
		for idx in gas_indices:
			y[idx, :] = np.clip(y[idx, :], 0.0, None)
	else:
		for idx in gas_indices:
			y[idx] = np.clip(y[idx], 0.0, None)
	return y

def surface_coverages(model, y: np.ndarray, dynamic_species: list[str], surfaces: list[str],
					  adsorbates_by_surface: dict[str, list[int]]) -> list[np.ndarray]:
	"""Compute surface coverages from adsorbate coverages."""
	values: list[np.ndarray] = []
	for surface_name in surfaces:
		coverage = np.ones(y.shape[1], dtype=float)
		for idx in adsorbates_by_surface[surface_name]:
			ads_name = dynamic_species[idx]
			coverage -= y[idx, :] * nsites(model, ads_name)
		values.append(np.clip(coverage, 0.0, 1.0))
	return values

def scale_dynamic_state(y, gas_indices, pressure_scale=float(pressure_std)):
	"""Convert physical state [Pa, coverage] to solver-scaled state."""
	y = np.asarray(y, dtype=float).copy()
	for idx in gas_indices:
		y[idx] /= pressure_scale
	return y

def unscale_dynamic_state(y, gas_indices, pressure_scale=float(pressure_std)):
	"""Convert solver-scaled state back to physical state [Pa, coverage]."""
	y = np.asarray(y, dtype=float).copy()
	if y.ndim == 1:
		for idx in gas_indices:
			y[idx] *= pressure_scale
	elif y.ndim == 2:
		for idx in gas_indices:
			y[idx, :] *= pressure_scale
	return y

def solve_ode_system(model, t_span: tuple[float, float], t_eval: np.ndarray,	dynamic_species: list[str],
					 surfaces: list[str], adsorbates_by_surface: dict[str, list[int]], gas_indices: list[int],
					 initial_conditions: list[float],	rhs_func, jac_func,	arguments: tuple[float, ...],
					 preferred_solver: str | None = None,):
	"""Solve the microkinetic ODE system for one set of external arguments."""
	y0_physical = np.asarray(initial_conditions, dtype=float)
	y0 = scale_dynamic_state(y0_physical, gas_indices,)     # escales the pressure so all species are in [0-1] regime
	gas_index_set = set(gas_indices)
	if y0.ndim != 1 or y0.size != len(dynamic_species):
		raise ValueError(f"Invalid initial state shape {y0.shape}; expected ({len(dynamic_species)},).")
	if not np.all(np.isfinite(y0)):
		raise ValueError(f"Initial state contains non-finite values: {y0}")
	if np.any(y0 < 0.0):
		invalid = [(name, value) for name, value in zip(dynamic_species, y0) if value < 0.0]
		raise ValueError(f"Negative initial conditions: {invalid}")
	for surface_name in surfaces:    # Check initial site balances.
		occupied = sum(nsites(model, dynamic_species[idx]) * y0[idx] for idx in adsorbates_by_surface[surface_name])
		if occupied > 1.0 + 1.0e-12:
			raise ValueError(f"Initial coverage exceeds the capacity of {surface_name}: occupied coverage = {occupied:.8e}")

	class SolverWorkLimit(RuntimeError):
		pass
	rhs_call_count = 0
	solver_limits = {"LSODA": 500000, "BDF": 15000, "Radau": 15000,}

	def ode_system(t_num, y_scaled, *args):
		nonlocal rhs_call_count
		rhs_call_count += 1
		y_scaled = np.asarray(y_scaled, dtype=float)
		# Convert solver variables back to physical variables.
		y_physical = unscale_dynamic_state(y_scaled, gas_indices,)
		# Evaluate the original physical microkinetic equations.
		dydt_physical = np.asarray(rhs_func(t_num, *args, *y_physical), dtype=float,).reshape(-1)

		if rhs_call_count % 5000 == 0:
			print(f"    T={args[0]:g} K  RHS={rhs_call_count}  t={t_num:.8e} s"
				  f"  min(y)={np.min(y_physical):.4e}  max(y)={np.max(y_physical):.4e}", flush=True,)
			for surface_name in surfaces:
				occupied = sum(nsites(model, dynamic_species[idx]) * y_physical[idx]
				               for idx in adsorbates_by_surface[surface_name])
				print(f"        {surface_name}: occupied={occupied:.8e}, free={1.0 - occupied:.8e}", flush=True,)
		if rhs_call_count > rhs_call_limit:
			raise SolverWorkLimit(f"Exceeded {rhs_call_limit} RHS evaluations at T={args[0]:g} K")
		if not np.all(np.isfinite(y_scaled)):
			raise FloatingPointError(f"Non-finite state at t={t_num:.16e}: {y_scaled}")

		if dydt_physical.shape != y_scaled.shape:
			raise ValueError(f"RHS returned shape {dydt_physical.shape}; expected {y_scaled.shape}.")
		if not np.all(np.isfinite(dydt_physical)):
			raise FloatingPointError(f"Non-finite RHS at t={t_num:.16e}; y={y_physical}; dydt={dydt_physical}")
		# Transform physical derivatives to solver coordinates.
		dydt_scaled = dydt_physical.copy()
		for idx in gas_indices:
			dydt_scaled[idx] /= float(pressure_std)
		return dydt_scaled

	def ode_jacobian(t_num, y_scaled, *args):
		y_scaled = np.asarray(y_scaled, dtype=float)
		y_physical = unscale_dynamic_state(y_scaled, gas_indices,)
		jac_physical = np.asarray(jac_func(t_num, *args, *y_physical), dtype=float,)
		expected_shape = (len(dynamic_species), len(dynamic_species),)
		if jac_physical.shape != expected_shape:
			raise ValueError(f"Jacobian returned shape {jac_physical.shape}; expected {expected_shape}.")
		# S diagonal:     gas variables -> pressure_std     coverages -> 1
		scales = np.ones(len(dynamic_species), dtype=float)
		scales[gas_indices] = float(pressure_std)
		jac_scaled = (jac_physical * scales[np.newaxis, :] / scales[:, np.newaxis])	    # J_scaled = S^-1 J_physical S
		if not np.all(np.isfinite(jac_scaled)):
			raise FloatingPointError(f"Non-finite scaled Jacobian at t={t_num:.16e}")
		return jac_scaled

	integration_options = {"fun": ode_system, "t_span": t_span, "y0": y0, "t_eval": t_eval, "args": arguments,
							   "rtol": 1.0e-7, "atol": 1.0e-10, "jac": ode_jacobian,}

	available_methods = ("BDF", "Radau", "LSODA")
	if preferred_solver in available_methods:
		methods_to_try = (preferred_solver, *(method for method in available_methods if method != preferred_solver),)
	else:
		methods_to_try = available_methods

	sol = None
	selected_method = None
	for method in methods_to_try:
		rhs_call_count = 0
		rhs_call_limit = solver_limits[method]
		#print(f"  Starting {method} at T={arguments[0]:g} K...", flush=True,)
		try:
			trial = solve_ivp(method=method,  **integration_options,)
		except SolverWorkLimit as err:
			print(f" WARNING: {method} abandoned at T={arguments[0]:g} K.\tReason: {err}", flush=True,)
			continue
		if trial.success:
			sol = trial
			selected_method = method
			print(f"  ODE at T={arguments[0]:g} K successfully solved with {method}. Function_calls={trial.nfev}",
				  flush=True,)
			break
		print(f" WARNING: {method} failed at T={arguments[0]:g} K.\n  Message           : {trial.message}\n"
			  f"  Last reached time : {trial.t[-1] if trial.t.size else t_span[0]:.16e}\n"
			  f"  Function calls    : {trial.nfev}\n  Jacobian calls    : {trial.njev}\n"
			  f"  Linear algebra operations : {trial.nlu}\n"
			  "  Moving to the next solver.", flush=True,)
		'''if trial.y.size:     # to print the species concentrations
			last_state = trial.y[:, -1]
			print("   Last solver state:")
			for name, value in zip(dynamic_species, last_state):
				print(f"\t{name:10s} = {value: .8e}")'''
	if sol is None:
		return None
	if len(sol.t) != len(t_eval):
		print(f" WARNING: Incomplete solution at T={arguments[0]:g} K: "
			  f"  returned {len(sol.t)} of {len(t_eval)} requested points.", flush=True,)
		return None

	y_physical = unscale_dynamic_state(sol.y, gas_indices,)
	y_dynamic = clip_dynamic_species(model, dynamic_species, surfaces, adsorbates_by_surface, gas_indices, y_physical)
	y_surfaces = surface_coverages(model, y_dynamic, dynamic_species, surfaces, adsorbates_by_surface)
	augmented = SimpleNamespace()
	augmented.t = sol.t
	augmented.y = np.vstack([y_dynamic, *y_surfaces])
	augmented.species = dynamic_species + surfaces
	rows = np.column_stack([np.tile(arguments, (len(sol.t), 1)), augmented.t, augmented.y.T,])
	return sol, rows, augmented, selected_method

def process_involves_gas(model, process,) -> bool:
	"""   Return True when any reactant or product of a process is a gas-phase molecule.    """
	for item in (*process.reactants, *process.products,):
		species = model.species[item.species]
		if is_molecule(species):
			return True
	return False


class Isothermal:
	"""Run constant-temperature/time-grid integrations and store results in memory."""
	def __init__(self, model):
		self.model = model
		self.result = self.build()
		ensure_experiment_results(model).isothermal = self.result

	def build(self) -> TimeSeriesResult:
		equations, _, _, _ = get_equation_tuple(self.model)
		initial_conditions, dynamic_species, surfaces, gas_indices, adsorbates_by_surface = initial_species(self.model)
		surface_expr = build_surface_expressions(self.model, dynamic_species, surfaces, adsorbates_by_surface,
												 smooth=True,)
		rhs = []
		for name in dynamic_species:
			expr = sp.Integer(0)
			for contribution in equations.get(name, []):
				expr += contribution
			expr = expr.subs(surface_expr).subs(constants)
			rhs.append(expr)

		dynamic_symbols = [species_symbol(name) for name in dynamic_species]
		conditions = [t, temp]
		rhs_func = sp.lambdify((*conditions, *dynamic_symbols), rhs, [{"exp": safe_exp}, "numpy"], cse=True)
		jacobian = sp.Matrix(rhs).jacobian(dynamic_symbols)
		jac_func = sp.lambdify((*conditions, *dynamic_symbols), jacobian, [{"exp": safe_exp}, "numpy"], cse=True)
		time_values = scan_values(self.model.conditions.time)
		if len(time_values) < 2:
			raise ValueError("Isothermal simulation requires at least two time points.")
		t_eval = np.asarray(time_values, dtype=float)
		t_span = (float(t_eval[0]), float(t_eval[-1]))
		temperature_values = scan_values(self.model.conditions.temperature, default=[300.0])

		rows = []
		expected_columns = (1                       # Temperature
							+ 1                     # time
							+ len(dynamic_species)
							+ len(surfaces))
		preferred_solver: str | None = None
		for temperature in temperature_values:
			results = solve_ode_system(self.model, t_span, t_eval, dynamic_species, surfaces, adsorbates_by_surface,
										  gas_indices, initial_conditions, rhs_func, jac_func, (float(temperature),),
									   preferred_solver=preferred_solver,)
			# Failed integration.
			if results is None:
				print(f"WARNING: all ODE solvers failed. No rows retained at T={float(temperature):g} K.", flush=True,)
				continue
			sol, data, augmented, selected_method = results
			preferred_solver = selected_method
			# Additional protection if the function returns a tuple of None values.
			if data is None:
				# print(f"WARNING: no rows retained at T={float(temperature):g} K.", flush=True,)
				continue
			data = np.asarray(data, dtype=float)
			if data.ndim != 2:
				# print(f"WARNING: skipping malformed result at T={float(temperature):g} K: "
				#      f"expected a 2D array, got shape {data.shape}.", flush=True,)
				continue
			if data.shape[0] == 0:
				# print(f"WARNING: skipping empty result at T={float(temperature):g} K.", flush=True,)
				continue
			if data.shape[1] != expected_columns:
				# print(f"WARNING: skipping malformed result at T={float(temperature):g} K: "
				#      f"expected {expected_columns} columns, got {data.shape[1]}.", flush=True,)
				continue
			rows.append(data)
			if rows:
				values = np.vstack(rows)
			else:
				values = np.empty((0, expected_columns), dtype=float,)
		headers = ["Temperature", "time", *dynamic_species, *surfaces]
		return TimeSeriesResult(headers=headers, values=values, species=dynamic_species, surfaces=surfaces,
								metadata={"adsorbates_by_surface": adsorbates_by_surface, "gas_indices": gas_indices,},)


class TPR:
	"""Run temperature-programmed simulations and store spectra in memory."""
	def __init__(self, model, *, heating_rates: list[float] | None = None, temperature_scan=None):
		self.model = model
		if heating_rates is not None:
			self.heating_rates = heating_rates
		elif (hasattr(model.conditions, "tpr") and model.conditions.tpr is not None):
			self.heating_rates = [model.conditions.tpr.heating_rate]
		else:
			self.heating_rates = [1.0]
		self.temperature_scan = temperature_scan
		self.results = self.build()
		ensure_experiment_results(model).tpr = self.results

	def build_surface_only_tpd_equations(self, *, tpd_equations: dict[str, list[sp.Expr],],
										 tpd_factors: dict[str, list[float],], adsorbates: list[str],) -> dict[str, list[sp.Expr],]:
		"""    Construct TPD equations containing only elementary processes that involve no gas-phase molecules.
			The TPD equation arrays contain only processes for which:            process.kind.upper() != "A"        """
		tpd_process_items = [(process_id, process,) for (process_id, process,) in self.model.processes.items()
							 if process.kind.upper() != "A"]
		number_of_tpd_processes = len(tpd_process_items)
		surface_only: dict[str, list[sp.Expr],] = {name: [] for name in adsorbates}

		# Verify that the equation arrays are aligned with the non-adsorption TPD process list.
		for adsorbate in adsorbates:
			contributions = tpd_equations.get(adsorbate, [],)
			factors = tpd_factors.get(adsorbate, [],)
			if len(contributions) != number_of_tpd_processes:
				raise ValueError(f"TPD contribution/process alignment failure for {adsorbate!r}: "
								 f"{len(contributions)} contributions but {number_of_tpd_processes} non-adsorption")
			if len(factors) != number_of_tpd_processes:
				raise ValueError(f"TPD factor/process alignment failure for {adsorbate!r}: {len(factors)} factors but "
								 f"{number_of_tpd_processes} non-adsorption processes.")

		# Retain only processes whose reactants and products contain no molecule species.
		#print("\nTPR: classifying preparation processes", flush=True,)
		active_process_count = 0
		disabled_process_count = 0
		for (tpd_index, (process_id, process,),)in enumerate(tpd_process_items):
			involves_gas = process_involves_gas(self.model, process,)
			if involves_gas:                        # BOTH as reactant and products in "A" and "D"
				disabled_process_count += 1
				#print(f"  {process_id}: DISABLED: gas exchange", flush=True,)
				continue
			active_process_count += 1
			#print(f"  {process_id}: ACTIVE: surface only", flush=True,)
			for adsorbate in adsorbates:
				contribution = tpd_equations[adsorbate][tpd_index]
				factor = float(tpd_factors[adsorbate][tpd_index])
				# The equation generator stores a zero expression and zero factor for species unaffected by this process.
				if factor == 0.0:
					continue
				surface_only[adsorbate].append(contribution)
		print(f"  TPR preparation mechanism: {active_process_count} active surface processes,"
			f" {disabled_process_count} gas processes disabled.", flush=True, )
		return surface_only

	def build(self) -> list[TPRResult]:
		_, _, tpd_equations, tpd_factors = get_equation_tuple(self.model)
		initial_conditions, dynamic_species, surfaces, gas_indices, adsorbates_by_surface = initial_species(self.model)
		gas_species = [dynamic_species[idx] for idx in gas_indices]
		adsorbates = [name for name in dynamic_species if is_adsorbate(self.model.species[name])]
		tpd_surface_equations = (self.build_surface_only_tpd_equations(tpd_equations=tpd_equations,
																	   tpd_factors=tpd_factors, adsorbates=adsorbates,))
		ads_symbols = [species_symbol(name) for name in adsorbates]
		surface_expr = build_surface_expressions(self.model, dynamic_species, surfaces, adsorbates_by_surface,
												 smooth=True, smoothing_epsilon=1.0e-16,)
		function_arguments = (t, temp, *ads_symbols,)

		# Helper for constructing one symbolic RHS vector.
		def construct_adsorbate_rhs(equation_contributions: dict[str, list[sp.Expr],], *, label: str,)-> list[sp.Expr]:
			""" Construct the ordered symbolic adsorbate RHS.
			 The order of the returned expressions matches ``adsorbates`` and therefore matches ``ads_symbols``.  """
			expressions: list[sp.Expr] = []
			#print(f"TPR: constructing {label} RHS expressions", flush=True,)
			for name in adsorbates:
				contributions = equation_contributions.get(name, [],)
				if contributions:
					expression = sp.Add(*contributions, evaluate=False,)
				else:
					expression = sp.Integer(0)
				expression = expression.subs(surface_expr)
				expression = expression.subs(constants)
				expressions.append(expression)
			if len(expressions) != len(ads_symbols):
				raise ValueError(f"The {label} RHS contains {len(expressions)} expressions, but "
								 f"{len(ads_symbols)} adsorbate symbols were defined.")
			return expressions

		rhs_ads = construct_adsorbate_rhs(tpd_equations, label="full",)
		rhs_matrix = sp.Matrix(rhs_ads)
		surface_rhs_ads = construct_adsorbate_rhs(tpd_surface_equations, label="surface-only",)
		surface_rhs_matrix = sp.Matrix(surface_rhs_ads)
		jacobian_matrix = rhs_matrix.jacobian(ads_symbols)
		modules = [{"exp": safe_exp}, "numpy",]
		rhs_func = sp.lambdify(function_arguments, rhs_matrix, modules, cse=True,)
		surface_rhs_func = sp.lambdify(function_arguments, surface_rhs_matrix, modules, cse=True,)
		jac_func = sp.lambdify(function_arguments, jacobian_matrix, modules, cse=True,)

		gas_rate_func: dict[str, Any] = {}
		for gas in gas_species:
			expr = sp.Integer(0)
			factors = tpd_factors.get(gas, [])
			contributions = tpd_equations.get(gas, [])
			for factor, contribution in zip(factors, contributions,):
				if factor > 0.0:
					expr += contribution
			expr = expr.subs(surface_expr).subs(constants)
			gas_rate_func[gas] = sp.lambdify(function_arguments, expr, modules, cse=True,)
		temperatures = self.temperature_values()
		results: list[TPRResult] = []
		for gas_name in gas_species:
			try:
				ads_name = self.desorbing_adsorbate(gas_name)
			except ValueError:
				continue
			if ads_name not in adsorbates:
				continue
			ads_idx = adsorbates.index(ads_name)
			for beta_min in self.heating_rates:
				result = self.run_one_tpr(gas_name=gas_name, adsorbates=adsorbates, ads_idx=ads_idx,
										  gas_species=gas_species, rhs_func=rhs_func, surface_rhs_func=surface_rhs_func,
										  jac_func=jac_func, gas_rate_func=gas_rate_func, temperatures=temperatures,
										  beta_min=float(beta_min),)
				if result is not None:
					results.append(result)
		return results

	def temperature_values(self) -> np.ndarray:
		return np.arange(150.0, 1273.0 + 5.0, 5.0)

	def pre_equilibrate_tpr_state(self, *, surface_rhs_func, y0: np.ndarray, temp_i: float, adsorbates: list[str],
								  site_occupancies: np.ndarray, maximum_relaxation_time: float = 100.0,
								  absolute_rate_tolerance: float = 1.0e-8, relative_rate_tolerance: float = 1.0e-6,
								  coverage_scale: float = 1.0e-10,  rtol: float = 1.0e-7, atol: float = 1.0e-10,
								  ) -> tuple[np.ndarray, dict[str, Any]]:
		"""  Pre-equilibrate a nominal TPD initial state at constant temperature.
				Only the reactions included in ``surface_rhs_func`` are evaluated.
			surface_rhs_func:            Lambdified surface-only RHS with the signature:
							surface_rhs_func(time, temperature, *coverages)
							It must contain no adsorption, desorption, or other elementary reaction involving a gas-phase.
			y0            Nominal initial adsorbate-coverage vector.
			temp_i            Initial/preparation temperature in K.
			adsorbates            Ordered adsorbate names corresponding to the state vector.
			site_occupancies            Number of sites occupied by one molecule of each adsorbate.
			maximum_relaxation_time             Maximum cumulative isothermal relaxation time in seconds.
			absolute_rate_tolerance            Absolute convergence threshold for max(abs(dtheta/dt)).
			relative_rate_tolerance            Scale-aware convergence threshold in s^-1.
			coverage_scale              Coverage floor used in the relative-rate calculation.
			rtol, atol            Solver tolerances used during pre-equilibration.   """
		y_initial = np.asarray(y0, dtype=float,).reshape(-1)
		y_current = y_initial.copy()
		site_occupancies = np.asarray(site_occupancies, dtype=float,).reshape(-1)
		n_adsorbates = len(adsorbates)
		if y_current.size != n_adsorbates:
			raise ValueError("The pre-equilibration state size does not match the number of adsorbates: "
							 f"state={y_current.size}, adsorbates={n_adsorbates}.")
		if site_occupancies.size != n_adsorbates:
			raise ValueError("The site-occupancy array size does not match the number of adsorbates: "
							 f"occupancies={site_occupancies.size}, adsorbates={n_adsorbates}.")
		if not np.all(np.isfinite(y_current)):
			raise ValueError("The nominal pre-equilibration state contains NaN or infinity.")
		if not np.all(np.isfinite(site_occupancies)):
			raise ValueError("The site-occupancy array contains NaN or infinity.")
		if np.any(site_occupancies <= 0.0):
			invalid_species = [name for name, occupancy in zip(adsorbates, site_occupancies,) if occupancy <= 0.0]
			raise ValueError(f"All site occupancies must be positive. Invalid species: {invalid_species}.")
		if np.any(y_current < 0.0):
			negative_species = [name for name, coverage in zip(adsorbates, y_current,) if coverage < 0.0]
			raise ValueError(f"The nominal pre-equilibration state contains negative coverages for: {negative_species}.")
		maximum_relaxation_time = float(maximum_relaxation_time)
		if (not np.isfinite(maximum_relaxation_time) or maximum_relaxation_time <= 0.0):
			raise ValueError("maximum_relaxation_time must be finite and positive.")

		occupied_fraction_initial = float(np.dot(y_current, site_occupancies,))
		if occupied_fraction_initial > 1.0 + 1.0e-8:
			raise ValueError(f"The nominal pre-equilibration state exceeds the "
							 f"available site capacity: {occupied_fraction_initial:.8e}.")

		rhs_call_count = 0

		def relaxation_rhs(relaxation_time: float, state: np.ndarray,) -> np.ndarray:
			nonlocal rhs_call_count
			rhs_call_count += 1
			state = np.asarray(state, dtype=float,).reshape(-1)
			if state.size != n_adsorbates:
				raise ValueError(f"The surface-only RHS received a state of the wrong "
								 f"size: expected {n_adsorbates}, got {state.size}.")
			if not np.all(np.isfinite(state)):
				raise FloatingPointError("The pre-equilibration solver generated a "
										 f"non-finite state at t={relaxation_time:.8e} s.")
			values = np.asarray(surface_rhs_func(float(relaxation_time), float(temp_i), *state,), dtype=float,).reshape(-1)
			if values.size != n_adsorbates:
				raise ValueError("The surface-only RHS returned the wrong number of "
								 f"components: expected {n_adsorbates}, got {values.size}.")
			if not np.all(np.isfinite(values)):
				invalid_indices = np.flatnonzero(~np.isfinite(values))
				invalid_species = [adsorbates[index] for index in invalid_indices]
				raise FloatingPointError("The surface-only RHS generated NaN or infinity at "
										 f"t={relaxation_time:.8e} s and T={temp_i:.8g} K. "
										 f"Invalid equations: {invalid_species}.")
			return values

		def evaluate_relaxation(state: np.ndarray, elapsed_time: float,) -> dict[str, Any]:
			derivative = relaxation_rhs(elapsed_time, state,)
			absolute_rates = np.abs(derivative)
			relative_rates = (absolute_rates / np.maximum(np.abs(state), coverage_scale,))
			maximum_absolute_rate = float(np.max(absolute_rates))
			maximum_relative_rate = float(np.max(relative_rates))
			absolute_converged = (maximum_absolute_rate <= absolute_rate_tolerance)
			relative_converged = (maximum_relative_rate <= relative_rate_tolerance)

			# - the absolute criterion handles species near zero;
			# - the relative criterion handles finite coverages.
			converged = bool(absolute_converged or relative_converged)
			return {"derivative": derivative, "maximum_absolute_rate": maximum_absolute_rate,
					"maximum_relative_rate": maximum_relative_rate, "absolute_converged": absolute_converged,
					"relative_converged": relative_converged, "converged": converged,}

		initial_diagnostics = evaluate_relaxation(y_current, 0.0,)

		'''print("\n Surface-only TPD pre-equilibration:"f"\n  temperature              = {temp_i:g} K"
			  f"\n  initial occupied fraction = {occupied_fraction_initial:.16e}"
			  f"\n  initial max |dtheta/dt|   = {initial_diagnostics['maximum_absolute_rate']:.6e}"
			  f"\n  initial max relative rate = {initial_diagnostics['maximum_relative_rate']:.6e} s^-1", flush=True,)
		print("\n Initial surface-only derivatives:", flush=True,)
		for (adsorbate_name, coverage, derivative,) in zip(adsorbates, y_current, initial_diagnostics["derivative"],):
			if (abs(coverage) > 1.0e-14 or abs(derivative) > 1.0e-14):
				print(f"  {adsorbate_name:18s} theta={coverage: .6e} dtheta/dt={derivative: .6e}", flush=True,)'''

		# Use successively larger cumulative relaxation times.
		# This avoids asking the solver to span immediately from a sub-picosecond transient to tens of seconds.

		candidate_end_times = np.asarray([1.0e-15, 1.0e-14, 1.0e-13, 1.0e-12, 1.0e-11, 1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7,
										  1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2,1.0e-1,1.0,10.0,100.0,],dtype=float,)
		relaxation_end_times = [float(value) for value in candidate_end_times if value <= maximum_relaxation_time]
		if (not relaxation_end_times or relaxation_end_times[-1] < maximum_relaxation_time):
			relaxation_end_times.append(maximum_relaxation_time)
		elapsed_time = 0.0
		converged = bool(initial_diagnostics["converged"])
		attempts: list[dict[str, Any]] = []
		final_diagnostics = (initial_diagnostics)
		#if converged:
		#    print("  The nominal state already satisfies the surface-relaxation criterion.", flush=True,)

		# Perform staged integrations.
		for target_time in relaxation_end_times:
			if converged:
				break
			interval_duration = (target_time - elapsed_time)
			if interval_duration <= 0.0:
				continue
			trial_solution = None
			selected_method = None
			for method in ("BDF", "LSODA", "Radau",):
				calls_before_attempt = (rhs_call_count)
				try:
					trial = solve_ivp(fun=relaxation_rhs, t_span=(elapsed_time, target_time,), y0=y_current,
									  method=method, rtol=rtol, atol=atol, max_step=max(interval_duration / 10.0,
																						np.finfo(float).tiny,),)
				except (FloatingPointError, ValueError, OverflowError,) as error:
					attempts.append({"method": method, "start_time_s": elapsed_time, "target_time_s": target_time,
									 "success": False, "message": (f"{type(error).__name__}: {error}"),
									 "nfev": (rhs_call_count - calls_before_attempt), "njev": 0, "nlu": 0,})
					print(f"  {method} pre-equilibration raised {type(error).__name__}: {error}", flush=True,)
					continue

				attempts.append({"method": method, "start_time_s": elapsed_time, "target_time_s": target_time,
								 "success": bool(trial.success), "message": str(trial.message),
								 "nfev": int(getattr(trial, "nfev", 0,)), "njev": int(getattr(trial, "njev", 0,)),
								 "nlu": int(getattr(trial, "nlu", 0,)),})
				if (trial.success and len(trial.t) > 0 and np.isclose(float(trial.t[-1]), target_time, rtol=0.0,
														atol=max(1.0e-18, 1.0e-10 * max(1.0,abs(target_time),),),)):
					trial_solution = trial
					selected_method = method
					break
				print(f"  {method} did not complete the pre-equilibration interval "
					  f"{elapsed_time:.3e}–{target_time:.3e} s: {trial.message}", flush=True,)
			if trial_solution is None:
				raise RuntimeError(f"Surface-only TPD pre-equilibration failed between {elapsed_time:.6e} and "
								   f"{target_time:.6e} s.")
			y_current = np.asarray(trial_solution.y[:, -1], dtype=float,).reshape(-1)
			elapsed_time = float(trial_solution.t[-1])
			if not np.all(np.isfinite(y_current)):
				raise FloatingPointError("Surface-only pre-equilibration produced NaN or infinity.")
			# Remove only negligible negative integration noise.
			tiny_negative = ((y_current < 0.0) & (y_current >= -1.0e-12))
			y_current[tiny_negative] = 0.0
			if np.any(y_current < -1.0e-12):
				negative_indices = np.flatnonzero(y_current < -1.0e-12)
				negative_species = [adsorbates[index] for index in negative_indices]
				minimum_coverage = float(np.min(y_current))
				raise RuntimeError("Surface-only pre-equilibration produced substantially negative coverages. "
								   f"Species: {negative_species}; minimum coverage={minimum_coverage:.6e}.")

			occupied_fraction = float(np.dot(y_current, site_occupancies,))
			occupancy_derivative = relaxation_rhs(elapsed_time, y_current,)
			occupancy_rate = float(np.dot(occupancy_derivative, site_occupancies,))
			capacity_repair_tolerance = 1.0e-5
			if occupied_fraction > 1.0:
				capacity_excess = (occupied_fraction - 1.0)
				if (capacity_excess <= capacity_repair_tolerance and occupancy_rate <= 0.0):
					# The solver has crossed the capacity boundary slightly, but the vector field points back into the
					# physical region. Project the state onto occupied_fraction = 1.
					y_current = (y_current / occupied_fraction)
					occupied_fraction = float(np.dot(y_current, site_occupancies,))
					print(f"  Corrected small surface-capacity overshoot: excess={capacity_excess:.6e}, "
						  f"d(occupied)/dt={occupancy_rate:.6e} s^-1, corrected occupied="
						  f"{occupied_fraction:.6e}", flush=True,)
				else:
					raise RuntimeError("Surface-only pre-equilibration exceeded the surface capacity at "
									   f"t={elapsed_time:.6e} s: occupied fraction={occupied_fraction:.16e}, "
									   f"excess={capacity_excess:.6e}, d(occupied)/dt={occupancy_rate:.16e}.")

			final_diagnostics = (evaluate_relaxation(y_current, elapsed_time,))
			occupied_fraction = float(np.dot(y_current, site_occupancies,))
			occupancy_rate = float(np.dot(final_diagnostics["derivative"], site_occupancies,))
			converged = bool(final_diagnostics["converged"])
			#print(f"  pre-equilibration: method={selected_method}, t={elapsed_time:.3e} s, max|dtheta/dt|="
			#      f"{final_diagnostics['maximum_absolute_rate']:.3e},"
			#      f" max relative rate={final_diagnostics['maximum_relative_rate']:.3e} s^-1,"
			#      f" occupied={occupied_fraction:.6e}", flush=True,)

		# Final validation and report.
		'''final_occupied_fraction = float(np.dot(y_current, site_occupancies,))
		print("\n Prepared TPD initial state:", flush=True,)
		for (adsorbate_name, before, after, derivative,) in zip(adsorbates, y_initial, y_current,
																final_diagnostics["derivative"],):
			if (abs(before) > 1.0e-14 or abs(after) > 1.0e-14 or abs(derivative) > 1.0e-14):
				print(f"  {adsorbate_name:18s} before={before: .6e} after={after: .6e} dtheta/dt={derivative: .6e}",
					  flush=True,)'''
		if not converged:
			print("WARNING: Surface-only pre-equilibration reached its maximum permitted time without satisfying the requested "
				  "rate tolerance. The final relaxed state will still be used for the TPD.", flush=True,)
		metadata = {"enabled": True, "surface_only": True, "gas_evolution_allowed": False, "temperature_K": float(temp_i),
					"converged": bool(converged), "elapsed_time_s": float(elapsed_time),
					"initial_occupied_fraction": (occupied_fraction_initial),
					"prepared_occupied_fraction": float(np.dot(y0, site_occupancies,)),
					"initial_maximum_absolute_rate": (initial_diagnostics["maximum_absolute_rate"]),
					"final_maximum_absolute_rate": (final_diagnostics["maximum_absolute_rate"]),
					"initial_maximum_relative_rate": (initial_diagnostics["maximum_relative_rate"]),
					"final_maximum_relative_rate": (final_diagnostics["maximum_relative_rate"]),
					"absolute_rate_tolerance": float(absolute_rate_tolerance),
					"relative_rate_tolerance": float(relative_rate_tolerance), "rtol": float(rtol), "atol": float(atol),
					"rhs_calls": int(rhs_call_count), "attempts": attempts,}
		return y_current, metadata

	def run_one_tpr(self, *, gas_name: str, adsorbates: list[str], ads_idx: int, gas_species: list[str], rhs_func,
					surface_rhs_func, jac_func, gas_rate_func: dict[str, Any], temperatures: np.ndarray,
					beta_min: float,) -> TPRResult | None:
		""" Run one temperature-programmed reaction/desorption simulation using the linear temperature ramp:
		T(t) = T_initial + beta * t
		gas_name:         Gas species associated with the initially populated adsorbate.
		adsorbates:       Ordered adsorbate names corresponding to the state vector.
		ads_idx:       Position of the initially populated adsorbate in the state vector.
		gas_species:        Gas-phase species for which TPR spectra are evaluated.
		rhs_func:        Lambdified kinetic RHS with arguments:   rhs_func(time, temperature, *coverages)
		jac_func:        Lambdified symbolic Jacobian. Accepted for compatibility with build(), but NONE here.
		gas_rate_func:        Dictionary mapping gas names to lambdified production-rate expressions.
		temperatures:        Requested output-temperature grid in K.
		beta_min:        Heating rate in K/min.    """

		#_ = jac_func
		temperatures = np.asarray(temperatures, dtype=float,).reshape(-1)
		if temperatures.size < 2:
			raise ValueError("TPR requires at least two temperature points.")
		if not np.all(np.isfinite(temperatures)):
			raise ValueError("TPR temperature grid contains NaN or infinity.")
		temperature_steps = np.diff(temperatures)
		if np.any(temperature_steps <= 0.0):
			raise ValueError("TPR temperatures (K) must be strictly increasing.")

		beta_min = float(beta_min)
		if not np.isfinite(beta_min) or beta_min <= 0.0:
			raise ValueError(f"TPR heating rate must be finite and positive; got {beta_min!r} K/min.")
		beta = beta_min / 60.0    # Convert K/min to K/s.

		temp_i = float(temperatures[0])
		temp_f = float(temperatures[-1])
		times = (temperatures - temp_i) / beta        # t(T) = (T - T_initial)/beta
		time_i = float(times[0])
		time_f = float(times[-1])
		time_span = (time_i, time_f,)

		# Validate and construct the initial state.
		n_adsorbates = len(adsorbates)
		if n_adsorbates == 0:
			raise ValueError("TPR requires at least one adsorbate.")
		if not 0 <= ads_idx < n_adsorbates:
			raise IndexError(f"ads_idx={ads_idx} is outside the valid range 0–{n_adsorbates - 1}.")
		site_occupancies = np.asarray([nsites(self.model, name,) for name in adsorbates], dtype=float,)
		if site_occupancies.shape != (n_adsorbates,):
			raise ValueError(f"Unexpected site-occupancy array shape: {site_occupancies.shape}; expected "
							 f"({n_adsorbates},).")
		if np.any(~np.isfinite(site_occupancies)):
			invalid_species = [name for name, occupancy in zip(adsorbates, site_occupancies,)
							   if not np.isfinite(occupancy)]
			raise ValueError(f"One or more adsorbates have non-finite nsites: {invalid_species}.")
		if np.any(site_occupancies <= 0.0):
			invalid_species = [name for name, occupancy in zip(adsorbates, site_occupancies,) if occupancy <= 0.0]
			raise ValueError(f"Adsorbate nsites must be positive. Invalid species: {invalid_species}.")
		initial_adsorbate = adsorbates[ads_idx]
		y0 = np.zeros(n_adsorbates, dtype=float,)

		# Initialise 95% of the available site capacity with the selected adsorbate.
		y0[ads_idx] = (0.95 / site_occupancies[ads_idx])
		total_initial_coverage = float(np.dot(y0, site_occupancies,))
		if total_initial_coverage > 1.0 + 1.0e-8:
			raise ValueError(f"Initial TPR coverage exceeds the available surface capacity: {total_initial_coverage:g}.")

		# Pre-equilibration
		y0, pre_equilibration_metadata = (self.pre_equilibrate_tpr_state(surface_rhs_func=surface_rhs_func, y0=y0,
																		 temp_i=temp_i, adsorbates=adsorbates,
																		 site_occupancies=site_occupancies,))
		if np.any(y0 < 0.0) or not np.all(np.isfinite(y0)):
			raise ValueError(f"Invalid initial TPR state: {y0}.")
		print("\n Starting TPR:"
			  f"\n  gas              = {gas_name}\n  initial adsorbate = {initial_adsorbate}"
			  f"\n  temperature       = {temp_i:g}–{temp_f:g} K\n  heating rate      = {beta_min:g} K/min"
			  f"\n  integration time  = {time_f:.3e} s\n  variables         = {n_adsorbates}", flush=True,)

		def temperature_profile(time_num: float | np.ndarray,) -> float | np.ndarray:
			""" Return the temperature corresponding to elapsed time. """
			values = (temp_i + beta * np.asarray(time_num, dtype=float,))
			if values.ndim == 0:
				return float(values)
			return values

		rhs_call_count = 0

		def ode_rhs(time_num: float, y: np.ndarray,) -> np.ndarray:
			""" Evaluate d(theta)/dt at one time and state. """
			nonlocal rhs_call_count
			rhs_call_count += 1
			time_num = float(time_num)
			temperature_num = float(temperature_profile(time_num))
			state = np.asarray(y, dtype=float,).reshape(-1)
			if state.size != n_adsorbates:
				raise ValueError("TPR solver supplied a state vector with the wrong "
								 f"size: expected {n_adsorbates}, got {state.size}.")
			if not np.all(np.isfinite(state)):
				raise FloatingPointError("TPR solver supplied a non-finite state at "
										 f"t={time_num:.8e} s and T={temperature_num:.8g} K.")
			values = np.asarray(rhs_func(time_num, temperature_num, *state,), dtype=float,).reshape(-1)
			if values.size != n_adsorbates:
				raise ValueError("TPR RHS returned the wrong number of components: "
								 f"expected {n_adsorbates}, got {values.size}.")
			if not np.all(np.isfinite(values)):
				invalid_indices = np.flatnonzero(~np.isfinite(values))
				invalid_names = [adsorbates[index] for index in invalid_indices]
				raise FloatingPointError(f"TPR RHS generated NaN or infinity at t={time_num:.8e} s and "
										 f"T={temperature_num:.8g} K. Invalid equations: {invalid_names}.")
			return values

		def ode_jacobian(time_num: float, state: np.ndarray,) -> np.ndarray:
			state = np.asarray(state, dtype=float,).reshape(-1)
			if state.size != n_adsorbates:
				raise ValueError(f"TPR Jacobian received a state of the wrong size: expected {n_adsorbates}, got {state.size}.")
			if not np.all(np.isfinite(state)):
				raise FloatingPointError(f"TPR Jacobian received a non-finite state at t={time_num:.8e} s.")
			temperature_num = float(temperature_profile(time_num))
			jacobian = np.asarray(jac_func(float(time_num), temperature_num, *state,), dtype=float,)
			expected_shape = (n_adsorbates, n_adsorbates,)
			if jacobian.shape != expected_shape:
				raise ValueError(f"TPR Jacobian returned the wrong shape: expected {expected_shape}, got {jacobian.shape}.")
			if not np.all(np.isfinite(jacobian)):
				raise FloatingPointError(f"TPR Jacobian generated NaN or infinity at t={time_num:.8e} s and "
										 f"T={temperature_num:.8g} K.")
			return jacobian

		# Report the initial derivatives before calling solve_ivp.
		initial_rhs = np.asarray(ode_rhs(time_i, y0,), dtype=float,).reshape(-1)
		print("\n Initial TPR derivatives:", flush=True,)
		for (adsorbate_name, initial_coverage, occupied_sites, derivative,) in zip(adsorbates, y0,
																				   site_occupancies, initial_rhs,):
			print(f"  {adsorbate_name:18s} theta={initial_coverage: .6e} nsites={occupied_sites: .6e}"
				  f" dtheta/dt={derivative: .6e}", flush=True,)
		print(f"  max |dtheta/dt| = {np.max(np.abs(initial_rhs)):.6e}", flush=True,)

		# Solver controls.
		# max_step is specified in seconds. Choose it so that one internal step corresponds to at most approximately 1 K:
		maximum_temperature_step = 1.0
		maximum_time_step = (maximum_temperature_step / beta)    #     delta_t = delta_T / beta
		absolute_tolerances = np.full(n_adsorbates, 1.0e-8, dtype=float,)
		absolute_tolerances[y0 > 1.0e-6] = 1.0e-10

		integration_options = {"fun": ode_rhs, "t_span": time_span, "y0": y0, "t_eval": times, "jac": ode_jacobian,
							   "rtol": 1.0e-5, "atol": absolute_tolerances, "max_step": maximum_time_step,}

		# BDF and Radau will estimate their Jacobians from ode_rhs.
		solution = None
		attempted_methods: list[dict[str, Any]] = []
		selected_method: str | None = None
		for method in ("BDF", "Radau", "LSODA"):
			calls_before_attempt = (rhs_call_count)
			try:

				trial = solve_ivp(method=method, **integration_options,)

			except (FloatingPointError, ValueError, OverflowError,) as error:
				attempted_methods.append({"method": method, "success": False, "complete": False,
										  "message": (f"{type(error).__name__}: {error}"),
										  "nfev": (rhs_call_count - calls_before_attempt),
										  "njev": 0, "nlu": 0, "returned_points": 0,})
				print(f"  {method}: raised "f"{type(error).__name__}: {error}", flush=True,)
				continue
			returned_points = len(trial.t)
			reached_final_time = (returned_points > 0 and
								  np.isclose(float(trial.t[-1]), time_f, rtol=0.0,
											 atol=max(1.0e-8, 1.0e-10 * max(1.0, abs(time_f),),),))
			complete = bool(trial.success and returned_points == times.size and reached_final_time)
			attempt_record = {"method": method, "success": bool(trial.success), "complete": complete,
							  "message": str(trial.message), "nfev": int(getattr(trial, "nfev", 0,)),
							  "njev": int(getattr(trial, "njev", 0,)),
							  "nlu": int(getattr(trial, "nlu", 0,)), "returned_points": (returned_points),}
			attempted_methods.append(attempt_record)
			print(f"\n  {method}: success={trial.success}, points={returned_points}/{times.size}, "
				  f"nfev={attempt_record['nfev']}, njev={attempt_record['njev']}, nlu={attempt_record['nlu']}\n",
				  flush=True,)
			if complete:
				solution = trial
				selected_method = method
				break
			if returned_points > 0:
				last_time = float(trial.t[-1])
				last_temperature = float(temperature_profile(last_time))
			else:
				last_time = time_i
				last_temperature = temp_i
			print(f"WARNING: TPR {method} failed or was incomplete for {gas_name}.\n"
				  f"  Last reached time        : {last_time:.8g} s\n"
				  f"  Last reached temperature : {last_temperature:.8g} K\n"
				  f"  Message                  : {trial.message}", flush=True,)

		if solution is None:
			print(f"WARNING: No TPR solver completed the {gas_name} calculation.", flush=True,)
			return None

		# Extract and validate the completed trajectory.
		solution_times = np.asarray(solution.t, dtype=float,).reshape(-1)
		result_temperatures = np.asarray(temperature_profile(solution_times), dtype=float,).reshape(-1)
		concentrations = np.asarray(solution.y.T, dtype=float,)
		expected_shape = (temperatures.size, n_adsorbates,)
		if concentrations.shape != expected_shape:
			raise ValueError(f"Unexpected TPR concentration array shape: {concentrations.shape}; expected "
							 f"{expected_shape}.")
		if solution_times.shape != times.shape:
			raise ValueError(f"Unexpected TPR time-array shape: {solution_times.shape}; expected "
							 f"{times.shape}.")
		if result_temperatures.shape != temperatures.shape:
			raise ValueError(f"Unexpected reconstructed temperature-array shape: {result_temperatures.shape}; expected "
							 f"{temperatures.shape}.")
		if not np.allclose(solution_times, times, rtol=0.0, atol=1.0e-8,):
			maximum_time_difference = float(np.max(np.abs(solution_times - times)))
			raise RuntimeError("TPR solver returned values at unexpected times. "
							   f"Maximum difference: {maximum_time_difference:.6e} s.")
		if not np.allclose(result_temperatures, temperatures, rtol=0.0, atol=1.0e-8,):
			maximum_temperature_difference = float(np.max(np.abs(result_temperatures - temperatures)))
			raise RuntimeError( "TPR solver times do not reproduce the requested "
								f"temperature grid. Maximum difference: {maximum_temperature_difference:.6e} K.")
		if not np.all(np.isfinite(concentrations)):
			raise FloatingPointError("Completed TPR trajectory contains NaN or infinity.")

		# Clean only small numerical violations after integration.
		negative_tolerance = 1.0e-10
		substantially_negative = (concentrations < -negative_tolerance)
		if np.any(substantially_negative):
			row_indices, column_indices = np.nonzero(substantially_negative)
			worst_position = int(np.argmin(concentrations[row_indices, column_indices,]))
			worst_row = int(row_indices[worst_position])
			worst_column = int(column_indices[worst_position])
			print("WARNING: Completed TPR trajectory contains negative coverages.\n"
				  f" Most negative species     : {adsorbates[worst_column]}\n"
				  f" Temperature               : {temperatures[worst_row]:.8g} K\n"
				  f" Coverage                  : {concentrations[worst_row, worst_column]:.6e}", flush=True,)
		concentrations = np.maximum(concentrations, 0.0,)

		# Enforce site capacity on the stored output; the solver's internal trajectory is not rescaled.
		occupied_fraction = (concentrations @ site_occupancies)
		overfilled = (occupied_fraction > 1.0)
		if np.any(overfilled):
			maximum_occupied_fraction = float(np.max(occupied_fraction[overfilled]))
			print(f"WARNING: Rescaling {np.count_nonzero(overfilled)} TPR output states that exceed the surface capacity."
				  f"Maximum occupied fraction before rescaling: {maximum_occupied_fraction:.6e}.", flush=True,)
			concentrations[overfilled, :,] /= occupied_fraction[overfilled, None,]

		# Evaluate the gas-production spectra.
		spectra = self.evaluate_tpr_spectra(gas_species=gas_species, gas_rate_func=gas_rate_func,
											temperatures=temperatures, concentrations=concentrations,
											temp_i=temp_i, beta=beta,)

		return TPRResult(gas=gas_name, heating_rate=beta_min, temperatures=temperatures.copy(),
						 adsorbates=list(adsorbates), gas_species=list(gas_species), concentrations=concentrations,
						 spectra=spectra, metadata={"initial_adsorbate": (initial_adsorbate), "beta_K_per_s": beta,
													"integration_variable": "time", "integration_time_s": time_f,
													"solver_method": selected_method, "solver_success": True,
													"solver_attempts": attempted_methods,
													"initial_occupied_fraction": (total_initial_coverage),
													"maximum_internal_step_s": (maximum_time_step),
													"maximum_internal_step_K": (maximum_temperature_step),
													"rtol": integration_options["rtol"],
													"atol": integration_options["atol"],},)

	def evaluate_tpr_spectra(self, *, gas_species: list[str], gas_rate_func: dict[str, Any], temperatures: np.ndarray,
							 concentrations: np.ndarray, temp_i: float, beta: float,) -> dict[str, np.ndarray]:
		"""Evaluate gas-production rates along a completed TPR trajectory."""
		temperatures = np.asarray(temperatures, dtype=float,)
		concentrations = np.asarray(concentrations, dtype=float,)
		if temperatures.ndim != 1:
			raise ValueError("TPR temperatures must be one-dimensional.")
		if concentrations.ndim != 2:
			raise ValueError("TPR concentrations must be a two-dimensional array.")
		if concentrations.shape[0] != temperatures.size:
			raise ValueError("TPR concentration and temperature arrays have "
							 f"different lengths: {concentrations.shape[0]} and {temperatures.size}.")
		if not np.isfinite(beta) or beta <= 0.0:
			raise ValueError(f"TPR heating rate must be positive, got {beta!r} K/s.")
		times = (temperatures - float(temp_i)) / float(beta)
		coverage_arguments = [concentrations[:, index] for index in range(concentrations.shape[1])]
		spectra: dict[str, np.ndarray] = {}
		for gas in gas_species:
			func = gas_rate_func.get(gas)
			if func is None:
				spectra[gas] = np.zeros_like(temperatures, dtype=float,)
				continue
			try:
				# Fast vectorised evaluation.
				rates = np.asarray(func(times, temperatures, *coverage_arguments,), dtype=float,)
				rates = np.squeeze(rates)
				if rates.ndim == 0:
					rates = np.full(temperatures.shape, float(rates), dtype=float,)
				else:
					rates = np.broadcast_to(rates, temperatures.shape,).astype(float, copy=True,)
			except (TypeError, ValueError):
				# Fallback for lambdified expressions that do not
				# support array-valued inputs.
				rates = np.empty_like(temperatures, dtype=float,)
				for index, (time_num, temperature_num,) in enumerate(zip(times, temperatures)):
					value = func(float(time_num), float(temperature_num), *concentrations[index],)
					value_array = np.asarray(value, dtype=float,).reshape(-1)
					if value_array.size != 1:
						raise ValueError(f"Gas-rate function for {gas!r} returned "
										 f"{value_array.size} values at one temperature; expected one scalar.")
					rates[index] = float(value_array[0])
			rates = np.nan_to_num(rates, nan=0.0, posinf=0.0, neginf=0.0,)
			# TPD/TPR spectra represent gas production. Negative
			# values correspond to gas consumption and are discarded.
			spectra[gas] = np.maximum(rates, 0.0,)
		return spectra

	def desorbing_adsorbate(self, gas_name: str) -> str:
		for process in self.model.processes.values():
			if process.kind.upper() != "D":
				continue
			product_names = [item.species for item in process.products]
			if gas_name not in product_names:
				continue
			for item in process.reactants:
				if is_adsorbate(self.model.species[item.species]):
					return item.species
		raise ValueError(f"No desorption process produces gas species {gas_name!r}.")
