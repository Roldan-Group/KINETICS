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

from Symbols_def import t, temp, constants


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
    raise ValueError("Rate equations are missing. Run REquations(model) before Experiments.")

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

def build_surface_expressions(model, dynamic_species: list[str], surfaces: list[str], adsorbates_by_surface: dict[str, list[int]], *, smooth: bool = False):
    """Build algebraic surface-coverage expressions."""
    expressions: dict[sp.Symbol, sp.Expr] = {}
    eps = sp.Float("1e-30")
    for surface_name in surfaces:
        coverage = sp.Float(1.0)
        for idx in adsorbates_by_surface[surface_name]:
            ads_name = dynamic_species[idx]
            coverage -= species_symbol(ads_name) * nsites(model, ads_name)
        if smooth:
            coverage = (coverage + sp.sqrt(coverage**2 + eps**2)) / 2
        expressions[species_symbol(surface_name)] = coverage
    return expressions

def clip_dynamic_species(model, dynamic_species: list[str], surfaces: list[str], adsorbates_by_surface: dict[str, list[int]], gas_indices: list[int], y: np.ndarray) -> np.ndarray:
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

def surface_coverages(model, y: np.ndarray, dynamic_species: list[str], surfaces: list[str], adsorbates_by_surface: dict[str, list[int]]) -> list[np.ndarray]:
    """Compute surface coverages from adsorbate coverages."""
    values: list[np.ndarray] = []
    for surface_name in surfaces:
        coverage = np.ones(y.shape[1], dtype=float)
        for idx in adsorbates_by_surface[surface_name]:
            ads_name = dynamic_species[idx]
            coverage -= y[idx, :] * nsites(model, ads_name)
        values.append(np.clip(coverage, 0.0, 1.0))
    return values

def solve_ode_system(model, t_span: tuple[float, float], t_eval: np.ndarray, dynamic_species: list[str],
                     surfaces: list[str], adsorbates_by_surface: dict[str, list[int]], gas_indices: list[int],
                     initial_conditions: list[float], rhs_func, jac_func, arguments: tuple[float, ...],):
    """Solve the microkinetic ODE system for one set of external arguments."""
    y0 = np.asarray(initial_conditions, dtype=float)

    def ode_system(t_num, y, *args):
        dydt = rhs_func(t_num, *args, *y)
        dydt = np.asarray(dydt, dtype=float)
        if not np.all(np.isfinite(dydt)):
            raise RuntimeError(f"Non-finite RHS at t={t_num}; y={y}; dydt={dydt}")
        return dydt

    def ode_jacobian(t_num, y, *args):
        jac = jac_func(t_num, *args, *y)
        return np.asarray(jac, dtype=float)

    sol = solve_ivp(ode_system, t_span, y0, t_eval=t_eval, args=arguments, method="LSODA", jac=ode_jacobian,
                    rtol=1e-6,   atol=1e-8,)

    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")
    y_dynamic = clip_dynamic_species(model, dynamic_species, surfaces, adsorbates_by_surface, gas_indices, sol.y)
    y_surfaces = surface_coverages(model, y_dynamic, dynamic_species, surfaces, adsorbates_by_surface)
    augmented = SimpleNamespace()
    augmented.t = sol.t
    augmented.y = np.vstack([y_dynamic, *y_surfaces])
    augmented.species = dynamic_species + surfaces
    rows = np.column_stack([np.tile(arguments, (len(sol.t), 1)), augmented.t, augmented.y.T,])
    return sol, rows, augmented


class Isothermal:
    """Run constant-temperature/time-grid integrations and store results in memory."""
    def __init__(self, model):
        self.model = model
        self.result = self.build()
        ensure_experiment_results(model).isothermal = self.result

    def build(self) -> TimeSeriesResult:
        equations, _, _, _ = get_equation_tuple(self.model)
        initial_conditions, dynamic_species, surfaces, gas_indices, adsorbates_by_surface = initial_species(self.model)
        surface_expr = build_surface_expressions(self.model, dynamic_species, surfaces, adsorbates_by_surface)
        rhs = []
        for name in dynamic_species:
            expr = sp.Integer(0)
            for contribution in equations.get(name, []):
                expr += contribution
            expr = expr.subs(surface_expr).subs(constants)
            rhs.append(expr)
        dynamic_symbols = [species_symbol(name) for name in dynamic_species]
        conditions = [t, temp]
        rhs_func = sp.lambdify((*conditions, *dynamic_symbols), rhs, [{"exp": safe_exp}, "numpy"])
        jacobian = sp.Matrix(rhs).jacobian(dynamic_symbols)
        jac_func = sp.lambdify((*conditions, *dynamic_symbols), jacobian, [{"exp": safe_exp}, "numpy"])
        time_values = scan_values(self.model.conditions.time)
        if len(time_values) < 2:
            raise ValueError("Isothermal simulation requires at least two time points.")
        t_eval = np.asarray(time_values, dtype=float)
        t_span = (float(t_eval[0]), float(t_eval[-1]))
        temperature_values = scan_values(self.model.conditions.temperature, default=[300.0])

        rows = []
        for temperature in temperature_values:
            _, data, _ = solve_ode_system(self.model, t_span, t_eval, dynamic_species, surfaces, adsorbates_by_surface,
                                          gas_indices, initial_conditions, rhs_func, jac_func, (float(temperature),),)
            rows.append(data)
        values = np.vstack(rows) if rows else np.empty((0, 0))
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

    def build(self) -> list[TPRResult]:
        _, _, tpd_equations, tpd_factors = get_equation_tuple(self.model)
        initial_conditions, dynamic_species, surfaces, gas_indices, adsorbates_by_surface = initial_species(self.model)
        gas_species = [dynamic_species[idx] for idx in gas_indices]
        adsorbates = [name for name in dynamic_species if is_adsorbate(self.model.species[name])]
        ads_symbols = [species_symbol(name) for name in adsorbates]
        surface_expr = build_surface_expressions(self.model, dynamic_species, surfaces, adsorbates_by_surface,
                                                 smooth=True,)
        rhs_ads = []
        max_rate = sp.Float(1e7)
        for name in adsorbates:
            expr = sp.Integer(0)
            for contribution in tpd_equations.get(name, []):
                expr += max_rate * sp.tanh(contribution / max_rate)
            expr = expr.subs(surface_expr).subs(constants)
            rhs_ads.append(expr)
        rhs_func = sp.lambdify((t, temp, *ads_symbols), rhs_ads, [{"exp": safe_exp}, "numpy"])
        gas_rate_expr: dict[str, sp.Expr] = {}
        for gas in gas_species:
            expr = sp.Integer(0)
            for factor, contribution in zip(tpd_factors.get(gas, []), tpd_equations.get(gas, [])):
                if factor > 0:
                    expr += contribution
            gas_rate_expr[gas] = expr.subs(surface_expr).subs(constants)

        gas_rate_func = {gas: sp.lambdify((t, temp, *ads_symbols), expr, [{"exp": safe_exp}, "numpy"])
                         for gas, expr in gas_rate_expr.items()}
        temperatures = self.temperature_values()
        if len(temperatures) < 2:
            raise ValueError("TPR requires at least two temperature points.")
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
                                          gas_species=gas_species, rhs_func=rhs_func, gas_rate_func=gas_rate_func,
                                          temperatures=temperatures, beta_min=float(beta_min),)
                results.append(result)
        return results

    def temperature_values(self) -> np.ndarray:
        if self.temperature_scan is not None:
            return np.asarray(scan_values(self.temperature_scan), dtype=float)
        if self.model.conditions.temperature is not None:
            values = scan_values(self.model.conditions.temperature)
            if len(values) >= 2:
                return np.asarray(values, dtype=float)
        return np.arange(10.0, 1273.0 + 1.0, 1.0)

    def run_one_tpr(self, *, gas_name: str, adsorbates: list[str], ads_idx: int, gas_species: list[str], rhs_func,
                    gas_rate_func: dict[str, Any], temperatures: np.ndarray, beta_min: float,) -> TPRResult:
        beta = beta_min / 60.0
        temp_i = float(temperatures[0])
        temp_f = float(temperatures[-1])
        t_eval = (temperatures - temp_i) / beta
        t_span = (float(t_eval[0]), float(t_eval[-1]))
        y0 = np.zeros(len(adsorbates), dtype=float)
        nmolsite = float(getattr(self.model.species[gas_name], "nmolsite", 1.0))
        y0[ads_idx] = 0.95 / nmolsite  # not 1 to give some wiggle

        def temperature_profile(t_num):
            return temp_i + beta * t_num

        def ode_rhs(t_num, y):
            temperature_num = temperature_profile(t_num)
            dydt = np.asarray(rhs_func(t_num, temperature_num, *y), dtype=float)
            if not np.all(np.isfinite(dydt)):
                raise RuntimeError(f"Non-finite TPR RHS at t={t_num}; T={temperature_num}; y={y}; dydt={dydt}")
            return dydt

        sol = solve_ivp(ode_rhs, t_span, y0, t_eval=t_eval, method="BDF", rtol=1e-3, atol=1e-8,)

        if not sol.success:
            raise RuntimeError(f"TPR solver failed for {gas_name}: {sol.message}")

        concentrations = np.zeros((len(temperatures), len(adsorbates)), dtype=float)
        if sol.y.size:
            concentrations[: len(sol.t), :] = sol.y.T
        concentrations = np.clip(concentrations, 0.0, None)
        for row in concentrations:
            total = sum(row[i] * nsites(self.model, ads) for i, ads in enumerate(adsorbates))
            if total > 1.0:
                row[:] /= total

        spectra: dict[str, np.ndarray] = {}
        for gas in gas_species:
            rates = np.zeros(len(temperatures), dtype=float)
            func = gas_rate_func.get(gas)
            if func is None:
                spectra[gas] = rates
                continue
            for i, (t_num, temperature_num) in enumerate(zip(t_eval, temperatures)):
                rates[i] = float(func(t_num, temperature_num, *concentrations[i]))
            spectra[gas] = np.maximum(rates, 0.0)

        return TPRResult(gas=gas_name, heating_rate=beta_min, temperatures=temperatures, adsorbates=adsorbates,
                         gas_species=gas_species, concentrations=concentrations, spectra=spectra,
                         metadata={"initial_adsorbate": adsorbates[ads_idx], "beta_K_per_s": beta,},)

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
