"""
Object-based diagnostics and reporting for the refactored microkinetic workflow.

Designed for Main_refactored.py, where stages are called as Routine(model).

This module has two responsibilities:
    1. Analyse in-memory isothermal and TPR/TPD results stored in model.experiments.
    2. Build printable/exportable tables for the model state, thermodynamics,
       kinetics, rate equations, profiles and experiment results.

No files are written and nothing is printed unless the user explicitly calls
`DiagnosticsPrinter.write_all(...)` or `DiagnosticsPrinter.print_all(...)`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any, Iterable

import csv
import math

import numpy as np
import sympy as sp
from scipy.optimize import least_squares
from scipy.signal import find_peaks

from Symbols_def import temp, constants



@dataclass(slots=True)
class TableResult:
    """Simple in-memory table."""
    headers: list[str]
    rows: list[list[Any]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

@dataclass(slots=True)
class ReactionPathwayResult:
    """Flux-weighted reaction-pathway edges for one temperature."""
    temperature: float
    edges: list[tuple[str, str, float]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

@dataclass(slots=True)
class IRSpectrumResult:
    """Calculated IR intensity map for a selected phase and temperature."""
    phase: str
    temperature: float
    wavenumbers: np.ndarray
    times: np.ndarray
    intensities: np.ndarray
    species: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

@dataclass(slots=True)
class TPRPeakResult:
    """Peak summary for a TPR/TPD spectrum."""
    source_gas: str
    signal_gas: str
    heating_rate: float
    peak_temperature: float
    peak_rate: float
    area: float
    metadata: dict[str, Any] = field(default_factory=dict)

@dataclass(slots=True)
class DiagnosticResults:
    """Main diagnostics container attached to model.diagnostics."""
    steady_state_net_rates: TableResult | None = None
    average_net_rates: TableResult | None = None
    equilibrium_control: TableResult | None = None
    drc_assessment: TableResult | None = None
    degree_of_rate_control: dict[str, TableResult] = field(default_factory=dict)
    degree_of_selectivity_control: dict[str, TableResult] = field(default_factory=dict)
    reaction_pathways: list[ReactionPathwayResult] = field(default_factory=list)
    ir_spectra: list[IRSpectrumResult] = field(default_factory=list)
    tpr_peaks: list[TPRPeakResult] = field(default_factory=list)
    model_tables: dict[str, TableResult] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)


# --- GENERIC
def ensure_diagnostics_results(model) -> DiagnosticResults:
    if not hasattr(model, "diagnostics") or model.diagnostics is None:
        model.diagnostics = DiagnosticResults()
    return model.diagnostics

def kind(species) -> str:
    return str(getattr(species, "kind", "")).lower()

def is_molecule(species) -> bool:
    return kind(species) == "molecule"


def is_surface(species) -> bool:
    return kind(species) == "surface"

def is_adsorbate(species) -> bool:
    return kind(species) in {"adsorbate", "transitionstate", "ts"}

def species_symbol(name: str) -> sp.Symbol:
    return sp.Symbol(name)

def nsites(model, name: str) -> float:
    value = getattr(model.species[name], "nsites", None)
    return float(value) if value is not None else 1.0

def safe_str(value: Any) -> str:
    if isinstance(value, np.ndarray):
        return np.array2string(value, precision=6, threshold=20)
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(safe_str(v) for v in value) + "]"
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)

def numeric_expr(expr, temperature: float) -> float:
    expr = sp.sympify(expr).subs(constants)
    func = sp.lambdify(temp, expr, ("numpy", "sympy"))
    return float(func(temperature))

def get_kinetics_results(model):
    if hasattr(model, "kinetics") and model.kinetics is not None:
        return model.kinetics
    try:
        from Kinetics import get_rate_constants  # type: ignore
    except Exception:
        return None
    return get_rate_constants(model)

def get_rate_equations(model):
    if hasattr(model, "equations") and model.equations is not None:
        equations = model.equations
        if isinstance(equations, tuple) and len(equations) == 4:
            return equations
        if all(hasattr(equations, attr) for attr in ("isothermal", "equation_factors", "tpd", "tpd_factors")):
            return (equations.isothermal, equations.equation_factors, equations.tpd, equations.tpd_factors,)
    try:
        from Kinetics import get_rate_equations as _get_rate_equations  # type: ignore
    except Exception:
        return None
    return _get_rate_equations(model)

def get_profile_data(model):
    if hasattr(model, "profile") and model.profile is not None:
        return model.profile
    try:
        from Kinetics import get_profile_data as _get_profile_data  # type: ignore
    except Exception:
        return None
    return _get_profile_data(model)

def dynamic_species_from_processes(model) -> list[str]:
    seen: list[str] = []
    for process in model.processes.values():
        for item in [*process.reactants, *process.products]:
            if item.species not in seen and not is_surface(model.species[item.species]):
                seen.append(item.species)
    return sorted(seen)

def surfaces_from_processes(model) -> list[str]:
    names: set[str] = set()
    for process in model.processes.values():
        for item in [*process.reactants, *process.products]:
            if is_surface(model.species[item.species]):
                names.add(item.species)
    return sorted(names)

def adsorbates_by_surface(model, dynamic_species: list[str], surfaces: list[str]) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for surface_name in surfaces:
        site = getattr(model.species[surface_name], "sites", None)
        result[surface_name] = [i for i, name in enumerate(dynamic_species) if is_adsorbate(model.species[name])
                                and getattr(model.species[name], "sites", None) == site]
    return result

def build_surface_expressions(model, dynamic_species: list[str], surfaces: list[str], ads_by_surface: dict[str, list[int]]):
    expressions: dict[sp.Symbol, sp.Expr] = {}
    for surface_name in surfaces:
        coverage = sp.Float(1.0)
        for idx in ads_by_surface[surface_name]:
            ads_name = dynamic_species[idx]
            coverage -= species_symbol(ads_name) * nsites(model, ads_name)
        expressions[species_symbol(surface_name)] = coverage
    return expressions

def clip_concentrations(model, dynamic_species: list[str], surfaces: list[str], ads_by_surface: dict[str, list[int]],
                        gas_indices: list[int], y: np.ndarray):
    y = np.asarray(y, dtype=float).copy()
    y =  np.nan_to_num(y, nan=0.0, neginf=0.0, posinf=1.0e30,)
    for surface_name in surfaces:
        indices = ads_by_surface[surface_name]
        total = sum(y[idx] * nsites(model, dynamic_species[idx]) for idx in indices)
        if total > 1.0:
            for idx in indices:
                y[idx] /= total
    for idx in gas_indices:
        y[idx] = max(y[idx], 0.0)
    return y


# --- DIAGNOSTICS
class Diagnostics:
    """
    Analyse stored results from Isothermal and TPR.
    The class stores all derived quantities in `model.diagnostics` and returns the same container in `self.results`.
    """
    def __init__(self, model, *, analyse_isothermal=True, analyse_tpr=True,):
        self.model = model
        self.results = ensure_diagnostics_results(model)
        if analyse_isothermal:
            self._safe(self.analyse_isothermal)
        if analyse_tpr:
            self._safe(self.analyse_tpr)

    def _safe(self, func):
        try:
            return func()
        except Exception as err:  # keep diagnostics non-fatal for workflow use
            self.results.errors.append(f"{func.__name__}: {err}")
            return None

    # --- Isothermal analysis
    def analyse_isothermal(self):
        experiments = getattr(self.model, "experiments", None)
        result = getattr(experiments, "isothermal", None)
        if result is None:
            return None
        equations = get_rate_equations(self.model)
        if equations is None:
            raise ValueError("Rate equations are unavailable. Run REquations(model) before Diagnostics(model).")
        context = self.build_context(result)
        functions, derivatives = self.get_functions(context)
        self.results.steady_state_net_rates, self.results.average_net_rates = self.rates(context, functions)
        self.results.equilibrium_control = self.equilibrium_control(context, functions, derivatives)
        self.results.drc_assessment = self.assess_drc_validity(context, functions, derivatives)
        self.kinetic_control(context, functions, derivatives)
        self.reaction_pathways(context, functions)
        self.ir_evolution(context)
        return self.results

    def build_context(self, isothermal_result):
        headers = list(isothermal_result.headers)
        values = np.asarray(isothermal_result.values, dtype=float)
        dynamic_species = list(isothermal_result.species)
        surfaces = list(isothermal_result.surfaces)
        gas_indices = list(isothermal_result.metadata.get("gas_indices", []))
        ads_by_surface = dict(isothermal_result.metadata.get("adsorbates_by_surface", {}))
        if not dynamic_species:
            dynamic_species = dynamic_species_from_processes(self.model)
        if not surfaces:
            surfaces = surfaces_from_processes(self.model)
        if not ads_by_surface:
            ads_by_surface = adsorbates_by_surface(self.model, dynamic_species, surfaces)
        if not gas_indices:
            gas_indices = [i for i, name in enumerate(dynamic_species) if is_molecule(self.model.species[name])]
        temp_idx = headers.index("Temperature")
        time_idx = headers.index("time") if "time" in headers else headers.index("Time")
        species_indices = [headers.index(name) for name in dynamic_species]
        temperatures = sorted(float(x) for x in np.unique(values[:, temp_idx]))
        times = sorted(float(x) for x in np.unique(values[:, time_idx]))
        concentration_by_temperature: dict[float, dict[float, np.ndarray]] = {temp_num: {} for temp_num in
                                                                              temperatures}
        for row in values:
            temp_num = float(row[temp_idx])
            time_value = float(row[time_idx])
            concentration_by_temperature[temp_num][time_value] = np.asarray(row[species_indices], dtype=float)
        surface_expr = build_surface_expressions(self.model, dynamic_species, surfaces, ads_by_surface)

        return {"headers": headers, "values": values, "dynamic_species": dynamic_species, "surfaces": surfaces,
                "gas_indices": gas_indices, "adsorbates_by_surface": ads_by_surface, "temperatures": temperatures,
                "times": times, "concentrations": concentration_by_temperature, "surface_expr": surface_expr,}

    def get_functions(self, context):
        equations = get_rate_equations(self.model)
        if equations is None:
            raise ValueError("Rate equations are unavailable.")
        _, equation_factors, _, _ = equations
        species_names = context["dynamic_species"]
        species_syms = {name: species_symbol(name) for name in species_names}
        process_ids = list(self.model.processes.keys())
        process_index = {pid: i for i, pid in enumerate(process_ids)}
        k_syms = {pid: sp.Symbol(f"k_{pid}") for pid in process_ids}
        k_sym_list = [k_syms[pid] for pid in process_ids]
        concentration_vec = sp.Matrix([species_syms[name] for name in species_names])
        all_variables = (temp, *concentration_vec, *k_sym_list)
        step_pairs = [(process_ids[i], process_ids[i + 1]) for i in range(0, len(process_ids) - 1, 2)]
        step_labels = [f"$r_{{{a}}} - r_{{{b}}}$" for a, b in step_pairs]
        rate_exprs: dict[str, sp.Expr] = {}
        for pid, process in self.model.processes.items():
            expr = k_syms[pid]
            for item in process.reactants:
                sym = species_syms.get(item.species)
                if sym is not None:
                    expr *= sym ** item.coefficient
                else:
                    expr *= context["surface_expr"][species_symbol(item.species)] ** item.coefficient
            rate_exprs[pid] = expr

        f_list = []
        for name in species_names:
            expr = sp.Integer(0)
            factors = equation_factors.get(name, [])
            for pid in process_ids:
                j = process_index[pid]
                factor = factors[j] if j < len(factors) else 0.0
                expr += factor * rate_exprs[pid]
            f_list.append(expr)
        f_vec = sp.Matrix(f_list)
        f_func = sp.lambdify(all_variables, f_vec, "numpy")
        products = [name for name, species in self.model.species.items() if is_molecule(species)]
                    # and float(getattr(species, "pressure0", 0.0)) == 0.0] Alberto 21/07/2026 test
        stoichi_pairs = np.array([[(equation_factors.get(product, [0.0] * len(process_ids))[process_index[forward]]
                                    - equation_factors.get(product, [0.0] * len(process_ids))[process_index[backward]])
                                   for forward, backward in step_pairs] for product in products], dtype=float)
        s_plus_pairs = np.array([[max(equation_factors.get(product, [0.0] * len(process_ids))[process_index[forward]],
                                      0.0) for forward, _ in step_pairs] for product in products], dtype=float)
        kinetics = get_kinetics_results(self.model)

        if kinetics is None:
            raise ValueError("Kinetic results are unavailable. Run RConstants(model) before Diagnostics(model).")
        k_func = {}
        for pid in process_ids:
            kin = kinetics.by_process[pid]
            k_func[pid] = sp.lambdify(temp, sp.sympify(kin.krate0).subs(constants), "numpy")
        r_forward_exprs = []
        r_backward_exprs = []
        for forward, backward in step_pairs:
            r_forward_exprs.append(rate_exprs[forward])
            r_backward_exprs.append(rate_exprs[backward])
        r_forward_vec = sp.Matrix(r_forward_exprs)
        r_backward_vec = sp.Matrix(r_backward_exprs)
        r_forward_func = sp.lambdify(all_variables, r_forward_vec, "numpy")
        r_backward_func = sp.lambdify(all_variables, r_backward_vec, "numpy")
        beta_vec = sp.Matrix([r_backward_exprs[i] / (r_forward_exprs[i] + sp.Float("1e-300"))
                              for i in range(len(step_pairs))])
        r_net_vec = sp.Matrix([r_forward_exprs[i] * (1 - beta_vec[i]) for i in range(len(step_pairs))])
        j_sym = f_vec.jacobian(concentration_vec)
        df_dk_sym = f_vec.jacobian(k_sym_list)
        dr_net_dk_sym = r_net_vec.jacobian(k_sym_list)
        dr_net_dtheta_sym = r_net_vec.jacobian(concentration_vec)
        dr_forward_dk_sym = r_forward_vec.jacobian(k_sym_list)
        dr_forward_dtheta_sym = r_forward_vec.jacobian(concentration_vec)
        derivatives = {"jacobian": sp.lambdify(all_variables, j_sym, "numpy"),
                       "df_dk": sp.lambdify(all_variables, df_dk_sym, "numpy"),
                       "dr_net_dk": sp.lambdify(all_variables, dr_net_dk_sym, "numpy"),
                       "dr_net_dtheta": sp.lambdify(all_variables, dr_net_dtheta_sym, "numpy"),
                       "dr_forward_dk": sp.lambdify(all_variables, dr_forward_dk_sym, "numpy"),
                       "dr_forward_dtheta": sp.lambdify(all_variables, dr_forward_dtheta_sym, "numpy"),}
        functions = {"process_ids": process_ids, "process_index": process_index, "k_func": k_func, "f_vec": f_vec,
                     "f_func": f_func, "rate_exprs": rate_exprs, "r_forward_func": r_forward_func,
                     "r_backward_func": r_backward_func, "step_pairs": step_pairs, "step_labels": step_labels,
                     "products": products, "stoichi_pairs": stoichi_pairs, "s_plus_pairs": s_plus_pairs,
                     "all_variables": all_variables,}
        return functions, derivatives

    def k_values(self, temperature_value: float, functions) -> np.ndarray:
        return np.asarray([functions["k_func"][pid](temperature_value) for pid in functions["process_ids"]], dtype=float)

    def rates_at(self, temperature_value: float, concentration: np.ndarray, functions):
        k_values = self.k_values(temperature_value, functions)
        all_values = (temperature_value, *concentration, *k_values)
        r_forward = np.asarray(functions["r_forward_func"](*all_values), dtype=float).flatten()
        r_backward = np.asarray(functions["r_backward_func"](*all_values), dtype=float).flatten()
        r_net = r_forward - r_backward
        return r_forward, r_backward, r_net, k_values, all_values

    def rates(self, context, functions) -> tuple[TableResult, TableResult]:
        times = context["times"]
        if len(times) < 2:
            raise ValueError("Rate diagnostics require at least two time points.")
        dt = times[-1] - times[0]
        if dt == 0.0:
            raise ValueError("Rate diagnostics require a non-zero time interval.")
        headers = ["Temperature", *functions["step_labels"]]
        ss_table = TableResult(headers=headers)
        avg_table = TableResult(headers=headers)

        for temp_num in [context["temperatures"][0], context["temperatures"][-1]]:
            r_net_t = []
            for time_value in times:
                concentration = context["concentrations"][temp_num][time_value]
                _, _, r_net, _, _ = self.rates_at(temp_num, concentration, functions)
                r_net_t.append(r_net)
            r_net_t = np.asarray(r_net_t, dtype=float)
            n_tail = max(1, int(0.1 * len(r_net_t)))
            ss = np.mean(r_net_t[-n_tail:], axis=0)
            avg = np.trapz(r_net_t, times, axis=0) / dt
            ss_table.rows.append([temp_num, *ss.tolist()])
            avg_table.rows.append([temp_num, *avg.tolist()])
        return ss_table, avg_table

    def steady_state_point(self, temp_num: float, context, functions):
        concentration_guess = np.asarray(context["concentrations"][temp_num][context["times"][-1]], dtype=float)
        # Validation test
        invalid = ( ~np.isfinite(concentration_guess) | (concentration_guess < 0.0))
        if np.any(invalid):
            species_names = self.dynamic_species_names()
            print("\nInvalid steady-state initial guess:", flush=True,)
            for index in np.flatnonzero(invalid):
                name = (species_names[index] if index < len(species_names) else f"index_{index}")
                print(f"    {name}: {concentration_guess[index]:.6e}", flush=True, )

        concentration_guess = np.asarray(concentration_guess, dtype=float,).copy()
        concentration_guess = np.nan_to_num(concentration_guess, nan=0.0, neginf=0.0, posinf=1.0e30,)
        concentration_guess = np.maximum(concentration_guess, 0.0,)
        k_values = self.k_values(temp_num, functions)
        def residual(x):
            return np.asarray(functions["f_func"](temp_num, *x, *k_values), dtype=float).ravel()
        solution = least_squares(residual, concentration_guess, bounds=(0, np.inf))
        gas_indices = context["gas_indices"]
        concentration = clip_concentrations(self.model, context["dynamic_species"], context["surfaces"],
                                            context["adsorbates_by_surface"], gas_indices, solution.x,)
        return concentration, k_values, (temp_num, *concentration, *k_values)

    def equilibrium_control(self, context, functions, derivatives, eps: float = 1e-30) -> TableResult:
        headers = ["Temperature", *functions["step_labels"]]
        table = TableResult(headers=headers)
        for temp_num in [context["temperatures"][0], context["temperatures"][-1]]:
            concentration, k_values, all_values = self.steady_state_point(temp_num, context, functions)
            r_forward = np.asarray(functions["r_forward_func"](*all_values), dtype=float).flatten()
            r_backward = np.asarray(functions["r_backward_func"](*all_values), dtype=float).flatten()
            r_net = r_forward - r_backward
            r_vec = functions["stoichi_pairs"] @ r_net
            jacobian = np.asarray(derivatives["jacobian"](*all_values), dtype=float)
            df_dk = np.asarray(derivatives["df_dk"](*all_values), dtype=float)
            dtheta_dk = -np.linalg.pinv(jacobian) @ df_dk
            dr_net_dk = np.asarray(derivatives["dr_net_dk"](*all_values), dtype=float)
            dr_net_dtheta = np.asarray(derivatives["dr_net_dtheta"](*all_values), dtype=float)
            dr_net_total = (functions["stoichi_pairs"] @ dr_net_dk) + (functions["stoichi_pairs"] @ dr_net_dtheta) @ dtheta_dk
            dec = np.zeros(len(functions["step_pairs"]), dtype=float)
            safe_r = np.where(np.abs(r_vec) < eps, np.nan, r_vec)
            for j, (forward, backward) in enumerate(functions["step_pairs"]):
                kf = k_values[functions["process_index"][forward]]
                kb = k_values[functions["process_index"][backward]]
                drc_f = np.nanmean((kf * dr_net_total[:, functions["process_index"][forward]]) / safe_r)
                drc_b = np.nanmean((kb * dr_net_total[:, functions["process_index"][backward]]) / safe_r)
                dec[j] = drc_f - drc_b
            table.rows.append([temp_num, *dec.tolist()])
        return table

    def assess_drc_validity(self, context, functions, derivatives, *, eq_threshold=1e-3, flux_threshold=1e-6,
                            rc_threshold=0.1, eps=1e-30) -> TableResult:
        table = TableResult(headers=["Temperature", "Step", "Classification", "FluxImportance",
                                     "DistanceFromEquilibrium", "GlobalDistance", "DRCRecommended"])
        for temp_num in [context["temperatures"][0], context["temperatures"][-1]]:
            concentration, k_values, all_values = self.steady_state_point(temp_num, context, functions)
            r_forward = np.asarray(functions["r_forward_func"](*all_values), dtype=float).flatten()
            r_backward = np.asarray(functions["r_backward_func"](*all_values), dtype=float).flatten()
            r_net = r_forward - r_backward
            distance_eq = np.abs(r_net) / (np.abs(r_forward) + eps)
            flux_importance = np.abs(r_net) / (np.sum(np.abs(r_net)) + eps)
            global_distance = float(np.sum(np.abs(r_net)) / (np.sum(np.abs(r_forward)) + eps))
            labels = []
            for i, label in enumerate(functions["step_labels"]):
                if flux_importance[i] < flux_threshold:
                    classification = "dormant"
                elif distance_eq[i] < eq_threshold:
                    classification = "quasi_equilibrated"
                elif flux_importance[i] > rc_threshold:
                    classification = "rate_controlling"
                else:
                    classification = "intermediate"
                labels.append(classification)
            recommended = not (global_distance < 1e-6 or labels.count("rate_controlling") == 0)
            for i, label in enumerate(functions["step_labels"]):
                table.rows.append([temp_num, label, labels[i], float(flux_importance[i]), float(distance_eq[i]),
                                   global_distance, recommended])
        return table

    def kinetic_control(self, context, functions, derivatives, eps: float = 1e-30):
        products = functions["products"]
        if not products:
            return
        headers = ["Temperature", *functions["step_labels"]]
        drc_data = {product: TableResult(headers=headers) for product in products}
        dsc_data = {product: TableResult(headers=headers) for product in products}
        for temp_num in [context["temperatures"][0], context["temperatures"][-1]]:
            concentration, k_values, all_values = self.steady_state_point(temp_num, context, functions)
            r_forward = np.asarray(functions["r_forward_func"](*all_values), dtype=float).flatten()
            r_backward = np.asarray(functions["r_backward_func"](*all_values), dtype=float).flatten()
            r_net = r_forward - r_backward
            r_vec = functions["stoichi_pairs"] @ r_net
            jacobian = np.asarray(derivatives["jacobian"](*all_values), dtype=float)
            df_dk = np.asarray(derivatives["df_dk"](*all_values), dtype=float)
            dtheta_dk = -np.linalg.pinv(jacobian) @ df_dk
            dr_net_dk = np.asarray(derivatives["dr_net_dk"](*all_values), dtype=float)
            dr_net_dtheta = np.asarray(derivatives["dr_net_dtheta"](*all_values), dtype=float)
            dr_forward_dk = np.asarray(derivatives["dr_forward_dk"](*all_values), dtype=float)
            dr_forward_dtheta = np.asarray(derivatives["dr_forward_dtheta"](*all_values), dtype=float)
            dr_net_total = ((functions["stoichi_pairs"] @ dr_net_dk)
                            + (functions["stoichi_pairs"] @ dr_net_dtheta) @ dtheta_dk)
            dr_forward_total = ((functions["s_plus_pairs"] @ dr_forward_dk)
                                + (functions["s_plus_pairs"] @ dr_forward_dtheta) @ dtheta_dk)
            alpha = np.abs(r_vec) / (np.abs(r_vec) + np.sum(np.abs(r_forward)) + eps)
            dr_total = alpha[:, None] * dr_net_total + (1 - alpha[:, None]) * dr_forward_total

            mask = np.abs(r_vec) > 1e-50
            if not np.any(mask):
                continue

            filtered_products = [p for i, p in enumerate(products) if mask[i]]
            r_vec_filtered = r_vec[mask]
            dr_total_filtered = dr_total[mask]
            filtered_products = [product for i, product in enumerate(products) if mask[i]]
            r_vec_filtered = r_vec[mask]
            dr_total_filtered = dr_total[mask]
            # Convert derivatives with respect to k into derivatives with respect to ln(k).
            scaled_derivatives = (dr_total_filtered * k_values[None, :])
            drc_process = scaled_derivatives / (r_vec_filtered[:, None] + eps)
            r_tot = np.sum(r_vec_filtered)
            if abs(r_tot) < 1.0e-20:
                dsc_process = np.zeros_like(drc_process)
            else:
                selectivity = r_vec_filtered / r_tot
                derivative_sum = np.sum(scaled_derivatives, axis=0,)
                dsc_process = ((scaled_derivatives - selectivity[:, None] * derivative_sum) / (
                        selectivity[:, None] * r_tot + eps))
            # Combine forward and backward process sensitivities into one value for each reversible elementary step.
            n_steps = len(functions["step_pairs"])
            drc = np.empty((len(filtered_products), n_steps), dtype=float,)
            dsc = np.empty_like(drc)
            for pair_index, (forward, backward) in enumerate(functions["step_pairs"]):
                forward_index = functions["process_index"][forward]
                backward_index = functions["process_index"][backward]
                drc[:, pair_index] = (drc_process[:, forward_index] + drc_process[:, backward_index])
                dsc[:, pair_index] = (dsc_process[:, forward_index] + dsc_process[:, backward_index])
            for i, product in enumerate(filtered_products):
                drc_data[product].rows.append([temp_num, *drc[i].tolist()])
                dsc_data[product].rows.append([temp_num, *dsc[i].tolist()])
        self.results.degree_of_rate_control = drc_data
        self.results.degree_of_selectivity_control = dsc_data

    def reaction_pathways(self, context, functions, flux_threshold: float = 1e-6):
        pathways: list[ReactionPathwayResult] = []
        for temp_num in [context["temperatures"][0], context["temperatures"][-1]]:
            concentration, _, all_values = self.steady_state_point(temp_num, context, functions)
            r_forward = np.asarray(functions["r_forward_func"](*all_values), dtype=float).flatten()
            r_backward = np.asarray(functions["r_backward_func"](*all_values), dtype=float).flatten()
            r_net = r_forward - r_backward
            edges: dict[tuple[str, str], float] = {}
            for i, (forward, backward) in enumerate(functions["step_pairs"]):
                flux = float(r_net[i])
                if abs(flux) < flux_threshold:
                    continue
                pid = forward if flux >= 0 else backward
                process = self.model.processes[pid]
                for reactant in process.reactants:
                    for product in process.products:
                        if reactant.species == product.species:
                            continue
                        key = (reactant.species, product.species)
                        edges[key] = edges.get(key, 0.0) + abs(flux)
            ordered_edges = [(a, b, w) for (a, b), w in sorted(edges.items(), key=lambda item: item[1], reverse=True)]
            pathways.append(ReactionPathwayResult(temperature=temp_num, edges=ordered_edges))
        self.results.reaction_pathways = pathways

    def ir_evolution(self, context, *, freq_min=500.0, freq_max=4000.0, points=3000, sigma=8.0):
        spectra: list[IRSpectrumResult] = []
        wavenumbers = np.linspace(freq_min, freq_max, points)
        temperatures = context["temperatures"]
        selected = [temperatures[0], temperatures[len(temperatures) // 2], temperatures[-1]]
        selected = sorted(set(selected))

        gases = [context["dynamic_species"][i] for i in context["gas_indices"]]
        adsorbates = [name for name in context["dynamic_species"] if is_adsorbate(self.model.species[name])]

        def gaussian(x, centre, amplitude):
            return amplitude * np.exp(-((x - centre) ** 2) / (2 * sigma**2))

        for temp_num in selected:
            times = np.asarray(context["times"], dtype=float)
            for phase, names in {"gas": gases, "adsorbates": adsorbates}.items():
                intensity = np.zeros((len(times), len(wavenumbers)), dtype=float)
                for i, time_value in enumerate(times):
                    concentration = context["concentrations"][temp_num][float(time_value)]
                    concentration_map = {name: concentration[j] for j, name in enumerate(context["dynamic_species"])}
                    for name in names:
                        species = self.model.species[name]
                        freqs = getattr(species, "frequencies", []) or []
                        ir = getattr(species, "ir_intensities", []) or []
                        if not freqs or not ir:
                            continue
                        amount = concentration_map.get(name, 0.0)
                        for nu, amp in zip(freqs, ir):
                            if freq_min <= nu <= freq_max:
                                intensity[i] += amount * gaussian(wavenumbers, float(nu), float(amp))
                spectra.append(IRSpectrumResult(phase=phase, temperature=temp_num, wavenumbers=wavenumbers, times=times,
                                                intensities=intensity, species=list(names)))
        self.results.ir_spectra = spectra

    # --- TPR analysis
    def analyse_tpr(self):
        experiments = getattr(self.model, "experiments", None)
        tpr_results = getattr(experiments, "tpr", None)
        if not tpr_results:
            return None

        peaks: list[TPRPeakResult] = []
        for result in tpr_results:
            temperatures = np.asarray(result.temperatures, dtype=float)
            for signal_gas, rates in result.spectra.items():
                rates = np.asarray(rates, dtype=float)
                if rates.size == 0 or np.max(rates) <= 0.0:
                    continue
                peak_indices, _ = find_peaks(rates, height=0.05 * np.max(rates))
                if len(peak_indices) == 0:
                    peak_indices = [int(np.argmax(rates))]
                area = float(np.trapz(rates, temperatures))
                for idx in peak_indices:
                    peaks.append(TPRPeakResult(source_gas=result.gas, signal_gas=signal_gas,
                                               heating_rate=float(result.heating_rate),
                                               peak_temperature=float(temperatures[idx]), peak_rate=float(rates[idx]),
                                               area=area,))
        self.results.tpr_peaks = peaks
        return peaks

