"""
Standalone model reporting and plotting utilities for the microkinetic workflow.
This script exports and visualises the results stored in model.
"""

from __future__ import annotations

import argparse
import itertools
from dataclasses import is_dataclass, fields
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import sympy as sp

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import SymmetricalLogLocator, MaxNLocator
from matplotlib import cm, colors
#plt.rcParams["path.simplify"] = False
from scipy.signal import find_peaks
from scipy.interpolate import PchipInterpolator

from Symbols_def import vext, ph, temp, constants, chem_label



icolour = ["b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki", "teal"]
iliner = ['-', '--', '-.', ':', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1))]
ipatterns = ["///", "...", "xx", "**", "\\", "|", "--", "++", "oo", "OO"]


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def safe_float(value: Any) -> float:
    try:
        value = complex(value)
        if abs(value.imag) > 1e-12:
            return float("nan")
        return float(value.real)
    except Exception:
        try:
            return float(value)
        except Exception:
            return float("nan")

def is_numeric_expression(value: Any) -> bool:
    return isinstance(value, (int, float, np.number, sp.Expr))

def evaluate_expression(expr: Any, temperature: float, ph_value: float | None = None,
                        potential_value: float | None = None) -> float:
    """Evaluate a scalar expression at the supplied conditions."""
    if expr is None:
        return float("nan")
    if isinstance(expr, (int, float, np.number)):
        return float(expr)
    if not isinstance(expr, sp.Expr):
        return float("nan")
    substitutions = dict(constants) if isinstance(constants, dict) else {}
    substitutions[temp] = temperature
    if ph_value is not None and np.isfinite(ph_value):
        substitutions[ph] = ph_value
    if potential_value is not None and np.isfinite(potential_value):
        substitutions[vext] = potential_value
    try:
        evaluated = expr.subs(substitutions)
        return safe_float(evaluated.evalf())
    except Exception:
        return float("nan")

def scan_values(scan: Any, default: list[float]) -> list[float]:
    if scan is None:
        return default
    if hasattr(scan, "values"):
        try:
            return [float(x) for x in scan.values()]
        except Exception:
            return default
    if isinstance(scan, (int, float)):
        return [float(scan)]
    if isinstance(scan, (list, tuple)):
        if len(scan) == 1:
            return [float(scan[0])]
        if len(scan) == 3:
            start, stop, step = [float(x) for x in scan]
            return list(np.arange(start, stop + 0.5 * step, step))
    return default

def condition_grid(model: Any) -> list[dict[str, float]]:
    conditions = getattr(model, "conditions", None)
    temperatures = scan_values(getattr(conditions, "temperature", None), default=[300.0])
    ph_values = scan_values(getattr(conditions, "ph", None), default=[float("nan")])
    potential_values = scan_values(getattr(conditions, "potential", None), default=[float("nan")])
    return [{"Potential[V]": float(vext), "pH": float(ph), "Temperature[K]": float(temp_num),}
        for temp_num, ph, vext in itertools.product(temperatures, ph_values, potential_values)]

def object_numeric_fields(obj: Any) -> dict[str, Any]:
    """Return scalar numeric/symbolic fields from an arbitrary object/dataclass."""
    if obj is None:
        return {}
    out: dict[str, Any] = {}
    if is_dataclass(obj):
        for f in fields(obj):
            value = getattr(obj, f.name, None)
            if is_numeric_expression(value):
                out[f.name] = value
        return out
    for name in dir(obj):
        if name.startswith("_"):
            continue
        try:
            value = getattr(obj, name)
        except Exception:
            continue
        if callable(value):
            continue
        if is_numeric_expression(value):
            out[name] = value
    return out

def finite(values: Iterable[float]) -> list[float]:
    return [float(v) for v in values if isinstance(v, (int, float, np.number)) and np.isfinite(v)]

def setup_signed_log_axis(ax, *value_arrays, padding=1.5, origin_gap=1.5,):
    """ Transform signed values into logarithmic-decade coordinates and configure the corresponding y-axis. """
    if not value_arrays:
        raise ValueError("At least one value array must be provided.")
    arrays = tuple(np.asarray(values, dtype=float) for values in value_arrays)
    all_values = np.concatenate([values.ravel() for values in arrays])
    finite_values = all_values[np.isfinite(all_values)]
    # Handle arrays containing only zeros, NaN, or infinity.
    if finite_values.size == 0:
        transformed = tuple(np.full_like(values, np.nan, dtype=float) for values in arrays)
        ax.set_yticks([0.0])
        ax.set_yticklabels([r"$0$"])
        ax.set_ylim(-1.0, 1.0)
        return transformed
    positive_values = finite_values[finite_values > 0.0]
    negative_values = finite_values[finite_values < 0.0]
    has_positive = positive_values.size > 0
    has_negative = negative_values.size > 0
    nonzero_values = finite_values[finite_values != 0.0]
    if nonzero_values.size == 0:     # All finite values are zero.
        transformed = tuple(np.where(np.isfinite(values), 0.0, np.nan) for values in arrays)
        ax.set_yticks([0.0])
        ax.set_yticklabels([r"$0$"])
        ax.set_ylim(-0.5, 0.5)
        return transformed
    magnitudes = np.abs(nonzero_values)
    min_exponent = int(np.floor(np.log10(np.min(magnitudes))))
    max_exponent = int(np.ceil(np.log10(np.max(magnitudes))))

    def transform(values):
        values = np.asarray(values, dtype=float)
        output = np.zeros_like(values, dtype=float)
        finite = np.isfinite(values)
        nonzero = finite & (values != 0.0)
        output[nonzero] = (np.sign(values[nonzero]) * (np.log10(np.abs(values[nonzero])) - min_exponent + origin_gap))
        output[~finite] = np.nan
        return output

    transformed_arrays = tuple(transform(values) for values in arrays)
    # Estimate how many exponent labels fit on the axis.
    available_tick_space = max(3, int(ax.yaxis.get_tick_space()),)
    # Mixed-sign plots share the available space between both sides.
    if has_positive and has_negative:
        target_ticks_per_side = max(3, (available_tick_space - 2) // 2,)
    else:
        target_ticks_per_side = max(3, available_tick_space - 2,)
    n_decades = max_exponent - min_exponent + 1
    # Use one constant exponent interval, so all displayed ticks are equidistant.
    decade_step = max(1, int(np.ceil(n_decades / target_ticks_per_side)),)

    exponents = np.arange(min_exponent, max_exponent + 1, decade_step, dtype=int,)
    # Use the same coordinate transform as the plotted values.
    positive_positions = (exponents - min_exponent + origin_gap).astype(float)
    # Remove the smallest decade adjacent to zero when both positive and negative values are present.
    if (has_positive and has_negative and origin_gap < 2.0 and len(exponents) > 1):
        exponents_to_plot = exponents[1:]
        positions_to_plot = positive_positions[1:]
    else:
        exponents_to_plot = exponents
        positions_to_plot = positive_positions

    tick_positions = []
    tick_labels = []
    if has_negative:
        tick_positions.extend((-positions_to_plot[::-1]).tolist())
        tick_labels.extend([rf"$-10^{{{exponent}}}$" for exponent in exponents_to_plot[::-1]])
    if has_positive:
        tick_positions.extend(positions_to_plot.tolist())
        tick_labels.extend([rf"$10^{{{exponent}}}$" for exponent in exponents_to_plot])
    # Include zero only when zero is relevant as the baseline.
    if has_negative and has_positive:
        tick_positions.append(0.0)
        tick_labels.append(r"$0$")

    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels)
    finite_transformed = np.concatenate([values[np.isfinite(values)] for values in transformed_arrays])
    transformed_min = float(np.min(finite_transformed))
    transformed_max = float(np.max(finite_transformed))

    if has_positive and not has_negative:
        # Positive-only plot: no unnecessary negative region.
        ax.set_ylim(0.0, transformed_max + padding,)
    elif has_negative and not has_positive:
        # Negative-only plot: no unnecessary positive region.
        ax.set_ylim(transformed_min - padding, 0.0,)
    else:
        # Mixed-sign plot: each side follows its actual data extent.
        ax.set_ylim(transformed_min - padding, transformed_max + padding,)
    return transformed_arrays

def plot_initial_final_bars(*, labels, initial_values, final_values, initial_temperature, final_temperature, output,
                            title, ylabel, xlabel="Process ID", logy=False,):
    """ Plot paired initial/final bar values.
        When logy=True, signed values are displayed using logarithmic-decade coordinates """
    labels = list(labels)
    if not labels:
        return
    # Ensure the arrays correspond directly to one bar per label.
    initial_values = np.asarray(initial_values, dtype=float).reshape(-1)
    final_values = np.asarray(final_values, dtype=float).reshape(-1)
    expected_size = len(labels)
    if not (expected_size == initial_values.size == final_values.size):
        raise ValueError("Labels and initial/final values must have equal lengths."
                         f"The {title} plot got labels={len(labels)}, initial_values={initial_values.size}, "
                         f"final_values={final_values.size}.")

    formatted_labels = [rf"{label}" for label in labels]
    x = np.arange(expected_size, dtype=float)
    width = 0.4
    fig, ax = plt.subplots(figsize=(max(10.0, 0.45 * expected_size), 6.0), clear=True,)
    if logy:
        (initial_plot_values, final_plot_values,) = setup_signed_log_axis(ax, initial_values, final_values,)
    else:
        initial_plot_values = initial_values
        final_plot_values = final_values
    ax.bar(x - width / 2.0, initial_plot_values, width, label=f"{initial_temperature:g} K", hatch=ipatterns[0],
           edgecolor="black",)
    ax.bar(x + width / 2.0, final_plot_values, width, label=f"{final_temperature:g} K", hatch=ipatterns[1],
           edgecolor="black",)
    ax.axhline(0.0, linewidth=0.8, color="black",)
    ax.set_xticks(x)
    ax.set_xticklabels(formatted_labels, rotation=45, ha="right", rotation_mode="anchor",)
    ax.tick_params(axis="x", labelsize=16, pad=2.5,)
    ax.tick_params(axis="y", labelsize=16,)
    ax.set_xlabel(xlabel, fontsize=18,)
    ax.set_ylabel(ylabel, fontsize=18,)
    ax.set_title(title, fontsize=18,)
    if not logy:
        ax.margins(y=0.15)
    if expected_size:
        ax.set_xlim(-0.75, expected_size - 0.25, )
    ax.legend(fontsize=16)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True,)
    fig.savefig(output, dpi=300, transparent=True, bbox_inches="tight",)
    plt.close(fig)


class ThermodynamicsReporter:
    def __init__(self, model, output_dir="./"):
        self.model = model
        self.output_dir = Path(output_dir)
        self.root = self.output_dir / "THERMODYNAMICS"

    def write_and_plot(self):
        for species in self.model.species.values():
            if not hasattr(species, "thermo") or species.thermo is None:
                continue
            species_dir = self.species_dir(species)
            self.write_partition_functions(species, species_dir)
            self.write_entropy(species, species_dir)
            self.write_enthalpy(species, species_dir)
            self.write_heat_capacity(species, species_dir)
            self.write_gibbs(species, species_dir)
            self.write_ir_spectrum(species, species_dir)
            self.plot_partition_functions(species, species_dir)
            self.plot_entropy(species, species_dir)
            self.plot_enthalpy(species, species_dir)
            self.plot_heat_capacity(species, species_dir)
            self.plot_gibbs(species, species_dir)
            self.plot_ir_spectrum(species, species_dir)

    def species_dir(self, species):
        path = self.root / species.name
        path.mkdir(parents=True, exist_ok=True)
        return path

    def temperatures(self):
        if self.model.conditions.temperature is None:
            return [300.0]
        return self.model.conditions.temperature.values()

    def evaluate(self, expr, temperature):
        if expr is None:
            return None
        try:
            expr = sp.sympify(expr).subs(constants)
            value = expr.subs(temp, temperature)
            return float(value)
        except Exception:
            try:
                return float(expr)
            except Exception:
                return None

    def table_rows(self, species, fields):
        rows = []
        thermo = species.thermo
        for temperature in self.temperatures():
            row = [temperature]
            for field in fields:
                expr = getattr(thermo, field, None)
                row.append(self.evaluate(expr, temperature))
            rows.append(row)
        return rows

    def write_table(self, filename, header, rows, title=None):
        Writer.write(filename=filename, header=header, rows=rows, title=title,)

    def write_partition_functions(self, species, species_dir):
        fields = ["qtrans3d", "qrot", "qelec", "qvib3d", "q3d", "qtrans2d", "qvib2d", "q2d"]
        self.write_table(species_dir / "Partition_functions.dat", ["Temperature(K)", *fields],
                         self.table_rows(species, fields),  title=f"Partition functions for {species.name}",)

    def write_entropy(self, species, species_dir):
        fields = ["strans3d", "srot", "selec", "svib3d", "sentropy3d", "strans2d", "svib2d", "sentropy2d",]
        self.write_table(species_dir / "Entropy.dat", ["Temperature(K)", *fields], self.table_rows(species, fields),
                         title=f"Entropy contributions for {species.name}",)

    def write_enthalpy(self, species, species_dir):
        fields = ["zpe3d", "htrans3d", "hvib3d", "hrot", "hthermal3d", "enthalpy3d",
                  "zpe2d", "htrans2d", "hvib2d", "hthermal2d", "enthalpy2d",]
        self.write_table(species_dir / "Enthalpy.dat", ["Temperature(K)", *fields], self.table_rows(species, fields),
                         title=f"Enthalpy contributions for {species.name}",)

    def write_heat_capacity(self, species, species_dir):
        fields = ["cp3d", "cp2d",]
        self.write_table(species_dir / "Heat_capacity.dat", ["Temperature(K)", *fields],
                         self.table_rows(species, fields), title=f"Heat capacity for {species.name}",)

    def write_gibbs(self, species, species_dir):
        fields = ["gibbs3d", "gibbs2d",]
        self.write_table(species_dir / "Gibbs.dat", ["Temperature(K)", *fields], self.table_rows(species, fields),
                         title=f"Gibbs free energy for {species.name}",)

    def write_ir_spectrum(self, species, species_dir):
        """Write calculated vibrational frequencies and IR intensities."""
        frequencies = getattr(species, "frequencies", []) or []
        intensities = getattr(species, "ir_intensities", []) or []
        if not frequencies or not intensities:
            return
        n_modes = min(len(frequencies), len(intensities),)
        frequencies = np.asarray(frequencies[:n_modes], dtype=float,)
        intensities = np.asarray(intensities[:n_modes], dtype=float,)
        # Remove invalid entries. Keep all valid vibrational modes;
        # the plotting routine can apply its own frequency window.
        mask = (np.isfinite(frequencies) & np.isfinite(intensities))
        frequencies = frequencies[mask]
        intensities = intensities[mask]
        if frequencies.size == 0:
            return
        output_file = species_dir / "IR_spectrum.dat"
        with open(output_file, "w") as output:
            output.write("# Wavenumber(cm^-1)\tIR_intensity(a.u.)\n")
            for frequency, intensity in zip(frequencies, intensities,):
                output.write(f"{frequency:15.3f}\t{intensity:15.3e}\n")

    def plot_fields(self, species, species_dir, fields, filename, ylabel):
        temperatures = np.asarray(self.temperatures(), dtype=float)
        fig, ax = plt.subplots(figsize=(8, 6), clear=True)
        for n,field in enumerate(fields):
            values = [self.evaluate(getattr(species.thermo, field, None), temp_num) for temp_num in temperatures]
            if all(value is None for value in values):
                continue
            values = np.asarray([np.nan if value is None else value for value in values], dtype=float,)
            ax.plot(temperatures, values, color=icolour[n % len(icolour)], linestyle=iliner[n % len(iliner)], label=field,)
        ax.tick_params(labelsize=16,)
        ax.set_xlabel('Temperature $(K)$', fontsize=18,)
        ax.set_ylabel(f'{ylabel}', fontsize=18, )
        ax.set_title(species.name, fontsize=18, )
        #ax.legend()
        fig.tight_layout()
        fig.savefig(species_dir / filename, dpi=300, transparent=True)
        plt.close(fig)

    def plot_partition_functions(self, species, species_dir):
        self.plot_fields(species, species_dir, ["q3d",], "Partition_functions.svg", "Partition function, $Q$",)

    def plot_entropy(self, species, species_dir):
        self.plot_fields(species, species_dir,["sentropy3d",], "Entropy.svg", "Entropy, $S$ $(eV/K)$",)

    def plot_enthalpy(self, species, species_dir):
        self.plot_fields(species, species_dir, ["enthalpy3d",], "Enthalpy.svg", "Enthalpy, $H$ $(eV)$",)

    def plot_heat_capacity(self, species, species_dir):
        self.plot_fields(species, species_dir, ["cp3d",], "Heat_capacity.svg", "Heat capacity, $C_{p}$ $(eV/K)$",)

    def plot_gibbs(self, species, species_dir):
        self.plot_fields(species, species_dir, ["gibbs3d"], "Gibbs.svg", "Gibbs free energy, $G$ $(eV)$",)

    def plot_ir_spectrum(self, species, species_dir, *, freq_min=500.0, freq_max=4500.0, points=5000, sigma=8.0,):
        """Plot the intrinsic IR spectrum of one species."""
        frequencies = getattr(species, "frequencies", []) or []
        intensities = getattr(species, "ir_intensities", []) or []
        if not frequencies or not intensities:
            return
        n_modes = min(len(frequencies), len(intensities),)
        frequencies = np.asarray(frequencies[:n_modes], dtype=float,)
        intensities = np.asarray(intensities[:n_modes], dtype=float,)
        # Keep only physical vibrational frequencies in the requested spectral window.
        mask = (np.isfinite(frequencies) & np.isfinite(intensities) & (frequencies >= freq_min) &
                (frequencies <= freq_max))
        frequencies = frequencies[mask]
        intensities = intensities[mask]
        if frequencies.size == 0:
            return
        wavenumbers = np.linspace(freq_min, freq_max, points,)
        spectrum = np.zeros_like(wavenumbers, dtype=float,)
        for frequency, intensity in zip(frequencies, intensities,):
            spectrum += intensity * np.exp(-((wavenumbers - frequency) ** 2 / (2.0 * sigma ** 2)))
        fig, ax = plt.subplots(figsize=(10, 6), clear=True,)
        ax.plot(wavenumbers, spectrum, linewidth=1.5,)
        # Optional sticks showing the actual calculated normal modes.
        #ax.vlines(frequencies, 0.0, intensities, linewidth=0.8, alpha=0.5,)
        # Conventional IR presentation: high wavenumber on the left.
        ax.set_xlim(freq_max, freq_min,)
        ax.set_xlabel(r"Wavenumber ($cm^{-1}$)", fontsize=18,)
        ax.set_ylabel("IR intensity ($a.u.$)", fontsize=18,)
        ax.set_title(f"IR spectrum of {chem_label(species.name)}", fontsize=18,)
        ax.tick_params(axis="both", labelsize=16,)
        ax.margins(y=0.10)
        fig.tight_layout()
        fig.savefig(species_dir / "IR_spectrum.svg", dpi=300, transparent=True, bbox_inches="tight",)
        plt.close(fig)

class KineticsReporter:
    def __init__(self, model, output_dir="./"):
        self.model = model
        self.output_dir = Path(output_dir)
        self.root = self.output_dir / "KINETICS"
        self._expression_cache = {}
        self._array_cache = {}

    def write_and_plot(self):
        self.root.mkdir(parents=True, exist_ok=True)
        self._array_cache = {}
        self.write_rate_constants()
        self.write_reaction_energies()
        self.write_symbolic_rate_equations()
        self.write_overall_stoichiometry()
        self.write_profile()
        self.plot_rate_constants()
        self.plot_activation_energies()
        self.plot_reaction_energies()

        self._array_cache.clear()

    def temperatures(self):
        if self.model.conditions.temperature is None:
            return [300.0]
        return self.model.conditions.temperature.values()

    def _compile_expression(self, expr):
        """ Compile a symbolic expression into a NumPy-compatible function.
        Each expression is compiled only once and stored in the cache. """
        if expr is None:
            return None
        key = id(expr)
        if key in self._expression_cache:
            return self._expression_cache[key]
        try:
            symbolic_expr = sp.sympify(expr).subs(constants)
            if temp in symbolic_expr.free_symbols:
                function = sp.lambdify(temp, symbolic_expr, modules="numpy",)
            else:
                constant_value = float(symbolic_expr)

                def function(temperatures, value=constant_value):
                    temperatures = np.asarray(temperatures, dtype=float,)
                    return np.full(temperatures.shape, value, dtype=float,)
            self._expression_cache[key] = function
            return function
        except (TypeError, ValueError, OverflowError, ZeroDivisionError, sp.SympifyError,):
            self._expression_cache[key] = None
            return None

    def evaluate_array(self, expr, temperatures):
        """ Evaluate an expression over a complete temperature array. """
        temperatures = np.asarray(temperatures, dtype=float,)
        if expr is None:
            return np.full(temperatures.shape, np.nan, dtype=float,)
        cache_key = (id(expr), temperatures.shape, temperatures.tobytes(),)
        if cache_key in self._array_cache:
            return self._array_cache[cache_key]
        function = self._compile_expression(expr)
        if function is None:
            values = np.full(temperatures.shape, np.nan, dtype=float,)
            self._array_cache[cache_key] = values
            return values
        try:
            with np.errstate(over="ignore", divide="ignore", invalid="ignore", under="ignore",):
                values = function(temperatures)
            values = np.asarray(values, dtype=float,)
            # Some lambdified expressions return one scalar even when an array of temperatures is supplied.
            if values.ndim == 0:
                values = np.full(temperatures.shape, float(values), dtype=float,)
            else:           # Ensure that the returned array has the expected shape.
                values = np.broadcast_to(values, temperatures.shape,).astype(float, copy=True,)
            values[~np.isfinite(values)] = np.nan
            return values
        except (TypeError, ValueError, OverflowError, ZeroDivisionError, sp.SympifyError,):
            values = np.full(temperatures.shape, np.nan, dtype=float,)
        self._array_cache[cache_key] = values
        return values

    def evaluate_kinetic_fields(self, fields):
        temperatures = np.asarray(self.temperatures(), dtype=float,)
        process_items = list(self.model.kinetics.by_process.items())
        evaluated = {pid: {field: self.evaluate_array(getattr(kinetics, field, None), temperatures,)
                           for field in fields} for pid, kinetics in process_items}
        return temperatures, process_items, evaluated

    def kinetic_rows(self, fields, temperatures=None, process_items=None, evaluated=None,):
        if temperatures is None or process_items is None or evaluated is None:
            temperatures, process_items, evaluated = (self.evaluate_kinetic_fields(fields))
        rows = []
        for index, temp_num in enumerate(temperatures):
            for pid, kinetics in process_items:
                process = self.model.processes[pid]
                rows.append([temp_num, pid, process.kind, *[evaluated[pid][field][index] for field in fields],])
        return rows

    def dynamic_species_names(self) -> list[str]:
        """   Return dynamic species in exactly the same order used by Kinetics.REquations.overall_stoichiometry().  """
        names: set[str] = set()
        for process in self.model.processes.values():
            names.update(item.species for item in process.reactants)
            names.update(item.species for item in process.products)
        return sorted(name for name in names if str(self.model.species[name].kind).lower() != "surface")

    @staticmethod
    def format_stoichiometric_number(value: sp.Expr | int | float,) -> str:
        """  Format an exact stoichiometric coefficient without unnecessary decimal notation.  """
        value = sp.nsimplify(value)
        if value == 1:
            return ""
        if value.is_Integer:
            return f"{int(value)} "
        return f"{value} "

    def write_overall_stoichiometry(self) -> Path | None:
        """  Write the overall gas-phase reaction stored in model.equations.overall_stoichiometry.
            Negative coefficients are reactants and positive coefficients are products.
            Adsorbates and bare surfaces must cancel from a valid overall reaction.     """
        equations = getattr(self.model, "equations", None,)
        if equations is None:
            print("WARNING: model.equations is unavailable; Overall_Reaction.dat was not written.", flush=True,)
            return None
        overall = list(getattr(equations, "overall_stoichiometry", [],))
        species_names = self.dynamic_species_names()
        if not overall:
            print("WARNING: no overall stoichiometry is available. Run REquations(model) before KineticsReporter.",
                  flush=True,)
            return None
        if len(overall) != len(species_names):
            raise ValueError(f"Overall-stoichiometry/species alignment failure: {len(overall)} coefficients but "
                             f"{len(species_names)} dynamic species.")
        coefficients = [sp.nsimplify(value) for value in overall]
        adsorbed_residuals: list[tuple[str, sp.Expr]] = []
        molecular_coefficients: list[tuple[str, sp.Expr]] = []
        for name, coefficient in zip(species_names, coefficients,):
            if coefficient == 0:
                continue
            species_kind = str(self.model.species[name].kind).lower()
            if species_kind == "molecule":
                molecular_coefficients.append((name, coefficient,))
            else:
                adsorbed_residuals.append((name, coefficient,))
        if adsorbed_residuals:
            residual_text = ", ".join(f"{name}={coefficient}" for name, coefficient in adsorbed_residuals)
            raise ValueError("The calculated overall reaction does not eliminate all adsorbed intermediates: "
                             f"{residual_text}.")
        if not molecular_coefficients:
            print("WARNING: the stored null-space solution produces no non-zero molecular reaction. The selected"
                  "null-space vector probably represents cancellation of a forward and reverse elementary process.",
                  flush=True,)
            return None
        reactants = [(name, -coefficient,) for name, coefficient in molecular_coefficients if coefficient < 0]
        products = [(name, coefficient,) for name, coefficient in molecular_coefficients if coefficient > 0]
        if not reactants or not products:
            raise ValueError("The overall stoichiometric vector does not contain both molecular reactants and "
                             f"molecular products: {molecular_coefficients}.")
        reactant_text = " + ".join(self.format_stoichiometric_number(coefficient) + name
                                   for name, coefficient in reactants)
        product_text = " + ".join(self.format_stoichiometric_number(coefficient) + name
                                  for name, coefficient in products)
        reaction_text = (f"{reactant_text} --> {product_text}")
        directory = self.root
        directory.mkdir(parents=True, exist_ok=True,)
        filename = (directory / "Overall_Reaction.dat")
        with open(filename, "w", encoding="utf-8",) as handle:
            handle.write("# Overall stoichiometry of the kinetic model\n")
            handle.write("# Negative coefficients: reactants\n# Positive coefficients: products\n")
            handle.write("# Adsorbed intermediates and surface sites must cancel.\n\n")
            handle.write(reaction_text + "\n\n")
            handle.write("# Complete dynamic-species vector\n# Species  Coefficient\n")
            for name, coefficient in zip(species_names, coefficients,):
                handle.write(f"{name:<20s} {str(coefficient):>12s}\n")
        #print(f"  Written overall reaction: {reaction_text}", flush=True,)
        return filename

    def write_profile(self) -> Path | None:
        """  Write the free-energy profile generated by Kinetics.Profile. """
        profile_results = getattr(self.model, "profile", None,)
        if profile_results is None:
            print("WARNING: model.profile is unavailable; KINETICS/Profile.txt was not written.", flush=True,)
            return None
        data = getattr(profile_results, "data", None,)
        if not data:
            print("WARNING: the kinetic profile is empty. Run Profile(model, temp_num=...) before KineticsReporter.",
                  flush=True,)
            return None
        if len(data) < 2:
            print("WARNING: the kinetic profile contains a header but no process rows.", flush=True,)
        header = [str(value) for value in data[0]]
        rows = [[str(value) for value in row] for row in data[1:]]
        expected_columns = 6
        if len(header) != expected_columns:
            raise ValueError(f"Unexpected profile-table structure: expected {expected_columns} columns, "
                             f"but the header contains {len(header)}.")
        for row_number, row in enumerate(rows, start=1,):
            if len(row) != expected_columns:
                raise ValueError(f"Malformed kinetic-profile row {row_number}: expected {expected_columns} "
                                 f"columns, got {len(row)}.")
        output_file = self.root / "Profile.txt"
        Writer.write(output_file, header, rows, title=("Free-energy profile of the elementary kinetic processes"),)
        #print(f"  Written kinetic profile: {output_file}", flush=True,)
        return output_file

    def write_rate_constants(self):
        fields = ["sticky", "arrhenius", "ktunneling", "krate0",]
        Writer.write(filename=self.root / "Rate_constants.dat", header=["Temperature(K)", "Process", "Kind", *fields,],
                     rows=self.kinetic_rows(fields), title="Rate constant contributions by process",)

    def write_reaction_energies(self):
        fields = ["activation", "reaction_energy",]
        Writer.write(filename=self.root / "Reaction_energies.dat",
                     header=["Temperature(K)", "Process", "Kind", "activation_energy(eV)", "reaction_energy(eV)",],
                     rows=self.kinetic_rows(fields), title="Reaction free energies by process",)

    @staticmethod
    def _coefficient_for(items, species_name: str) -> float:
        """ Return the total stoichiometric coefficient of one species. """
        return sum(float(item.coefficient) for item in items if item.species == species_name)
    def _dynamic_species_names(self) -> list[str]:
        """ Return all integrated species. Surfaces are excluded because their coverages are reconstructed from
        the site balance rather than integrated as independent variables."""
        names: set[str] = set()
        for process in self.model.processes.values():
            names.update(item.species for item in process.reactants)
            names.update(item.species for item in process.products)
        return sorted(name for name in names if str(self.model.species[name].kind).lower() != "surface")

    def _species_rate_symbol(self, species_name: str) -> str:
        """ Use P_species for gas-phase molecules and the species name for adsorbates and surfaces. """
        species = self.model.species[species_name]
        kind = str(species.kind).lower()
        if kind == "molecule":
            return f"P_{species_name}"
        return species_name

    @staticmethod
    def _format_number(value: float) -> str:
        """ Format an integer-like stoichiometric coefficient without a decimal. """
        if np.isclose(value, round(value)):
            return str(int(round(value)))
        return f"{value:g}"

    def _format_reactant_factor(self, item) -> str:
        """ Format one reactant using its stoichiometric coefficient as the elementary reaction order. """
        symbol = self._species_rate_symbol(item.species)
        order = float(item.coefficient)
        if np.isclose(order, 1.0):
            return symbol
        return f"{symbol}^{self._format_number(order)}"

    def _format_elementary_rate(self, process_id, process) -> str:
        """ Construct: k_i reactant_1 reactant_2^order ... """
        factors = [f"k_{process_id}"]
        factors.extend(self._format_reactant_factor(item) for item in process.reactants)
        return " ".join(factors)

    def _symbolic_rhs_contributions(self, species_name: str,) -> list[tuple[float, str]]:
        """ Return every elementary contribution to d[species]/dt. """
        contributions: list[tuple[float, str]] = []
        for process_id, process in self.model.processes.items():
            product_coefficient = self._coefficient_for(process.products, species_name,)
            reactant_coefficient = self._coefficient_for(process.reactants, species_name,)
            net_factor = product_coefficient - reactant_coefficient
            if np.isclose(net_factor, 0.0):
                continue
            rate_expression = self._format_elementary_rate(process_id, process,)
            contributions.append((net_factor, rate_expression))
        return contributions

    def _format_species_rhs(self, species_name: str) -> str:
        """ Format a complete species rate equation. """
        contributions = self._symbolic_rhs_contributions(species_name)
        if not contributions:
            return f"d[{species_name}]/dt = 0"
        terms: list[str] = []
        for index, (net_factor, rate_expression) in enumerate(contributions):
            sign = "+" if net_factor > 0.0 else "-"
            magnitude = abs(net_factor)
            if np.isclose(magnitude, 1.0):
                coefficient = ""
            else:
                coefficient = f"{self._format_number(magnitude)} "
            term = f"{sign} {coefficient}{rate_expression}"
            if index == 0:
                terms.append(f"d[{species_name}]/dt = {term}")
            else:
                indentation = " " * (len(species_name) + 9)
                terms.append(f"{indentation}{term}")
        return "\n".join(terms)

    def symbolic_rate_equations_text(self) -> str:
        """ Return all symbolic species rate equations as plain text. """
        title = "Symbolic species rate equations"
        lines = [title,"=" * len(title),"",]
        for species_name in self._dynamic_species_names():
            lines.append(self._format_species_rhs(species_name))
            lines.append("")
        return "\n".join(lines).rstrip() + "\n"

    def write_symbolic_rate_equations(self) -> None:
        """ Write the symbolic RHS equations and also print them to the terminal. """
        filename = self.root / "Equations.dat"
        filename.write_text(self.symbolic_rate_equations_text(), encoding="utf-8")

    def plot_process_field(self, field, filename, ylabel, logy=False,):
        temperatures = np.asarray(self.temperatures(), dtype=float,)
        if temperatures.size == 0:
            return
        temp_initial = float(temperatures[0])
        temp_final = float(temperatures[-1])
        process_ids = list(self.model.kinetics.by_process.keys())
        if not process_ids:
            return
        process_ids = []
        initial_values = []
        final_values = []
        for pid, kinetics in self.model.kinetics.by_process.items():
            values = self.evaluate_array(getattr(kinetics, field, None), temperatures,)
            process_ids.append(pid)
            initial_values.append(values[0])
            final_values.append(values[-1])
        plot_initial_final_bars(labels=process_ids, initial_values=initial_values, final_values=final_values,
                                initial_temperature=temp_initial, final_temperature=temp_final,
                                output=self.root / f"{filename}.svg", title=f"{field}",
                                ylabel=ylabel, xlabel="Process ID", logy=logy,)

    def plot_rate_constants(self):
        self.plot_process_field(field="krate0", filename="Rate_constants",
                                ylabel="Rate constant, $K$ $(M^{\\dagger}/s)$", logy=True,)

    def plot_activation_energies(self):
        self.plot_process_field(field="activation", filename="Activation_energies",
                                ylabel="Activation energy, $G_{A}$ $(eV)$", logy=False,)

    def plot_reaction_energies(self):
        self.plot_process_field(field="reaction_energy", filename="Reaction_energies",
                                ylabel="Reaction free energy, $G_{R}$ $(eV)$", logy=False,)


class ExperimentReporter:
    def __init__(self, model, output_dir="./"):
        self.model = model
        self.output_dir = Path(output_dir)
        self.root = self.output_dir / "EXPERIMENTS"

    def write(self):
        self.root.mkdir(parents=True, exist_ok=True)
        if hasattr(self.model, "experiments"):
            if self.model.experiments.isothermal is not None:
                self.write_isothermal()
            if self.model.experiments.tpr:
                self.write_tpr()

    def write_isothermal(self):
        result = self.model.experiments.isothermal
        Writer.write(filename=self.root / "Isothermal.dat", header=result.headers, rows=result.values.tolist(),
                     title="Isothermal microkinetic simulation",)

    def write_tpr(self):
        tpr_dir = self.root / "TPR"
        tpr_dir.mkdir(parents=True, exist_ok=True)
        for result in self.model.experiments.tpr:
            label = (f"{result.gas}_{result.heating_rate:g}Kmin")
            gas_dir = tpr_dir / label
            gas_dir.mkdir(parents=True, exist_ok=True)
            self.write_tpr_concentrations(result, gas_dir / f"{label}_concentrations.dat",)
            self.write_tpr_spectra(result, gas_dir / f"{label}_spectra.dat",)
            #self.write_tpr_metadata(result, gas_dir / f"{label}_metadata.dat",)

    def write_tpr_concentrations(self, result, filename):
        rows = []
        for i, temperature in enumerate(result.temperatures):
            row = [temperature]
            for j in range(len(result.adsorbates)):
                row.append(result.concentrations[i, j])
            rows.append(row)
        Writer.write(filename=filename, header=["Temperature(K)", *result.adsorbates,], rows=rows,
                     title=(f"TPR adsorbate concentrations for {result.gas} at {result.heating_rate:g} K/min"),)

    def write_tpr_spectra(self, result, filename):
        rows = []
        gas_species = list(result.spectra.keys())
        for i, temperature in enumerate(result.temperatures):
            row = [temperature]
            for gas in gas_species:
                row.append(result.spectra[gas][i])
            rows.append(row)
        Writer.write(filename=filename, header=["Temperature(K)", *gas_species,], rows=rows,
                     title=(f"TPR gas production rates for {result.gas} at {result.heating_rate:g} K/min"),)

    def write_tpr_metadata(self, result, filename):
        rows = [["Gas", result.gas], ["Heating_rate(K/min)", result.heating_rate],
                ["Beta(K/s)", result.metadata.get("beta_K_per_s")],
                ["Initial_adsorbate", result.metadata.get("initial_adsorbate")],]
        Writer.write(filename=filename, header=["Property", "Value",], rows=rows,
                     title=f"TPR metadata for {result.gas}",)


class DiagnosticsReporter:
    def __init__(self, model, output_dir="./"):
        self.model = model
        self.output_dir = Path(output_dir)
        self.root = self.output_dir / "EXPERIMENTS"
        self.net_rates_root = self.root / "Rates"
        self.drc_root = self.root / "DRC"

    def write_and_plot(self):
        self.net_rates_root.mkdir(parents=True, exist_ok=True,)
        self.drc_root.mkdir(parents=True, exist_ok=True)
        diagnostics = getattr(self.model, "diagnostics", None)
        if diagnostics is None:
            print("WARNING: model.diagnostics is None.", flush=True,)
            return

        self.write_table_result(diagnostics.steady_state_net_rates, self.net_rates_root / "Steady_state_net_rates.dat",)
        self.write_table_result(diagnostics.average_net_rates, self.net_rates_root / "Average_net_rates.dat",)
        self.write_table_result(diagnostics.equilibrium_control, self.drc_root / "Equilibrium_control.dat",)
        self.write_table_result(diagnostics.drc_assessment, self.drc_root / "DRC_assessment.dat",)
        self.write_control_tables(diagnostics.degree_of_rate_control, self.drc_root / "Degree_of_rate_control",)
        self.write_control_tables(diagnostics.degree_of_selectivity_control, self.drc_root/"Degree_of_selectivity_control",)
        self.write_reaction_pathways()
        self.write_tpr_peaks()
        self.write_ir_spectra()

        self.plot_table_bars(diagnostics.steady_state_net_rates, self.net_rates_root / "Steady_state_net_rates.svg",
                             "Steady-state net rates",  "Net rate, $R_{N}^{SS}$ $(s^{-1})$", logy=True,)
        self.plot_table_bars(diagnostics.average_net_rates, self.net_rates_root / "Average_net_rates.svg",
                             "Average net rates", "Average net rate, $R_{N}^{Av}$ $(s^{-1})$", logy=True,)

        self.plot_table_bars(diagnostics.equilibrium_control, self.drc_root / "Equilibrium_control.svg",
                             "Equilibrium control", "DEC",)
        self.plot_control_tables(diagnostics.degree_of_rate_control, self.drc_root / "Degree_of_rate_control", ylabel="DRC")
        self.plot_control_tables(diagnostics.degree_of_selectivity_control, self.drc_root / "Degree_of_selectivity_control",
                                 ylabel="DSC",)
        self.plot_reaction_pathways()
        self.plot_reaction_pathway_networks()
        #self.plot_tpr_peaks()  # no plot an bar image with TPR peacks
        self.plot_tpr_spectra()
        self.plot_ir_spectra3d()
        self.plot_isothermal_species()

    def write_control_tables(self, tables, directory):  # WITH several tables
        if not tables:
            return
        directory.mkdir(parents=True, exist_ok=True,)
        for name, table in tables.items():
            if table is None:
                continue
            self.write_table_result(table, directory / f"{name}.dat",)

    def write_table_result(self, table, filename):
        if table is None:
            print(f"WARNING: table for {filename.name} is None; file not written.", flush=True,)
            return
        if not table.rows:
            print(f"WARNING: table for {filename.name} has no rows; file not written.", flush=True,)
            return
        Writer.write(filename=filename, header=table.headers, rows=table.rows, title=filename.stem.replace("_", " "),)

    def write_reaction_pathways(self):
        diagnostics = self.model.diagnostics
        directory = self.root / "Reaction_pathways"
        directory.mkdir(parents=True, exist_ok=True)
        rows = []
        for pathway in diagnostics.reaction_pathways:
            for reactant, product, flux in pathway.edges:
                rows.append([pathway.temperature, reactant, product, flux,])
        Writer.write(filename=directory / "Reaction_pathways.dat", header=["Temperature($K$)", "Reactant", "Product",
                                                                           "NetFlux($s^{-1}$) ",],
                     rows=rows, title="Flux-weighted reaction pathways",)

    def write_tpr_peaks(self):
        diagnostics = self.model.diagnostics
        directory = self.root / "TPR"
        directory.mkdir(parents=True, exist_ok=True)
        rows = []
        for peak in diagnostics.tpr_peaks:
            rows.append([peak.source_gas, peak.signal_gas, peak.heating_rate, peak.peak_temperature,
                         peak.peak_rate, peak.area,])
        Writer.write(filename=directory / "TPR_peaks.dat", header=["SourceGas", "SignalGas", "HeatingRate(K/min)",
                                                                   "PeakTemperature(K)", "PeakRate", "Area",],
                     rows=rows, title="TPR peak summary",)

    def write_ir_spectra(self):
        diagnostics = self.model.diagnostics
        directory = self.root / "IR_spectra"
        directory.mkdir(parents=True, exist_ok=True)
        for spectrum in diagnostics.ir_spectra:
            filename = (directory / f"{spectrum.phase}_{spectrum.temperature:g}K.dat")
            rows = []
            for i, time_value in enumerate(spectrum.times):
                for j, wavenumber in enumerate(spectrum.wavenumbers):
                    rows.append([spectrum.temperature, time_value, wavenumber, spectrum.intensities[i, j],])
            Writer.write(filename=filename, header=["Temperature(K)", "Time(s)", "Wavenumber(cm-1)", "Intensity",],
                         rows=rows, title=(f"IR evolution: {spectrum.phase}, {spectrum.temperature:g} K"),)

    def plot_table_bars(self,table, filename, title, ylabel, *, logy=False,):
        if table is None or not table.rows:
            return
        rows = sorted(table.rows, key=lambda row: float(row[0]),)
        initial_row = rows[0]
        final_row = rows[-1]
        initial_temperature = float(initial_row[0])
        final_temperature = float(final_row[0])
        labels = list(table.headers[1:])
        initial_values = [np.nan if value is None else float(value) for value in initial_row[1:]]
        final_values = [np.nan if value is None else float(value) for value in final_row[1:]]
        plot_initial_final_bars(labels=labels, initial_values=initial_values, final_values=final_values,
                                initial_temperature=initial_temperature, final_temperature=final_temperature,
                                output=filename, title=title, ylabel=ylabel, xlabel="Process ID", logy=logy,)

    def plot_control_tables(self, control_tables, directory, ylabel,):
        directory.mkdir(parents=True, exist_ok=True,)
        if not control_tables:
            print("No control tables available to plot.", flush=True,)
            return
        for name, table in control_tables.items():
            if table is None:
                print(f"Skipping empty control table: {table_name}", flush=True,)
                continue
            headers = getattr(table, "headers", None,)
            rows = getattr(table, "rows", None,)
            if headers is None:
                raise ValueError(f"Control table {name!r} has no 'headers' attribute.")
            if rows is None:
                raise ValueError(f"Control table {name!r} has no 'rows' attribute.")
            expected_columns = len(headers)
            for row_index, row in enumerate(rows):
                if len(row) != expected_columns:
                    raise ValueError(f"{name}: malformed control table. Headers={expected_columns}, "
                                     f"row {row_index}={len(row)}. Expected one value per header.")
            initial_row = rows[0]
            final_row = rows[-1]
            initial_temperature = float(initial_row[0])
            final_temperature = float(final_row[0])
            labels = list(table.headers[1:])
            initial_values = [np.nan if value is None else float(value) for value in initial_row[1:]]
            final_values = [np.nan if value is None else float(value) for value in final_row[1:]]
            plot_initial_final_bars(labels=labels, initial_values=initial_values, final_values=final_values,
                                    initial_temperature=initial_temperature, final_temperature=final_temperature,
                                    output=(directory / f"{name}_initial_final.svg"), title=f"{name}: {ylabel}",
                                    ylabel=ylabel, xlabel="Process ID", logy=False,)

    def plot_reaction_pathways(self):
        diagnostics = self.model.diagnostics
        directory = self.root / "Reaction_pathways"
        directory.mkdir(parents=True, exist_ok=True)
        for pathway in diagnostics.reaction_pathways:
            if not pathway.edges:
                continue
            labels = [(f"{chem_label(reactant)}$\\rightarrow$"
                       f"{chem_label(product)}") for reactant, product, _ in pathway.edges]
            values = [flux for _, _, flux in pathway.edges]
            fig, ax = plt.subplots(figsize=(10, 6), clear=True)
            y = np.arange(len(labels))
            ax.barh(y, values)
            ax.tick_params(axis="x", labelsize=16, pad=2.5,)
            ax.tick_params(axis="y",)
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=18,)
            ax.set_xlabel("Net Flux ($s^{-1}$)", fontsize=18,)
            ax.set_title(f"Reaction pathway at {pathway.temperature:g} K", fontsize=18,)
            fig.tight_layout()
            fig.savefig(directory / f"pathway_{pathway.temperature:g}K.svg", dpi=300, transparent=True,)
            plt.close(fig)

    def plot_reaction_pathway_networks(self):
        """Plot flux-weighted reaction pathways as directed NetworkX graphs."""
        try:
            import networkx as nx
        except ImportError:
            print("WARNING: networkx is not installed; skipping reaction-pathway network plots.", flush=True,)
            return
        diagnostics = self.model.diagnostics
        directory = self.root / "Reaction_pathways"
        directory.mkdir(parents=True, exist_ok=True)
        for pathway in diagnostics.reaction_pathways:
            if not pathway.edges:
                continue
            graph = nx.DiGraph()
            gas_species = {name for name, species in self.model.species.items()
                           if str(species.kind).lower() == "molecule"}
            adsorbates = {name for name, species in self.model.species.items()
                          if str(species.kind).lower() not in {"molecule", "surface",}}
            allowed_species = (gas_species | adsorbates)
            # pathway.edges already contains the net-flux direction determined by Diagnostics.reaction_pathways().
            for reactant, product, flux in pathway.edges:
                if (reactant not in allowed_species or product not in allowed_species):
                    continue
                flux = abs(float(flux))
                if graph.has_edge(reactant, product):
                    graph[reactant][product]["weight"] += flux
                else:
                    graph.add_edge(reactant, product, weight=flux,)
            if graph.number_of_edges() == 0:
                continue
            labels = {name: chem_label(name) for name in graph.nodes()}
            # Reproducible positioning.
            pos = nx.spring_layout(graph, seed=7, k=4,)
            edges = list(graph.edges(data=True))
            weights = np.asarray([data["weight"] for _, _, data in edges], dtype=float,)
            max_weight = np.max(weights)
            if max_weight <= 0.0:
                continue
            # Width shows relative net flux.
            #widths = [0.5 + 5.0 * weight / max_weight for weight in weights] # linear
            positive_weights = weights[weights > 0.0]
            min_weight = np.min(positive_weights)
            max_weight = np.max(positive_weights)
            if max_weight > min_weight: # log weight
                log_weights = np.log10(weights)
                widths = (0.5 + 5.0 * (log_weights - np.log10(min_weight))
                          / (np.log10(max_weight) - np.log10(min_weight)))
            else:
                widths = np.full_like(weights, 2.0,)

            fig, ax = plt.subplots(figsize=(12, 8), clear=True,)
            node_size = 5000
            nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=node_size, node_color="lightblue", alpha=0.75,
                                   edgecolors="grey", linewidths=1.0,)
            nx.draw_networkx_labels(graph, pos, labels=labels, ax=ax, font_size=16,)
            nx.draw_networkx_edges(graph, pos, ax=ax, width=widths, node_size=node_size, arrows=True,
                                   arrowstyle="->", arrowsize=20, connectionstyle="arc3, rad=0.05",
                                   min_source_margin=5, min_target_margin=5,)
            edge_labels = {(u, v): f"{data['weight']:.2e}" for u, v, data in edges}    # Numerical net flux on edges
            #nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, ax=ax, font_size=9, rotate=False,)
            ax.set_title(f"Reaction pathway at {pathway.temperature:g} K", fontsize=18,)
            ax.text(0.01, 0.01, r"Edge labels: Net Flux (s$^{-1}$)", transform=ax.transAxes, fontsize=14,)
            ax.axis("off")
            fig.tight_layout()
            fig.savefig(directory / f"pathway_network_{pathway.temperature:g}K.svg", dpi=300, transparent=True,
                        bbox_inches="tight",)
            plt.close(fig)

    def plot_tpr_peaks(self):
        diagnostics = self.model.diagnostics
        if not diagnostics.tpr_peaks:
            return
        labels = [f"${peak.source_gas}$ - ${peak.signal_gas}$" for peak in diagnostics.tpr_peaks]
        temperatures = [peak.peak_temperature for peak in diagnostics.tpr_peaks]
        fig, ax = plt.subplots(figsize=(10, 6), clear=True)
        x = np.arange(len(labels))
        ax.bar(x, temperatures)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor")
        ax.tick_params(axis="x", labelsize=16, pad=2.5,)
        ax.tick_params(axis="y", labelsize=16,)
        ax.set_ylabel("Peak temperature (K)", fontsize=18)
        ax.set_title("TPR peak temperatures", fontsize=18)
        ax.margins(y=0.15)
        ax.legend(fontsize=16)
        fig.tight_layout()
        fig.savefig(self.root / "TPR_peak_temperatures.svg", dpi=300, transparent=True, bbox_inches="tight",)
        plt.close(fig)

    def plot_ir_spectra2d(self):      # 2D
        diagnostics = self.model.diagnostics
        directory = self.root / "IR_spectra"
        directory.mkdir(parents=True, exist_ok=True)
        for spectrum in diagnostics.ir_spectra:
            if spectrum.intensities.size == 0:
                continue
            fig, ax = plt.subplots(figsize=(10, 6), clear=True)
            mesh = ax.pcolormesh(spectrum.wavenumbers, spectrum.times, spectrum.intensities, shading="auto",)
            fig.colorbar(mesh, ax=ax, label="Intensity (a.u.)",)
            ax.set_xlabel("Wavenumber (cm$^{-1}$)", fontsize=18)
            ax.set_ylabel("Time (s)", fontsize=18)
            ax.tick_params(axis="x", labelsize=16, pad=2.5,)
            ax.tick_params(axis="y", labelsize=16,)
            ax.set_title(f"IR evolution: {spectrum.phase}, {spectrum.temperature:g} K", fontsize=18)
            ax.margins(y=0.15)
            ax.legend(fontsize=16)
            fig.tight_layout()
            fig.savefig(directory / f"{spectrum.phase}_{spectrum.temperature:g}K.svg", dpi=300, transparent=True,)
            plt.close(fig)

    def plot_ir_spectra3d(self):
        """Plot time-resolved IR spectra as clean 3D waterfall plots."""
        diagnostics = self.model.diagnostics
        directory = self.root / "IR_spectra"
        directory.mkdir(parents=True, exist_ok=True)
        max_time_slices = 25
        max_frequency_points = 2000
        for spectrum in diagnostics.ir_spectra:
            if spectrum.intensities.size == 0:
                continue
            wavenumbers = np.asarray(spectrum.wavenumbers, dtype=float,)
            times = np.asarray(spectrum.times, dtype=float,)
            intensities = np.asarray(spectrum.intensities, dtype=float,)
            if intensities.ndim != 2:
                continue
            if intensities.shape != (len(times), len(wavenumbers),):
                continue
            if len(wavenumbers) > max_frequency_points:        # Frequency subsampling
                freq_indices = np.linspace(0, len(wavenumbers) - 1, max_frequency_points, dtype=int,)
            else:
                freq_indices = np.arange(len(wavenumbers), dtype=int,)
            x = wavenumbers[freq_indices]
            if len(times) > max_time_slices:        # Time sampling: many spectra early, progressively fewer later.
                t0 = float(times[0])
                tf = float(times[-1])
                if tf > t0:
                    scaled = np.geomspace(1.0, 1.0 + (tf - t0), max_time_slices,) - 1.0
                    target_times = t0 + scaled
                    time_indices = np.asarray([np.argmin(np.abs(times - target_time)) for target_time in target_times],
                                              dtype=int,)
                    time_indices = np.unique(time_indices)
                else:
                    time_indices = np.array([0], dtype=int,)
            else:
                time_indices = np.arange(len(times), dtype=int,)

            # Always retain the first and final states.
            time_indices = np.unique(np.concatenate(([0], time_indices, [len(times) - 1],)))
            # Normalise all spectra using ONE maximum for this complete time-resolved spectrum.
            plotted_intensities = intensities[:, freq_indices]
            maximum_intensity = np.nanmax(plotted_intensities)
            if (not np.isfinite(maximum_intensity) or maximum_intensity <= 0.0):
                continue
            normalised_intensities = (plotted_intensities / maximum_intensity)
            # Colour represents spectral intensity.
            colour_norm = colors.Normalize(vmin=0.0, vmax=1.0,)
            cmap = cm.viridis
            fig = plt.figure(figsize=(13, 8), clear=True,)
            ax = fig.add_subplot(111, projection="3d",)
            ax.set_proj_type("ortho")
            # Draw each spectrum as short coloured line segments so the colour varies continuously with local IR intensity.
            for time_index in time_indices:
                y = np.full(x.shape, times[time_index], dtype=float,)
                z = normalised_intensities[time_index, :]
                for j in range(len(x) - 1):
                    local_intensity = (0.5 * (z[j] + z[j + 1]))
                    ax.plot(x[j:j + 2], y[j:j + 2], z[j:j + 2], color=cmap(colour_norm(local_intensity)), linewidth=1.2)

            ax.set_xlim(np.max(x), np.min(x),)
            ax.set_ylim(np.min(times), np.max(times),)
            ax.set_zlim(0.0, 1.05,)
            ax.set_xlabel(r"Wavenumber ($cm^{-1}$)", fontsize=18, labelpad=15,)
            ax.set_ylabel("Time ($s$)", fontsize=18, labelpad=15,)
            #ax.text2D(0.03, 0.2, "Time ($s$)", transform=ax.transAxes, fontsize=18,rotation=90,va="center",ha="center",)
            ax.zaxis.set_rotate_label(False)
            ax.set_zlabel("IR intensity", fontsize=18, rotation=90, labelpad=-5)
            ax.set_title(f"IR evolution: {spectrum.phase}, {spectrum.temperature:g} K", fontsize=18,)
            ax.tick_params(axis="x", labelsize=16, pad=2)
            ax.tick_params(axis="y", labelsize=16, pad=5)
            ax.set_zticklabels([])
            #ax.tick_params(axis="z", labelsize=16,)
            time_min = float(np.min(times))
            time_max = float(np.max(times))
            ax.set_ylim(time_min, time_max,)
            time_ticks = np.linspace(time_min, time_max, 6,)
            set_yticks(time_ticks)
            ax.grid(False)
            # Grid only on the x-z plane at the back of the waterfall.
            y_grid = np.max(times)
            #x_ticks = ax.get_xticks()
            z_ticks = np.linspace(0.0, 1.0, 6)
            # Vertical grid lines: constant x, varying z.
            #for x_tick in x_ticks:
            #    ax.plot([x_tick, x_tick], [y_grid, y_grid], [0.0, 1.0], linewidth=0.5,alpha=0.3,color="grey",zorder=0,)
            # Horizontal grid lines: varying x, constant z.
            for z_tick in z_ticks:
                ax.plot([np.min(x), np.max(x)], [y_grid, y_grid], [z_tick, z_tick], linewidth=0.5,
                        alpha=0.3, color="grey", zorder=0,)
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            #ax.xaxis.pane.set_edgecolor("none")
            #ax.yaxis.pane.set_edgecolor("none")
            #ax.zaxis.pane.set_edgecolor("none")
            ax.view_init(elev=20, azim=-91,)    # before 18, -90
            ax.set_box_aspect((3.2, 2.6, 1.6))
            scalar_map = cm.ScalarMappable(norm=colour_norm, cmap=cmap,)
            scalar_map.set_array([])
            colourbar = fig.colorbar(scalar_map, ax=ax, fraction=0.035, pad=0.01, shrink=0.7, aspect=25)
            colourbar.set_label("Normalised IR intensity (a.u.)", fontsize=18, labelpad=5)
            colourbar.ax.yaxis.set_label_position("left")
            colourbar.ax.tick_params(labelsize=16,)
            fig.subplots_adjust(left=0.12, right=0.90, bottom=0.12, top=0.95,)
            #fig.tight_layout()
            fig.savefig(directory / f"{spectrum.phase}_{spectrum.temperature:g}K.svg", dpi=300, transparent=True,
                        bbox_inches=None)
            plt.close(fig)

    def plot_tpr_spectra(self):
        experiments = getattr(self.model, "experiments", None)
        if experiments is None:
            return
        tpr_results = getattr(experiments, "tpr", None)
        if not tpr_results:
            return
        tpr_dir = self.root / "TPR"
        tpr_dir.mkdir(parents=True, exist_ok=True)
        for result in tpr_results:
            self.plot_single_tpr_result(result, tpr_dir)
        self.plot_combined_tpr_spectra(tpr_results, tpr_dir)

    def plot_single_tpr_result(self, result, tpr_dir):
        source_dir = (tpr_dir / f"{result.gas}_{result.heating_rate:g}Kmin")
        source_dir.mkdir(parents=True, exist_ok=True)
        temperatures = np.asarray(result.temperatures, dtype=float,)
        rates = np.asarray(result.spectra.get(result.gas, []), dtype=float,)
        if rates.size == 0:
            return
        if np.nanmax(rates) <= 0.0:
            return
        fig, ax = plt.subplots(figsize=(10, 6), clear=True,)
        ax.plot(temperatures, rates, label=chem_label(result.gas), color=icolour[0], linestyle=iliner[0],)
        self.annotate_tpr_peaks(ax, temperatures, rates, result.gas, 0,)
        ax.set_xlabel("Temperature (K)", fontsize=16,)
        ax.set_ylabel("Desorption Rate", fontsize=16,)
        ax.set_yticklabels([])
        ax.tick_params(axis="x", labelsize=16, pad=2.5,)
        ax.tick_params(axis="y", labelsize=16,)
        ax.set_title(f"TPR spectrum from {result.gas} at {result.heating_rate:g} K/min", fontsize=18,)
        ax.margins(y=0.15)
        ax.set_ylim(bottom=0.0, )  # Pressures cannot physically be negative.
        fig.tight_layout()
        fig.savefig(source_dir / "TPR_single.svg", dpi=300, transparent=True, bbox_inches="tight",)
        plt.close(fig)

    def plot_combined_tpr_spectra(self, tpr_results, tpr_dir,):
        combined_dir = tpr_dir / "Combined"
        combined_dir.mkdir(parents=True, exist_ok=True,)
        fig, ax = plt.subplots(figsize=(10, 6), clear=True,)
        plotted = False
        for i, result in enumerate(tpr_results):
            temperatures = np.asarray(result.temperatures, dtype=float,)
            rates = np.asarray(result.spectra.get(result.gas, []), dtype=float,)
            if rates.size == 0:
                continue
            if np.nanmax(rates) <= 0.0:
                continue
            colour = icolour[i % len(icolour)]
            linestyle = iliner[i % len(iliner)]
            ax.plot(temperatures, rates, color=colour, linestyle=linestyle, label=chem_label(result.gas),)
            self.annotate_tpr_peaks(ax, temperatures, rates, result.gas, i,)
            plotted = True
        if not plotted:
            plt.close(fig)
            return
        ax.set_xlabel("Temperature ($K$)", fontsize=16,)
        ax.set_ylabel("Desorption Rate", fontsize=16,)
        ax.set_yticklabels([])
        ax.tick_params(axis="x", labelsize=16, pad=2.5,)
        ax.tick_params(axis="y", labelsize=16,)
        heating_rates = {float(result.heating_rate) for result in tpr_results}
        if len(heating_rates) == 1:
            heating_rate = next(iter(heating_rates))
            title = (f"Combined TPR spectra at {heating_rate:g} $K/min$")
        else:
            title = "Combined TPR spectra"
        ax.set_title(title, fontsize=18,)
        ax.margins(y=0.15)
        ax.set_ylim(bottom=0.0,)        # Pressures cannot physically be negative.
        ax.legend(fontsize=16,)
        fig.tight_layout()
        fig.savefig(combined_dir / "TPR_combined.svg", dpi=300, transparent=True, bbox_inches="tight",)
        plt.close(fig)

    def annotate_tpr_peaks(self, ax, temperatures, rates, gas, i):
        rates = np.asarray(rates, dtype=float)
        if rates.size == 0:
            return
        maximum = np.nanmax(rates)
        if maximum <= 0.0:
            return
        peak_indices, _ = find_peaks(rates, height=0.25 * maximum,)
        if len(peak_indices) == 0:
            peak_indices = [int(np.nanargmax(rates))]
        for idx in peak_indices:
            peak_temperature = temperatures[idx]
            peak_rate = rates[idx]
            ax.scatter([peak_temperature], [peak_rate],)
            ax.annotate(f"{chem_label(gas)}: {peak_temperature:.0f} K", xy=(peak_temperature, peak_rate),
                        xytext=(5, 5), textcoords="offset points", fontsize=14, color=icolour[i % len(icolour)],)

    def plot_isothermal_species(self):
        """  Plot the time evolution of: (a) all gas-phase species pressures (b) all adsorbate coverages
         at three representative temperatures: lowest, midpoint, and highest available temperature
        Time is displayed on a symmetric-logarithmic axis.   """
        experiments = getattr(self.model, "experiments", None,)
        if experiments is None:
            return
        result = getattr(experiments, "isothermal", None,)
        if result is None:
            return
        values = np.asarray(result.values, dtype=float,)
        if values.ndim != 2 or values.shape[0] == 0:
            return
        headers = list(result.headers)
        if "Temperature" not in headers:
            print("WARNING: isothermal data contain no Temperature column.", flush=True,)
            return
        if "time" not in headers:
            print("WARNING: isothermal data contain no time column.", flush=True,)
            return
        temperature_column = headers.index("Temperature")
        time_column = headers.index("time")

        # result.species contains the dynamic variables: gas-phase molecules + adsorbates.
        # Bare-surface coverages in result.surfaces are intentionally not included in the adsorbate figure.
        dynamic_species = list(result.species)
        gas_indices = list(result.metadata.get("gas_indices", [],))
        adsorbates_by_surface = dict(result.metadata.get("adsorbates_by_surface", {},))
        gas_species = [dynamic_species[index] for index in gas_indices]
        adsorbate_indices = []
        for indices in adsorbates_by_surface.values():
            for index in indices:
                if index not in adsorbate_indices:
                    adsorbate_indices.append(index)
        adsorbates = [dynamic_species[index] for index in adsorbate_indices]
        if not gas_species and not adsorbates:
            return
        # Determine the available temperatures.
        available_temperatures = np.unique(values[:, temperature_column])
        available_temperatures = (available_temperatures[np.isfinite(available_temperatures)])
        if available_temperatures.size == 0:
            return
        available_temperatures.sort()
        minimum_temperature = float(available_temperatures[0])
        maximum_temperature = float(available_temperatures[-1])
        if available_temperatures.size == 1:
            selected_temperatures = [minimum_temperature]
        else:
            target_mid_temperature = (minimum_temperature + maximum_temperature) / 2.0
            middle_index = int(np.argmin(np.abs(available_temperatures - target_mid_temperature)))
            middle_temperature = float(available_temperatures[middle_index])
            selected_temperatures = [minimum_temperature, middle_temperature, maximum_temperature,]
        # Remove duplicates while retaining the selected order.
        selected_temperatures = list(dict.fromkeys(selected_temperatures))

        directory = (self.root / "Isothermal")
        directory.mkdir(parents=True, exist_ok=True,)
        # Determine a sensible linear region for the symlog time axis.
        # symlog is preferable to log because the integration normally contains t = 0.  The smallest strictly positive
        # stored time is used to define the approximately linear region around zero.
        for temperature in selected_temperatures:
            temperature_mask = np.isclose(values[:, temperature_column], temperature, rtol=0.0, atol=1.0e-8,)
            data = values[temperature_mask]
            if data.shape[0] == 0:
                continue
            time_order = np.argsort(data[:, time_column])   # Ensure increasing time even if the table ordering changes.
            data = data[time_order]
            times = np.asarray(data[:, time_column], dtype=float,)
            positive_times = times[np.isfinite(times) & (times > 0.0)]
            if positive_times.size > 0:
                time_linthresh = float(np.min(positive_times))
            else:
                time_linthresh = 1.0e-12
                time_linthresh = max(time_linthresh, np.finfo(float).tiny,)

            if gas_species:
                fig, ax = plt.subplots(figsize=(10, 6), clear=True,)
                plotted = False
                for i, name in enumerate(gas_species):
                    if name not in headers:
                        continue
                    column = headers.index(name)
                    pressures = np.asarray(data[:, column], dtype=float,)
                    if not np.any(np.isfinite(pressures)):
                        continue
                    plot_times, plot_pressures = (self.interpolate_for_symlog_plot(times, pressures,
                                                                                points_per_decade=50,))
                    ax.plot(plot_times, plot_pressures, linewidth=2, color=icolour[ i % len(icolour)],
                            linestyle=iliner[i % len(iliner)], label=chem_label(name),)
                    plotted = True
                if plotted:
                    ax.set_xscale("symlog", linthresh=time_linthresh,)
                    ax.xaxis.set_major_locator(SymmetricalLogLocator(base=10, linthresh=time_linthresh,
                                                                     subs=np.arange(1, 10),))
                    ax.set_yscale("symlog", linthresh=1.0e-6,)
                    ax.set_xlabel("Time ($s$)", fontsize=18,)
                    ax.set_ylabel("Gas pressure ($Pa$)", fontsize=18,)
                    ax.set_title(f"Gas-phase species at {temperature:g} K", fontsize=18,)
                    ax.tick_params(axis="both", labelsize=16,)
                    ax.set_ylim(bottom=0.0,)        # Pressures cannot physically be negative.
                    ax.legend(fontsize=16,)
                    ax.margins(x=0.02, y=0.08,)
                    fig.tight_layout()
                    fig.savefig(directory / f"Gas_{temperature:g}K.svg", dpi=300, transparent=True, bbox_inches="tight")
                plt.close(fig)

            if adsorbates:
                fig, ax = plt.subplots(figsize=(10, 6), clear=True,)
                plotted = False
                for i, name in enumerate(adsorbates):
                    if name not in headers:
                        continue
                    column = headers.index(name)
                    coverages = np.asarray(data[:, column], dtype=float,)
                    if not np.any(np.isfinite(coverages)):
                        continue
                    plot_times, plot_coverages = (self.interpolate_for_symlog_plot(times, coverages,
                                                                                points_per_decade=50,))
                    ax.plot(plot_times, plot_coverages, linewidth=2, color=icolour[i % len(icolour)],
                            linestyle=iliner[i % len(iliner)], label=chem_label(name),)
                    plotted = True
                if plotted:
                    ax.set_xscale("symlog", linthresh=time_linthresh,)
                    ax.xaxis.set_major_locator(SymmetricalLogLocator(base=10, linthresh=time_linthresh,
                                                                     subs=np.arange(1, 10),))
                    ax.set_xlabel("Time ($s$)", fontsize=18)
                    ax.set_ylabel("Coverage ($ML$)", fontsize=18,)
                    ax.set_title(f"Adsorbed species at {temperature:g} K",fontsize=18,)
                    ax.tick_params(axis="both", labelsize=16,)
                    ax.set_ylim(bottom=0.0,)        # Coverages cannot physically be negative.
                    ax.legend(fontsize=16,)
                    ax.margins(x=0.02, y=0.08,)
                    fig.tight_layout()
                    fig.savefig(directory / f"Adsorbates_{temperature:g}K.svg", dpi=300, transparent=True,
                                bbox_inches="tight",)
                plt.close(fig)

    def interpolate_for_symlog_plot(self, times, values, points_per_decade=50,):      # create a dense plotting grid
        """Interpolate a trajectory for smooth display on a symlog time axis."""
        times = np.asarray(times, dtype=float,)
        values = np.asarray(values, dtype=float,)
        finite = (np.isfinite(times) & np.isfinite(values))
        times = times[finite]
        values = values[finite]
        if times.size == 0:
            return times, values
        plot_times = []
        plot_values = []
        # Preserve the calculated value at t = 0.
        zero = times == 0.0
        if np.any(zero):
            plot_times.append(0.0)
            plot_values.append(values[np.flatnonzero(zero)[0]])
        # Positive times are interpolated in log10(time).
        positive = times > 0.0
        positive_times = times[positive]
        positive_values = values[positive]
        if positive_times.size < 3:
            plot_times.extend(positive_times)
            plot_values.extend(positive_values)
            return (np.asarray(plot_times), np.asarray(plot_values),)
        log_times = np.log10(positive_times)
        decades = (log_times[-1] - log_times[0])
        number_of_points = max(int(np.ceil(decades * points_per_decade)) + 1, len(positive_times),)
        dense_log_times = np.linspace(log_times[0], log_times[-1], number_of_points,)
        interpolation = PchipInterpolator(log_times, positive_values, extrapolate=False,)
        dense_values = interpolation(dense_log_times)
        dense_times = (10.0 ** dense_log_times)
        plot_times.extend(dense_times)
        plot_values.extend(dense_values)
        return (np.asarray(plot_times), np.asarray(plot_values),)

class Writer:
    """
    Write human-readable fixed-width ASCII tables.
    """
    float_format = "{:.3e}"
    @classmethod
    def format_value(cls, value: Any) -> str:
        if value is None:
            return "---"
        if isinstance(value, bool):
            return str(value)
        if isinstance(value, int):
            return str(value)
        if isinstance(value, float):
            if abs(value) < 1e-50:
                value = 0.0
            return cls.float_format.format(value)
        return str(value)

    @classmethod
    def write(cls, filename: str | Path, header: list[str], rows: list[list[Any]],title: str | None = None,) -> None:
        filename = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)
        header_fmt = [str(item) for item in header]
        rows_fmt = [[cls.format_value(value) for value in row] for row in rows]
        ncols = len(header_fmt)
        widths = []
        for i in range(len(header_fmt)):
            column_widths = [len(header_fmt[i])]
            for row in rows_fmt:
                column_widths.append(len(row[i]))
            widths.append(max(column_widths))

        def line(values):
            return "  ".join(str(value).rjust(width) for value, width in zip(values, widths))

        with open(filename, "w", encoding="utf-8") as handle:
            if title:
                handle.write(f"# {title}\n")
            handle.write("# " + line(header_fmt) + "\n")
            handle.write("# " + line(["-" * w for w in widths]) + "\n")
            for row in rows_fmt:
                handle.write("  " + line(row) + "\n")

def build_model(input_file: str, run_experiments: bool = True) -> Any:
    if read_input is None:
        raise RuntimeError("Could not import Main_refactored.read_input")
    model = read_input(input_file)
    for routine in [PartitionFunctions, Energy, Profile, RConstants, REquations]:
        if routine is None:
            continue
        try:
            routine(model)
        except TypeError:
            continue
    if run_experiments:
        for routine in [Isothermal, TPR]:
            if routine is None:
                continue
            try:
                routine(model)
            except TypeError:
                continue
            except Exception:
                continue
    return model

def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate and plot all stored thermodynamic, kinetic and experimental model results.")
    parser.add_argument("input", help="Input .mk.in file")
    parser.add_argument("--out", default="MODEL_REPORT", help="Output directory")
    parser.add_argument("--skip-experiments", action="store_true", help="Do not run Isothermal/TPR before printing.")
    args = parser.parse_args()
    model = build_model(args.input, run_experiments=not args.skip_experiments)
    ModelPrinter(model, args.out).write_all()

if __name__ == "__main__":
    main()
