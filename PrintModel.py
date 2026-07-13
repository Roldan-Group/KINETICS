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
from scipy.signal import find_peaks

from Symbols_def import vext, ph, temp, constants


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
            self.plot_partition_functions(species, species_dir)
            self.plot_entropy(species, species_dir)
            self.plot_enthalpy(species, species_dir)
            self.plot_heat_capacity(species, species_dir)
            self.plot_gibbs(species, species_dir)

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

    def plot_fields(self, species, species_dir, fields, filename, ylabel):
        temperatures = np.asarray(self.temperatures(), dtype=float)
        fig, ax = plt.subplots(figsize=(8, 6), clear=True)
        for n,field in enumerate(fields):
            values = [self.evaluate(getattr(species.thermo, field, None), temp_num) for temp_num in temperatures]
            if all(value is None for value in values):
                continue
            values = np.asarray([np.nan if value is None else value for value in values], dtype=float,)
            ax.plot(temperatures, values, color=icolour[n % len(icolour)], linestyle=iliner[n % len(iliner)], label=field,)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(ylabel)
        ax.set_title(species.name)
        ax.legend()
        fig.tight_layout()
        fig.savefig(species_dir / filename, dpi=300, transparent=True)
        plt.close(fig)

    def plot_partition_functions(self, species, species_dir):
        self.plot_fields(species, species_dir, ["qtrans3d", "qrot", "qelec", "qvib3d", "q3d",
                                                "qtrans2d", "qvib2d", "q2d",],
                         "Partition_functions.svg", "Partition function",)

    def plot_entropy(self, species, species_dir):
        self.plot_fields(species, species_dir,["strans3d", "selec", "srot", "svib3d", "sentropy3d",
                                               "strans2d", "svib2d", "sentropy2d",], "Entropy.svg", "Entropy (eV/K)",)

    def plot_enthalpy(self, species, species_dir):
        self.plot_fields(species, species_dir, ["zpe3d", "hthermal3d", "enthalpy3d",
                                                "zpe2d", "hthermal2d", "enthalpy2d",], "Enthalpy.svg", "Enthalpy (eV)",)

    def plot_heat_capacity(self, species, species_dir):
        self.plot_fields(species, species_dir, ["cp3d", "cp2d",], "Heat_capacity.svg", "Heat capacity (eV/K)",)

    def plot_gibbs(self, species, species_dir):
        self.plot_fields(species, species_dir, ["gibbs3d", "gibbs2d"], "Gibbs.svg", "Gibbs free energy (eV)",)


class KineticsReporter:
    def __init__(self, model, output_dir="./"):
        self.model = model
        self.output_dir = Path(output_dir)
        self.root = self.output_dir / "KINETICS"

    def write_and_plot(self):
        self.root.mkdir(parents=True, exist_ok=True)
        self.write_rate_constants()
        self.write_reaction_energies()
        self.plot_rate_constants()
        self.plot_activation_energies()
        self.plot_reaction_energies()

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

    def kinetic_rows(self, fields):
        rows = []
        for temp_num in self.temperatures():
            for pid, kinetics in self.model.kinetics.by_process.items():
                process = self.model.processes[pid]
                row = [temp_num, pid, process.kind,]
                for field in fields:
                    expr = getattr(kinetics, field, None)
                    row.append(self.evaluate(expr, temp_num))
                rows.append(row)
        return rows

    def write_rate_constants(self):
        fields = ["sticky", "arrhenius", "ktunneling", "krate0",]
        Writer.write(filename=self.root / "Rate_constants.dat", header=["Temperature(K)", "Process", "Kind", *fields,],
                     rows=self.kinetic_rows(fields), title="Rate constant contributions by process",)

    def write_reaction_energies(self):
        fields = ["activation", "reaction_energy",]
        Writer.write(filename=self.root / "Reaction_energies.dat", header=["Temperature(K)", "Process", "Kind",
                                                                        "activation_energy(eV)", "reaction_energy(eV)",],
                     rows=self.kinetic_rows(fields), title="Reaction free energies by process",)

    def plot_process_field(self, field, filename, ylabel, logy=False):
        temperatures = self.temperatures()
        process_ids = list(self.model.kinetics.by_process.keys())

        for temp_num in temperatures:
            values = []
            for pid in process_ids:
                kinetics = self.model.kinetics.by_process[pid]
                value = self.evaluate(getattr(kinetics, field, None), temp_num,)
                if value is None:
                    value = np.nan
                values.append(value)
            fig, ax = plt.subplots(figsize=(10, 6), clear=True)
            x = np.arange(len(process_ids))
            ax.bar(x, values)
            ax.set_xticks(x)
            ax.set_xticklabels(process_ids, rotation=45, ha="right")
            ax.set_xlabel("Process ID")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{field} at T = {temp_num:g} K")
            if logy:
                positive = [value for value in values if value is not None and np.isfinite(value) and value > 0]
                if positive:
                    ax.set_yscale("log")
            fig.tight_layout()
            fig.savefig(self.root / f"{filename}_{temp_num:g}K.svg", dpi=300, transparent=True,)
            plt.close(fig)

    def plot_rate_constants(self):
        self.plot_process_field(field="krate0", filename="Rate_constants", ylabel="Rate constant", logy=True,)

    def plot_activation_energies(self):
        self.plot_process_field(field="activation", filename="Activation_energies", ylabel="Activation energy (eV)",
                                logy=False,)

    def plot_reaction_energies(self):
        self.plot_process_field(field="reaction_energy", filename="Reaction_energies", ylabel="Reaction free energy (eV)",
                                logy=False,)


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
            label = (f"{result.gas}_{result.heating_rate:g}K_min")
            self.write_tpr_concentrations(result, tpr_dir / f"{label}_concentrations.dat",)
            self.write_tpr_spectra(result, tpr_dir / f"{label}_spectra.dat",)
            self.write_tpr_metadata(result, tpr_dir / f"{label}_metadata.dat",)

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
        self.root = self.output_dir / "DIAGNOSTICS"

    def write_and_plot(self):
        self.root.mkdir(parents=True, exist_ok=True)
        diagnostics = getattr(self.model, "diagnostics", None)
        if diagnostics is None:
            return
        self.write_table_result(diagnostics.steady_state_net_rates, self.root / "Steady_state_net_rates.dat",)
        self.write_table_result(diagnostics.average_net_rates, self.root / "Average_net_rates.dat",)
        self.write_table_result(diagnostics.equilibrium_control, self.root / "Equilibrium_control.dat",)
        self.write_table_result(diagnostics.drc_assessment, self.root / "DRC_assessment.dat",)
        self.write_control_tables(diagnostics.degree_of_rate_control, self.root / "Degree_of_rate_control",)
        self.write_control_tables(diagnostics.degree_of_selectivity_control, self.root/"Degree_of_selectivity_control",)
        self.write_reaction_pathways()
        self.write_tpr_peaks()
        self.write_ir_spectra()
        self.plot_table_bars(diagnostics.steady_state_net_rates, "Steady_state_net_rates.svg",
                             "Steady-state net rates", "Net rate",)
        self.plot_table_bars(diagnostics.average_net_rates, "Average_net_rates.svg",
                             "Average net rates", "Average net rate",)
        self.plot_table_bars(diagnostics.equilibrium_control, "Equilibrium_control.svg",
                             "Equilibrium control", "DEC",)
        self.plot_control_tables(diagnostics.degree_of_rate_control, self.root /"Degree_of_rate_control", ylabel="DRC",)
        self.plot_control_tables(diagnostics.degree_of_selectivity_control, self.root / "Degree_of_selectivity_control",
                                 ylabel="DSC",)
        self.plot_reaction_pathways()
        self.plot_tpr_peaks()
        self.plot_tpr_spectra()
        self.plot_ir_spectra()

    def write_table_result(self, table, filename):
        if table is None:
            return
        Writer.write(filename=filename, header=table.headers, rows=table.rows, title=filename.stem.replace("_", " "),)

    def write_control_tables(self, tables, directory):
        directory.mkdir(parents=True, exist_ok=True)
        for name, table in tables.items():
            self.write_table_result(table, directory / f"{name}.dat",)

    def write_reaction_pathways(self):
        diagnostics = self.model.diagnostics
        filename = self.root / "Reaction_pathways.dat"
        rows = []
        for pathway in diagnostics.reaction_pathways:
            for reactant, product, flux in pathway.edges:
                rows.append([pathway.temperature, reactant, product, flux,])
        Writer.write(filename=filename, header=["Temperature(K)", "Reactant", "Product", "Flux",],
                     rows=rows, title="Flux-weighted reaction pathways",)

    def write_tpr_peaks(self):
        diagnostics = self.model.diagnostics
        filename = self.root / "TPR_peaks.dat"
        rows = []
        for peak in diagnostics.tpr_peaks:
            rows.append([peak.source_gas, peak.signal_gas, peak.heating_rate, peak.peak_temperature,
                         peak.peak_rate, peak.area,])
        Writer.write(filename=filename, header=["SourceGas", "SignalGas", "HeatingRate(K/min)", "PeakTemperature(K)",
                                                "PeakRate", "Area",], rows=rows, title="TPR peak summary",)

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

    def plot_table_bars(self, table, filename, title, ylabel):
        if table is None or not table.rows:
            return
        headers = table.headers
        xlabels = headers[1:]
        for row in table.rows:
            temperature = row[0]
            values = np.asarray(row[1:], dtype=float)
            fig, ax = plt.subplots(figsize=(10, 6), clear=True)
            x = np.arange(len(xlabels))
            ax.bar(x, values)
            ax.set_xticks(x)
            ax.set_xticklabels(xlabels, rotation=45, ha="right")
            ax.set_xlabel("Reaction step")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{title} at {temperature:g} K")
            fig.tight_layout()
            output = (self.root / filename.replace(".svg", f"_{temperature:g}K.svg",))
            fig.savefig(output, dpi=300, transparent=True,)
            plt.close(fig)

    def plot_control_tables(self, tables, directory, ylabel):
        directory.mkdir(parents=True, exist_ok=True)
        for name, table in tables.items():
            if table is None or not table.rows:
                continue
            headers = table.headers
            xlabels = headers[1:]
            for row in table.rows:
                temperature = row[0]
                values = np.asarray(row[1:], dtype=float)
                fig, ax = plt.subplots(figsize=(10, 6), clear=True)
                x = np.arange(len(xlabels))
                ax.bar(x, values)
                ax.set_xticks(x)
                ax.set_xticklabels(xlabels, rotation=45, ha="right")
                ax.set_xlabel("Reaction step")
                ax.set_ylabel(ylabel)
                ax.set_title(f"{name} {ylabel} at {temperature:g} K")
                fig.tight_layout()
                fig.savefig(directory / f"{name}_{temperature:g}K.svg", dpi=300, transparent=True,)
                plt.close(fig)

    def plot_reaction_pathways(self):
        diagnostics = self.model.diagnostics
        directory = self.root / "Reaction_pathways"
        directory.mkdir(parents=True, exist_ok=True)
        for pathway in diagnostics.reaction_pathways:
            if not pathway.edges:
                continue
            labels = [f"{reactant}_to_{product}" for reactant, product, _ in pathway.edges]
            values = [flux for _, _, flux in pathway.edges]
            fig, ax = plt.subplots(figsize=(10, 6), clear=True)
            y = np.arange(len(labels))
            ax.barh(y, values)
            ax.set_yticks(y)
            ax.set_yticklabels(labels)
            ax.set_xlabel("Flux")
            ax.set_title(f"Reaction pathway at {pathway.temperature:g} K")
            fig.tight_layout()
            fig.savefig(directory / f"pathway_{pathway.temperature:g}K.svg", dpi=300, transparent=True,)
            plt.close(fig)

    def plot_tpr_peaks(self):
        diagnostics = self.model.diagnostics
        if not diagnostics.tpr_peaks:
            return
        labels = [f"{peak.source_gas}/{peak.signal_gas}" for peak in diagnostics.tpr_peaks]
        temperatures = [peak.peak_temperature for peak in diagnostics.tpr_peaks]
        fig, ax = plt.subplots(figsize=(10, 6), clear=True)
        x = np.arange(len(labels))
        ax.bar(x, temperatures)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel("Peak temperature (K)")
        ax.set_title("TPR peak temperatures")
        fig.tight_layout()
        fig.savefig(self.root / "TPR_peak_temperatures.svg", dpi=300, transparent=True,)
        plt.close(fig)

    def plot_ir_spectra(self):
        diagnostics = self.model.diagnostics
        directory = self.root / "IR_spectra"
        directory.mkdir(parents=True, exist_ok=True)
        for spectrum in diagnostics.ir_spectra:
            if spectrum.intensities.size == 0:
                continue
            fig, ax = plt.subplots(figsize=(10, 6), clear=True)
            mesh = ax.pcolormesh(spectrum.wavenumbers, spectrum.times, spectrum.intensities, shading="auto",)
            fig.colorbar(mesh, ax=ax, label="Intensity (a.u.)",)
            ax.set_xlabel("Wavenumber (cm$^{-1}$)")
            ax.set_ylabel("Time (s)")
            ax.set_title(f"IR evolution: {spectrum.phase}, {spectrum.temperature:g} K")
            fig.tight_layout()
            fig.savefig(directory / f"{spectrum.phase}_{spectrum.temperature:g}K.svg", dpi=300, transparent=True,)
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
        temperatures = np.asarray(result.temperatures, dtype=float)
        fig, ax = plt.subplots(figsize=(10, 6), clear=True)
        for gas, rates in result.spectra.items():
            rates = np.asarray(rates, dtype=float)
            if rates.size == 0:
                continue
            if np.nanmax(rates) <= 0.0:
                continue
            ax.plot(temperatures, rates, label=gas,)
            self.annotate_tpr_peaks(ax, temperatures, rates, gas,)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Desorption Rate")
        ax.set_title(f"TPR spectrum from {result.gas} at {result.heating_rate:g} K/min")
        handles, labels = ax.get_legend_handles_labels()
        if not handles:
            plt.close(fig)  #prevents saving the figure when all spectra is zero
            return
        ax.legend()
        fig.tight_layout()
        fig.savefig(source_dir / "TPR_single.svg", dpi=300, transparent=True,)
        plt.close(fig)

    def plot_combined_tpr_spectra(self, tpr_results, tpr_dir):
        combined_dir = tpr_dir / "Combined"
        combined_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10, 6), clear=True)
        plotted = False
        for result in tpr_results:
            temperatures = np.asarray(result.temperatures, dtype=float)
            for gas, rates in result.spectra.items():
                rates = np.asarray(rates, dtype=float)
                if rates.size == 0:
                    continue
                if np.nanmax(rates) <= 0.0:
                    continue
                label = f"{result.gas}_to_{gas} ({result.heating_rate:g} K/min)"
                ax.plot(temperatures, rates, label=label,)
                plotted = True
        if not plotted:
            plt.close(fig)
            return
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Desorption Rate")
        ax.set_title("Combined TPR spectra")
        ax.legend()
        fig.tight_layout()
        fig.savefig(combined_dir / "TPR_combined.svg", dpi=300, transparent=True,)
        plt.close(fig)

    def annotate_tpr_peaks(self, ax, temperatures, rates, gas,):
        rates = np.asarray(rates, dtype=float)
        if rates.size == 0:
            return
        maximum = np.nanmax(rates)
        if maximum <= 0.0:
            return
        peak_indices, _ = find_peaks(rates, height=0.05 * maximum,)
        if len(peak_indices) == 0:
            peak_indices = [int(np.nanargmax(rates))]
        for idx in peak_indices:
            peak_temperature = temperatures[idx]
            peak_rate = rates[idx]
            ax.scatter([peak_temperature], [peak_rate],)
            ax.annotate(f"{gas}: {peak_temperature:.0f} K", xy=(peak_temperature, peak_rate),
                        xytext=(5, 5), textcoords="offset points", fontsize=10,)


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
