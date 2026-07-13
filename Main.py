"""
Main driver for the microkinetic workflow.

Workflow:
    Atomistic outputs -> Input2mk.py -> file.mk.in -> Main_refactored.py
    -> Thermodynamics.py -> Kinetics.py -> Experiments.py

by Alberto Roldan
"""

from __future__ import annotations
import argparse
import re
import time
from dataclasses import dataclass, field
from enum import Enum
import sympy as sp

from Thermodynamics import Energy, PartitionFunctions
from Kinetics import Profile, RConstants, REquations
from Experiments import Isothermal, TPR
from Kinetics import KineticsResults
from Experiments import ExperimentResults
from Diagnostics import DiagnosticResults, Diagnostics
from PrintModel import ThermodynamicsReporter, KineticsReporter, ExperimentReporter, DiagnosticsReporter


class SpeciesType(str, Enum):
    molecule = "Molecule"
    surface = "Surface"
    adsrobate = "Adsorbate"
    ts = "TransitionState"
    unknown = "Unknown"

@dataclass
class CoveragePoint:    # defines the coverage
    coverage: float
    energy0: float = 0.0
    frequencies: list[float] = field(default_factory=list)

@dataclass
class Species:      # SYSTEMS in the reaction (not in coverage)
    name: str
    energy0: float = 0.0
    kind: str = SpeciesType.unknown.value
    frequencies: list[float] = field(default_factory=list)
    ir_intensities: list[float] = field(default_factory=list)
    frequencies_2d: list[float] = field(default_factory=list)
    degeneracy: int = 1
    mass: float | None = None
    symmetry_factor: float = 1.0
    pressure0: float = 0.0
    coverage0: float = 0.0
    linear: bool = False
    inertia: list[float] = field(default_factory=list)
    nmolsites: float = 1.0
    molsite: str | None = None
    volume: sp.Expr | None = None
    pressure: sp.Symbol | None = None
    area: float | None = None  # surface properties
    sites: str | None = None
    nsites: int | None = None
    coverage_model: list[CoveragePoint] = field(default_factory=list)

@dataclass
class SpeciesCoefficient:   # reaction stoichiometry
    species: str
    coefficient: float = 1.0

@dataclass(slots=True)
class KineticProperties:
    activation: sp.Expr | float = 0.0
    sticky: sp.Expr | float = 1.0
    arrhenius: sp.Expr | float = 1.0
    ktunneling: sp.Expr | float = 1.0
    krate0: sp.Expr | float = 0.0

@dataclass(slots=True)
class Process:
    id: int
    kind: str
    reactants: list[SpeciesCoefficient]
    products: list[SpeciesCoefficient]
    ts: list[SpeciesCoefficient] = field(default_factory=list)
    kinetics: KineticProperties = field(default_factory=KineticProperties)
    @property
    def reactant_names(self) -> list[str]:
        return [item.species for item in self.reactants]
    @property
    def product_names(self) -> list[str]:
        return [item.species for item in self.products]
    @property
    def ts_names(self) -> list[str]:
        return [item.species for item in self.ts]
    @property
    def reactant_stoich(self) -> list[float]:
        return [item.coefficient for item in self.reactants]
    @property
    def product_stoich(self) -> list[float]:
        return [item.coefficient for item in self.products]

@dataclass(slots=True)
class Scan:
    start: float
    stop: float | None = None
    step: float | None = None
    def values(self) -> list[float]:
        if self.stop is None:
            return [self.start]
        if self.step is None or self.step == 0:
            raise ValueError("Scan requires a non-zero step when stop is provided.")
        values: list[float] = []
        x = self.start
        if self.step > 0:
            while x <= self.stop + 1e-12:
                values.append(x)
                x += self.step
        else:
            while x >= self.stop - 1e-12:
                values.append(x)
                x += self.step
        return values

@dataclass(slots=True)
class TPRConditions:
    temperature_initial: float = 50.0
    temperature_final: float = 1273.0
    temperature_step: float = 1.0
    heating_rate: float = 5.0      # K/min

@dataclass(slots=True)
class Conditions:
    temperature: Scan | None = None
    time: Scan | None = None
    potential: Scan | None = None
    ph: Scan | None = None
    tpr: TPRConditions = field(default_factory=TPRConditions)

@dataclass(slots=True)
class EquationResults:
    isothermal: dict[str, list[sp.Expr]] = field(default_factory=dict)
    equation_factors: dict[str, list[float]] = field(default_factory=dict)

    tpd: dict[str, list[sp.Expr]] = field(default_factory=dict)
    tpd_factors: dict[str, list[float]] = field(default_factory=dict)

    overall_stoichiometry: list[sp.Expr] = field(default_factory=list)

@dataclass(slots=True)
class ProfileResults:
    data: list[list[str]] = field(default_factory=list)

@dataclass(slots=True)
class MKModel:
    conditions: Conditions = field(default_factory=Conditions)
    species: dict[str, Species] = field(default_factory=dict)
    processes: dict[str, Process] = field(default_factory=dict)
    equations: EquationResults = field(default_factory=EquationResults)
    kinetics: KineticsResults = field(default_factory=KineticsResults)
    experiments: ExperimentResults = field(default_factory=ExperimentResults)
    profile: ProfileResults = field(default_factory=ProfileResults)
    diagnostics: DiagnosticResults | None = None

# --- PARSERS
def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()

def split_keyword(line: str) -> tuple[str, str]:
    clean = strip_comment(line)
    if "=" not in clean:
        raise ValueError(f"Expected KEY = VALUE, got: {line!r}")
    key, value = clean.split("=", 1)
    return key.strip().upper(), value.strip()

def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "yes", "y", "1", ".true."}

def parse_scan(value: str) -> Scan:
    values = [float(x) for x in value.split()]
    if len(values) == 1:
        return Scan(start=values[0])
    elif len(values) == 2:
        return Scan(start=0.0, stop=values[0], step=values[1])
    elif len(values) == 3:
        return Scan(start=values[0], stop=values[1], step=values[2])
    raise ValueError( f"Scan expects 1, 2, or 3 values, got {len(values)}: {value}")

def parse_species_process(text: str) -> list[SpeciesCoefficient]:
    items: list[SpeciesCoefficient] = []
    stoich_re = re.compile(r"^([0-9]*\.?[0-9]+)\s*([A-Za-z_][A-Za-z0-9_]*)$")
    for raw_term in text.split("+"):
        term = raw_term.strip()
        if not term:
            continue
        match = stoich_re.match(term)
        if match:
            coeff = float(match.group(1))
            name = match.group(2)
        else:
            coeff = 1.0
            name = term
        items.append(SpeciesCoefficient(species=name, coefficient=coeff,))
    return items

def parse_process(line: str, number: int) -> Process:
    _, right = line.split("=", 1)
    fields = right.strip().split(maxsplit=1)
    if len(fields) != 2:
        raise ValueError(f"Invalid PROCESS line: {line}")
    kind = fields[0]
    reaction = fields[1]
    pieces = [piece.strip() for piece in reaction.split(">")]
    if len(pieces) == 2:
        reactants = parse_species_process(pieces[0])
        ts = []
        products = parse_species_process(pieces[1])
    elif len(pieces) == 3:
        reactants = parse_species_process(pieces[0])
        if pieces[1] in {"", "1", "none", "None"}:
            ts = []
        else:
            ts = parse_species_process(pieces[1])
        products = parse_species_process(pieces[2])
    else:
        raise ValueError(f"Cannot parse process:\n{line}")
    return Process(id=number, kind=kind, reactants=reactants, ts=ts, products=products,)

def parse_system(lines, i):
    name = lines[i].split("=", 1)[1].strip()
    species = Species(name=name)
    i += 1
    while not lines[i].startswith("END"):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue
        if line.startswith("MODEL_COVERAGE"):
            i += 1
            while not lines[i].startswith("END_MODEL_COVERAGE"):
                key, value = split_keyword(lines[i])
                if key == "POINT":
                    fields = value.split()
                    coverage = float(fields[0])
                    energy0 = float(fields[1])
                    frequencies = []
                    if len(fields) > 2:
                        frequencies = [float(x) for x in fields[2:]]
                    species.coverage_model.append(CoveragePoint(coverage=coverage, energy0=energy0,
                                                                frequencies=frequencies,))
                i += 1
        else:
            key, value = split_keyword(line)
            if key == "TYPE":
                species.kind = value
            elif key in {"E0", "ENERGY", "ENERGY0"}:
                species.energy0 = float(value)
            elif key in {"DEGENERATION", "DEGENERACY"}:
                species.degeneracy = int(float(value))
            elif key == "DEGENERACY":
                species.degeneracy = int(value)
            elif key == "FREQ":
                species.frequencies = [float(x) for x in value.split()]
            elif key == "IRINTENSITY":
                species.ir_intensities = [float(x) for x in value.split()]
            elif key == "FREQ2D":
                species.frequencies_2d = [float(x) for x in value.split()]
            elif key == "MASS":
                species.mass = float(value)
            elif key == "SYMFACTOR":
                species.symmetry_factor = float(value)
            elif key == "LINEAR":
                species.linear = parse_bool(value)
            elif key == "INERTIA":
                species.inertia = [float(x) for x in value.split()]
            elif key == "IPRESSURE":
                species.pressure0 = float(value)
            elif key == "ICOVERAGE":
                species.coverage0 = float(value)
            elif key == "MOLSITE":
                fields = value.split()
                if len(fields) == 1:
                    species.nmolsites = 1.0
                    species.molsite = fields[0]
                elif len(fields) == 2:
                    species.nmolsites = float(fields[0])
                    species.molsite = fields[1]
                else:
                    raise ValueError(f"{species.name}: invalid MOLSITE definition: {value}")
            elif key == "IACAT":
                species.area = float(value)
            elif key == "ISITES":
                fields = value.split()
                if len(fields) == 1:
                    species.nsites = 1
                    species.sites = fields[0]
                elif len(fields) == 2:
                    species.nsites = int(float(fields[0]))
                    species.sites = fields[1]
                else:
                    raise ValueError(f"{species.name}: invalid ISITES definition: {value}")

        i += 1
    return species, i + 1

def read_input(filename: str) -> MKModel:
    model = MKModel()
    with open(filename) as fh:
        lines = [line.strip() for line in fh]
    i = 0
    while i < len(lines):
        line = strip_comment(lines[i])
        if not line:
            i += 1
            continue
        if line.startswith("TEMPERATURE"):
            _, value = split_keyword(line)
            model.conditions.temperature = parse_scan(value)
        elif line.startswith("TIME"):
            _, value = split_keyword(line)
            model.conditions.time = parse_scan(value)
        elif line.startswith("PH"):
            _, value = split_keyword(line)
            model.conditions.ph = parse_scan(value)
        elif line.startswith("POTENTIAL") or line.startswith("ELECTROCHEMICAL"):
            _, value = split_keyword(line)
            model.conditions.potential = parse_scan(value)
        elif line.startswith("TPR_TEMPERATURE"):
            _, value = split_keyword(line)
            values = [float(x) for x in value.split()]
            model.conditions.tpr.temperature_initial = values[0]
            model.conditions.tpr.temperature_final = values[1]
            model.conditions.tpr.temperature_step = values[2]
        elif line.startswith("TPR_HEATING_RATE"):
            _, value = split_keyword(line)
            model.conditions.tpr.heating_rate = float(value)
        elif line.startswith("PROCESS"):
            process_id = len(model.processes) + 1
            process = parse_process(line, process_id)
            model.processes[str(process.id)] = process
        elif line.startswith("SYSTEM"):
            species, i = parse_system(lines, i)
            model.species[species.name] = species
            continue
        else:
            raise ValueError(f"Unknown top-level instruction: {line}")
        i += 1
    return model

# --- OTHERS
def run_stage(label, func, model):
    start = time.time()
    print(f"{label}...")
    result = func(model)
    elapsed = time.time() - start
    if elapsed < 30:
        print(f"\t{elapsed:.3f} seconds")
    else:
        print(f"\t{elapsed/60:.3f} minutes")
    return result

def main():
    start_total = time.time()
    print("Reading input...")
    start = time.time()
    model = read_input(args.input)
    print(f"\t{time.time() - start:.3f} seconds")

    run_stage("Building partition functions", PartitionFunctions, model)
    run_stage("Building thermodynamics", Energy, model)
    ThermodynamicsReporter(model, output_dir="./").write_and_plot()
    run_stage("Building rate constants", RConstants, model)
    KineticsReporter(model, output_dir="./").write_and_plot()
    run_stage("Building rate equations", REquations, model)
    run_stage("Building energy profile", Profile, model)
    run_stage("Running microkinetics", Isothermal, model)
    ExperimentReporter(model, output_dir="./").write()
    run_stage("Running TPR/D", TPR, model)
    ExperimentReporter(model, output_dir="./").write()
    run_stage("Diagnosis", Diagnostics, model)
    DiagnosticsReporter(model, output_dir="./").write_and_plot()

    total = time.time() - start_total
    if total < 60:
        print(f"\nFinished in {total:.3f} seconds")
    else:
        print(f"\nFinished in {total/60:.3f} minutes")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    args = parser.parse_args()
    main()