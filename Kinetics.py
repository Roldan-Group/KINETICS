"""
Side-effect-free kinetics module for the object-based microkinetic workflow.

This module is designed to work with Main.py, where each stage is called as Routine(model).

Expected model structure:
    model.species[name]      -> Species object
    model.processes[id]     -> Process object
    species.thermo          -> ThermodynamicProperties object

Main outputs kept in memory:
    RConstants(model).results
    REquations(model).all_equations
    Profile(model).data
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import reduce
from math import gcd
from typing import Iterable
import numpy as np
import sympy as sp
from Symbols_def import temp, h, kb, hc, JtoeV, constants


@dataclass(slots=True)
class KineticProperties:
    activation: sp.Expr | float = 0.0
    sticky: sp.Expr | float = 1.0
    arrhenius: sp.Expr | float = 1.0
    ktunneling: sp.Expr | float = 1.0
    krate0: sp.Expr | float = 0.0
    reaction_energy: sp.Expr | float = 0.0

@dataclass(slots=True)
class KineticsResults:
    by_process: dict[str, KineticProperties] = field(default_factory=dict)

# --- DEFINITIONS
def _model_key(model) -> int:
    return id(model)

def kinetic_results(model) -> KineticsResults:
    """Return the in-memory kinetic results associated with a model."""
    key = _model_key(model)
    return model.kinetics

def process_key(process) -> str:
    return str(process.id)

def process_kinetics(model, process) -> KineticProperties:
    results = kinetic_results(model)
    key = process_key(process)
    if key not in results.by_process:
        results.by_process[key] = KineticProperties()
    kin = results.by_process[key]
    if hasattr(process, "kinetics"):
        try:
            if process.kinetics is None:
                process.kinetics = kin
            else:
                process.kinetics.activation = kin.activation
                process.kinetics.sticky = kin.sticky
                process.kinetics.arrhenius = kin.arrhenius
                process.kinetics.ktunneling = kin.ktunneling
                process.kinetics.krate0 = kin.krate0
        except (AttributeError, TypeError):
            pass
    return kin

def get_thermo(species):
    if not hasattr(species, "thermo"):
        raise AttributeError(f"{species.name}: missing `thermo`; run Thermodynamics before Kinetics.")
    return species.thermo

def species_kind(species) -> str:
    return str(species.kind).lower()

def is_molecule(species) -> bool:
    return species_kind(species) == "molecule"

def is_surface(species) -> bool:
    return species_kind(species) == "surface"

def is_adsorbed_like(species) -> bool:
    return species_kind(species) in {"adsorbate", "transitionstate", "ts"}

def first_existing_attr(obj, names: Iterable[str], default=None):
    for name in names:
        if hasattr(obj, name):
            value = getattr(obj, name)
            if value is not None:
                return value
    return default


class RConstants:
    def __init__(self, model):
        self.model = model
        self.build()

    def build(self):
        for pid, process in self.model.processes.items():
            kinetics = KineticProperties()
            kinetics.activation = self.activation(process)
            kinetics.reaction_energy = self.reaction_energy(process)
            kinetics.arrhenius = self.arrhenius(process)
            kinetics.ktunneling = self.tunneling(process)
            if process.kind.upper() == "A":
                kinetics.sticky = self.sticky(process, kinetics.activation,)
                kinetics.krate0 = sp.simplify(kinetics.sticky * kinetics.arrhenius * kinetics.ktunneling)
            else:
                kinetics.sticky = sp.Integer(1)
                kinetics.krate0 = sp.simplify(kinetics.arrhenius * sp.exp(-kinetics.activation / (kb * temp * JtoeV))
                                                * kinetics.ktunneling)
            self.model.kinetics.by_process[pid] = kinetics
        return self.model

    def species(self, name):
        return self.model.species[name]

    @staticmethod
    def is_molecule(species):
        return species.kind.lower() == "molecule"

    def gibbs3d(self, name):
        return self.species(name).thermo.gibbs3d

    def gibbs2d(self, name):
        return self.species(name).thermo.gibbs2d

    def q3d(self, name):
        return self.species(name).thermo.q3d

    def state_gibbs3d(self, items):
        total = sp.Integer(0)
        for item in items:
            total += item.coefficient * self.gibbs3d(item.species)
        return sp.simplify(total)

    def reaction_energy(self, process):
        return sp.simplify(self.state_gibbs3d(process.products) - self.state_gibbs3d(process.reactants))

    def activation(self, process):
        ''' the activation energy in adsorption and desorption processes is considered as the difference between
        a state in which the molecule has only two degrees of freedom (being the third degree the reaction coordinate)
        and the reactants '''
        reactants_g = self.state_gibbs3d(process.reactants)
        if process.ts:
            ts_g = self.state_gibbs3d(process.ts)
        elif any(self.is_molecule(self.species(item.species)) for item in process.reactants):
            ts_g = self.adsorption_like_ts_gibbs(process.reactants)
        elif any(self.is_molecule(self.species(item.species)) for item in process.products):
            ts_g = self.adsorption_like_ts_gibbs(process.products)
        else:
            ts_g = self.state_gibbs3d(process.products)
        return sp.Max(sp.simplify(ts_g - reactants_g), sp.Integer(0),)

    def adsorption_like_ts_gibbs(self, items):
        ''' the activation energy in adsorption and desorption processes is considered as the difference between
        a state in which the molecule has only two degrees of freedom (being the third degree the reaction coordinate)
        and the reactants '''
        total = sp.Integer(0)
        for item in items:
            species = self.species(item.species)
            if self.is_molecule(species):
                total += item.coefficient * species.thermo.gibbs2d
            else:
                total += item.coefficient * species.thermo.gibbs3d
        return sp.simplify(total)

    def state_q3d(self, items):
        total = sp.Integer(1)
        for item in items:
            species = self.species(item.species)
            total *= species.thermo.q3d ** item.coefficient
        return sp.simplify(total)

    def sticky(self, process, activation):
        ''' the sticky coefficient is evaluate as the reduction of degrees of freedom, i.e. from a 3D free molecule
        to a 2D trapped molecule moving parallel to the surface (being the third degree the reaction coordinate) '''
        qr = self.state_q3d(process.reactants)
        if process.ts:
            qts = self.state_q3d(process.ts)
            return sp.simplify(qts / qr)
        qts = sp.Integer(1)
        for item in process.reactants:
            species = self.species(item.species)
            if self.is_molecule(species):
                qts *= self.adsorption_qts(species, item.coefficient, activation,)
            else:
                qts *= species.thermo.q3d ** item.coefficient
        return sp.simplify(qts / qr)

    def adsorption_qts(self, species, coefficient, activation):
        thermo = species.thermo
        expr1 = thermo.q3d ** coefficient
        expr2 = (thermo.qrot * thermo.qelec * thermo.qtrans3d * thermo.qvib2d) ** coefficient
        expr3 = (thermo.qrot * thermo.qelec * thermo.qtrans2d * thermo.qvib2d) ** coefficient
        expr4 = (thermo.qrot * thermo.qelec * thermo.qvib2d) ** coefficient
        expr5 = (thermo.qelec * thermo.qvib2d) ** coefficient
        w1 = 1 - self.smooth_step(activation, 0.01, 0.01)
        w2 = self.smooth_step(activation, 0.01, 0.01) * (1 - self.smooth_step(activation, 0.25, 0.02))
        w3 = self.smooth_step(activation, 0.25, 0.02) * (1 - self.smooth_step(activation, 0.70, 0.05))
        w4 = self.smooth_step(activation, 0.70, 0.05) * (1 - self.smooth_step(activation, 1.00, 0.05))
        w5 = self.smooth_step(activation, 1.00, 0.05)
        return sp.simplify(w1 * expr1 + w2 * expr2 + w3 * expr3 + w4 * expr4 + w5 * expr5 *
                           sp.exp(-activation / (kb * temp * JtoeV)))

    @staticmethod
    def smooth_step(x, x0, width):
        return 1 / (1 + sp.exp(-(x - x0) / width))

    def arrhenius(self, process):
        if process.kind.upper() == "A":
            molecule = self.adsorbing_molecule(process)
            mass = molecule.mass
            area = self.adsorption_area(molecule)
            if mass is None:
                raise ValueError(f"{molecule.name}: molecular mass is required for adsorption.")
            if area is None:
                raise ValueError(f"{molecule.name}: adsorption area could not be resolved from the surface site information.")
            return sp.simplify(area / sp.sqrt(2 * sp.pi * mass * kb * temp))
        qts = self.transition_state_partition_function(process)
        qr = self.state_q3d(process.reactants)
        return sp.simplify(kb * temp / h * qts / qr)

    def transition_state_partition_function(self, process):
        if process.ts:
            return self.state_q3d(process.ts)
        qts = sp.Integer(1)
        for item in process.products:
            species = self.species(item.species)
            if self.is_molecule(species):
                qts *= species.thermo.q2d ** item.coefficient
            else:
                qts *= species.thermo.q3d ** item.coefficient
        return sp.simplify(qts)

    def adsorbing_molecule(self, process):
        for item in process.reactants:
            species = self.species(item.species)
            if self.is_molecule(species):
                return species
        raise ValueError(f"Process {process.id}: adsorption process has no molecular reactant.")

    def adsorption_area(self, molecule):
        if molecule.area is not None:
            return molecule.area
        molsite = molecule.molsite
        if isinstance(molsite, list):
            molsite = molsite[0] if molsite else None
        if molsite is None:
            return None
        for species in self.model.species.values():
            if species.kind.lower() != "surface":
                continue
            if species.sites == molsite:
                return species.area
        return None

    def tunneling(self, process):
        ''' Second order harmonic Wigner approach to shallow quantum tunneling valid for
        vast numbers of reaction including surface-catalysed --> DOI: 10.1039/C4CP03235G '''
        correction = sp.Integer(1)
        for item in process.ts:
            species = self.species(item.species)
            if hasattr(species, "ifreq") and species.ifreq is not None:
                correction *= (1 + sp.Rational(1, 24) * (hc * species.ifreq / (2 * sp.pi * kb * temp)) ** 2)
        return sp.simplify(correction)


class REquations:
    def __init__(self, model):
        self.model = model
        self.build()

    def build(self):
        constant_temperature: dict[str, list[sp.Expr]] = {}
        equation_factors: dict[str, list[float]] = {}
        tpd: dict[str, list[sp.Expr]] = {}
        tpd_factors: dict[str, list[float]] = {}
        for name in self.dynamic_species_names():
            constant_temperature[name] = []
            equation_factors[name] = []
            tpd[name] = []
            tpd_factors[name] = []
            for process in self.model.processes.values():
                equation, factor = self.equation(process, name)
                constant_temperature[name].append(equation)
                equation_factors[name].append(factor)
                if process.kind.upper() != "A":
                    tpd[name].append(equation)
                    tpd_factors[name].append(factor)
        eq = self.model.equations
        eq.isothermal = constant_temperature
        eq.equation_factors = equation_factors
        eq.tpd = tpd
        eq.tpd_factors = tpd_factors
        eq.overall_stoichiometry = self.overall_stoichiometry()
        return self.model

    def dynamic_species_names(self) -> list[str]:
        names = set()
        for process in self.model.processes.values():
            names.update(item.species for item in process.reactants)
            names.update(item.species for item in process.products)
        return sorted(name for name in names if not is_surface(self.model.species[name]))

    @staticmethod
    def coefficient_for(items, name: str) -> float:
        return sum(float(item.coefficient) for item in items if item.species == name)

    def equation(self, process, name: str) -> tuple[sp.Expr, float]:
        product_coeff = self.coefficient_for(process.products, name)
        reactant_coeff = self.coefficient_for(process.reactants, name)
        net_factor = product_coeff - reactant_coeff
        if net_factor == 0.0:
            return sp.Integer(0), 0.0
        kin = process_kinetics(self.model, process)
        if kin.krate0 == 0.0:
            raise ValueError(f"Process {process.id}: missing krate0. Run RConstants(model) before REquations(model).")
        rate = kin.krate0
        for item in process.reactants:
            rate *= sp.Symbol(item.species) ** item.coefficient
        return sp.simplify(net_factor * rate), net_factor

    def overall_stoichiometry(self) -> list[sp.Expr]:
        """ Return one overall stoichiometric vector that eliminates adsorbed species.

        No files are written. The returned vector follows the order of the dynamic
        species names used in the rate equations. """
        names = self.dynamic_species_names()
        if not names:
            return []
        factors = self.model.equations.equation_factors
        stoich_matrix = sp.Matrix([factors[name] for name in names])
        ads_matrix = sp.Matrix([factors[name] for name in names if is_adsorbed_like(self.model.species[name])])
        if ads_matrix.rows == 0:
            return [sp.Integer(0) for _ in names]
        nullspace = ads_matrix.nullspace()
        if not nullspace:
            return []
        v_int = self.normalise_integer_vector(nullspace[0])
        if any(value < 0 for value in v_int):
            v_int = -v_int
        return list(stoich_matrix * v_int)

    @staticmethod
    def normalise_integer_vector(vector):
        values = [sp.nsimplify(value) for value in vector]
        denominators = [value.as_numer_denom()[1] for value in values]
        def lcm(a, b):
            return abs(int(a) * int(b)) // gcd(int(a), int(b))
        lcm_denominator = reduce(lcm, denominators, 1)
        integers = [int(value * lcm_denominator) for value in values]
        nonzero = [abs(value) for value in integers if value != 0]
        if nonzero:
            common = reduce(gcd, nonzero)
            integers = [value // common for value in integers]
        return sp.Matrix(integers)


class Profile:
    """Build an in-memory free-energy profile table."""
    def __init__(self, model, temp_num: float = 300.0):
        self.model = model
        self.temp_num = temp_num
        self.cache: dict[tuple, float] = {}
        self.data = self.build()
        self.model.profile = self.data

    def build(self) -> list[list[str]]:
        data: list[list[str]] = [["Process", "Kind", "Reaction", f"G(T={self.temp_num} K) [eV]", "Ea [eV]", "Er [eV]"]]
        for process in self.model.processes.values():
            react_state = self.state_key(process.reactants)
            product_state = self.state_key(process.products)
            g_react = self.state_energy(react_state)
            g_product = self.state_energy(product_state)
            if process.ts:
                ts_state = self.state_key(process.ts)
                ts_label = self.state_label(ts_state)
                g_ts = self.state_energy(ts_state)
            else:
                ts_label = ""
                if process.kind.upper() == "A":
                    g_ts = self.state_energy_adsorption_ts(react_state)
                elif process.kind.upper() == "D":
                    g_ts = self.state_energy_adsorption_ts(product_state)
                else:
                    g_ts = None
            ea = g_ts - g_react if g_ts is not None else None
            er = g_product - g_react
            reaction = " > ".join([self.state_label(react_state), ts_label, self.state_label(product_state),])
            energy_path = " > ".join([f"{g_react:.2f}", f"{g_ts:.2f}" if g_ts is not None else "", f"{g_product:.2f}",])
            data.append([str(process.id), process.kind, reaction, energy_path,
                         f"{ea:.2f}" if ea is not None else "---", f"{er:.2f}",])
        return data

    def species_gibbs(self, name: str):
        species = self.model.species[name]
        thermo = get_thermo(species)
        if is_molecule(species):
            return thermo.gibbs3d
        return first_existing_attr(thermo, ["gibbs_theta", "gibbs3d"])

    def species_gibbs_adsorption_ts(self, name: str):
        species = self.model.species[name]
        thermo = get_thermo(species)
        if is_molecule(species):
            return thermo.gibbs2d
        return first_existing_attr(thermo, ["gibbs_theta", "gibbs3d"])

    def state_energy(self, state) -> float:
        if state in self.cache:
            return self.cache[state]
        value = 0.0
        for name, coefficient in state:
            expr = self.species_gibbs(name).subs(constants)
            value += coefficient * self.evaluate_at_temperature(expr)
        self.cache[state] = value
        return value

    def state_energy_adsorption_ts(self, state) -> float:
        value = 0.0
        for name, coefficient in state:
            expr = self.species_gibbs_adsorption_ts(name).subs(constants)
            value += coefficient * self.evaluate_at_temperature(expr)
        return value

    def evaluate_at_temperature(self, expr) -> float:
        func = sp.lambdify(temp, expr, ("numpy", "sympy"))
        return float(func(self.temp_num))

    @staticmethod
    def state_key(items):
        return tuple(sorted((item.species, float(item.coefficient)) for item in items))

    @staticmethod
    def state_label(state) -> str:
        pieces = []
        for name, coefficient in state:
            if float(coefficient).is_integer():
                coefficient = int(coefficient)
            pieces.append(f"{coefficient}.{name}")
        return " + ".join(pieces)

