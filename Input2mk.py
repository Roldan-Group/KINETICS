'''

	by Alberto Roldan:	04/2020
    Updated A. Roldan:  06/2026
	Reads information from outputs to generate a input file for microkinetics

- dataclass architecture
- VASP + FHI-aims support via ASE
- preserves MOLSITE, ISITES, IPRESSURE, ICOVERAGE, IACAT
- global keyword support (TEMPERATURE, TIME, PH, ELECTROCHEMICAL)
- reaction/system containers
- dipole parsing hooks and IR intensity support
'''

from __future__ import annotations
import sys
import re
import numpy as np
from dataclasses import dataclass, field
from collections import defaultdict, Counter
from pathlib import Path
from ase import Atoms
from ase.io import read
from ase.io import vasp


# --- Dataclasses
@dataclass
class GlobalSettings:
    temperature: float | None = None
    time: float | None = None
    ph: float | None = None
    electrochemical: float | None = None


@dataclass
class Reaction:
    process_type: str      # A, D, R
    expression: str


@dataclass
class SystemDefinition:
    name: str
    syspath: str | None = None
    freqpath: str | None = None
    symfactor: float | None = None
    molsite: tuple[float, str] | None = None
    ipressure: list[float] | None = None
    icoverage: list[float] | None = None
    nsites: float | None = None
    sites: str | None = None
    system_type: str | None = None    # classification


@dataclass
class InputModel:
    globals: GlobalSettings
    reactions: list[Reaction]
    systems: list[SystemDefinition]


@dataclass
class ParsedSystem:
    energy: float    # thermochemistry
    frequencies: list[float] = field(default_factory=list)    # vibrational data
    frequencies_2d: list[float] = field(default_factory=list)
    ir_intensities: list[float] = field(default_factory=list)    # IR data
    degeneracy: int = 1
    symfactor: int | None = None
    linear: str | None = None
    inertia: list[float] = field(default_factory=list)
    masses: list[float] = field(default_factory=list)
    natoms_per_mass: list[int] = field(default_factory=list)
    area: float | None = None    # surface properties
    system_type: str | None = None


@dataclass
class DisplacedGeometry:
    positions: np.ndarray
    dipole: np.ndarray | None = None
    energy: float | None = None


# --- Parser
def parse_global(line, globals_):
    key, value = [x.strip() for x in line.split("=",1)]
    value = remove_comment(value)
    values = [float(v) for v in value.split()]
    if key == "TEMPERATURE":
        globals_.temperature = values
    elif key == "TIME":
        globals_.time = values
    elif key == "PH":
        globals_.ph = values
    elif key == "ELECTROCHEMICAL":
        globals_.electrochemical = values

def parse_process(line):
    rhs = line.split("=",1)[1].strip()
    process_type = rhs.split()[0]
    expression = rhs[len(process_type):].strip()
    return Reaction(process_type=process_type, expression=expression)

def parse_system_keyword(system, line):
    line = remove_comment(line)
    if "=" not in line:
        return
    key, value = [x.strip() for x in line.split("=", 1)]
    if key == "SYSPATH":
        system.syspath = value
    elif key == "FREQPATH":
        system.freqpath = value
    elif key == "SYMFACTOR":
        system.symfactor = float(value)
    elif key == "MOLSITE":
        words = value.split()
        if len(words) == 1:
            system.molsite = (1.0, words[0])
        else:
            system.molsite = (float(words[0]), words[1])
    elif key == "ISITES":
        words = value.split()
        site_words = []
        system.nsites = None
        for word in words:
            try:
                system.nsites = float(word)
            except ValueError:
                site_words.append(word)
        system.sites = " ".join(site_words) if site_words else None
    elif key == "IPRESSURE":
        system.ipressure = [float(v) for v in value.split()]
    elif key == "ICOVERAGE":
        system.icoverage = [float(v) for v in value.split()]

    if system.nsites is None and system.molsite is None:
        system.system_type = "Surface"
    elif isinstance(system.nsites, (float, int)):
        system.system_type = "Catalyst"
    else:
        system.system_type = "Molecule"

def parse_input(path):
    globals_ = GlobalSettings()
    reactions = []
    systems = []
    current = None
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            line = line.replace("##","").strip()
            if line.startswith("TEMPERATURE"):
                parse_global(line, globals_)
            elif line.startswith("TIME"):
                parse_global(line, globals_)
            elif line.startswith("PH"):
                parse_global(line, globals_)
            elif line.startswith("ELECTROCHEMICAL"):
                parse_global(line, globals_)
            elif line.startswith("PROCESS"):
                reactions.append(parse_process(line))
            elif line.startswith("SYSTEM"):
                name = line.split("=",1)[1].strip()
                current = SystemDefinition(name=name)
            elif line.upper() == "END":
                systems.append(current)
                current = None
            elif current is not None:
                parse_system_keyword(current,line)
    return InputModel(globals=globals_, reactions=reactions, systems=systems)


# --- Readers
def get_reader(path):
    text = open(path).read(10000)
    if "FHI-aims" in text:
        return FHIAimsReader()
    if "vasp" in text.lower():
        return VaspReader()
    raise RuntimeError(f"Unknown code: {path}")

class BaseReader:
    def detect(self, path: str) -> str:
        txt = Path(path).read_text(errors="ignore")[:10000]
        if "FHI-aims" in txt:
            return "FHI-aims"
        if "vasp" in txt.lower() or "OUTCAR" in txt:
            return "VASP"
        raise ValueError(f"Unknown code for {path}")


class VaspReader(BaseReader):
    def read_system(self, path: str):
        return vasp.read_vasp_out(path, index=-1)

    def read_geometries(self, outcar):
        geometries = []
        with open(outcar, errors="ignore") as fh:
            lines = fh.readlines()
        i = 0
        while i < len(lines):
            if "POSITION" in lines[i] and "TOTAL-FORCE" in lines[i]:
                i += 2
                positions = []
                while i < len(lines):
                    line = lines[i]
                    if not line.strip():
                        break
                    words = line.split()
                    if len(words) < 3:
                        break
                    try:
                        positions.append([float(words[0]), float(words[1]), float(words[2]),])
                    except ValueError:
                        break
                    i += 1
                geometries.append(np.asarray(positions))
            i += 1
        return geometries

    def extract_dipoles(self, path: str):
        dipoles = []
        current_dipole = None
        with open(path, errors="ignore") as f:
            for line in f:
                if "dipolmoment" in line:
                    tail = line.lower().split("dipolmoment", 1)[1]
                    nums = re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", tail)
                    if len(nums) >= 3:
                        current_dipole = np.array([float(nums[0]), float(nums[1]), float(nums[2])])
                elif "free  energy   TOTEN" in line:
                    if current_dipole is not None:
                        dipoles.append(current_dipole)
                        current_dipole = None
        return dipoles

    def read_displaced_geometries(self, outcar):
        geometries = self.read_geometries(outcar)
        dipoles = self.extract_dipoles(outcar)
        n = min(len(geometries), len(dipoles))
        result = []
        for geom, dip in zip(geometries[:n], dipoles[:n]):
            result.append(DisplacedGeometry(positions=geom, dipole=dip))
        return result

    def read_frequencies(self, freqfile):
        freq = []
        with open(freqfile) as fh:
            lines = fh.readlines()
        nline = 0
        while nline < len(lines):
            words = lines[nline].split()
            nline += 1
            if len(words) >7 and words[1] == "f" and "cm-1" in words:
                freq.append(float(words[7]))
            elif len(words) >7 and words[1] == "f/i=" and "cm-1" in words:
                if float(words[6]) < 200:
                    freq.append(float(words[6]))
                else:
                    freq.append(-float(words[6]))
        #freq.sort(reverse=True)    # Not sorted, making sure these are align with the intensities
        return freq

    def read_eigenvectors(self, outcar):
        with open(outcar, errors="ignore") as fh:
            lines = fh.readlines()
        eigenvectors = []
        i = 0
        while i < len(lines):
            line = lines[i]
            if ("cm-1" in line and ("f" in line or "f/i" in line)):
                vec = []
                i += 2
                while i < len(lines):
                    words = lines[i].split()
                    if len(words) < 6:
                        break
                    try:
                        vec.extend([float(words[-3]), float(words[-2]), float(words[-1]),])
                    except ValueError:
                        break
                    i += 1
                eigenvectors.append(np.asarray(vec))
            i += 1
        return np.asarray(eigenvectors)

    def read_displacement_step(self, outcar):
        with open(outcar, errors="ignore") as fh:
            for line in fh:
                if "POTIM" in line:
                    m = re.search(r"POTIM\s*=\s*([-\d.]+)", line)
                    if m:
                        return float(m.group(1))
        return 0.05


class FHIAimsReader(BaseReader):
    def read_system(self, path: str):
        return read(path)

    def extract_dipole(self, path: str):
        text = Path(path).read_text(errors="ignore")
        m = re.findall(r"Total dipole moment\s*:\s*([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)", text,)
        if not m:
            return None
        return np.asarray(m[-1], dtype=float)

    def read_frequencies(self, freqfile):
        frequencies = []
        with open(freqfile) as fh:
            for line in fh:
                words = line.split()
                if (len(words) > 4 and words[0] == "Mode"):
                    freq = words[4]
                    if freq.endswith("i"):
                        freq = -float(freq[:-1])
                    else:
                        freq = float(freq)
                    frequencies.append(freq)
        frequencies.sort(reverse=True)
        return frequencies

# --- Writer
def write_globals(fh, globals_):
    if globals_.temperature is not None:
        fh.write("TEMPERATURE = " + " ".join(map(str, globals_.temperature)) + "\n")
    if globals_.time is not None:
        fh.write("TIME = " + " ".join(map(str, globals_.time)) + "\n")
    if globals_.electrochemical is not None:
        fh.write("ELECTROCHEMICAL = " + " ".join(map(str, globals_.electrochemical)) + "\n")
    if globals_.ph is not None:
        fh.write("PH = " + " ".join(map(str, globals_.ph)) + "\n")
    fh.write("\n")

def write_reactions(fh, reactions):
    fh.write("# Reaction processes\n")
    for r in reactions:
        fh.write(f"PROCESS = {r.process_type} {r.expression}\n")
    fh.write("\n")

def write_system(fh, system: SystemDefinition, parsed: ParsedSystem):
    fh.write(f"SYSTEM = {system.name}\n")
    fh.write(f" TYPE = {system.system_type}\n")
    fh.write(f" SYSPATH = {system.syspath}\n")
    if system.freqpath:
        fh.write(f" FREQPATH = {system.freqpath}\n")
    if parsed.energy is not None:
        fh.write(f" E0 = {parsed.energy:.6f}\t# eV\n")
    if parsed.frequencies:
        fh.write(" FREQ =")
        for freq in parsed.frequencies:
            fh.write(f" {freq:.1f}")
        fh.write("\t# cm^-1\n")
    if parsed.ir_intensities:
        fh.write(" IRINTENSITY =")
        for i in parsed.ir_intensities:
            fh.write(f" {i:.5f}")
        fh.write("\n")
    if parsed.degeneracy is not None:
        fh.write(f" DEGENERATION = {parsed.degeneracy}\n")
    if system.system_type == "Molecule":    # ONLY for Moleccules
        if len(parsed.frequencies_2d) > 0:
            fh.write(" FREQ2D =")
            for freq in parsed.frequencies_2d:
                fh.write(f" {freq:.1f}")
            fh.write("\t# cm^-1\n")
        else:
            fh.write(" FREQ2D = 0.0\t# no 2D freqs\n")
        if system.symfactor:
            fh.write(f" SYMFACTOR = {system.symfactor:g}\n")
        if parsed.masses:
            fh.write(" IMASS =")
            for m in parsed.masses:
                fh.write(f" {m:.3f}")
            fh.write("\n")
        if parsed.natoms_per_mass:
            fh.write(" INATOMS =")
            for n in parsed.natoms_per_mass:
                fh.write(f" {n}")
            fh.write("\n")
        if parsed.linear:
            fh.write(f" LINEAR = {parsed.linear}\n")
        if parsed.inertia:
            fh.write(" INERTIA =")
            for i in parsed.inertia:
                fh.write(f" {i:.5e}")
            fh.write("\t# kg m^2\n")
        if system.molsite:
            nsites, site = system.molsite
            fh.write(f" MOLSITE = {nsites} {site}\n")
        if system.ipressure is not None:
            fh.write(" IPRESSURE = " + " ".join(map(str, system.ipressure)) + "\t# Pa\n")
        else:
            fh.write(" IPRESSURE = 0.0\t# Pa\n")
    elif system.system_type == "Surface":   # ONLY for Surfaces
        if system.sites:
            fh.write(f" ISITES = {system.sites}\n")
        if parsed.area is not None:
            fh.write(f" IACAT = {parsed.area:.6e}\t# m^2\n")
    elif system.system_type == "Catalyst":  # ONLY for Catalysts
        if system.nsites:
            fh.write(f" ISITES = {system.nsites} {system.sites}\n")
        if system.icoverage is not None:
            fh.write(" ICOVERAGE = " + " ".join(map(str, system.icoverage)) + "\t# ML\n")
        else:
            fh.write(" ICOVERAGE = 0.0\t# ML\n")
    fh.write("END\n\n")

def write_mk(outfile, model, results):
    with open(outfile, "w") as fh:
        write_globals(fh, model.globals)
        write_reactions(fh, model.reactions)
        for system in model.systems:
            parsed = results[system.name]
            write_system(fh, system, parsed)

# --- Frequencies and Others
def remove_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()

def molecular_data(atoms, tol=1e-3):
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)
    masses = []
    natoms = []
    for el, n in sorted(counts.items()):
        idx = symbols.index(el)
        masses.append(atoms[idx].mass)
        natoms.append(n)
    if len(atoms) < 3:
        linear = 'yes'
    pos = atoms.get_positions()
    r0 = pos[1] - pos[0]
    for i in range(2, len(pos)):
        ri = pos[i] - pos[0]
        if np.linalg.norm(np.cross(r0, ri)) > tol:
            linear = 'no'
    inertia = atoms.get_moments_of_inertia() * 1.66053906660e-47    # in kg m^2
    return masses, natoms, linear, inertia.tolist()

def surface_area(atoms):
    cell = atoms.get_cell()
    area_ang2 = np.linalg.norm(np.cross(cell[0], cell[1]))
    area_m2 = area_ang2 * 1.0e-20
    return area_m2

def find_equilibrium(displacements):
    return displacements[0]

def displacement_signature(positions, reference, tol=1e-6):
    delta = positions - reference
    idx = np.argwhere(np.abs(delta) > tol)
    if len(idx) != 1:
        return None
    atom, coord = idx[0]
    sign = np.sign(delta[atom, coord])
    return atom, coord, int(sign)

def pair_displacements(displaced_geometries, equilibrium):
    pairs = defaultdict(dict)
    ref = equilibrium.positions
    for geom in displaced_geometries:
        sig = displacement_signature(geom.positions, ref)
        if sig is None:
            continue
        atom, coord, sign = sig
        key = (atom, coord)
        pairs[key][sign] = geom
    return pairs

def cartesian_dipole_derivatives(pairs, delta):
    derivatives = {}
    for key, pair in pairs.items():
        if +1 not in pair:
            continue
        if -1 not in pair:
            continue
        mu_plus = pair[+1].dipole
        mu_minus = pair[-1].dipole
        derivatives[key] = (mu_plus - mu_minus) / (2.0 * delta)
    return derivatives

def derivative_matrix(derivatives, natoms):
    result = []
    for atom in range(natoms):
        for coord in range(3):
            result.append(derivatives.get((atom,coord), np.zeros(3)))
    return np.array(result)

def mode_dipole_derivatives(derivative_matrix, eigenvectors):
    mode_derivs = []
    for mode in eigenvectors:
        dmu = np.zeros(3)
        for i in range(len(mode)):
            dmu += (mode[i] * derivative_matrix[i])
        mode_derivs.append(dmu)
    return np.array(mode_derivs)

def ir_intensities(derivative_matrix, eigenvectors):
    mode_derivs = mode_dipole_derivatives(derivative_matrix, eigenvectors)
    intensities = []
    for dmu in mode_derivs:
        intensities.append(float(np.dot(dmu,dmu)))
    return intensities

def normalize_ir(intensities):
    intensities = np.asarray(intensities, dtype=float)
    if len(intensities) == 0:
        return []
    imax = np.max(intensities)
    if imax > 0:
        intensities /= imax
    return intensities.tolist()

def compute_ir_intensities(freqfile, reader, eigenvectors):
    displaced = reader.read_displaced_geometries(freqfile)
    equilibrium = displaced[0]
    delta = reader.read_displacement_step(freqfile)
    pairs = pair_displacements(displaced[1:], equilibrium)
    derivs = cartesian_dipole_derivatives(pairs, delta)
    dmu_dx = derivative_matrix(derivs, len(equilibrium.positions))
    intens = ir_intensities(dmu_dx, eigenvectors)
    return intens

def get_frequencies_2d(frequencies, eigenvectors, threshold=0.25):
    """ Classify modes according to their displacement direction.
    threshold:     fraction of displacement in x/y required to be considered 2D. """
    frequencies_2d = []
    for i in range(min(len(frequencies), len(eigenvectors))):
        freq = float(frequencies[i])
        mode = np.asarray(eigenvectors[i]).reshape(-1, 3)
        xy = np.sum(mode[:, 0]**2 + mode[:, 1]**2)
        z  = np.sum(mode[:, 2]**2)
        frac_xy = xy / (xy + z)
        if frac_xy > threshold:
            frequencies_2d.append(freq)
    return frequencies_2d

# --- Main
def main(inputfile):
    model = parse_input(inputfile)
    results = {}
    for system in model.systems:
        if system.syspath is None:
            continue
        reader = get_reader(system.syspath)
        atoms = reader.read_system(system.syspath)
        freqfile = system.freqpath or system.syspath
        frequencies = reader.read_frequencies(freqfile)
        eigenvectors = reader.read_eigenvectors(freqfile)
        ir_intensities = compute_ir_intensities(freqfile, reader, eigenvectors)

        if system.system_type == "Molecule":
            masses, natoms, linear, inertia = molecular_data(atoms)
            if linear == 'yes':
                nfreq = 3 * sum(natoms) - 5
            else:
                nfreq = 3 * sum(natoms) - 6
            frequencies = frequencies[:nfreq]
            ir_intensities = normalize_ir(ir_intensities[:nfreq])
            frequencies_2d = get_frequencies_2d(frequencies, eigenvectors[:nfreq])

            results[system.name] = ParsedSystem(energy=atoms.get_total_energy(),
                                                frequencies=frequencies,
                                                frequencies_2d=frequencies_2d,
                                                ir_intensities=ir_intensities,
                                                masses=masses,
                                                natoms_per_mass=natoms,
                                                linear=linear,
                                                inertia=inertia,
                                                )
        else:
            area = None
            if system.system_type == "Surface":
                area = surface_area(atoms)
            ir_intensities = normalize_ir(ir_intensities)
            results[system.name] = ParsedSystem(energy=atoms.get_total_energy(),
                                                frequencies=frequencies,
                                                ir_intensities=ir_intensities,
                                                area=area,
                                                )
    output = inputfile + ".mk.in"
    write_mk(output, model, results)
    print(f"Wrote {output}")

if __name__ == "__main__":
    main(sys.argv[1])
