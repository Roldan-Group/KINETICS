'''
Workflow:  Partition Functions  |   Energy          |     Interpolation
									 ├── Entropy
									 ├── Enthalpy
									 └── Gibbs
'''

from dataclasses import dataclass, field
import sympy as sp
import numpy as np
from Symbols_def import temp, h, kb, hc, cov, JtoeV, constants


@dataclass
class ThermodynamicProperties:
	# Partition functions
	qtrans3d: sp.Expr | float = 1.0
	qtrans2d: sp.Expr | float = 1.0
	qrot: sp.Expr | float = 1.0
	qelec: sp.Expr | float = 1.0
	qvib3d: sp.Expr | float = 1.0
	qvib2d: sp.Expr | float = 1.0
	# Total Partition functions
	q3d: sp.Expr | float = 1.0
	q2d: sp.Expr | float = 1.0
	# Entropy contributions
	strans3d: sp.Expr | float = 0.0
	strans2d: sp.Expr | float = 0.0
	srot: sp.Expr | float = 0.0
	svib3d: sp.Expr | float = 0.0
	svib2d: sp.Expr | float = 0.0
	selec: sp.Expr | float = 0.0
	# Total entropy
	sentropy3d: sp.Expr | float = 0.0
	sentropy2d: sp.Expr | float = 0.0
	# Enthalpy contributions
	zpe3d: sp.Expr | float = 0.0
	zpe2d: sp.Expr | float = 0.0
	htrans3d: sp.Expr | float = 0.0
	htrans2d: sp.Expr | float = 0.0
	hrot: sp.Expr | float = 0.0
	hvib3d: sp.Expr | float = 0.0
	hvib2d: sp.Expr | float = 0.0
	# Heat capacities
	cp3d: sp.Expr | float = 0.0
	cp2d: sp.Expr | float = 0.0
	hthermal3d: sp.Expr | float = 0.0
	hthermal2d: sp.Expr | float = 0.0
	# Thermodynamic functions
	enthalpy3d: sp.Expr | float = 0.0
	enthalpy2d: sp.Expr | float = 0.0
	# Gibbs Energy
	gibbs3d: sp.Expr | float = 0.0
	gibbs2d: sp.Expr | float = 0.0
	# Coverage interpolation
	q3d_theta: sp.Expr | None = None
	gibbs3d_theta: sp.Expr | None = None


class PartitionFunctions:
	def __init__(self, model):
		self.model = model
		self.build()

	def build(self):
		for species in self.model.species.values():
			self.build_species(species)
		return self.model

	def build_species(self, species):
		if not hasattr(species, "thermo") or species.thermo is None:
			species.thermo = ThermodynamicProperties()
		thermo = species.thermo
		kind = species.kind.lower()
		thermo.qelec = self.qelec(species)
		if kind == "molecule":
			thermo.qtrans3d = self.qtrans3d(species)
			thermo.qtrans2d = self.qtrans2d(species)
			thermo.qrot = self.qrot(species)
			thermo.qvib3d = self.qvib(species.frequencies)
			thermo.qvib2d = self.qvib(species.frequencies_2d)
			thermo.q3d = thermo.qtrans3d * thermo.qrot * thermo.qelec * thermo.qvib3d
			thermo.q2d = thermo.qtrans2d * thermo.qrot * thermo.qelec * thermo.qvib2d
		elif kind in {"surface", "adsorbate"}:
			thermo.qtrans3d = sp.Integer(1)
			thermo.qrot = sp.Integer(1)
			thermo.qvib3d = self.qvib(species.frequencies)
			thermo.q2d = 0
			thermo.q3d = thermo.qelec * thermo.qvib3d
		else:
			raise ValueError(f"{species.name}: unknown species kind {species.kind}")

	@staticmethod
	def qtrans3d(species):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		Adsorption Journal Of The International Adsorption Society
		(Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
		page 88		 PER m^3'''
		return ((2 * sp.pi * species.mass * kb * temp) ** sp.Rational(3, 2)/ h**3)

	@staticmethod
	def qtrans2d(species):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		Adsorption Journal Of The International Adsorption Society
		(Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
		page 88		  PER m^2 '''
		return (2 * sp.pi * species.mass * kb * temp / h**2)

	@staticmethod
	def qrot(species):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		Adsorption Journal Of The International Adsorption Society
		(Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.		 page 90
		approximation that omits nuclear spin and assumes many states are occupied, i.e. no nuclear spin included'''
		inertia = np.asarray(species.inertia, dtype=float)
		if species.linear:
			moment = max(inertia)
			return (8 * sp.pi**2 * moment * kb * temp / (species.symmetry_factor * h**2))
		return (sp.sqrt(sp.pi) / species.symmetry_factor * ((8 * sp.pi**2 * kb * temp) / h**2) ** sp.Rational(3, 2)
				* sp.sqrt(np.prod(inertia)))

	@staticmethod
	def qvib(freqs):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		Adsorption Journal Of The International Adsorption Society
		(Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.		 page 89
		---> The equation used is respect the lowest occupied state, not the bottom of the potential energy curve.
		It also consider small frequencies, when h*v ~ kb*T (v=frequencies).
		The Zero Point Energy should be added to the energy as this qvib is to calculate the entropy (S),
		specific head (Cp), and pre-exponential factor of Arrhenius (A). '''
		positive_freqs = np.asarray([freq for freq in freqs if freq > 0.0], dtype=float,)
		if positive_freqs.size == 0:
			return sp.Integer(1)
		x_values = hc * positive_freqs / (kb * temp)
		return sp.prod(1 / (1 - sp.exp(-x)) for x in x_values)

	@staticmethod
	def qelec(species):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		Adsorption Journal Of The International Adsorption Society
		(Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.		 page 92
		>>> It is consider that the contribution of excited states is negligible at the working temperatures '''
		return sp.Integer(species.degeneracy)


class Energy:
	def __init__(self, model):
		self.model = model
		self.build()

	def build(self):
		print("\t... Generating entropic contribution ...")
		Entropy(self.model)
		print("\t... Generating enthalpic contribution ...")
		Enthalpy(self.model)
		print("\t... Generating Gibbs free energies ...")
		self.build_gibbs()
		print("\t... Generating coverage interpolation ...")
		Interpolation(self.model)
		return self.model

	def build_gibbs(self):
		for species in self.model.species.values():
			thermo = species.thermo
			thermo.gibbs3d = thermo.enthalpy3d - temp * thermo.sentropy3d
			thermo.gibbs2d = thermo.enthalpy2d - temp * thermo.sentropy2d


class Entropy:
	''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
	Adsorption Journal Of The International Adsorption Society
	(Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
	page 88 :: S = kb*ln(Q)+kb*T*diff(ln(Q), T)'''
	def __init__(self, model):
		self.model = model
		self.build()
	def build(self):
		for species in self.model.species.values():
			thermo = species.thermo
			kind = species.kind.lower()
			thermo.selec = self.electronic_entropy(species)
			if kind == "molecule":
				thermo.strans3d = self.translational_entropy_3d(species)
				thermo.strans2d = self.translational_entropy_2d(species)
				thermo.srot = self.rotational_entropy(species)
				thermo.svib3d = self.vibrational_entropy(species.frequencies)
				thermo.svib2d = self.vibrational_entropy(species.frequencies_2d)
			elif kind in {"surface", "adsorbate"}:
				thermo.strans3d = sp.Integer(0)
				thermo.strans2d = sp.Integer(0)
				thermo.srot = sp.Integer(0)
				thermo.svib3d = self.vibrational_entropy(species.frequencies)
				thermo.svib2d = self.vibrational_entropy(species.frequencies_2d)
			else:
				raise ValueError(f"{species.name}: unknown species kind {species.kind}")
			thermo.sentropy3d = thermo.strans3d + thermo.srot + thermo.svib3d + thermo.selec
			thermo.sentropy2d = thermo.strans2d + thermo.srot + thermo.svib2d + thermo.selec

	@staticmethod
	def positive_frequencies(freqs):
		return np.asarray([freq for freq in freqs if freq > 0.0], dtype=float,)

	@staticmethod
	def vibrational_entropy(freqs):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
		 harmonic oscillator '''
		freqs = Entropy.positive_frequencies(freqs)
		if freqs.size == 0:
			return sp.Integer(0)
		x_values = hc * freqs / (kb * temp)
		return kb * sp.Add(*[x / (sp.exp(x) - 1) - sp.log(1 - sp.exp(-x)) for x in x_values]) * JtoeV

	@staticmethod
	def translational_entropy_3d(species):
		'''re-Formulation from explicit derivatives::
		https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
		(Note that the translational component also includes components from the Stirling approximation)
		the Sackur-Tetrode equation for ideal gases'''
		q = species.thermo.qtrans3d
		return kb * (sp.log(q) + sp.Rational(5, 2)) * JtoeV

	@staticmethod
	def translational_entropy_2d(species):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
		 (Note that the translational component also includes components from the Stirling approximation)'''
		q = species.thermo.qtrans2d
		return kb * (sp.log(q) + 2) * JtoeV

	@staticmethod
	def rotational_entropy(species):
		'''re-Formulation from explicit derivatives::
		https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
		qrot = species.thermo.qrot
		return kb * (sp.log(qrot) + 1) * JtoeV

	@staticmethod
	def electronic_entropy(species):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
		return kb * sp.log(species.degeneracy) * JtoeV


class Enthalpy:
	""" C.J. Cramer. Essentials of Computational Chemistry, Second Edition. Wiley, 2004.
	and Raymand Chang, "PHYSICAL CHEMISTRY for Chemical and Biological Science", ISBN: 1-891389-06-8,
		page 91 :: H=[dCp/dT](T1,T2)    	H =  E + ZPE + integral(Cp, 0 --> T) """
	def __init__(self, model):
		self.model = model
		self.build()

	def build(self):
		for species in self.model.species.values():
			self.build_species(species)

	def build_species(self, species):
		thermo = species.thermo
		kind = species.kind.lower()
		thermo.zpe3d = self.zero_point_energy(species.frequencies)
		thermo.zpe2d = self.zero_point_energy(species.frequencies_2d)
		thermo.cp3d = self.heat_capacity_3d(species)
		thermo.cp2d = self.heat_capacity_2d(species)
		thermo.hthermal3d = self.hthermal3d(species)
		thermo.hthermal2d = self.hthermal2d(species)
		thermo.enthalpy3d = species.energy0 + thermo.zpe3d + thermo.hthermal3d
		thermo.enthalpy2d = species.energy0 + thermo.zpe2d + thermo.hthermal2d

	@staticmethod
	def zero_point_energy(freqs):
		positive_freqs = np.asarray([freq for freq in freqs if freq > 0.0], dtype=float,)
		if positive_freqs.size == 0:
			return sp.Integer(0)
		return (sp.Rational(1, 2) * hc * np.sum(positive_freqs) * JtoeV)

	def hthermal3d(self, species):
		kind = species.kind.lower()
		if kind == "molecule":
			return (self.htrans3d()	+ self.hrot(species) + self.hvib(species.frequencies))
		if kind in {"surface", "adsorbate"}:
			return self.hvib(species.frequencies)
		return sp.Integer(0)

	def hthermal2d(self, species):
		kind = species.kind.lower()
		if kind == "molecule":
			return (self.htrans2d()	+ self.hrot(species) + self.hvib(species.frequencies_2d))
		return sp.Integer(0)

	@staticmethod
	def hvib(freqs):
		total = sp.Integer(0)
		for freq in freqs:
			if freq > 0.0:
				x = hc * freq / (kb * temp)
				total += hc * freq * sp.exp(-x) / (1 - sp.exp(-x))
		return total * JtoeV

	@staticmethod
	def htrans3d():
		return sp.Rational(5, 2) * kb * temp * JtoeV

	@staticmethod
	def htrans2d():
		return sp.Rational(3, 2) * kb * temp * JtoeV

	@staticmethod
	def hrot(species):
		if species.linear:
			return kb * temp * JtoeV
		return sp.Rational(3, 2) * kb * temp * JtoeV

	@staticmethod
	def heat_capacity_3d(species):
		'''Hans Kuhn, Horst-Dieter Försterling, David Hennessey Waldeck, "Principles of Physical Chemistry"
		ISBN: 9780470089644, page 551 :: Cp=T*[dS/dT](N,P)
		Also, Cp = Cp_trans + Cp_rot + Cp_vib
		Later on to calculate H[T], Cp has to be integrated, which is very demanding in resources.
		For that reason, it is done by steps'''
		kind = species.kind.lower()
		cp = sp.Integer(0)
		if kind == "molecule":
			cp += Enthalpy.translational_cp_3d()
			cp += Enthalpy.rotational_cp(species)
			cp += Enthalpy.vibrational_cp(species.frequencies)
		elif kind in {"surface", "adsorbate"}:
			cp += Enthalpy.vibrational_cp(species.frequencies)
		else:
			raise ValueError(f"{species.name}: unknown species kind '{species.kind}'")
		return cp

	@staticmethod
	def heat_capacity_2d(species):
		'''Hans Kuhn, Horst-Dieter Försterling, David Hennessey Waldeck, "Principles of Physical Chemistry"
		ISBN: 9780470089644, page 551 :: Cp=T*[dS/dT](N,P)
		Also, Cp = Cp_trans + Cp_rot + Cp_vib
		Later on to calculate H[T], Cp has to be integrated, which is very demanding in resources.
		For that reason, it is done by steps'''
		kind = species.kind.lower()
		cp = sp.Integer(0)
		if kind == "molecule":
			cp += Enthalpy.translational_cp_2d()
			cp += Enthalpy.rotational_cp(species)
			cp += Enthalpy.vibrational_cp(species.frequencies_2d)
		elif kind in {"surface", "adsorbate"}:
			cp += Enthalpy.vibrational_cp(species.frequencies_2d)
		else:
			raise ValueError(f"{species.name}: unknown species kind '{species.kind}'" )
		return cp

	@staticmethod
	def translational_cp_3d():
		return sp.Rational(5, 2) * kb * JtoeV

	@staticmethod
	def translational_cp_2d():
		return kb * JtoeV

	@staticmethod
	def rotational_cp(species):
		if species.linear:
			return kb
		return sp.Rational(3, 2) * kb * JtoeV

	@staticmethod
	def vibrational_cp(freqs):
		positive_freqs = np.asarray([freq for freq in freqs if freq > 0.0], dtype=float,)
		if positive_freqs.size == 0:
			return sp.Integer(0)
		x_values = hc * positive_freqs / (kb * temp)
		return kb * sp.Add(*[x**2 * sp.exp(x) / (sp.exp(x) - 1)**2 for x in x_values ]) * JtoeV


class Interpolation:
	def __init__(self, model):
		self.model = model
		self.build()

	def build(self):
		for species in self.model.species.values():
			self.build_species(species)

	def build_species(self, species):
		thermo = species.thermo
		if not species.coverage_model:
			thermo.q3d_theta = thermo.q3d
			thermo.gibbs3d_theta = thermo.gibbs3d
			return
		points = sorted(species.coverage_model, key=lambda point: point.coverage,)
		theta_values = []
		q3d_values = []
		gibbs3d_values = []
		for point in points:
			theta_values.append(point.coverage)
			q3d, gibbs3d = self.build_point(species, point)
			q3d_values.append(q3d)
			gibbs3d_values.append(gibbs3d)
		thermo.q3d_theta = self.interpolate(theta_values, q3d_values)
		thermo.gibbs3d_theta = self.interpolate(theta_values, gibbs3d_values)

	def build_point(self, point):
		qvib = PartitionFunctions.qvib(point.frequencies)
		zpe = Enthalpy.zero_point_energy(point.frequencies)
		cpvib = Enthalpy.vibrational_cp(point.frequencies)
		hvib = sp.integrate(cpvib, temp)
		svib = Entropy.vibrational_entropy(point.frequencies)
		q = qvib
		enthalpy = point.energy0 + zpe + hvib
		gibbs = enthalpy - temp * svib
		return q, gibbs

	@staticmethod
	def interpolate(theta_values, values):
		if len(theta_values) == 1:
			return values[0]
		return sp.interpolate(list(zip(theta_values, values)), cov,)