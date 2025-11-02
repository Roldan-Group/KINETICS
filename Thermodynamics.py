"""
	This script builds on the perl version by A.Roldan.

"""

import pathlib
import time
import sympy as sp
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import ticker


def printdata(rconditions, name, nadsorbates, properties, datalabel, dataname):
	data = getdata(rconditions, properties, datalabel)
	maxlen = [max([len(f"{data[r][c]}")+1 for r in range(len(data))]) for c in range(len(data[0]))] # max length per
	# column
	folder = './THERMODYNAMICS/DATA/'+ name + "/" + nadsorbates
	outputfile = folder + "/" + str(dataname) + ".dat"
	if not pathlib.Path(folder).exists():
		pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
	output = open(outputfile, "w")
	for i in range(len(data[0])):
		output.write(" {val:>{wid}s}".format(wid=maxlen[i], val=data[0][i]))  # headings
	output.write("\n")
	for row in data[1:]:
		for i in range(len(row)):
			output.write(" {val:>{wid}.3{c}}".format(wid=maxlen[i], val=row[i],
													 c='e' if row[i] > 1e3 or np.abs(row[i]) < 1e-2 else 'f'))
		output.write("\n")
	output.close()


def getdata(rconditions, properties, datalabel):
	data = []
	headings = ["# Temperature[K]"]
	for i in datalabel:
		units = ''
		if i.startswith('enthalpy') or i.startswith('zpe') or i.startswith('energ'):
			units = '[eV]'  # energy and ZPE
		elif i.startswith('sentropy'):
			units = '[eV K^-1]'  # entropy or specific heat as "defined here"
		headings.append(i + units)
	data.append(headings)
	temp = sp.symbols("temperature", positive=True, real=True)
	if isinstance(rconditions["temperature"], float):
		row = [rconditions["temperature"]]
		for i in datalabel:
			row.append(round(float(sp.lambdify(temp, properties[str(i)])(float(rconditions["temperature"]))), 6))
		data.append(row)
	else:
		ramp = [float(i) for i in rconditions["temperature"]]
		for t in np.arange(ramp[0], ramp[1], ramp[2]):
			row = [t]
			for i in datalabel:
				row.append(round(float(sp.lambdify(temp, properties[str(i)], "numpy")(t)), 6))
			data.append(row)
	return data


def interpolate(rconditions, systems, name, restricted_arg, ykey):
	startint = time.time()
	folder = './THERMODYNAMICS/DATA/'+ name + "/"
	''' Reaction conditions are set as symbols using SYMPY '''
	temp, cov = sp.symbols("temperature coverage", positive=True, real=True)
	x0 = [] # list of nadsorbates
	x = []  # independent coordinate, e.g. coverage
	y = []  # dependent coordinate, e.g. energy
	for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
		if nadsorbates not in restricted_arg:       # only for nadsorbates
			x0.append(nadsorbates)
			x.append(float(nadsorbates))
			y.append(systems[name][nadsorbates][ykey])
	function = sp.interpolate(list(zip(x, y)), cov)    # generates the interpolation function with the x and y.
	''' generate a plot for the interpolated data saving it as name.key.svg '''
	fig, ax1 = plt.subplots(figsize=(8, 6), clear=True)
	''' generate a loop for the temperatures and the key'''
	#flambdified = sp.lambdify((temp, cov), function, modules='numpy')
	xdata = np.linspace(0, 1, 301)
	if isinstance(rconditions["temperature"], float):
		t = float(rconditions["temperature"])
		y = [function.subs({temp: t, cov: i}).evalf() for i in x]
		ax1.plot(x, y, marker='o', color='k', fillstyle='none', linestyle='none',
				 label='original at T='+str(float(rconditions["temperature"]))+' K')
		ydata = [function.subs({temp: t, cov: i}).evalf() for i in xdata]
		ax1.plot(xdata, ydata, color='b', linestyle='-',
				 label='interpolated at T='+str(float(rconditions["temperature"]))+' K')
	else:
		ramp = [int(i) for i in rconditions["temperature"]]
		y = [function.subs({temp: ramp[0], cov: i}).evalf() for i in x]
		ax1.plot(x, y, marker='o', color='k', fillstyle='none', linestyle='none',
				 label='original at T='+str(ramp[0])+' K')
		for t in range(ramp[0], ramp[1]+ramp[2], ramp[2]):
			ydata = [function.subs({temp: t, cov: i}).evalf() for i in xdata]
			ax1.plot(xdata, ydata, color='b', linestyle='-', alpha=1-(t-ramp[0])/(ramp[1]-ramp[0]),
				 label='interpolated at T='+str(t)+' K')
	ax1.set_xlabel('coverage (ML)', fontsize=18)
	ax1.set_xlim([0, 1])
	ax1.tick_params(axis='both', rotation=0, labelsize=16)
	ax1.yaxis.set_major_formatter(ticker.FuncFormatter('{:.2f}'.format))
	if ykey.startswith("i"):
		ax1.set_ylabel(ykey+" $(cm^{-1})$", fontsize=18)
	elif ykey.startswith("e"):
		ax1.set_ylabel(ykey+" $(eV)$", fontsize=18)
	else:
		ax1.set_ylabel(ykey, fontsize=18)
	legend = ax1.legend(loc="best", fontsize=14)
	fig.tight_layout()
	plt.ion()
	plt.show()
	plt.savefig(folder+ykey+"_coverage.svg", dpi=300, orientation='landscape', transparent=True)
	print("\t... Interpolating", name, ykey,"...", round(time.time()-startint, 3), " seconds")
	return function


class PartitionFunctions:
	def __init__(self, rconditions, systems, constants, restricted_arg, tstart):
		''' Reaction conditions are set as symbols using SYMPY '''
		temp = sp.symbols("temperature", positive=True, real=True)
		for name in systems.keys():     # species
			if systems[name]["kind"] == "molecule":
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						systems[name][nadsorbates]["qtrans3d"] = self.qtrans3d(systems[name][nadsorbates], constants)
						systems[name][nadsorbates]["qtrans2d"] = self.qtrans2d(systems[name][nadsorbates], constants)
						systems[name][nadsorbates]["qrot"] = self.qrot(systems[name][nadsorbates], constants)
						systems[name][nadsorbates]["qelec"] = self.qelec(systems[name][nadsorbates])
						systems[name][nadsorbates]["qvib3d"] = self.qvib(systems[name][nadsorbates]['freq3d'], constants)
						systems[name][nadsorbates]["qvib2d"] = self.qvib(systems[name][nadsorbates]['freq2d'], constants)
						datalabel3d = ["qrot", "qelec", "qtrans3d", "qvib3d"]
						q3d = 1.0
						for i in datalabel3d:
							q3d *= systems[name][nadsorbates][str(i)]
						systems[name][nadsorbates]["q3d"] = q3d
						datalabel3d.append("q3d")
						datalabel2d = ["qrot", "qelec", "qtrans2d", "qvib2d"]
						q2d = 1.0
						for i in datalabel2d:
							q2d *= systems[name][nadsorbates][str(i)]
						systems[name][nadsorbates]["q2d"] = q2d
						datalabel2d.append("q2d")
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates],
							  datalabel3d + ["qtrans2d", "qvib2d", "q2d"], "PartitionFunctions")
			else:
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						systems[name][nadsorbates]["qtrans3d"] = 1
						systems[name][nadsorbates]["qrot"] = 1
						systems[name][nadsorbates]["qelec"] = self.qelec(systems[name][nadsorbates])
						if 'freq3d' in systems[name][nadsorbates].keys():
							systems[name][nadsorbates]["qvib3d"] = self.qvib(systems[name][nadsorbates]['freq3d'], constants)
						else:
							systems[name][nadsorbates]["qvib3d"] = 1.0
						datalabel = ["qrot", "qelec", "qtrans3d", "qvib3d"]
						q3d = 1.0
						for i in datalabel:
							q3d *= systems[name][nadsorbates][str(i)]
						datalabel.append("q3d")
						systems[name][nadsorbates]["q3d"] = q3d
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates],
								  datalabel, "PartitionFunctions")
		''' Partition function interpolation between nadsorbates systems of the same name '''
		for name in systems.keys():     # species
			if len([adsorbate for adsorbate in systems[name] if adsorbate not in restricted_arg]) > 1:
				if systems[name]["kind"] == "molecule":
					systems[name]["q3d"] = interpolate(rconditions, systems, name, restricted_arg, "q3d")
					systems[name]["q2d"] = interpolate(rconditions, systems, name, restricted_arg, "q2d")
				else:
					systems[name]["q3d"] = interpolate(rconditions, systems, name, restricted_arg, "q3d")
			else:
				adsorbate =  [i for i in systems[name] if i not in restricted_arg][0]
				if systems[name]["kind"] == "molecule":
					systems[name]["q3d"] = systems[name][adsorbate]["q3d"]
					systems[name]["q2d"] = systems[name][adsorbate]["q2d"]
				else:
					systems[name]["q3d"] = systems[name][adsorbate]["q3d"]

		self.systems = systems
		print("\t\t\t\t", round((time.time() - tstart) / 60, 3), " minutes")


	@staticmethod
	def qtrans3d(properties, constants):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		 Adsorption Journal Of The International Adsorption Society
		 (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
		 page 88
		 huge numbers are expected, e.g. for CO: 6.8*1010 m-1 at 500 K in one dimension'''
		temp = sp.symbols("temperature", positive=True, real=True)
		return properties["volume"]*((2*sp.pi*properties["mass"]*constants["kb"]*temp)**(3/2))/(constants["h"]**3)

	@staticmethod
	def qtrans2d(properties, constants):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		 Adsorption Journal Of The International Adsorption Society
		 (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
		 page 88 '''
		temp = sp.symbols("temperature", positive=True, real=True)
		return properties["marea"]*(2*sp.pi*properties["mass"]*constants["kb"]*temp)/(constants["h"]**2)

	@staticmethod
	def qrot(properties, constants):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		 Adsorption Journal Of The International Adsorption Society
		 (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
		 page 90
		 large values are expected, e.g. for CO: 180 at 500 K'''
		temp = sp.symbols("temperature", positive=True, real=True)
		prod_inertia = 1.0
		for i in properties["inertia"]:
			prod_inertia *= i
		if properties["linear"] == "yes":
			qrot = (8*sp.pi**2*prod_inertia*constants["kb"]*temp)/(properties["symfactor"]*constants["h"]**2)
		else:
			qrot = ((sp.sqrt(sp.pi*prod_inertia)/properties["symfactor"]) *
					(8*sp.pi**2*constants["kb"]*temp/(constants["h"]**2))**(3/2))
		return qrot

	@staticmethod
	def qelec(properties):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		 Adsorption Journal Of The International Adsorption Society
		 (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
		 page 92 '''
		''' It is consider that the contribution of excited states is negligible
		at the working temperatures '''
		return properties["degeneration"]

	@staticmethod
	def qvib(freqs, constants):
		''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
		 Adsorption Journal Of The International Adsorption Society
		 (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
		 page 89 '''
		''' The equation used is respect the lowest occupied state, not the bottom of the potential energy curve.
		It also consider small frequencies, when h*v ~ kb*T (v=frequencies). 
		The Zero Point Energy should be added to the energy as this qvib is to calculate the entropy (S),
		specific head (Cp), and pre-exponential factor of Arrhenius (A). '''
		temp = sp.symbols("temperature", positive=True, real=True)
		qvib = 1
		for freq in freqs:
			if freq > 0.0:
				qvib *= 1/(1-sp.exp((-constants["hc"]*freq)/(constants["kb"]*temp)))
				qvib = sp.powsimp(qvib, force=True)
		return qvib


class Energy:        # Gibbs free energy in eV
	def __init__(self, rconditions, processes, systems, constants, restricted_arg):
		''' Reaction conditions are set as symbols using SYMPY '''
		temp = sp.symbols("temperature", positive=True, real=True)
		''' Entropy is required to calculate the specific heat (Cp), which contributes 
		to the enthalpy at specific temperatures '''
		start = time.time()
		print("\t... Generating Entropic Contribution ...")
		self.systems = Entropy(rconditions, systems, constants, restricted_arg).systems         # in eV
		print("\t\t\t\t", round((time.time()-start)/60, 3), " minutes")
		start = time.time()
		print("\t... Generating Enthalpic Contribution ...")
		self.systems = Enthalpy(rconditions, systems, constants, restricted_arg).systems        # in eV
		print("\t\t\t\t", round((time.time()-start)/60, 3), " minutes")
		start = time.time()
		print("\t... Generating Gibbs Free Energy ...")
		for name in systems.keys():     # species
			if systems[name]["kind"] == "molecule":
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						systems[name][nadsorbates]["energy3d"] = (systems[name][nadsorbates]["enthalpy3d"] -
																   temp * systems[name][nadsorbates]["sentropy3d"])
						systems[name][nadsorbates]["energy2d"] = (systems[name][nadsorbates]["enthalpy2d"] -
																   temp * systems[name][nadsorbates]["sentropy2d"])
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates],
							  ["energy3d", "energy2d"], "GibbsFreeEnergy")
			else:
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						systems[name][nadsorbates]["energy3d"] = (systems[name][nadsorbates]["enthalpy3d"] -
																   temp * systems[name][nadsorbates]["sentropy3d"])
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates],
							  ["energy3d"], "GibbsFreeEnergy")
		''' Energy interpolation between nadsorbates of the same name '''
		tss = [processes[pr]['ts'][i] for pr in processes.keys() for i in range(len(processes[pr]['ts']))]
		for name in systems.keys():     # species
			if len([adsorbate for adsorbate in systems[name] if adsorbate not in restricted_arg]) > 1:
				if systems[name]["kind"] == "molecule":
					systems[name]["energy3d"] = interpolate(rconditions, systems, name, restricted_arg, "energy3d")
					systems[name]["energy2d"] = interpolate(rconditions, systems, name, restricted_arg, "energy2d")
				else:
					systems[name]["energy3d"] = interpolate(rconditions, systems, name, restricted_arg, "energy3d")
					if name in tss:
						systems[name]["ifreq"] = interpolate(rconditions, systems, name, restricted_arg, "ifreq")
			else:
				adsorbate =  [i for i in systems[name] if i not in restricted_arg][0]
				if systems[name]["kind"] == "molecule":
					systems[name]["energy3d"] = systems[name][adsorbate]["energy3d"]
					systems[name]["energy2d"] = systems[name][adsorbate]["energy2d"]
				else:
					systems[name]["energy3d"] = systems[name][adsorbate]["energy3d"]
					if name in tss:
						systems[name]["ifreq"] = systems[name][adsorbate]["ifreq"]

		self.systems = systems
		print("\t\t\t\t", round((time.time()-start)/60, 3), " minutes")


class Entropy:
	''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
	Adsorption Journal Of The International Adsorption Society
	(Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
	page 88 :: S = kb*ln(Q)+kb*T*diff(ln(Q), T)'''
	def __init__(self, rconditions, systems, constants, restricted_arg):
		''' Reaction conditions are set as symbols using SYMPY '''
		for name in systems.keys():     # species
			if systems[name]["kind"] == "molecule":
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						systems[name][nadsorbates]["strans3d"] = self.strans3d(systems[name][nadsorbates], constants)
						systems[name][nadsorbates]["strans2d"] = self.strans2d(systems[name][nadsorbates], constants)
						systems[name][nadsorbates]["srot"] = self.srot(systems[name][nadsorbates], constants)
						systems[name][nadsorbates]["selec"] = self.selec(systems[name][nadsorbates], constants)
						systems[name][nadsorbates]["svib3d"] = self.svib(systems[name][nadsorbates]['freq3d'], constants)
						systems[name][nadsorbates]["svib2d"] = self.svib(systems[name][nadsorbates]['freq2d'], constants)
						datalabel3d = ["srot", "selec", "strans3d", "svib3d"]

						entropy = 0.0
						for i in datalabel3d:
							entropy += systems[name][nadsorbates][str(i)]

						#print(name, entropy)

						systems[name][nadsorbates]["sentropy3d"] = sp.factor(sp.logcombine(sp.powsimp(entropy,
																							force=True)), modulus=None)

						datalabel3d.append("sentropy3d")
						datalabel2d = ["srot", "selec", "strans2d", "svib2d"]
						entropy = 0.0
						for i in datalabel2d:
							entropy += systems[name][nadsorbates][str(i)]
						systems[name][nadsorbates]["sentropy2d"] = sp.logcombine(sp.powsimp(entropy, force=True))
						datalabel2d.append("sentropy2d")
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates],
								  datalabel3d + ["strans2d", "svib2d", "sentropy2d"], "Entropy")
			else:
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						systems[name][nadsorbates]["strans3d"] = 1e-15      # epsilon to avoid nan or inf
						systems[name][nadsorbates]["srot"] = 1e-15          # epsilon to avoid nan or inf
						systems[name][nadsorbates]["selec"] = self.selec(systems[name][nadsorbates], constants)
						if 'freq3d' in systems[name][nadsorbates].keys():
							systems[name][nadsorbates]["svib3d"] = self.svib(systems[name][nadsorbates]['freq3d'], constants)
						else:
							systems[name][nadsorbates]["svib3d"] = 1e-15    # epsilon to avoid nan or inf
						datalabel3d = ["strans3d", "srot", "selec", "svib3d"]
						entropy = 0.0
						for i in datalabel3d:
							entropy += systems[name][nadsorbates][str(i)]
						systems[name][nadsorbates]["sentropy3d"] = sp.logcombine(sp.powsimp(entropy, force=True))
						datalabel3d.append("sentropy3d")
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates],
								  datalabel3d, "Entropy")
		self.systems = systems

	@staticmethod
	def strans3d(properties, constants):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
		 (Note that the translational component also includes components from the Stirling approximation)'''
		temp = sp.symbols("temperature", positive=True, real=True)
		return (constants["kb"]*(sp.log(((2*sp.pi*properties["mass"]*constants["kb"]*temp)/constants["h"]**2)**(3/2) *
										(constants["kb"]*temp)/1 ) + 5/2) * constants["JtoeV"])

	@staticmethod
	def strans2d(properties, constants):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
		 (Note that the translational component also includes components from the Stirling approximation)'''
		temp = sp.symbols("temperature", positive=True, real=True)
		return (constants["kb"]*(sp.log(((2*sp.pi*properties["mass"]*constants["kb"]*temp)/constants["h"]**2)**(2/2) *
										(constants["kb"]*temp)/1 ) + 3/2)*constants["JtoeV"])

	@staticmethod
	def srot(properties, constants):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
		temp = sp.symbols("temperature", positive=True, real=True)
		prod_inertia = 1.0
		for i in properties["inertia"]:
			prod_inertia *= i
		if properties["linear"] == "yes":
			srot = constants["kb"]*(sp.log((1/properties["symfactor"])*(8*sp.pi**2*constants["kb"]*temp*
													prod_inertia/(constants["h"]**2))) + 1)*constants["JtoeV"]
		else:
			srot = constants["kb"]*(sp.log((sp.sqrt(sp.pi*prod_inertia)/properties["symfactor"]) *
								  (8*sp.pi**2*constants["kb"]*temp/(constants["h"]**2))**(3/2)) +3/2)*constants["JtoeV"]
		return srot

	@staticmethod
	def selec(properties, constants):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
		return constants["kb"]*(sp.log(2*properties["degeneration"] + 1)) * constants["JtoeV"]

	@staticmethod
	def svib(freqs, constants):
		'''re-Formulation from explicit derivatives::
		 https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
		 harmonic oscillator '''
		temp = sp.symbols("temperature", positive=True, real=True)
		qvib = 1e-15    # to avoid nan or inf
		for freq in freqs:
			if freq > 0.0:
				x = (constants['hc'] * freq) / (constants['kb'] * temp)
				qvib += (x / (sp.exp(x) - 1)) - (sp.log(1 - sp.exp(-x)))
				qvib = sp.powsimp(qvib, force=True)  # combine powers/exp patterns
		return constants["kb"] * qvib * constants["JtoeV"]


class Enthalpy:
	""" C.J. Cramer. Essentials of Computational Chemistry, Second Edition. Wiley, 2004.
	and Raymand Chang, "PHYSICAL CHEMISTRY for Chemical and Biological Science", ISBN: 1-891389-06-8,
		page 91 :: H=[dCp/dT](T1,T2)
		H =  E + ZPE + integral(Cp, 0 --> T) """
	def __init__(self, rconditions, systems, constants, restricted_arg):
		''' Reaction conditions are set as symbols using SYMPY '''
		temp = sp.symbols("temperature", positive=True, real=True)
		(temp0, temp1, tempstep) = rconditions['temperature']
		for name in systems.keys():     # species
			if systems[name]["kind"] == "molecule":
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						systems[name][nadsorbates]["zpe3d"] = self.zpe(systems[name][nadsorbates]['freq3d'], constants) # in eV
						systems[name][nadsorbates]["zpe2d"] = self.zpe(systems[name][nadsorbates]['freq2d'], constants)  # in eV
						systems[name][nadsorbates]["cp3d"], cp3d_integral = self.cp3d(systems[name][nadsorbates][
																					   'freq3d'], constants)
						systems[name][nadsorbates]["cp2d"] = self.cp2d(systems[name][nadsorbates]['freq2d'], constants)

						print(name, systems[name][nadsorbates]["cp3d"])

						enthalpy3d = (systems[name][nadsorbates]["energy0"] + systems[name][nadsorbates]["zpe3d"] +
									  cp3d_integral)
						systems[name][nadsorbates]["enthalpy3d"] = enthalpy3d
						enthalpy2d = (systems[name][nadsorbates]["energy0"] + systems[name][nadsorbates]["zpe2d"] +
									  sp.integrate(systems[name][nadsorbates]["cp2d"], temp, meijerg=True))
						systems[name][nadsorbates]["enthalpy2d"] = enthalpy2d
						datalabel = ["zpe3d", "cp3d", "enthalpy3d", "cp2d", "zpe2d", "enthalpy2d"]
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates],
							  datalabel, "Enthalpy")
			else:
				for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
					if nadsorbates not in restricted_arg:       # only for nadsorbates
						if 'freq3d' in systems[name][nadsorbates].kesy():
							systems[name][nadsorbates]["zpe3d"] = self.zpe(systems[name][nadsorbates]['freq3d'], constants) # in eV
						else:
							systems[name][nadsorbates]["zpe3d"] = 1e-15 # epsilon to avoid nan or inf
						systems[name][nadsorbates]["cp3d"] = self.cp3d(systems[name][nadsorbates]['freq3d'], constants)
						enthalpy3d = (systems[name][nadsorbates]["energy0"] +  systems[name][nadsorbates]["zpe3d"] +
									  sp.integrate(systems[name][nadsorbates]["cp3d"], temp, meijerg=True))
						systems[name][nadsorbates]["enthalpy3d"] = enthalpy3d
						datalabel = ["zpe3d", "cp3d", "enthalpy3d"]
						printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], datalabel, "Enthalpy")
		self.systems = systems


	@staticmethod
	def zpe(freqs, constants):
		''' Temperature dependent quantum vibrational energy (exact for harmonic oscillator),
		 useful if you want the thermal contribution on top of ZPE '''
		temp = sp.symbols("temperature", positive=True, real=True)
		zpe = 1/2*constants["hc"]
		freq_sum = 0
		uvib = 0
		for freq in freqs:
			if freq > 0.0:
				freq_sum += freq
				uvib += constants["hc"]*freq * (1/2 + 1/(sp.exp(constants["hc"]*freq/(constants["kb"]*temp)) - 1))
				uvib = sp.powsimp(uvib, force=True)
		return (zpe*freq_sum + uvib) * constants["JtoeV"]

	@staticmethod
	def cp3d(freqs, constants):
		'''Hans Kuhn, Horst-Dieter Försterling, David Hennessey Waldeck, "Principles of Physical Chemistry"
		ISBN: 9780470089644, page 551 :: Cp=T*[dS/dT](N,P)
		Also, Cp = Cp_trans + Cp_rot + Cp_vib
		Later on to calculate H[T], Cp has to be integrated, which is very demanding in resources.
		 For that reason, Cp analytical expression is simplidied with a chain of simplidiers'''
		temp = sp.symbols("temperature", positive=True, real=True)
		#cp = temp * sp.diff(entropy, temp)  # already in eV
		cp_trans = 5/2 * constants['kb']
		cp_rot = constants['kb']
		cp_vib = 0.
		cp_integral = sp.integrate(sp.together(sp.cancel(sp.factor(cp_trans + cp_rot))), temp)
		for freq in freqs:
			x = (constants['hc']*freq) / (constants['kb']*temp)
			cp_vib += constants['kb'] * (x / (2*sp.sinh(x/2)))**2
			''' Best targeted simplifiers '''
			cp_vib = sp.powsimp(cp_vib)  # simplify powers and exponents
			cp_vib = sp.cancel(cp_vib)  # cancel common factors
			cp_vib = sp.radsimp(cp_vib)  # rational simplification
			cp_vib = sp.factor_terms(cp_vib)  # factor by common subexpressions
			cp_integral += sp.integrate(sp.together(sp.cancel(sp.factor(cp_vib))), temp)
		cp = sp.together(sp.cancel(sp.factor(cp_trans + cp_rot + cp_vib)))
		return constants['JtoeV'] * cp,  constants['JtoeV'] * cp_integral

	@staticmethod
	def cp2d(freqs, constants):
		'''Hans Kuhn, Horst-Dieter Försterling, David Hennessey Waldeck, "Principles of Physical Chemistry"
		ISBN: 9780470089644, page 551 :: Cp=T*[dS/dT](N,P)
		Also, Cp = Cp_trans + Cp_rot + Cp_vib
		Later on to calculate H[T], Cp has to be integrated, which is very demanding in resources.
		 For that reason, Cp analytical expression is simplidied with a chain of simplidiers'''
		temp = sp.symbols("temperature", positive=True, real=True)
		#cp = temp * sp.diff(entropy, temp)  # already in eV
		cp_trans = 2 * constants['kb']  # in 2 dimensions
		cp_rot = constants['kb']
		cp_vib = 0.
		for freq in freqs:
			x = (constants['hc'] * freq) / (constants['kb'] * temp)
			cp_vib += (x / (2 * sp.sinh(x / 2))) ** 2
			''' Best targeted simplifiers '''
			cp_vib = sp.powsimp(cp_vib)  # simplify powers and exponents
			cp_vib = sp.cancel(cp_vib)  # cancel common factors
			cp_vib = sp.radsimp(cp_vib)  # rational simplification
			cp_vib = sp.factor_terms(cp_vib)  # factor by common subexpressions
		cp = sp.together(sp.cancel(sp.factor(cp_trans + cp_rot + constants['kb']*cp_vib)))
		return constants['JtoeV'] * cp

