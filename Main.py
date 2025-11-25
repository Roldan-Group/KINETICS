"""

make first Input2mk.py


	This script reads and executes the packages for KINETICS

		by Alberto Roldan
"""
import sys, re
import time
import numpy as np
import sympy as sp
from Thermodynamics import PartitionFunctions, Energy
from Kinetics import RConstants, REquations
from Experiments import ConsTemperature, TPR
from Symbols_def import temp, kb




''' Read input file with the systems and kinetic simulation conditions'''
def mkread(inputfile, restricted_arg):
	rconditions = {}    # reaction conditions
	processes = {}      # reaction processes
	process = 0         # process number, key of processes starting from 1
	systems = {}        # species
	print(inputfile)
	for line in open(inputfile, 'r'):
		if not re.match(r'^\s*$', line) and line.startswith("#") is False:
			line = line.split("=")
			head = line[0].strip()
			tail = []
			for i in line[1].split():
				if i != "#" :
					tail.append(i)
				else:
					break
			''' Reaction conditions are stored in a dictionary (rconditions), including:
			to include the last point dx should be include, e.g. 100, 200, 10 -> 100, 210, 10
				- External potential (vext): constant or ramp (initial, final, step)
				- pH (ph): constant or ramp (initial, final, step) 
				- Temperature (temperature): constant or ramp (initial, final, step)
				- time (time): from 0 to time in a timestep
				*** conditions to be added
			'''
			if head == "ELECTRO" or head == "ELECTROPOTENTIAL":
				if len(tail) == 1:
					rconditions["vext"] = float(tail[0])    # constant
				else:
					rconditions["vext"] = [float(i) for i in tail] # ramp
					rconditions["vext"][1] = rconditions['vext'][1] + rconditions['vext'][2]
			if head == "PH" or head == "pH":
				if len(tail) == 1:
					rconditions["ph"] = float(tail[0])    # constant
				else:
					rconditions["ph"] = [float(i) for i in tail[:-1]] + [int(tail[-1])] # ramp
					rconditions["ph"][1] = rconditions['phvext'][1] + rconditions['ph'][2]
			if head == "TEMP" or head == "TEMPERATURE":
				if len(tail) == 1:
					rconditions["temperature"] = float(tail[0])     # constant
				else:
					rconditions["temperature"] = [float(i) for i in tail[:-1]] + [int(tail[-1])]   # ramp
					rconditions["temperature"][1] = rconditions['temperature'][1] + rconditions['temperature'][2]

			if head == "TIME" or head == "Time":
				rconditions["time"] = [0, float(tail[0]), int(tail[-1])]  # initial time is 0
				rconditions["time"][1] = rconditions['time'][1] + rconditions['time'][2]

			''' Processes (Adsorption, Reaction, Desorption) in a dictionary (processes)
			with key = number of the process, including:
				- kind of process (kind = a, r, d)
				- [reactant species] ('reactants')
				- [transition states] ('ts')
				- [product species] ('products')
				- [reactants stoichiometry] ('rstoichio')
				- [products stoichiometry] ('pstoichio')
				*** conditions to be added
			'''
			if head == "PROCESS":
				process += 1
				processes[str(process)] = {}    # dictionary of items for process
				processes[str(process)]["kind"] = str(tail[0][0])  # kind of process
				reaction = ''.join(tail[1:]).split(">")
				for i in range(3):      # [reactants, TSs, products]
					species = []
					stoichio = []
					for j in reaction[i].split("+"):
						k = j[0:len(j.rsplit(r'[0-9]'))]
						try:
							stoichio.append(float(k))
							species.append(str(j[len(j.rsplit(r'[0-9]')):]))
						except ValueError:
							if len(j) > 0:
								stoichio.append(1.0)
								species.append(str(j))
					if i == 0:
						processes[str(process)]["reactants"] = species  # reactants
						processes[str(process)]["rstoichio"] = stoichio
					elif i == 1:
						processes[str(process)]["ts"] = species
						processes[str(process)]["tsstoichio"] = stoichio
					else:
						processes[str(process)]["products"] = species  # products
						processes[str(process)]["pstoichio"] = stoichio
			''' Systems, i.e. the species involved, in a nested dictionary (systems) with key = name,
			including:
				- number of adsorbates (nadsorbates), accounting for the systems' coverage.
				- kind (surface, molecule, adsorbed)
				- path the input, i.e. QM data, (syspath)
				- path to frequencies (freqpath)
				- DFT energy (energy0)
				- frequencies (freq3d & freq2d) either 3D or only along the plane x and y axis (2D)
				- different adsorption sites (site), only for naked systems, e.g. surface
				- different areas related to the sites (area) in m^2, only for naked systems, e.g. surface
				- molecule's mass in kg (mass)
				- number of atoms in the molecule (natoms)
				- molecule's symmetry factor (symfactor)
				- molecule's inertia moment(s) in Kg.m^-2 (inertia)
				- molecule's adsorption sites (molsite), which must be in sites
				- electronic multiplicity (degeneration)
				- initial molecular partial pressure in Pa, constant or ramp (initial, final, step)
				- initial molecular coverage in ML, constant or ramp (initial, final, step)
				
				?? species = e from electrons when  "vext" in rconditions
				*** conditions to be added
			'''
			if head == "SYSTEM":
				name = str(tail[0])         # species name
				systems[name] = {}
				try:
					nadsorbates = str(tail[1])    # number of adsorbates in system
				except:
					nadsorbates = str(1)  # number of adsorbates in system by default
					pass
				systems[name][nadsorbates] = {}
			if head == "SYSPATH":
				systems[name][nadsorbates]["syspath"] = tail[0]
			if head == "FREQPATH":
				systems[name][nadsorbates]["freqpath"] = str(tail[0])
			if head == "E0":             # species energy
				systems[name][nadsorbates]["energy0"] = float(tail[0])
			if head == "ISITES":
				nsites = []
				sites = []                  # catalyst adsorption sites
				for i in range(int(np.ceil(len(tail)/2))):
					try:
						nsites.append(float(tail[i]))
					except ValueError:
						nsites.append(1.0)
						sites.append(str(tail[i]))
					else:
						sites.append(str(tail[i+1]))
				systems[name][nadsorbates]["nsites"] = nsites   # number of sites occupied by adsorbate
				systems[name][nadsorbates]["sites"] = sites   # catalyst adsorption sites
				systems[name]["sites"] = sites
			if head == "IACAT":
				areas = []                   # adsorption areas
				for i in tail:
					try:
						areas.append(float(i))
					except ValueError:
						pass
				systems[name][nadsorbates]["area"] = areas   # adsorption areas
			if head == "FREQ":
				freq = []                   # species frequencies 3D
				for i in tail:
					try:
						freq.append(float(i))
					except ValueError:
						pass
				systems[name][nadsorbates]["freq3d"] = sorted(freq, reverse=True)   # species frequencies 3D
			if head == "FREQ2D":            # species frequencies only considering x and y displacements
				freq2d = []                 # i.e. displacement of their center of mass < 0.1 on the Z-axis.
				for i in tail:
					try:
						freq2d.append(float(i))
					except ValueError:
						pass
				systems[name][nadsorbates]["freq2d"] = sorted(freq2d, reverse=True)   # species frequencies only considering x and y displacements
			if head == "IMASS":
				systems[name][nadsorbates]["imass"] = list([float(i) for i in tail])
			if head == "INATOMS":
				systems[name][nadsorbates]["natoms"] = list([int(i) for i in tail])
			if head == "SYMFACTOR":
				systems[name][nadsorbates]["symfactor"] = int(tail[0])
			if head == "INERTIA":
				inertia = []          # molecular inertia moment(s)
				for i in tail:
					try:
						inertia.append(np.abs(float(i)))
					except ValueError:
						pass
				if len(inertia) == 1:
					systems[name][nadsorbates]["inertia"] = inertia[0]   # molecule's inertia moment(s) in Kg/m^2
				else:
					systems[name][nadsorbates]["inertia"] = inertia   # molecule's inertia moment(s) in Kg/m^2
				if isinstance(inertia, float):      # molecular linearity
					systems[name][nadsorbates]["linear"] = "yes"
				elif inertia[0] == inertia[1] or inertia[0] == inertia[2] or inertia[1] == inertia[2]:
					systems[name][nadsorbates]["linear"] = "yes"
				else:
					systems[name][nadsorbates]["linear"] = "no"
			'''A molecule will adsorb on one a site with a particular area (marea). If the molecules has more than 
			site to adsorbed, differente systems needs to be described'''
			if head == "MOLSITE":
				equivalent = None
				for i in tail:
					try:
						equivalent = float(i)   # molecular area equivalent of molsite
					except ValueError:
						molsite = str(i)   # molecular adsorption site)
				''' It may be the case equivalent = 1 is neglected.
				Then, it is check that the len(equivalents) is the same than
				the number of molsites or 1 is added to equivalent.'''
				if equivalent == None:
					equivalent = float(1)
				systems[name][nadsorbates]["molsite"] = molsite   # molecular adsorption site
				systems[name][nadsorbates]["nmolsite"] = equivalent   # molecular area equivalent of molsite
			if head == "DEGENERATION":
				systems[name][nadsorbates]["degeneration"] = int(tail[0])
			if head == "IPRESSURE":
				if len(tail) == 1:
					systems[name]["pressure0"] = float(tail[0])     # constant
				else:
					systems[name]["pressure0"] = [float(i) for i in tail]   # ramp
			if head == "ICOVERAGE":
				if len(tail) == 1:
					systems[name]["coverage0"] = float(tail[0])     # constant
				else:
					systems[name]["coverage0"] = [float(i) for i in tail]   # ramp
	''' Checks in systems to assign the kind of system and check:
		- the kind, 
		- the path for the frequencies,
		- the sites,
		- the degeneration,
		- the initial pressures/coverages,
		- the 3N-6(5) frequencies,
		- The number of adsorbates in an adsorbed system. '''
	surf = []
	for name in systems.keys():
		if "area" in systems[name][list(systems[name].keys())[0]].keys():
			systems[name]["kind"] = "surface"
			surf.append(name)       # temporal list to later remove the surfaces name from systems
		elif "imass" in systems[name][list(systems[name].keys())[0]].keys():
			systems[name]["kind"] = "molecule"
		else:
			systems[name]["kind"] = "adsorbate"
	''' Swaping the surfaces from dict(systems) for the site '''
	for name in surf:
		for site in systems[name][list(systems[name].keys())[0]]['sites']:
			systems[site] = systems[name]
		del systems[name]

	for name in systems.keys():
		if systems[name]["kind"] == "surface":
			for key in systems[name].keys():
				if key not in restricted_arg:       # only for nadsorbates
					if len(systems[name][str(key)]['area']) != len(systems[name][str(key)]["sites"]):
						print("   ERROR: the number of sites and areas in {}{} is not the same.".format(name,key))
						exit()
					if ("freq3d" not in systems[name][key].keys() or len(systems[name][key]["freq3d"]) == 0):
						print("   NOTE: frequencies for {}{} are not provided".format(name, key))
		if systems[name]["kind"] == "molecule":
			if "pressure0" not in systems[name].keys():
				systems[name]["pressure0"] = 0.
			for key in systems[name].keys():
				if key not in restricted_arg:       # only for nadsorbates
					mass = 0
					for i in range(len(systems[name][key]["imass"])):
						mass += systems[name][key]['imass'][i] * systems[name][key]['natoms'][i]
					systems[name][key]["mass"] = mass / (6.02214076e23 * 1000)    # in kg
					if "freqpath" not in systems[name][key].keys():
						systems[name][key]["freqpath"] = systems[name][key]["syspath"]
					if "degeneration" not in systems[name][key].keys():
						systems[name][key]["degeneration"] = 1
					if "freq3d" not in systems[name][key].keys() or len(systems[name][key]["freq3d"]) == 0:
						print("   ERROR: frequencies for {}{} are not provided".format(name, key))
						exit()
					if "freq2d" not in systems[name][key].keys() or len(systems[name][key]["freq2d"]) == 0:
						print("   ERROR: 2D frequencies for {}{} are not provided".format(name, key))
						exit()
		if systems[name]["kind"] == "adsorbate":
			if "coverage0" not in systems[name].keys():
				systems[name]["coverage0"] = 0.
			for key in systems[name].keys():
				if key not in restricted_arg and isinstance(systems[name][key], dict):       # only for nadsorbates
					if ("freq3d" not in systems[name][key].keys() or len(systems[name][key]["freq3d"]) == 0):
						print("   NOTE: frequencies for {}{} are not provided".format(name, key))
						exit()

	''' Once the systems has defined as molecule or surface, 
	looking for the area occupied for each molecule (marea)'''
	sites = {}
	for name in systems.keys():
		if systems[name]["kind"] == "surface":
			for key in systems[name].keys():
				if key not in restricted_arg:       # only for nadsorbates
					for i in range(len(systems[name][key]["sites"])):
						sites[str(systems[name][key]["sites"][i])] = float(systems[name][key]["area"][i])
	'''A molecule will adsorb on one site with a particular area (marea). If the molecules has more than site to 
	adsorbed, differente systems needs to be described'''
	for name in systems.keys():
		if systems[name]["kind"] == "molecule":
			for key in systems[name].keys():
				if key not in restricted_arg:       # only for nadsorbates
					if systems[name][key]["molsite"] in sites:
						molsite = systems[name][key]["molsite"]
						systems[name][key]["marea"] = (sites[str(systems[name][key]["molsite"])] /
										 systems[name][key]["nmolsite"])    # should be the same as the stoichiometry in processes
						nmolsite = systems[name][key]["nmolsite"]
					pressure = 101325       # Pa == kg⋅m^−1⋅s^−2
					systems[name][key]["volume"] = kb*temp/pressure        # Assuming ideal behaviour of gases
			systems[name]['molsite'] = molsite
			systems[name]['nmolsite'] = nmolsite
		elif systems[name]["kind"] == "adsorbate":
			adsorbate_site = []     # there is ONLY one kind of site per adsorbate-name
			for key in systems[name].keys():
				if key not in restricted_arg:       # only for nadsorbates
					if systems[name][key]["sites"][0] in sites: # there in ONLY one kind of site per adsorbade-name
						adsorbate_site.append(systems[name][key]['nsites'][0] / float(key))
			systems[name]["nsites"] = sum(adsorbate_site) / len(adsorbate_site)     # nsites per single adsorbate

	for name in systems.keys():
		freq = []
		ts = [i for n in processes.keys() for i in processes[n]["ts"]]
		for key in systems[name].keys():
			if key not in restricted_arg and "freq3d" in systems[name][key].keys(): # only for adsorbates
				if name in ts:
					for f in range(len(systems[name][key]["freq3d"])-1):
						i = systems[name][key]['freq3d'][f]
						if i < -100:
							print("   ALERT: {}_{} has more than one significant imaginary frequency (abs({}))".
								  format(name, key, i))
							freq.append(np.abs(i))
						else:
							freq.append(i)
					systems[name][key]["ifreq"] = systems[name][key]["freq3d"][-1]
				elif name not in ts:  # only for nadsorbates
					for i in systems[name][key]["freq3d"]:
						if i < -100:
							print("   ALERT: {}_{} has a significant imaginary frequency (abs({}))".format(name, key, i))
							freq.append(np.abs(i))
						else:
							freq.append(i)
				systems[name][key]["freq3d"] = freq

	''' Check the species in processes and convert the "surfaces" in kind of 
	sites with the corresponding stoichiometry'''
	def take_molecule(names, systems):
		for i in range(len(names)):
			try:
				molsite = systems[names[i]]['molsite']
				site_stoi = systems[names[i]]['nmolsite']
				molecule = i
			except:
				continue
		return molsite, site_stoi, [i for i in range(len(names)) if i != molecule][0]

	species = []
	for pr in processes.keys():
		if processes[str(pr)]['kind'] == "A":
			molsite, site_stoi, i = take_molecule(processes[str(pr)]['reactants'], systems)
			processes[str(pr)]['reactants'][i] = molsite   # in case the Surf name is changed by the site
			processes[str(pr)]['rstoichio'][i] = site_stoi
		elif processes[str(pr)]['kind'] == "D":
			molsite, site_stoi, i = take_molecule(processes[str(pr)]['products'], systems)
			processes[str(pr)]['products'][i] = molsite   # in case the Surf name is changed by the site
			processes[str(pr)]['pstoichio'][i] = site_stoi
		species.extend(processes[str(pr)]['reactants'])
		species.extend(processes[str(pr)]['ts'])
		species.extend(processes[str(pr)]['products'])
	species = list(set(species))
	for i in species:
		if i not in systems.keys():
				print("   ERROR: {} is not defined in dict(systems).".format(str(i)))
				exit()

	return rconditions, processes, systems

''' list of restricted argunments in systems[name] containing the interpolated functions'''
restricted_arg = ["kind", "pressure0", "coverage0", "sites", "nsites", 'molsite', 'nmolsite',
				  'q3d', 'q2d', 'energy3d', 'energy2d', 'ifreq']

start0 = time.time()
rconditions, processes, systems = mkread(str(sys.argv[1]), list(restricted_arg))
print("... Reading ...\t\t\t", round(time.time()-start0, 3), " seconds")
start = time.time()
print("... Generating Partition Functions ...")
systems = PartitionFunctions(dict(rconditions), dict(systems), list(restricted_arg), start).systems
print("... Generating Thermodynamics ...")
systems = Energy(dict(rconditions), dict(processes), dict(systems), list(restricted_arg)).systems
start = time.time()
print("... Generating Reaction Constants ...")
processes = RConstants(dict(rconditions), dict(systems), dict(processes), list(restricted_arg)).processes
print("\t\t\t\t", round((time.time()-start)/60, 3), " minutes")
start = time.time()
print("... Generating Rate Equations ...")
constemperature = REquations(dict(processes), dict(systems)).constemperature
#surf_equations = REquations(dict(processes), dict(systems)).surfequations
tpd = REquations(dict(processes), dict(systems)).tpd
print("\t\t\t\t", round((time.time()-start), 3), " seconds")
print("... Computing Microkinetics ...")
ConsTemperature(dict(rconditions), dict(systems), dict(processes), dict(constemperature))
TPR(dict(rconditions), dict(systems), dict(processes), dict(tpd))
print("... Microkinetics Completed ...")
print("\t\t\t\tTotal time:", round((time.time()-start0)/60, 3), " minutes")
