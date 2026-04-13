"""
	This script builds on the perl version by A.Roldan.

"""
import os, pathlib
import sympy as sp
import numpy as np
from Symbols_def import t, temp, h, kb, hc, JtoeV, constants
from sympy import Max, Piecewise
from math import gcd
from functools import reduce
from itertools import product
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from collections import defaultdict



class RConstants:
	def __init__(self, rconditions, systems, processes, restricted_arg):
		''' Reaction conditions are set as symbols using SYMPY '''
		ramp = [float(i) for i in rconditions["temperature"]]
		constants_data = [[temp_num for temp_num in np.arange(ramp[0], ramp[1], ramp[2])]]
		for process in processes:
			processes[process]["activation"] = self.activation(processes[process], systems)
			if processes[process]['kind'] == 'A':
				processes[process]["sticky"] = self.sticky(processes[process], systems)
				processes[process]["arrhenius"] = self.arrhenius(processes[process], systems, restricted_arg)
				processes[process]["ktunneling"] = self.tunneling(processes[process], systems)
				processes[process]['krate0'] = (processes[process]["sticky"] * processes[process]["arrhenius"] *
												processes[process]["ktunneling"])
				# units of m*kg^-1*s^-1 |when multiplied by Pa = s^-1
				datalabel = ["activation", "sticky", "arrhenius", 'ktunneling', "krate0"]
				#self.printdata(rconditions, processes[process], datalabel, "Process"+str(process))
			else:
				processes[process]["arrhenius"] = self.arrhenius(processes[process], systems, restricted_arg)
				processes[process]["ktunneling"] = self.tunneling(processes[process], systems, )   # NO units
				processes[process]['krate0'] = (processes[process]["arrhenius"] *
												sp.exp(-processes[process]['activation']/ (kb*temp*JtoeV)) *
												processes[process]["ktunneling"])   # units s^-1
				datalabel = ["activation", "arrhenius", 'ktunneling', "krate0"]
			data = self.getdata(rconditions, processes[process], datalabel)
			self.printdata(rconditions, processes[process], "Process"+str(process), data)
			constants_data.append([row[-1] for row in data[1:]])
		stack = np.column_stack([i for i in constants_data])
		constants_data = [["# Temperature", *list(processes.keys())]]
		constants_data.append(stack)
		self.barplot("Elementary Steps", "Reaction Constants $(M^{ \dagger } \cdot s^{-1})$", constants_data, 0.5)
		self.processes = processes

	@staticmethod
	def activation(process, systems):
		''' the activation energy in adsorption and desorption processes is considered as the difference between
		a state in which the molecule has only two degrees of freedom (being the third degree the reaction coordinate)
		and the reactants '''
		if len(process['ts']) > 0:
			ets = 0  # total energy for transition states
			for i in range(len(process['ts'])):
				ets += process['tsstoichio'][i] * systems[process['ts'][i]]['energy3d']
		elif 'molecule' in [systems[i]['kind'] for i in process['reactants']]:  # for adsorption processes
			''' This elif considers molecules in reactants as in adsorption processes '''
			ets = 0  # total energy for transition states
			for i in range(len(process['reactants'])):
				if 'energy2d' in systems[process['reactants'][i]]:
					ets += process['rstoichio'][i] * systems[process['reactants'][i]]['energy2d']
				else:
					ets += process['rstoichio'][i] * systems[process['reactants'][i]]['energy3d']
		elif 'molecule' in [systems[i]['kind'] for i in process['products']]:   # for desorption processes
			''' This elif considers molecules in products as in desorption processes '''
			ets = 0  # total energy for transition states
			for i in range(len(process['products'])):
				if 'energy2d' in systems[process['products'][i]]:
					ets += process['pstoichio'][i] * systems[process['products'][i]]['energy2d']
				else:
					ets += process['pstoichio'][i] * systems[process['products'][i]]['energy3d']
		else:
			''' In the rare case that there is two transtion states and none of the reactants is a molecule the 
			energy of the transition states will be the energy of the final state.'''
			ets = 0  # total energy for transition states
			for i in range(len(process['products'])):
				ets += process['pstoichio'][i] * systems[process['products'][i]]['energy3d']

		er = 0      # total energy for reactants
		for i in range(len(process['reactants'])):
			er += process['rstoichio'][i] * systems[process['reactants'][i]]['energy3d']
		e_activation = sp.Max(ets - er, 0.0)    # ensures that the activation energy is never below 0
		return e_activation

	@staticmethod
	def sticky(process, systems):
		''' the sticky coefficient is evaluate as the reduction of degrees of freedom, i.e. from a 3D free molecule
		to a 2D trapped molecule moving parallel to the surface (being the third degree the reaction coordinate) '''
		''' Reaction conditions are set as symbols using SYMPY '''
		qr = 1     # total partition function for reactants
		for i in range(len(process['reactants'])):
			qr *= systems[process['reactants'][i]]['q3d'] ** process['rstoichio'][i]

		def build_qts(name, q_list):	# to calculate the qts of molecular Adsorptions
			expr = 1
			for q in q_list:
				expr *= systems[name][q] ** process['rstoichio'][i]
			return expr

		qts = 1     # total partition function for transition states
		if len(process['ts']) > 0:
			for i in range(len(process['ts'])):
				qts *= systems[process['ts'][i]]['q3d']**process['tsstoichio'][i]
			for i in range(len(process['reactants'])):
				qr *= systems[process['reactants'][i]]['q3d'] ** process['rstoichio'][i]
		else:
			for i in range(len(process['reactants'])):
				name = process['reactants'][i]
				if systems[name]["kind"] == "molecule":
					''' TS can be mobile or immobile (Chorkendorff, I. & Niemantsverdriet, 
					J. W. "Concepts of Modern Catalysis and Kinetics." doi:10.1002/3527602658. page 119-121 '''
					e_a = process['activation']
					expr1 = build_qts(name, ["q3d"])	# minimal distortion
					expr2 = build_qts(name, ["qrot", "qelec", "qtrans3d", "qvib2d"])	# mobile TS
					expr3 = build_qts(name, ["qrot", "qelec", "qtrans2d", "qvib2d"])	# mobile TS in 2D
					expr4 = build_qts(name, ["qrot", "qelec", "qvib2d"])	# partially immobile TS (Direct adsorption)
					expr5 = build_qts(name, ["qelec", "qvib2d"])	# immobile TS
					qts *= Piecewise((expr1, (e_a <= 0.01)),
									(expr2, (0.01 < e_a) & (e_a <= 0.25)),
									(expr3, (0.25 < e_a) & (e_a <= 0.7)),
									(expr4, (0.7 < e_a) & (e_a <= 1.)),
									(expr5 * sp.exp(-e_a / (kb * temp * JtoeV)), (e_a > 1.)))
				else:
					qts *= systems[name]['q3d']**process['rstoichio'][i]
		return qts/qr

	@staticmethod
	def arrhenius(process, systems, restricted_arg):
		''' Reaction conditions are set as symbols using SYMPY '''
		if process['kind'] == 'A':
			area = 0    # molecular area (marea)
			mass = 0    # Molecular mass
			for i in process['reactants']:
				if systems[i]['kind'] == 'molecule':
					''' In principle, the area of a molecule will be practically the same independently of 
					the working coverages (mean-field is not applicable at high coverages), for that reason it takes
					the area of the first nadsorbate. '''
					area = float([systems[i][j]['marea'] for j in systems[i].keys() if j not in restricted_arg][0])
					''' Same reasoning is applied for the molecular mass'''
					mass = float([systems[i][j]['mass'] for j in systems[i].keys() if j not in restricted_arg][0])
			arrhenius = area * 1/sp.sqrt(2*sp.pi*mass*kb*temp)
			# units of m*kg^-1*s^-1 |when multiplied by Pa = s^-1
		else:
			qts = 1  # total partition function for transition states
			qr = 1  # total partition function for reactants
			if len(process['ts']) > 0:
				for i in range(len(process['ts'])):
					qts *= systems[process['ts'][i]]['q3d'] ** process['tsstoichio'][i]
				for i in range(len(process['reactants'])):
					qr *= systems[process['reactants'][i]]['q3d'] ** process['rstoichio'][i]
			else:   # e.g. for desorption processes
				for i in range(len(process['products'])):
					if 'q2d' in systems[process['products'][i]]:
						qts *= systems[process['products'][i]]['q2d'] ** process['pstoichio'][i]
					elif 'q3d' in systems[process['products'][i]]:
						qts *= systems[process['products'][i]]['q3d'] ** process['pstoichio'][i]
				for i in range(len(process['reactants'])):
					qr *= systems[process['reactants'][i]]['q3d'] ** process['rstoichio'][i]
			arrhenius = kb*temp/h * qts/qr    # units of s^-1
		return arrhenius

	@staticmethod
	def tunneling(process, systems, ):
		''' Second order harmonic Wigner approach to shallow quantum tunneling valid for
		vast numbers of reaction including surface-catalysed --> DOI: 10.1039/C4CP03235G '''
		''' Reaction conditions are set as symbols using SYMPY '''
		k = 1
		for i in range(len(process['ts'])):     # systems[process['ts'][i] is the system's name
			try:
				k = 1 + 1/24 * (hc * systems[process['ts'][i]]['ifreq'] /(2*sp.pi * kb*temp))**2
			except:
				pass
		return k

	@staticmethod
	def ext_pressure(process, systems, ):
		''' Reaction conditions are set as symbols using SYMPY '''
		k = 1
		return k

	@staticmethod
	def electric(process, systems, ):
		''' Phys. Rev. Lett. 2007, 99, 126101                             DOI:https://doi.org/10.1103/PhysRevLett.99.126101
			J. Phys. Chem. C, 2010, 114 (42), pp 18182–18197              DOI: 10.1021/jp1048887
			The hydrogen coverage will be dependent on the potential via the reaction:
								H+ + e- + *  --> H*
			At standard conditions (298 K, pH 0, 1 bar H2) and U = 0 V vs NHE,
			the left-hand side is in equilibrium with hydrogen gas.
			At finite bias, U, the chemical potential of the electron will be linearly dependent on the bias.
			The reaction free energy can be written as:
								ΔGH* = AG + AG(U)= AG + −eU
			defines the chemical potential of H*.'''
		''' Reaction conditions are set as symbols using SYMPY '''
		k = 1
		return k

	@staticmethod
	def ph(process, systems, ):
		''' J. Phys. Chem. B, 2004, 108 (46), pp 17886–17892    DOI: 10.1021/jp047349j
			At a pH different from 0, we can correct the free energy of H+ ions by the concentration dependence of the entropy:
						G = H -TS + kT ln(Products/Reactants)   ;    pH = -log[H3O+]
					   G(pH) = −kT ln[H+]= kT ln (10) × pH.
		'''
		''' Reaction conditions are set as symbols using SYMPY '''
		k = 1
		return k

	@staticmethod
	def printdata(rconditions, process, dataname, data):
		maxlen = [max([len(f"{data[r][c]}")+1 for r in range(len(data))]) for c in range(len(data[0]))] # max length per column
		folder = './KINETICS/PROCESSES'
		outputfile = folder + "/" + str(dataname) + ".dat"
		if not pathlib.Path(folder).exists():
			pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
		''' details on conditrions and raction for such process '''
		output = open(outputfile, "w")
		output.write("#")
		for i in rconditions.keys():
			output.write("\t {val:>{wid}s}:".format(wid=len(i), val=i))
			for value in rconditions[i]:
				output.write(" {val:>5.3f}".format(val=value))
		output.write("\n# {:}\t".format(process['kind']))
		species = process['reactants'] + ['>'] + process['ts'] + ['>'] + process['products']
		stoi = process['rstoichio'] + [' '] + process['tsstoichio'] + [' '] + process['pstoichio']
		for i in range(len(species)):
			output.write(" {stoi:}{val:>{wid}s}".format(stoi=stoi[i], wid=len(species[i]) + 1, val=species[i]))
		output.write("\n")
		for row in data:
			for i in range(len(row)):
				if isinstance(row[i], (float, int)):
					output.write("{val:>{wid}.3{c}}".format(val=row[i], wid=maxlen[i],
														c='f' if 1e-5 < np.abs(row[i]) < 1e3 or row[i] == 0. else 'e'))
				else:
					output.write("{val:>{wid}}".format(val=row[i], wid=maxlen[i]))
			output.write("\n")
		output.close()

	@staticmethod
	def getdata(rconditions, process,  datalabel):
		''' substitute the symbolic constants, e.g. h, kb, hc, ... by numeric values'''
		equations = []
		for i in datalabel:
			if not isinstance(process[str(i)], (int, float)):
				equations.append(process[str(i)].subs(constants))
			else:
				equations.append(float(process[str(i)]))
		''' lambdify the equation and substitute rconditions'''
		data = [["# Temperature[K]", *datalabel]]
		if isinstance(rconditions["temperature"], float):
			row = [rconditions["temperature"]]
			for eq in equations:
				row.append(float(sp.lambdify(temp, eq, ['numpy', 'sympy'])(rconditions["temperature"])))
			data.append(row)
		else:
			ramp = [float(i) for i in rconditions["temperature"]]
			for temp_num in np.arange(ramp[0], ramp[1], ramp[2]):
				row = [temp_num]
				for eq in equations:
					a = float(sp.lambdify(temp, eq, ['numpy', 'sympy'])(temp_num))
					value = "{val:>.3{c}}".format(val=a, c='f' if 1e-3 < np.abs(a) < 1e3 or np.abs(a) == 0. else 'e')
					row.append(value)
				data.append(row)
		return data

	@staticmethod
	def barplot(experiment, y_label, data, bar_width):
		icolour = ["b", "r", "c", "g", "m", "y", "grey", "olive", "brown", "pink", "darkgreen", "seagreen", "khaki",
		   "teal"]
		ipatterns = ["///", "...", "xx", "**", "\\", "|", "--", "++", "oo", "OO"]
		gap = 0.7   # gap between group of columns, e.g. processes
		labels = data[0][1:].copy()
		x_label = experiment
		x = np.arange(0, len(labels)) * (1 + gap)
		temps = [float(data[1][0][0]), float(data[1][-1][0])]  # temperatures; first row is for labels
		y = np.array([data[1][0][1:], data[1][-1][1:]], dtype=float)		# elements after 1 but only T0 and T-1
		y[y <= 0] = 1e-20  # to prevent log crash

		fig, ax1 = plt.subplots(figsize=(10, 6), clear=True)
		ax1.set_ylim(bottom=np.min(y) * 0.5, top=np.max(y) * 5)  # Prevents from hitting zero during autoscaling
		for n in range(len(y)):    # processes | the first column is temperature
			offset = (n - len(temps)/2) * bar_width + bar_width/2
			ax1.bar(x + offset, y[n], log=True,
					width=bar_width, hatch=ipatterns[n], color=icolour[n], label=f"{temps[n]} K", alpha=0.7)

		x_limit = [ax1.get_xlim()[0], ax1.get_xlim()[1]]
		#ax1.plot(x_limit, [1e-20, 1e-20], "k-", lw=1.5)
		ax1.set_xlim(x_limit)
		ax1.set_xlabel(x_label, fontsize=16)
		ax1.set_xticks(x)
		ax1.tick_params(axis='x', rotation=0, labelsize=14)
		ax1.set_xticklabels(labels, rotation=0, ha="center")

		#ax1.set_ylim([1e-10, ax1.get_ylim()[1]])
		ax1.set_ylabel(y_label, fontsize=18)
		ax1.tick_params(axis='y', rotation=0, labelsize=16)
		leg_lines, leg_labels = ax1.get_legend_handles_labels()
		legend = ax1.legend(leg_lines[0:len(temps)], leg_labels[0:len(temps)], loc='best', fontsize=16)
		fig.tight_layout()
		plt.ion()
		plt.savefig('./KINETICS/PROCESSES/' + "_".join(experiment.split()) + ".svg", dpi=300, orientation='landscape',
					transparent=True)


class REquations:
	def __init__(self, processes, systems):
		constemperature = {}  # dictionary of equations, e.g. equation[A] = - K1[A][B] + K2[C]
		equation_factors = {}	# dictionary with the factors for rate eqations, i.e. 2 * k_1[A][B]
		print_constemperature = {}
		tpd = {}  # dictionary of equations WITHOUT adsorptions
		tpd_factors = {}
		print_tpd = {}
		'''surfequations = {}
		for name in systems.keys():
			if systems[name]['kind'] == 'surface':
				surfequations[name] = 1'''
		no_ts = []
		for process in processes:
			no_ts.extend(processes[process]['reactants'])
			no_ts.extend(processes[process]['products'])
		for name in no_ts:  # species without TSs
			if systems[name]['kind'] != "surface": # excluding surfaces to avoid DAE
				constemperature[name] = []
				equation_factors[name] = []
				print_constemperature[name] = []
				tpd[name] = []
				tpd_factors[name] = []
				print_tpd[name] = []
				for process in processes:
					eq, rfactor, peq = self.equation(processes[process], name, str(process))
					constemperature[name].append(eq) #if rfactor != 0 else 0
					equation_factors[name].append(rfactor) #if rfactor != 0 else None
					print_constemperature[name] += peq
					if processes[process]["kind"] != "A":
						tpd[name].append(eq)
						tpd_factors[name].append(rfactor)
						print_tpd[name] += peq
				if sum(equation_factors[name]) != 0:
					raise ValueError(f"The species {name} does not considers equilibrium")

		self.all_equations = (constemperature, equation_factors, tpd, tpd_factors)
		self.printdata(print_constemperature, 'Cons_Temperature')
		self.printdata(print_tpd, 'TPD')
		self.overall_stoichiometry(systems, equation_factors)

	@staticmethod
	def equation(process, name, i):  # process is processes[process]; i indicates the process number
		rfactor = 0
		requation = 0
		pequation = []      # list of equations to print
		if name in process['products']:
			for r in range(len(process['products'])):
				if name == process['products'][r]:
					rfactor = float(process['pstoichio'][r])
			requation = process['krate0']
			pequation.append("+" + str(rfactor) + "*K_" + str(i))
			for r in range(len(process['reactants'])):
				requation *= sp.symbols(f"{process['reactants'][r]}") ** process['rstoichio'][r]
				pequation.append("[" + process['reactants'][r] + "]^" + str(process['rstoichio'][r]))
		elif name in process['reactants']:
			for r in range(len(process['reactants'])):
				if name == process['reactants'][r]:
					rfactor = -1*float(process['rstoichio'][r])	# positive number as the "-" is in equation
			requation = process['krate0']
			pequation.append(str(rfactor) + "*K_" + str(i))
			for r in range(len(process['reactants'])):
				requation *= sp.symbols(f"{process['reactants'][r]}") ** process['rstoichio'][r]
				pequation.append("[" + process['reactants'][r] + "]^" + str(process['rstoichio'][r]))
		return requation, rfactor, pequation

	@staticmethod
	def overall_stoichiometry(systems, factors):
		coeff_range = range(1, 2)  # range for vectors combination. Grows very fast so keep the range small.
		names = list(factors.keys())
		stoich_matrix = np.array([factors[name] for name in names])
		ads_matrix = sp.Matrix([factors[name] for name in names if systems[name]['kind'] == 'adsorbate'])

		def lcm(a, b):
			return abs(int(a) * int(b)) // gcd(int(a), int(b))
		def lcm_list(lst):
			return reduce(lcm, lst)
		def normalize_integer_vector(v):
			v = [sp.nsimplify(x) for x in v]	# Convert to rationals
			denoms = [x.as_numer_denom()[1] for x in v]    # Get denominators
			lcm_denom = lcm_list(denoms)
			v_int = [int(x * lcm_denom) for x in v]    # Scale to integers
			nonzero = [abs(x) for x in v_int if x != 0]    # Remove common gcd
			if nonzero:
				g = reduce(gcd, nonzero)
				v_int = [x // g for x in v_int]
			return sp.Matrix(v_int)
		def is_nonzero_reaction(overall_reaction):
			return any(x != 0 for x in overall_reaction)

		ns = ads_matrix.nullspace()    # Nullspace: eliminate surface species
		if not ns:
			raise ValueError("No solution eliminating adsorbate species")
		solutions = []
		for i, vec in enumerate(ns):
			v_int = normalize_integer_vector(vec)
			if any(x < 0 for x in v_int):    # Optional: enforce positive direction
				v_int = -v_int
		overall_reaction = stoich_matrix * v_int
		solutions.append((v_int, overall_reaction))
		filtered = [sol for sol in solutions if is_nonzero_reaction(sol[1])]
		react = []
		prod = []
		for i, name in enumerate(names):
			f =  filtered[0][1][i]
			if f != 0 and systems[name]["kind"] == 'molecule' and systems[name]["pressure0"] > 0:
				react.append(f"{int(f) * -1}.{name}")
			elif f !=0 and systems[name]["kind"] == 'molecule' and systems[name]["pressure0"] == 0:
				prod.append(f"{int(f)}.{name}")
		folder = './KINETICS/EQUATIONS'
		outputfile = folder + "/Overall_Reaction.dat"
		output = open(outputfile, "w")
		output.write("\n\t")
		for i in range(len(react)-1):
			output.write("{val} + ".format(val=react[i]))
		output.write("{val} --> ".format(val=react[-1]))
		for i in range(len(prod)-1):
			output.write("{val} + ".format(val=prod[i]))
		output.write("{val}\n\n".format(val=prod[-1]))
		return list(filtered[0][1])

	@staticmethod
	def printdata(equations, experiment):
		folder = './KINETICS/EQUATIONS'
		outputfile = folder + "/" + experiment + ".dat"
		if not pathlib.Path(folder).exists():
			pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
			os.chmod(folder, 0o755)
		output = open(outputfile, "w")
		output.write("# species/dt\t reactions\n")
		for name in equations.keys():
			output.write(" [{val:>{wid}s}]\dt = ".format(wid=len(name), val=name))
			output.write("{val:>s}".format(val=str.join(" ", equations[name])))
			output.write("\n")
		'''
		for name in surf_equations.keys():
			output.write(" [{val:>{wid}s}] = ".format(wid=len(name), val=name))
			output.write("{val:>s}".format(val=str.join(" ", surf_equations[name])))
			output.write("\n")
		'''
		output.close()


class Profile:
	''' generates a Processes.txt with all the reactions and energies.
		generates a Profile.svg image ONLY with the ODD reaction steps simulating a forward process. '''
	def __init__(self, processes, systems, temp_num=300.0):
		''' Generate an energy Profile from the species in processes '''
		graph = self.build_graph(processes)  # graph_forward only includes ODD steps in processes
		p0 = list(processes.values())[0]
		reference = self.state_key(p0["reactants"], p0["rstoichio"])
		cache = {}
		self.printprofile(processes, systems, temp_num, cache)
		''' there are many REQUISITES to generate a sensible profile 
		self.plot_shared_reactant(graph, reference, systems, temp_num, cache) '''

	@staticmethod
	def build_graph(processes):
		graph = defaultdict(list)        # ALL the steps, i.e. back and forwards.
		for pr in processes.values():
			r = Profile.state_key(pr["reactants"], pr["rstoichio"])
			p = Profile.state_key(pr["products"], pr["pstoichio"])
			ts = None
			if pr["ts"]:
				ts = Profile.state_key(pr["ts"], pr["tsstoichio"])
			graph[r].append({"ts": ts, "products": p})
		return graph

	@staticmethod
	def printprofile(processes, systems, temp_num, cache):
		data = [["#Process", "Kind", "Reaction", f"G(T={temp_num} K)[eV]", "Ea[eV]", "Er[eV]"]]
		for n, pr in enumerate(processes.values()):
			react = Profile.state_key(pr['reactants'], pr['rstoichio'])
			r = Profile.state_label(react)
			e0 = Profile.state_energy(react, systems, temp_num, cache)
			if pr['ts']:
				ts = Profile.state_key(pr['ts'], pr["tsstoichio"])
				t = Profile.state_label(ts)
				ets = Profile.state_energy(ts, systems, temp_num, cache)
				ea = ets - e0
			else:
				t = ''
				if pr['kind'] == "A":
					ts = Profile.state_key(pr['reactants'], pr['rstoichio'])
					ets = Profile.state_energy_ts(ts, systems, temp_num)
					ea = ets - e0
				elif pr['kind'] == "D":
					ts = Profile.state_key(pr['products'], pr['pstoichio'])
					ets = Profile.state_energy_ts(ts, systems, temp_num)
					ea = ets - e0
				else:
					ets = None
			pro = Profile.state_key(pr["products"], pr["pstoichio"])
			p = Profile.state_label(pro)
			ep = Profile.state_energy(pro, systems, temp_num, cache)
			er = ep - e0
			reaction = " > ".join([r, t, p])
			ereaction = " > ".join([f"{e0:.2f}", f"{ets:.2f}" if ets is not None else "", f"{ep:.2f}"])
			data.append([f"{n+1}", f"{pr['kind']}", f"{reaction}", f"{ereaction}",
						 f"{ea:5.2f}" if ets is not None else "  ---", f"{er:.2f}"])
		maxlen = [max([len(f"{data[r][c]}")+1 for r in range(len(data))]) for c in range(len(data[0]))] # max length per column
		outputfile = "Processes.txt"
		output = open(outputfile, "w")
		for i in range(len(data)):
			for j in range(len(data[i])):
				output.write("{val:{wid}s} ".format(wid=maxlen[j], val=data[i][j]))
			output.write("\n")
		output.close()

	@staticmethod
	def state_key(names, stoich):
		return tuple(sorted(zip(names, stoich)))

	@staticmethod
	def state_energy(state, systems, temp_num, cache):
		if state in cache:
			return cache[state]
		e = 0.0
		if state:
			for name, stoi in state:
				expr = systems[name]["energy3d"].subs(constants)
				f = sp.lambdify(temp, expr, ("numpy", "sympy"))
				e += stoi * float(f(temp_num))
		cache[state] = e
		return e

	@staticmethod
	def state_energy_ts(state, systems, temp_num):
		e = 0.0
		if state:
			for name, stoi in state:
				if systems[name]['kind'] == "molecule":
					expr = systems[name]["energy2d"].subs(constants)
				else:
					expr = systems[name]["energy3d"].subs(constants)
				f = sp.lambdify(temp, expr, ("numpy", "sympy"))
				e += stoi * float(f(temp_num))
		return e

	@staticmethod
	def state_label(state):
		item = []
		for name, stoi in state:
			try:
				stoi = int(stoi)
			except:
				pass
			item.append(f"{stoi}.{name}")
		return " + ".join(item)

	@staticmethod
	def plot_shared_reactant(graph, reference, systems, temp_num, cache):
		def note(text, xy):
			ax.annotate(f"{text}", xy=(xy[0], xy[1]), xytext=(0, -5),  # 15 points below
						textcoords="offset points", ha="center", fontsize=12, rotation=90, va='top',
						bbox=dict(boxstyle="round, pad=0.1", fc="white", ec="white", lw=1, alpha=0.8))

		e0 = Profile.state_energy(reference, systems, temp_num, cache)
		xtick = 0.0
		xticks = [xtick]
		xlabels = [Profile.state_label(reference)]
		step_length = 1.0/2 #    length of every minima

		fig, ax = plt.subplots(figsize=(10, 6), clear=True)  # prepares a figure

		ax.plot([xtick - step_length, xtick + step_length], [0.0, 0.0], linestyle='-', lw=5, color="k", alpha=1)
		#ax.text(xtick, -0.2, Profile.state_label(reference), ha="center", va="top", fontsize=14)
		note(Profile.state_label(reference), [xtick, 0.0])

		for key in graph.keys():
			branches = graph[key]
			n = len(branches)
			offsets = [0.0] #np.linspace(-0.8, 0.8, n) if n > 1 else [0.0]     # horizontal spacing for branches

			x_react = xtick  # x positions
			x_ts = xtick + 1
			x_prod = xtick + 2
			y_react = Profile.state_energy(key, systems, temp_num, cache) - e0

			for offset, step in zip(offsets, branches):
				prod = step["products"]
				y_prod = Profile.state_energy(prod, systems, temp_num, cache) - e0
				ax.plot([x_prod - step_length, x_prod + step_length], [y_prod, y_prod], lw=2.5, color="k")
				# ax.text(x_prod, y_prod - 0.2, Profile.state_label(prod), ha="center", va="top", fontsize=14)
				note(Profile.state_label(prod), [x_prod, y_prod])
				xticks.append(x_prod)
				xlabels.append(Profile.state_label(prod))

				if step["ts"]:
					ts = step["ts"]
					y_ts = Profile.state_energy(ts, systems, temp_num, cache) - e0
					ea = y_ts - y_react
					# spline through TS
					x0 = [x_react +step_length, x_ts, x_prod -step_length]
					y0 = [y_react, y_ts, y_prod]
					spl = splrep(x0, y0, k=2)
					xt = np.linspace(min(x0), max(x0), 60)
					ax.plot(xt, splev(xt, spl), linestyle='--', lw=1.0, color="k", alpha=1)
					#ax.text(x_ts, y_ts + 0.2, f"$E_A$ = {ea:.2f} eV", ha="center", fontsize=14)
					ax.annotate(f"$E_{{A}}$ = {ea:.2f} eV", xy=(x_ts, y_ts), xytext=(0, 15),  # 15 points below
								textcoords="offset points", ha="center", fontsize=12)
				else:   # linking reactants and products
					ax.plot([x_react +step_length, x_prod -step_length], [y_react, y_prod], linestyle='--',
							lw=1.0, color="k", alpha=1)
			xtick += 2  # connecting with the previous process

		ax.axhline(0, lw=1, ls=":", color="k")
		ax.tick_params(labelsize=14)
		ax.set_xticks([])
		#ax.set_xticks(xticks)
		#ax.set_xticklabels(xlabels, rotation=25, ha="right")  # rotation=0, ha="center")
		ax.set_ylabel(f"$\Delta G_{{T={temp_num}\,K}}$ (eV)", fontsize=16)
		fig.tight_layout()
		plt.ion()
		plt.savefig("Profile.svg", dpi=300, orientation='landscape', transparent=True)
