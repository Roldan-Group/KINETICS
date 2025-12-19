"""
	This script builds on the perl version by A.Roldan.

"""
import os, pathlib
import sympy as sp
import numpy as np
from Symbols_def import t, temp, h, kb, hc, JtoeV, constants
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from collections import defaultdict



class RConstants:
	def __init__(self, rconditions, systems, processes, restricted_arg):
		''' Reaction conditions are set as symbols using SYMPY '''
		for process in processes:
			processes[process]["activation"] = self.activation(processes[process], systems)
			if processes[process]['kind'] == 'A':
				processes[process]["sticky"] = self.sticky(processes[process], systems, )
				processes[process]["arrhenius"] = self.arrhenius(processes[process], systems, restricted_arg)
				processes[process]["ktunneling"] = self.tunneling(processes[process], systems, )
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
		self.processes = processes

	@staticmethod
	def activation(process, systems):
		''' the activation energy in adsorption and desorption processes is considered as the difference between
		a state in which the molecule has only two degrees of freedom (being the third degree the reaction coordinate)
		and the reactants '''
		ets = 0     # total energy for transition states
		if len(process['ts']) > 0:
			for i in range(len(process['ts'])):
				ets += process['tsstoichio'][i] * systems[process['ts'][i]]['energy3d']
		elif 'molecule' in [systems[i]['kind'] for i in process['reactants']]:  # for adsorption processes
			''' This elif considers molecules in reactants as in adsorption processes '''
			for i in range(len(process['reactants'])):
				if 'energy2d' in systems[process['reactants'][i]]:
					ets += process['rstoichio'][i] * systems[process['reactants'][i]]['energy2d']
				else:
					ets += process['rstoichio'][i] * systems[process['reactants'][i]]['energy3d']
		elif 'molecule' in [systems[i]['kind'] for i in process['products']]:   # for desorption processes
			''' This elif considers molecules in products as in desorption processes '''
			for i in range(len(process['products'])):
				if 'energy2d' in systems[process['products'][i]]:
					ets += process['pstoichio'][i] * systems[process['products'][i]]['energy2d']
				else:
					ets += process['pstoichio'][i] * systems[process['products'][i]]['energy3d']
		else:
			''' In the rare case that there is to transtion state and none of the reactants is a molecule the energy of 
			the transition states will be the energy of the final state.'''
			for i in range(len(process['products'])):
				ets += process['pstoichio'][i] * systems[process['products'][i]]['energy3d']
		er = 0      # total energy for reactants
		for i in range(len(process['reactants'])):
			er += process['rstoichio'][i] * systems[process['reactants'][i]]['energy3d']
		return ets - er

	@staticmethod
	def sticky(process, systems, ):
		''' the sticky coefficient is evaluate as the reduction of degrees of freedom, i.e. from a 3D free molecule
		to a 2D trapped molecule moving parallel to the surface (being the third degree the reaction coordinate) '''
		''' Reaction conditions are set as symbols using SYMPY '''
		qts = 1     # total partition function for transition states
		qr = 1     # total partition function for reactants
		if len(process['ts']) > 0:
			for i in range(len(process['ts'])):
				qts *= systems[process['ts'][i]]['q3d']**process['tsstoichio'][i]
			for i in range(len(process['reactants'])):
				qr *= systems[process['reactants'][i]]['q3d'] ** process['rstoichio'][i]
		else:
			for i in range(len(process['reactants'])):
				if 'q2d' in systems[process['reactants'][i]].keys():
					qts *=  systems[process['reactants'][i]]['q2d']**process['rstoichio'][i]
				else:
					qts *=  systems[process['reactants'][i]]['q3d']**process['rstoichio'][i]
			for i in range(len(process['reactants'])):
				qr *=  systems[process['reactants'][i]]['q3d']**process['rstoichio'][i]
		return (qts/qr)*sp.exp(-process['activation']/(kb*temp*JtoeV))

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
			else:
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
				output.write(" {val:>{wid}}".format(wid=maxlen[i], val=row[i]))
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
				value = float(sp.lambdify(temp, eq, ['numpy', 'sympy'])(rconditions["temperature"]))
				c = 'e' if value > 1e3 or np.abs(value) < 1e-2 else 'f'
				row.append(f"{value:.3{c}}")
			data.append(row)
		else:
			ramp = [float(i) for i in rconditions["temperature"]]
			for t in np.arange(ramp[0], ramp[1], ramp[2]):
				row = [t]
				for eq in equations:
					value = float(sp.lambdify(temp, eq, ['numpy', 'sympy'])(t))
					c = 'e' if value > 1e3 or np.abs(value) < 1e-2 else 'f'
					row.append(f"{value:.3{c}}")
				data.append(row)
		return data


class REquations:
	def __init__(self, processes, systems):
		constemperature = {}  # dictionary of equations, e.g. equation[A] = - K1[A][B] + K2[C]
		print_constemperature = {}
		tpd = {}  # dictionary of equations WITHOUT adsorptions
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
				constemperature[name] = 0
				print_constemperature[name] = []
				tpd[name] = 0
				print_tpd[name] = []
				for process in processes:
					eq, peq = self.equation(processes[process], name, str(process))
					constemperature[name] += eq
					print_constemperature[name] += peq
					if processes[process]["kind"] != "A":
						tpd[name] += eq
						print_tpd[name] += peq
				'''
				if systems[name]['kind'] == 'adsorbate':
					for s in surfequations:
						if systems[s]['sites'] == systems[name]['sites']:
							\''' the coverage is multiplied by the number of sites the adsorbate occupies because the
							site is defined as the "area" of a single site and a molecule may occupy more than 
							one, e.g. O2 occupies to O sites.\'''
							surfequations[s] -= sp.symbols(f"{name}") * systems[name]['nsites']'''

		self.constemperature = constemperature
		#self.surfequations = surfequations
		self.tpd = tpd

		self.printdata(print_constemperature, 'Cons_Temperature')
		self.printdata(print_tpd, 'TPD')

	@staticmethod
	def equation(process, name, i):  # process is processes[process]; i indicates the process number
		equation = 0
		pequation = []      # list of equations to print
		if name in process['products']:
			for r in range(len(process['products'])):
				if name == process['products'][r]:
					rfactor = process['pstoichio'][r]
					prfactor = process['pstoichio'][r]
			pequation.append("   +" + str(prfactor) + " * K_" + str(i))
			equation += rfactor * process['krate0']
			for r in range(len(process['reactants'])):
				equation *= sp.symbols(f"{process['reactants'][r]}") ** process['rstoichio'][r]
				pequation.append("[" + process['reactants'][r] + "]^" + str(process['rstoichio'][r]))
		elif name in process['reactants']:
			for r in range(len(process['reactants'])):
				if name == process['reactants'][r]:
					rfactor = process['rstoichio'][r]
					prfactor = process['rstoichio'][r]
			pequation.append("   -" + str(prfactor) + " * K_" + str(i))
			equation -= rfactor * process['krate0']
			for r in range(len(process['reactants'])):
				equation *= sp.symbols(f"{process['reactants'][r]}") ** process['rstoichio'][r]
				pequation.append("[" + process['reactants'][r] + "]^" + str(process['rstoichio'][r]))
		return equation, pequation

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
				ets = None
				ea = None
			pro = Profile.state_key(pr["products"], pr["pstoichio"])
			p = Profile.state_label(pro)
			ep = Profile.state_energy(pro, systems, temp_num, cache)
			er = ep - e0
			reaction = " > ".join([r, t, p])
			ereaction = " > ".join([f"{e0:.2f}", f"{ets:.2f}" if ets is not None else "", f"{ep:.2f}"])
			data.append([f"{n+1}", f"{pr['kind']}", f"{reaction}", f"{ereaction}",
						 f"{ea:.2f}" if ets is not None else "---", f"{er:.2f}"])
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
