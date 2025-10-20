"""
	This script builds on the perl version by A.Roldan.

"""

import os, pathlib
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp


class ConsTemperature:
	def __init__(self, rconditions, systems, equations):
		t , temp = sp.symbols('time temperature')
		ics, species = self.initial_species(systems)
		rhs = []
		for name in species:        # lists of names with the order of ics
			if name in equations.keys():
				rhs.append(equations[name])     # rhs for SciPy

		conditions = [t, temp]
		t_span = (rconditions["time"][0], rconditions["time"][1])   # time grid
		t_eval = np.arange(*t_span, rconditions["time"][2])
		temp_span = (rconditions["temperature"][0], rconditions["temperature"][1])      # temperature grid
		temp_eval = np.arange(*temp_span, rconditions["temperature"][2])

		f_ode = sp.lambdify((*conditions, *species), rhs,"numpy")      # Convert symbolic RHS into numerical function

		# define ode function compatible with solve_ivp
		def ode_system(t, species, temp):
			return f_ode(t, temp, *species)  # ics in concentrations at t=0

		for temp in temp_eval:     # integrate at different temperatures
			sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval, args=(temp, ), method='BDF', rtol=1e-8, atol=1e-10)
			data = list(zip(sol.t, *sol.y))
			self.printdata(dict(rconditions), species, str("Cons_Temperature"), [temp], data)
			'''
			import matplotlib.pyplot as plt
			for i,name in enumerate(species):
				plt.plot(sol.t, sol.y[i], label=f'{name}')
			plt.legend()
			plt.xlabel('Time (s)')
			plt.ylabel('Concentration / Coverage')
			plt.show()
			'''


	# tpd
	# Reaction control
	# Rate
	# Selectivity control

	@staticmethod
	def initial_species(systems):  # process is processes[process]
		ics = []    # list of initial concentrations in the order of systems[names]
		species = []    # list of species in the order on systems
		surf0 = {}
		sites_name = []
		for name in systems.keys():
			if systems[name]['kind'] == 'surface':
				sites_name.extend(systems[name]['sites'])
		for s in list(set(sites_name)):
			surf0[s] = 1
		for name in systems.keys():
			if systems[name]['kind'] == "molecule":
				ics.append(systems[name]["pressure0"])
				species.append(name)
			elif systems[name]['kind'] == "adsorbate":
				ics.append(systems[name]["coverage0"])
				''' adsorbates have only one kind of adsorption site per system '''
				surf0[systems[name]['sites'][0]] -= systems[name]["coverage0"] * systems[name]["nsites"]
				species.append(name)
		for s in surf0.keys():
			ics.append(surf0[s])
			species.append(s)
		return ics, species

	@staticmethod
	def printdata(rconditions, species, experiment, conditions, data):
		folder = './KINETICS/DATA'
		outputfile = folder + "/" + str(experiment) + ".dat"
		headings = []
		if not pathlib.Path(folder).exists():
			pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
		output = open(outputfile, "a")
		# heading
		output.write("#")
		maxlen = max(len(f"{v:.3f}") for row in data for v in row)+2
		for key in rconditions.keys():      # includes Vext, pH, Temp, time, ...
			wid = len(key)+1    # for the space
			if wid < maxlen:
				wid = maxlen
			output.write(" {val:>{wid}}".format(wid=wid, val=key))
			headings.append(wid)
		for i, name in enumerate(species):
			output.write(" {val:>{wid}}".format(wid=maxlen, val=name))
			if isinstance(i/3, int):
				output.write("\t")
		output.write("\n")

		for i, row in enumerate(data):
			# conditions
			for j, con in enumerate(conditions):    # includes Vext, pH, Temp, ... BUT not time!
				output.write(" {val:>{wid}.1f}".format(wid=headings[j]+1, val=con))
			# data
			for n, value in enumerate(row):
				output.write(" {val:>{wid}.3{c}}".format(wid=maxlen, val=value, c='e' if 1e-3 < value > 1e3 else 'f'))
				if isinstance(n/3, int):
					output.write("\t")
			output.write("\n")
		output.write("\n")
		output.close()
