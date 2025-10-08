"""
    This script builds on the perl version by A.Roldan.

"""

import os, pathlib
import sympy as sp

from Main import rconditions

{'temperature': [300.0, 320.0, 10.0], 'time': [0, 100.0, 10.0]}


class Experiments:
    def __init__(self, rconditions, systems, equations):
        t , temp = sp.symbols('time temperature')
        ics = self.initial_species(systems)
        species = [sp.Function(name)(t) for name in systems.keys()]
        ode = []
        rhs = []
        for name in species:
            ode.append(sp.Eq(sp.diff(name, t), equations[name]))
            rhs.append(ode[-1].rhs)     # rhs for SciPy

        conditions = [t, temp]

        # Convert symbolic RHS into numerical function
        f_ode = sp.lambdify((*conditions, *species), rhs,"numpy")
        y = [name for name in systems.keys()]
        # Time grid
        t_span = (rconditions["time"][0], rconditions["time"][1])
        t_eval = np.linspace(*t_span, rconditions["time"][2])
        for temp in range(rconditions["temperature"]):
            sol = solve_ivp(   lambda t, temp, y: f_ode(t, temp, *y), t_span, ics, t_eval=t_eval, rtol=1e-8, atol=1e-10)
            print(sol)
            #self.printdata("Constant_Temperature", temp, sol)

    # Reaction control
    # Rate
    # Selectivity control

    @staticmethod
    def initial_species(self, systems):  # process is processes[process]
        ics = {}
        surf0 = {}
        sites_name = list(set([systems[name]['site'] for name in systems.keys()]))
        for s in sites_name:
            surf0[s] = 1
        for name in systems.keys():
            if systems[name]['kind'] == "molecule":
                ics[systems[name].subs(t, 0)] = systems[name]["pressure0"]
            elif systems[name]['kind'] == "adsorbate":
                ics[systems[name].subs(t, 0)] = systems[name]["coverage0"]
                surf0[systems[name]['site']] -= systems[name]["coverage0"] * systems[name]["nsites2"]
        for s in surf0.keys():
            ics[s.subs(t, 0)] = surf0[s]
        return ics

    @staticmethod
    def printdata(self, experiment, temp, sol):
        folder = './KINETICS/DATA'
        outputfile = folder + "/" + str(experiment) + ".dat"
        if not pathlib.Path(folder).exists():
            pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
        if not pathlib.Path(outputfile).exists():
            output = open(outputfile, "w")
            # heading
            output.write("# ")
            for key in rconditions.keys():
                output.write(key)
            for name in systems.keys():
                output.write(name)
            output.write("\n")
        # data
        for value in sol:
            output.write(" {val:>5.5{c}}".format(val=value, c='e' if 1e-3 < value > 1e3 else 'f'))

