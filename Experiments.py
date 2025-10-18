"""
    This script builds on the perl version by A.Roldan.

"""

import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp


class ConsTemperature:
    def __init__(self, rconditions, systems, equations):
        t , temp = sp.symbols('time temperature')
        ics, species = self.initial_species(systems)
        rhs = []
        for name in species:        # lists of names with the order of ics
            rhs.append(equations[name])     # rhs for SciPy
        ''' ODEs do not contain algrebraic equations such as those for coverage.
        To solve Differention-Algebraic Equations (DAE) use CasADi or Assimulo
        for name in surf:
            sp_species.append(sp.Function(name)(t))
            ode.append(sp.Eq(name, surf_equations[name]))
            rhs.append(ode[-1].rhs)     # rhs for SciPy
        '''
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
            sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval, args=(temp, ), rtol=1e-8, atol=1e-10)
            self.printdata(rconditions, species, "Cons_Temperature", [temp], sol)

            import matplotlib.pyplot as plt
            for i,name in enumerate(species):
                plt.plot(sol.t, sol.y[i], label=f'{name}')
            plt.legend()
            plt.xlabel('Time (s)')
            plt.ylabel('Concentration / Coverage')
            plt.show()


    # tpd
    # Reaction control
    # Rate
    # Selectivity control

    @staticmethod
    def initial_species(systems):  # process is processes[process]
        t = sp.symbols('time')
        ics = []    # list of initial concentrations in the order of systems[names]
        species = []    # list of species in the order on systems
        ics_surf = []   # list of initial free sites
        surf = []   # list of surfaces in the order of ics
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
    def printdata(self, rconditions, species, experiment, conditions, sol):
        folder = './KINETICS/DATA'
        outputfile = folder + "/" + str(experiment) + ".dat"
        if not pathlib.Path(folder).exists():
            pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
        if not outputfile.exist():
            output = open(outputfile, "w")
            # heading
            output.write("# ")
            for key in rconditions.keys():
                output.write(key)
            for name in species:
                output.write(name)
            output.write("\n")
            output.close()
        output = open(outputfile, "a")
        for i in range(len(sol.t)):
            # conditions
            for c in conditions:
                output.write(" {val:>5.3f}".format(val=c))
            output.write(" {val:>5.3{c}".format(val=sol.t[i], c='e' if 1e-3 < sol.t[i] > 1e3 else 'f'))
            # data
            for n in range(len(species)):
                output.write(" {val:>5.5{c}}".format(val=sol.y[i][n], c='e' if 1e-3 < sol.y[i][n] > 1e3 else 'f'))
            output.write("\n")
        output.write("\n")
        output.close()
