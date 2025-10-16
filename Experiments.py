"""
    This script builds on the perl version by A.Roldan.

"""

import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp


class ConsTemperature:
    def __init__(self, rconditions, systems, equations, surf_equations):
        t , temp = sp.symbols('time temperature')
        ics, species, ics_surf, surf = self.initial_species(systems)
        ode = []
        rhs = []
        sp_species = []
        for name in species:        # lists of names with the order of ics
            sp_species.append(sp.Function(name)(t))
            ode.append(sp.Eq(sp.diff(name, t), equations[name]))
            rhs.append(ode[-1].rhs)     # rhs for SciPy
        ''' ODEs do not contain algrebraic equations such as those for coverage.
        To solve Differention-Algebraic Equations (DAE) use CasADi or Assimulo
        for name in surf:
            sp_species.append(sp.Function(name)(t))
            ode.append(sp.Eq(name, surf_equations[name]))
            rhs.append(ode[-1].rhs)     # rhs for SciPy
        '''
        print(species)
        print(surf_equations)

        conditions = [t, temp]
        # Convert symbolic RHS into numerical function
        f_ode = sp.lambdify((*conditions, *species), rhs,"numpy")
        y = [name for name in systems.keys()]   #__ especies

        # Time grid
        t_span = (rconditions["time"][0], rconditions["time"][1])
        t_eval = np.linspace(*t_span, rconditions["time"][2])
        # Temperature grid
        temp_span = (rconditions["temperature"][0], rconditions["temperature"][1])
        temp_eval = np.linspace(*temp_span, rconditions["temperature"][2])

        # define ode function compatible with solve_ivp
        def ode_system(t, species, temperature):
            return f_ode(t, temperature, *species)  # ics in concentrations at t=0

        # integrate at different temperatures
        for temp in temp_eval:
            sol = solve_ivp(ode_system, t_span, ics, t_eval=t_eval, args=(temp, ), rtol=1e-8, atol=1e-10)
            print(sol)
            exit()
            #self.printdata("Constant_Temperature", temp, sol)

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
            ics_surf.append(surf0[s])
            surf.append(s)
        return ics, species, ics_surf, surf

    '''
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
        '''
