"""
    This script builds on the perl version by A.Roldan.

"""

import os, pathlib
import sympy as sp


def printdata(rconditions, process, constants, datalabel, dataname):
    folder = './KINETICS/PROCESSES'
    outputfile = folder + "/" + str(dataname) + ".dat"
    if not pathlib.Path(folder).exists():
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
        os.chmod(folder, 0o755)
    if not pathlib.Path(outputfile).exists():
        output = open(outputfile, "w")
        output.write("# Temperature[K]")
        for i in datalabel:
            output.write(" {val:>{wid}s}".format(wid=len(i)+3, val=i))
        output.write("\n")
        temp = sp.symbols("temperature")
        if isinstance(rconditions["temperature"], float):
            output.write(" {val:>{wid}.1f}".format(wid=len("Temperature[K]"), val=rconditions["temperature"]))
            for i in datalabel:
                value = sp.lambdify(temp, process[str(i)])(float(rconditions["temperature"]))
                output.write(" {val:>{wid}.3{c}}".format(wid=len(i)+3, val=value, c='e' if value > 1e3 else 'f'))
        else:
            ramp = [int(i) for i in rconditions["temperature"]]
            for t in range(ramp[0], ramp[1], ramp[2]):
                output.write(" {val:>{wid}.1f}".format(wid=len("Temperature[K]"), val=t))
                for i in datalabel:
                    value = sp.lambdify(temp, process[str(i)])(t)
                    output.write(" {val:>{wid}.3{c}}".format(wid=len(i)+3, val=value, c='e' if value > 1e3 else 'f'))
                output.write("\n")
        output.close()

class RConstants:
    def __init__(self, rconditions, systems, constants, processes, restricted_arg):
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
        for process in processes:
            processes[process]["activation"] = self.activation(processes[process], systems)
            if processes[process]['kind'] == 'A':
                processes[process]["sticky"] = self.sticky(processes[process], systems, constants)
                processes[process]["arrhenius"] = self.arrhenius(processes[process], systems, constants)
                processes[process]["ktunneling"] = self.tunneling(processes[process], systems, constants)
                processes[process]['krate0'] = (processes[process]["sticky"] * processes[process]["arrhenius"] *
                                                processes[process]["ktunneling"])
                # units of m*kg^-1*s^-1 |when multiplied by Pa = s^-1
                datalabel = ["activation", "sticky", "arrhenius", 'ktunneling', "krate0"]
                printdata(rconditions, processes[process], constants, datalabel, "Process"+str(process))
            else:
                processes[process]["arrhenius"] = self.arrhenius(processes[process], systems, constants)
                processes[process]["ktunneling"] = self.tunneling(processes[process], systems, constants)   # NO units
                processes[process]['krate0'] = (processes[process]["arrhenius"] *
                                                sp.exp(-processes[process]['activation']/
                                                       (constants['kb']*temp*constants['JtoeV'])) *
                                                processes[process]["ktunneling"])   # units s^-1
                datalabel = ["activation", "arrhenius", 'ktunneling', "krate0"]
                printdata(rconditions, processes[process], constants, datalabel, "Process"+str(process))
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
    def sticky(process, systems, constants):
        ''' the sticky coefficient is evaluate as the reduction of degrees of freedom, i.e. from a 3D free molecule
        to a 2D trapped molecule moving parallel to the surface (being the third degree the reaction coordinate) '''
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
        qts = 1     # total partition function for transition states
        qr = 1     # total partition function for reactants
        if len(process['ts']) > 0:
            for i in range(len(process['ts'])):
                qts *= systems[process['ts'][i]]['q3d']**process['tsstoichio'][i]
            for i in range(len(process['reactants'])):
                qr *= systems[process['reactants'][i]]['q3d'] ** process['rstoichio'][i]
        else:
            for i in range(len(process['reactants'])):
                if systems[process['reactants'][i]]['q2d']:
                    qts *=  systems[process['reactants'][i]]['q2d']**process['rstoichio'][i]
                else:
                    qts *=  systems[process['reactants'][i]]['q3d']**process['rstoichio'][i]
            for i in range(len(process['reactants'])):
                qr *=  systems[process['reactants'][i]]['q3d']**process['rstoichio'][i]
        return (qts/qr)*sp.exp(-process['activation']/(constants['kb']*temp*constants['JtoeV']))

    @staticmethod
    def arrhenius(process, systems, constants, restricted_arg):
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
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
            arrhenius = area * 1/sp.sqrt(2*sp.pi*mass*constants['kb']*temp)
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
            arrhenius = constants['kb']*temp/constants['h'] * qts/qr    # units of s^-1
        return arrhenius

    @staticmethod
    def tunneling(process, systems, constants):
        ''' Second order harmonic Wigner approach to shallow quantum tunneling valid for
        vast numbers of reaction including surface-catalysed --> DOI: 10.1039/C4CP03235G '''
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
        k = 1
        for i in range(len(process['ts'])):
            for f in systems[process['ts'][i]]['ifreq']:
                k = 1 + 1/24 * (constants['hc']* f /(2*sp.pi * constants['kb']*temp))**2
        return k

    @staticmethod
    def electric(process, systems, constants):
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
        temp = sp.symbols("temperature")
        k = 1
        return k

    @staticmethod
    def ph(process, systems, constants):
        ''' J. Phys. Chem. B, 2004, 108 (46), pp 17886–17892    DOI: 10.1021/jp047349j
            At a pH different from 0, we can correct the free energy of H+ ions by the concentration dependence of the entropy:
        		        G = H -TS + kT ln(Products/Reactants)   ;    pH = -log[H3O+]
                       G(pH) = −kT ln[H+]= kT ln (10) × pH.
        '''
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
        k = 1
        return k


class ODE:
    def __init__(self, processes, systems, restricted_arg):
        ''' Coverage (cov_i) and Pressures (pro_i) at time t are set as symbols using SYMPY '''
        for name in systems.keys():     # species
            if kind == 'molecule':
                pretemp = sp.symbols("temperature")
        ode = {}
        for name in systems.keys():     # species
            ode[name] = 0
            for process in processes:
                species = 0
                for r in range(len(processes[process]['reactants'])):
                    species =



                for i, species in enumerate(processes[process]['products']):
                    if name == species:
                        ode[name] +=  XX processes[process]['krate0']*





        self.ode = ode



