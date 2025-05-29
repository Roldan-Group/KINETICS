"""
    This script reads and executes the packages for KINETICS

        by Alberto Roldan
"""

import os, sys
from sympy as sp
from Thermodynamics import PartitionFunctions


constants = {"h": 6.62607588705515e-34,     # kg m^2 s^-1 == J s
             "kb": 1.380658045430573e-23,   # J K^-1
             "c": 299792458,      # m s^-1
             "hc": 1.98644586e-23,      #  J cm (to balance the frequencies units of cm^-1)
             "R": 8.31446261815324,     # J⋅K−1⋅mol−1 ----------needed?
             "Av": 6.022139922973909e23,    # mols^-1
             "Fa": 96485.3321233100184,     # C⋅mol^−1
             "qelectron": -1.60217648740e-19,   # C
             "JtoeV": 6.24150974e18     # 1 J = 6.24..e18 eV
             }

''' Read input file with the systems and kinetic simulation conditions'''
def mkread(inputfile):
    rconditions = {}    # reaction conditions
    processes = {}      # reaction processes
    process = 0         # process number, key of processes starting from 1
    systems = {}        # species
    for line in open(inputfile).readlines():
        if len(line) >= 1 and line.startswith("#") is False:
            line = line.split("=")
            head = line[0].strip()
            tail = []
            i = 0
            while line[1].split()[i][0] != "#" and i < len(line[1].split()):
                tail.append(line[1].split()[i].strip())
                i += 1
            ''' Reaction conditions are stored in a dictionary (rconditions), including:
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
                    rconditions["vext"] = [float(i) for i in tail]  # ramp
            if head == "PH" or head == "pH":
                if len(tail) == 1:
                    rconditions["ph"] = float(tail[0])    # constant
                else:
                    rconditions["ph"] = [float(i) for i in tail]  # ramp
            if head == "TEMP" or head == "TEMPERATURE":
                if len(tail) == 1:
                    rconditions["temperature"] = float(tail[0])     # constant
                else:
                    rconditions["temperature"] = [float(i) for i in tail]   # ramp
            if head == "TIME" or head == "Time":
                rconditions["time"] = [0] + [float(i) for i in tail]  # initial time is 0
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
                processes[str(process + 1)] = {}    # dictionary of items for process
                processes[str(process + 1)]["kind"] = str(tail[0][0])  # kind of process
                # getting reactants
                i = 1
                reactants = []
                rstoichio = []
                while tail[i] != ">" or tail[i] != ">>":
                    try:
                        rstoichio.append(float(tail[i]))
                    except ValueError:
                        if tail[i] != "+":
                            reactants.append(str(tail[i]))
                    i += 1
                ''' It may be the case stoichiometry = 1 is neglected.
                Then, it is check that the number of species is the same than
                the number of stoichiometries or it is added a stoichio = 1.'''
                while len(reactants) > len(rstoichio):
                    rstoichio.append(float(1))
                processes[str(process + 1)]["reactants"] = reactants  # reactants
                processes[str(process + 1)]["rstoichio"] = rstoichio
                # getting ts
                if tail[i] == ">":
                    i += 1  # added due to the ">"
                    ts = []
                    while tail[i] != ">":
                        if tail[i] != "+":
                            ts.append(str(tail[i]))
                        i += 1
                    processes[str(process + 1)]["ts"] = ts  # ts
                # getting products
                i += 1  # added due to the ">" or ">>"
                products = []
                pstoichio = []
                for j in range(len(tail[i:])):
                    try:
                        pstoichio.append(float(tail[j]))
                    except ValueError:
                        if tail[j] != "+":
                            products.append(str(tail[j]))
                ''' It may be the case that stoichiometry = 1 is neglected.
                Then, it is check that the number of species is the same than
                the number of stoichiometries or 1 it is added to stoichio.'''
                while len(products) > len(pstoichio):
                    pstoichio.append(float(1))
                processes[str(process + 1)]["products"] = products  # products
                processes[str(process + 1)]["pstoichio"] = pstoichio
            ''' Systems, i.e. the species involved, in a nested dictionary (systems) with key = name,
            including:
                - number of adsorbates (nadsorbates), accounting for the systems' coverage.
                - kind (surface, molecule, adsorbed)
                - path the input, i.e. QM data, (syspath)
                - path to frequencies (freqpath)
                - energy (energy)
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
                nadsorbates = 1
                try:
                    nadsorbates = int(tail[1])    # number of adsorbates in system
                except ValueError:
                        pass
                systems[name][str(nadsorbates)] = {}
            if head == "SYSPATH":
                systems[name][nadsorbates]["syspath"] = str(tail[0])
            if head == "FREQPATH":
                systems[name][nadsorbates]["freqpath"] = str(tail[0])
            if head == "E0":             # species energy
                systems[name][nadsorbates]["energy0"] = float(tail[0])
            if head == "ISITES":
                sites = []                  # catalyst adsorption sites
                for i in tail:
                    try:
                        a = float(i)
                    except ValueError:
                        sites.append(str(i))
                systems[name][nadsorbates]["site"] = sites   # catalyst adsorption sites
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

                print(freq)

                systems[name][nadsorbates]["freq3d"] = sorted(freq, reverse=True)   # species frequencies 3D
            if head == "FREQ2D":
                freq2d = []                   # species frequencies only considering x and y displacements
                for i in tail:
                    try:
                        freq2d.append(float(i))
                    except ValueError:
                        pass
                systems[name][nadsorbates]["freq2d"] = sorted(freq2d, reverse=True)   # species frequencies only considering x and y displacements
            if head == "IMASS":
                systems[name][nadsorbates]["mass"] = float(tail[0])
            if head == "INATOMS":
                systems[name][nadsorbates]["natoms"] = int(tail[0])
            if head == "SYMFACTOR":
                systems[name][nadsorbates]["symfactor"] = int(tail[0])
            if head == "INERTIA":
                inertia = []          # molecular inertia moment(s)
                for i in tail:
                    try:
                        inertia.append(float(i))
                    except ValueError:
                        pass
                systems[name][nadsorbates]["inertia"] = inertia   # molecule's inertia moment(s) in Kg/m^2
                if isinstance(inertia, float):      # molecular linearity
                    systems[name][nadsorbates]["linear"] = "yes"
                elif inertia[0] ==inertia[1] or inertia[0] == inertia[2] or inertia[1] == inertia[2]:
                    systems[name][nadsorbates]["linear"] = "yes"
                else:
                    systems[name][nadsorbates]["linear"] = "no"
            if head == "MOLSITE":
                equivalents = []        # molecular area equivalent of molsite
                molsites = []           # molecular adsorption sites
                for i in tail:
                    try:
                        equivalents.append(float(i))
                    except ValueError:
                        molsites.append(str(i))
                ''' It may be the case equivalent = 1 is neglected.
                Then, it is check that the len(equivalents) is the same than
                the number of molsites or 1 is added to equivalent.'''
                while len(molsites) > len(equivalents):
                    equivalents.append(float(1))
                systems[name][nadsorbates]["molsite"] = molsites   # molecular adsorption sites
                systems[name][nadsorbates]["nmolsite"] = equivalents   # molecular area equivalent of molsite
            if head == "DEGENERATION":
                systems[name][nadsorbates]["degeneration"] = int(tail[0])
            if head == "IPRESSURE":
                if len(tail) == 1:
                    rconditions[name]["pressure0"] = float(tail[0])     # constant
                else:
                    rconditions[name]["pressure0"] = [float(i) for i in tail]   # ramp
            if head == "ICOVERAGE":
                if len(tail) == 1:
                    rconditions[name]["coverage0"] = float(tail[0])     # constant
                else:
                    rconditions[name]["coverage0"] = [float(i) for i in tail]   # ramp
    ''' Checks in systems to assign the kind of system and check:
        - the kind, 
        - the path for the frequencies,
        - the sites,
        - the degeneration,
        - the initial pressures/coverages,
        - the 3N-6(5) frequencies,
        - The number of adsorbates in an adsorbed system. '''
    for name in systems.keys():
        if "area" in systems[name][list(systems[name].keys())[0]].keys():
            systems[name]["kind"] = "surface"
        elif "mass" in systems[name][list(systems[name].keys())[0]].keys():
            systems[name]["kind"] = "molecule"
        else:
            systems[name]["kind"] = "adsorbed"
        if systems[name]["kind"] == "surface":
            for nadsorbates in systems[name].keys():
                if len(systems[name][nadsorbates]["area"]) != len(systems[name][nadsorbates]["site"]):
                    print("   ERROR: the number of sites and areas in {}{} is\ "
                          "not the same.".format(name,nadsorbates))
                    exit(0)
        if systems[name]["kind"] == "molecule" and "pressure0" not in systems[name].keys():
            systems[name]["pressure0"] = 0.
        if systems[name]["kind"] == "adsorbed" and "coverage0" not in systems[nadsorbates].keys():
            systems[name]["coverage0"] = 0.
        for nadsorbates in systems[name].keys():
            if "freqpath" not in systems[name][nadsorbates].keys():
                systems[name][nadsorbates]["freqpath"] = systems[name][nadsorbates]["syspath"]
            if "degeneration" not in systems[name][nadsorbates].keys():
                systems[name][nadsorbates]["degeneration"] = 1
            if systems[name][nadsorbates]["kind"] != "surface" and\
                    len(systems[name][nadsorbates]["freq3d"]) == 0:
                print("   ERROR: frequencies for {}{} are not provided".format(name,nadsorbates))
                exit(0)
            if systems[name][nadsorbates]["kind"] == "molecule" and\
                    len(systems[name][nadsorbates]["freq2d"]) == 0:
                print("   ERROR: 2D frequencies for {}{} are not provided".format(name,nadsorbates))
                exit(0)
            ''' It may be the case that the number of adsorbates in an adsorbed 
                system has not been neglected or provided (nadsorbates = 0). 
                To maintain the physical meaning, nadsorbates is redefined as 1.'''
            if systems[name][nadsorbates]["kind"] == "adsorbed" and nadsorbates == str(0):
                if len(systems[name].keys()) == 1:
                    systems[name][str(1)] = systems[name][str(0)]
                    del systems[name][str(0)]
                else:
                    print("   ERROR: the number of adsorbates in {} is conflicting or not provided".format(name))
                    exit(0)
    ''' Once the systems has defined as molecule or surface, 
    looking for the area occupied for each molecule (marea)'''
    sites = {}
    for name in systems.keys():
        if systems[name]["kind"] == "surface":
            for nadsorbates in systems[name].keys():
                sites[str(systems[name][nadsorbates]["site"])] = float(systems[name][nadsorbates]["area"])
    for name in systems.keys():
        if systems[name]["kind"] == "molecule":
            for nadsorbates in systems[name].keys():
                if systems[name][nadsorbates]["molsite"] in sites:
                    systems[name][nadsorbates]["marea"] = (sites[str(systems[name][nadsorbates]["molsite"])] /
                                                           int(systems[name][nadsorbates]["nmolsite"]))
                temp = sp.symbol("temperature")
                pressure = 101325       # Pa == kg⋅m^−1⋅s^−2
                systems[name][nadsorbates]["volume"] = constants["kb"]*temp/pressure        # Assuming ideal behaviour of gases


    for key in processes.keys():
        freq = []
        if name in [i for i in processes[key]["ts"]]:
            for i in systems[name][nadsorbates]["freq3d"][:-2]:
                if i < -100:
                    print("   ERROR: {}{} has more than one significant imaginary frequency".format(name,nadsorbates))
                    exit(0)
                else:
                    freq.append(systems[name][nadsorbates]["freq3d"][i])
            freq.append(systems[name][nadsorbates]["freq3d"][-1])
            systems[name][nadsorbates]["freq3d"] = freq
        else:
            for i in systems[name][nadsorbates]["freq3d"]:
                if i < -100:
                    print("   ERROR: {}{} has a significant imaginary frequency".format(name,nadsorbates))
                    exit(0)
                else:
                    freq.append(systems[name][nadsorbates]["freq3d"][i])
            systems[name][nadsorbates]["freq3d"] = freq

    return rconditions, processes, systems


rconditions, processes, systems = mkread(str(sys.argv[1]))

systems = PartitionFunctions(dict(rconditions), dict(systems), dict(constants))
