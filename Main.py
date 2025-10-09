"""

make first Input2mk.py


    This script reads and executes the packages for KINETICS

        by Alberto Roldan
"""
import sys, re
import time
import sympy as sp
from Thermodynamics import PartitionFunctions , Energy
from Kinetics import RConstants, REquations



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
                for i in tail:
                    try:
                        nsites.append(float(i))
                    except ValueError:
                        sites.append(str(i))
                if len(nsites) < len(sites):
                    systems[name][nadsorbates]["nsite"] = [1 for i in range(len(sites))]
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
                        inertia.append(float(i))
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
    for name in systems.keys():
        if "area" in systems[name][list(systems[name].keys())[0]].keys():
            systems[name]["kind"] = "surface"
        elif "imass" in systems[name][list(systems[name].keys())[0]].keys():
            systems[name]["kind"] = "molecule"
        else:
            systems[name]["kind"] = "adsorbate"
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
                if key not in ["kind", "pressure0"]:       # only for nadsorbates
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
    '''A molecule will adsorb on one a site with a particular area (marea). If the molecules has more than site to 
    adsorbed, differente systems needs to be described'''
    for name in systems.keys():
        if systems[name]["kind"] == "molecule":
            for key in systems[name].keys():
                if key not in ["kind", "pressure0"]:       # only for nadsorbates
                    if systems[name][key]["molsite"] in sites:
                        systems[name][key]["marea"] = (sites[str(systems[name][key]["molsite"])] /
                                         int(systems[name][key]["nmolsite"]))
                    temp = sp.symbols("temperature")
                    pressure = 101325       # Pa == kg⋅m^−1⋅s^−2
                    systems[name][key]["volume"] = constants["kb"]*temp/pressure        # Assuming ideal behaviour of gases
    for nproc in processes.keys():
        for name in systems.keys():
            freq = []
            if name in [i for i in processes[nproc]["ts"]]:
                for key in systems[name].keys():
                    if key not in ["kind", "pressure0", "coverage0"]:  # only for nadsorbates
                        for i in systems[name][key]["freq3d"][:-2]:
                            if i < -100:
                                print("   ALERT: {}{} has more than one significant imaginary frequency".format(
                                    name, key))
                            else:
                                freq.append(i)
                        freq.append(systems[name][key]["freq3d"][-1])
                        systems[name][key]["freq3d"] = freq
            else:
                for key in systems[name].keys():
                    if key not in restricted_arg and "freq3d" in systems[name][key].keys():  # only for nadsorbates
                        for i in systems[name][key]["freq3d"]:
                            if i < -100:
                                print("   ALERT: {}{} has a significant imaginary frequency".format(name, key))
                            else:
                                freq.append(i)
                        systems[name][key]["freq3d"] = freq
    return rconditions, processes, systems

''' list of restricted argunments in systems[name] containing the interpolated functions'''
restricted_arg = ["kind", "pressure0", "coverage0", "sites", 'q3d', 'q2d', 'energy3d', 'energy2d', 'ifreq']

start0 = time.time()
rconditions, processes, systems = mkread(str(sys.argv[1]), list(restricted_arg))
print("... Reading ...", round(time.time()-start0, 3), " seconds")
start = time.time()
systems = PartitionFunctions(dict(rconditions), dict(systems), dict(constants), list(restricted_arg)).systems
print("... Generating Partition Functions ...", round(time.time()-start, 3), " seconds")
start = time.time()
systems = Energy(dict(rconditions), dict(systems), dict(constants), list(restricted_arg)).systems

''' INTEGRATION for the calculation of ENTHALPY has been shitched OFF '''

print("... Generating Thermodynamics ...", round(time.time()-start, 3), " seconds")
start = time.time()
processes = RConstants(dict(rconditions), dict(systems), dict(constants), dict(processes), list(restricted_arg)).processes
print("... Generating Reaction Constants ...", round(time.time()-start, 3), " seconds")
start = time.time()
constemperature_equations = REquations(dict(processes), dict(systems)).constemperature
tpd_equations = REquations(dict(processes), dict(systems)).tpd
print("... Generating Rate Equations ...", round(time.time()-start, 3), " seconds")
start = time.time()
print("... Microkinetics Completed ...", round((time.time()-start0)/60, 3), " minutes")
