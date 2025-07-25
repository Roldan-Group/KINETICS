"""
    This script builds on the perl version by A.Roldan.

"""

import os, pathlib
import time
import sympy as sp
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from withMatlab.Input2mk import interpolate


def printdata(rconditions, name, nadsorbates, properties, constants, datalabel, dataname):
    folder = './THERMODYNAMICS/DATA/'+ name + "/" + nadsorbates
    outputfile = folder + "/" + str(dataname) + ".dat"
    if not pathlib.Path(folder).exists():
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)
        os.chmod(folder, 0o755)
    if not pathlib.Path(outputfile).exists():
        output = open(outputfile, "w")
        output.write("# Temperature[K]")
        units = ''
        for i in datalabel:
            if i.startswith('e') or i.startswith('z'):
                units = '[eV]'     # energy and ZPE
            elif i.startswith('s') or i.startswith('c'):
                units = '[eV K^-1]' # entropy or specific heat as "defined here"
            output.write(" {val:>{wid}s}".format(wid=len(i+units)+3, val=i+units))
        output.write("\n")
        temp = sp.symbols("temperature")
        if isinstance(rconditions["temperature"], float):
            output.write(" {val:>{wid}.1f}".format(wid=len("Temperature[K]"), val=rconditions["temperature"]))
            for i in datalabel:
                value = sp.lambdify(temp, properties[str(i)])(float(rconditions["temperature"]))
                output.write(" {val:>{wid}.3{c}}".format(wid=len(i)+3, val=value, c='e' if value > 1e3 else 'f'))
        else:
            ramp = [int(i) for i in rconditions["temperature"]]
            for t in range(ramp[0], ramp[1], ramp[2]):
                output.write(" {val:>{wid}.1f}".format(wid=len("Temperature[K]"), val=t))
                for i in datalabel:
                    value = sp.lambdify(temp, properties[str(i)])(t)
                    output.write(" {val:>{wid}.3{c}}".format(wid=len(i)+3, val=value, c='e' if value > 1e3 else 'f'))
                output.write("\n")
        output.close()

def interpolate(rconditions, systems, name, restricted_arg, ykey):
    startint = time.time()
    folder = './THERMODYNAMICS/DATA/'+ name + "/"
    ''' Reaction conditions are set as symbols using SYMPY '''
    temp, cov = sp.symbols("temperature coverage")
    x0 = [] # list of nadsorbates
    x = []  # independent coordinate, e.g. coverage
    y = []  # dependent coordinate, e.g. energy
    for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
        if nadsorbates not in restricted_arg:       # only for nadsorbates
            x0.append(nadsorbates)
            x.append(float(nadsorbates))
            y.append(systems[name][nadsorbates][ykey])
    function = sp.interpolate(list(zip(x, y)), cov)    # generates the interpolation function with the x and y.
    ''' generate a plot for the interpolated data saving it as name.key.svg '''
    fig, ax1 = plt.subplots(figsize=(8, 6), clear=True)
    ''' generate a loop for the temperatures and the key'''
    #flambdified = sp.lambdify((temp, cov), function, modules='numpy')
    xdata = np.linspace(0, 1, 301)
    if isinstance(rconditions["temperature"], float):
        t = float(rconditions["temperature"])
        y = [function.subs({temp: t, cov: i}).evalf() for i in x]
        ax1.plot(x, y, marker='o', color='k', fillstyle='none', linestyle='none',
                 label='original at T='+str(float(rconditions["temperature"]))+' K')
        ydata = [function.subs({temp: t, cov: i}).evalf() for i in xdata]
        ax1.plot(xdata, ydata, color='b', linestyle='-',
                 label='interpolated at T='+str(float(rconditions["temperature"]))+' K')
    else:
        ramp = [int(i) for i in rconditions["temperature"]]
        y = [function.subs({temp: ramp[0], cov: i}).evalf() for i in x]
        ax1.plot(x, y, marker='o', color='k', fillstyle='none', linestyle='none',
                 label='original at T='+str(ramp[0])+' K')
        for t in range(ramp[0], ramp[1]+ramp[2], ramp[2]):
            ydata = [function.subs({temp: t, cov: i}).evalf() for i in xdata]
            ax1.plot(xdata, ydata, color='b', linestyle='-', alpha=1-(t-ramp[0])/(ramp[1]-ramp[0]),
                 label='interpolated at T='+str(t)+' K')
    ax1.set_xlabel('coverage (ML)', fontsize=18)
    ax1.set_xlim([0, 1])
    ax1.tick_params(axis='both', rotation=0, labelsize=16)
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter('{:.2f}'.format))
    if ykey.startswith("i"):
        ax1.set_ylabel(ykey+" $(cm^{-1})$", fontsize=18)
    elif ykey.startswith("e"):
        ax1.set_ylabel(ykey+" $(eV)$", fontsize=18)
    else:
        ax1.set_ylabel(ykey, fontsize=18)
    legend = ax1.legend(loc="best", fontsize=14)
    fig.tight_layout()
    plt.ion()
    plt.show()
    plt.savefig(folder+ykey+"_coverage.svg", dpi=300, orientation='landscape', transparent=True)
    print("\t... Interpolating", name, ykey,"...", round(time.time()-startint, 3), " seconds")
    return function


class PartitionFunctions:
    def __init__(self, rconditions, systems, constants, restricted_arg):
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
        for name in systems.keys():     # species
            if systems[name]["kind"] == "molecule":
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["qtrans3d"] = self.qtrans3d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["qtrans2d"] = self.qtrans2d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["qrot"] = self.qrot(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["qelec"] = self.qelec(systems[name][nadsorbates])
                        systems[name][nadsorbates]["qvib3d"] = self.qvib3d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["qvib2d"] = self.qvib2d(systems[name][nadsorbates], constants)
                        datalabel3d = ["qrot", "qelec", "qtrans3d", "qvib3d"]
                        q3d = 1.0
                        for i in datalabel3d:
                            q3d *= systems[name][nadsorbates][str(i)]
                        systems[name][nadsorbates]["q3d"] = q3d
                        datalabel3d.append("q3d")
                        datalabel2d = ["qrot", "qelec", "qtrans2d", "qvib2d"]
                        q2d = 1.0
                        for i in datalabel2d:
                            q2d *= systems[name][nadsorbates][str(i)]
                        systems[name][nadsorbates]["q2d"] = q2d
                        datalabel2d.append("q2d")
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                              datalabel3d + datalabel2d, "PartitionFunctions")
            else:
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["qtrans3d"] = 1
                        systems[name][nadsorbates]["qrot"] = 1
                        systems[name][nadsorbates]["qelec"] = self.qelec(systems[name][nadsorbates])
                        systems[name][nadsorbates]["qvib3d"] = self.qvib3d(systems[name][nadsorbates], constants)
                        datalabel = ["qtrans3d", "qrot", "qelec", "qvib3d"]
                        q3d = 1.0
                        for i in datalabel:
                            q3d *= systems[name][nadsorbates][str(i)]
                        datalabel.append("q3d")
                        systems[name][nadsorbates]["q3d"] = q3d
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                                  datalabel, "PartitionFunctions")
        ''' Partition function interpolation between nadsorbates of the same name '''
        for name in systems.keys():     # species
            if systems[name]["kind"] == "molecule":
                systems[name]["q3d"] = interpolate(rconditions, systems, name, restricted_arg, "q3d")
                systems[name]["q2d"] = interpolate(rconditions, systems, name, restricted_arg, "q2d")
            else:
                systems[name]["q3d"] = interpolate(rconditions, systems, name, restricted_arg, "q3d")
        self.systems = systems

    @staticmethod
    def qtrans3d(properties, constants):
        ''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
         Adsorption Journal Of The International Adsorption Society
         (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
         page 88 '''
        temp = sp.symbols("temperature")
        return properties["volume"]*((2*sp.pi*properties["mass"]*constants["kb"]*temp)**(3/2))/(constants["h"]**3)

    @staticmethod
    def qtrans2d(properties, constants):
        ''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
         Adsorption Journal Of The International Adsorption Society
         (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
         page 88 '''
        temp = sp.symbols("temperature")
        return properties["marea"]*(2*sp.pi*properties["mass"]*constants["kb"]*temp)/(constants["h"]**2)

    @staticmethod
    def qrot(properties, constants):
        ''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
         Adsorption Journal Of The International Adsorption Society
         (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
         page 90 '''
        temp = sp.symbols("temperature")
        if properties["linear"] == "yes":
            qrot = (8*sp.pi**2*properties["inertia"]*constants["kb"]*temp)/(properties["symfactor"]*constants["h"]**2)
        else:
            prod_inertia = 1.0
            for i in properties["inertia"]:
                prod_inertia *= i
            qrot = ((sp.sqrt(sp.pi*prod_inertia)/properties["symfactor"]) *
                    (8*sp.pi**2*constants["kb"]*temp/(constants["h"]**2))**(3/2))
        return qrot

    @staticmethod
    def qelec(properties):
        ''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
         Adsorption Journal Of The International Adsorption Society
         (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
         page 92 '''
        ''' It is consider that the contribution of excited states is negligible
        at the working temperatures '''
        return properties["degeneration"]

    @staticmethod
    def qvib3d(properties, constants):
        ''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
         Adsorption Journal Of The International Adsorption Society
         (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
         page 89 '''
        ''' The equation used is respect the lowest occupied state, not the bottom of the potential energy curve.
        It also consider small frequencies, when h*v ~ kb*T (v=frequencies). 
        The Zero Point Energy should be added to the energy as this qvib is to calculate the entropy (S),
        specific head (Cp), and pre-exponential factor of Arrhenius (A). '''
        temp = sp.symbols("temperature")
        qvib = 1
        if "freq3d" in properties:
            for freq in properties["freq3d"]:
                if freq > 0.0:
                    #qvib *= ((sp.exp((-1/2*constants["hc"]*freq)/(constants["kb"]*temp))) /
                    #         (1-sp.exp((-constants["hc"]*freq)/(constants["kb"]*temp))))
                    qvib *= 1/(1-sp.exp((-constants["hc"]*freq)/(constants["kb"]*temp)))
        return qvib

    @staticmethod
    def qvib2d(properties, constants):
        ''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
         Adsorption Journal Of The International Adsorption Society
         (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
         page 89 '''
        ''' The equation used is respect the lowest occupied state, not the bottom of the potential energy curve.
        It also consider small frequencies, when h*v ~ kb*T (v=frequencies). 
        The Zero Point Energy should be added to the energy as this qvib is to calculate the entropy (S),
        specific head (Cp), and pre-exponential factor of Arrhenius (A). '''
        temp = sp.symbols("temperature")
        qvib = 1
        for freq in properties["freq2d"]:
            if freq > 0.0:
                qvib *= 1/(1-sp.exp((-constants["hc"]*freq)/(constants["kb"]*temp)))
        return qvib


class Energy:        # Gibbs free energy in eV
    def __init__(self, rconditions, systems, constants, restricted_arg):
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
        ''' The ZPE is needed first as it is required to calculate the vibrational 
        contribution to the entropy '''
        for name in systems.keys():     # species
            if systems[name]["kind"] == "molecule":
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        zpe3d, ifreq = self.zpe3d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["zpe3d"] = zpe3d     # in eV
                        systems[name][nadsorbates]["ifreq"] = ifreq     # in cm^-1
                        systems[name][nadsorbates]["zpe2d"] = self.zpe2d(systems[name][nadsorbates], constants) # in eV
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                              ["zpe3d", "zpe2d"], "ZeroPointEnergy")
            else:
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        zpe3d, ifreq = self.zpe3d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["zpe3d"] = zpe3d     # in eV
                        systems[name][nadsorbates]["ifreq"] = ifreq     # in cm-1
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                              ["zpe3d"], "ZeroPointEnergy")

        ''' Entropy is required to calculate the specific heat (Cp), which contributes 
        to the enthalpy at specific temperatures '''
        self.systems = Entropy(rconditions, systems, constants, restricted_arg).systems         # in eV
        self.systems = Enthalpy(rconditions, systems, constants, restricted_arg).systems        # in eV

        for name in systems.keys():     # species
            if systems[name]["kind"] == "molecule":
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["energy3d"] = (systems[name][nadsorbates]["enthalpy3d"] -
                                                                   temp * systems[name][nadsorbates]["sentropy3d"])
                        systems[name][nadsorbates]["energy2d"] = (systems[name][nadsorbates]["enthalpy2d"] -
                                                                   temp * systems[name][nadsorbates]["sentropy2d"])
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                              ["energy3d", "energy2d"], "GibbsFreeEnergy")
            else:
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["energy3d"] = (systems[name][nadsorbates]["enthalpy3d"] -
                                                                   temp * systems[name][nadsorbates]["sentropy3d"])
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                              ["energy3d"], "GibbsFreeEnergy")
        ''' Energy interpolation between nadsorbates of the same name '''
        for name in systems.keys():     # species
            if systems[name]["kind"] == "molecule":
                systems[name]["energy3d"] = interpolate(rconditions, systems, name, restricted_arg, "energy3d")
                systems[name]["energy2d"] = interpolate(rconditions, systems, name, restricted_arg, "energy2d")

            else:
                systems[name]["ifreq"] = interpolate(rconditions, systems, name, restricted_arg, "ifreq")
                systems[name]["energy3d"] = 1
                # Activate for production
                #systems[name]["energy3d"] = interpolate(rconditions, systems, name, restricted_arg, "energy3d")

        print("\n  >> integrals in Enthalpy deactivated <<\n",
              "   >> interpolate 3d deactivated <<\n")
        self.systems = systems

    @staticmethod
    def zpe3d(properties, constants):
        ''' Quantum ZPE corrected with the Wigner's harmonic oscillator approach  (DOI: 10.1063/1.216119).
            Harmonic approach is also applied to include quantum tunneling in the calculation of reaction constants.'''
        temp = sp.symbols("temperature")
        zpe = 0
        ifreq = 0
        if "freq3d" in properties:
            for freq in properties["freq3d"]:
                if freq > 0.0:
                    # Quantum-mechanical Zero-Point-Energy
                    zpe += (1/2*constants["hc"]*freq) + (constants["hc"]*freq /
                                                         (sp.exp(constants["hc"]*freq/(constants["kb"]*temp)) - 1))
                else:
                    ifreq = freq
        return zpe * constants["JtoeV"], ifreq

    @staticmethod
    def zpe2d(properties, constants):
        ''' Quantum ZPE corrected with the Wigner's harmonic oscillator approach  (DOI: 10.1063/1.216119).
            Harmonic approach is also applied to include quantum tunneling in the calculation of reaction constants.'''
        temp = sp.symbols("temperature")
        zpe = 0
        if "freq2d" in properties:
            for freq in properties["freq2d"]:
                if freq > 0.0:
                    # Quantum-mechanical Zero-Point-Energy
                    zpe += (1/2*constants["hc"]*freq) + (constants["hc"]*freq /
                                                         (sp.exp(constants["hc"]*freq/(constants["kb"]*temp)) - 1))
        return zpe * constants["JtoeV"]


class Entropy:
    ''' Chorkendorff, I. & Niemantsverdriet, J. W. "Concepts of Modern Catalysis and Kinetics."
    Adsorption Journal Of The International Adsorption Society
    (Wiley, Weinheim, FRG, 2003). doi:10.1002/3527602658.
    page 88 :: S = kb*ln(Q)+kb*T*diff(ln(Q), T)'''
    def __init__(self, rconditions, systems, constants, restricted_arg):
        ''' Reaction conditions are set as symbols using SYMPY '''
        for name in systems.keys():     # species
            if systems[name]["kind"] == "molecule":
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["strans3d"] = self.strans3d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["strans2d"] = self.strans2d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["srot"] = self.srot(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["selec"] = self.selec(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["svib3d"] = self.svib3d(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["svib2d"] = self.svib2d(systems[name][nadsorbates], constants)
                        datalabel3d = ["strans3d", "srot", "selec", "svib3d"]
                        entropy = 0.0
                        for i in datalabel3d:
                            entropy =+ systems[name][nadsorbates][str(i)]
                        systems[name][nadsorbates]["sentropy3d"] = entropy
                        datalabel3d.append("sentropy3d")
                        datalabel2d = ["strans2d", "srot", "selec", "svib2d"]
                        entropy = 0.0
                        for i in datalabel2d:
                            entropy =+ systems[name][nadsorbates][str(i)]
                        systems[name][nadsorbates]["sentropy2d"] = entropy
                        datalabel2d.append("sentropy2d")
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                                  datalabel3d + datalabel2d, "Entropy")
            else:
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["strans3d"] = 0
                        systems[name][nadsorbates]["srot"] = 0
                        systems[name][nadsorbates]["selec"] = self.selec(systems[name][nadsorbates], constants)
                        systems[name][nadsorbates]["svib3d"] = self.svib3d(systems[name][nadsorbates], constants)
                        datalabel3d = ["strans3d", "srot", "selec", "svib3d"]
                        entropy = 0.0
                        for i in datalabel3d:
                            entropy =+ systems[name][nadsorbates][str(i)]
                        systems[name][nadsorbates]["sentropy3d"] = entropy
                        datalabel3d.append("sentropy3d")
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                                  datalabel3d, "Entropy")
        self.systems = systems

    @staticmethod
    def strans3d(properties, constants):
        '''re-Formulation from explicit derivatives::
         https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
        temp = sp.symbols("temperature")
        return constants["kb"]*(sp.log(((2*sp.pi*properties["mass"]*constants["kb"]*temp)/constants["h"]**2)**(3/2) *
                                       properties["volume"]) + 3/2)*constants["JtoeV"]

    @staticmethod
    def strans2d(properties, constants):
        '''re-Formulation from explicit derivatives::
         https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
        temp = sp.symbols("temperature")
        return constants["kb"]*(sp.log(((2*sp.pi*properties["mass"]*constants["kb"]*temp)/constants["h"]**2) *
                                       properties["marea"]) + 2/2)*constants["JtoeV"]

    @staticmethod
    def srot(properties, constants):
        '''re-Formulation from explicit derivatives::
         https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
        temp = sp.symbols("temperature")
        if properties["linear"] == "yes":
            srot = constants["kb"]*(sp.log((1/properties["symfactor"])*(8*sp.pi**2*constants["kb"]*temp*
                                                    properties["inertia"]/(constants["h"]**2)))+1)*constants["JtoeV"]
        else:
            prod_inertia = 1.0
            for i in properties["inertia"]:
                prod_inertia *= i
            srot = constants["kb"]*(sp.log((sp.sqrt(sp.pi*prod_inertia)/properties["symfactor"])*
                                  (8*sp.pi**2*constants["kb"]*temp/(constants["h"]**2))**(3/2))+5/2)*constants["JtoeV"]
        return srot

    @staticmethod
    def selec(properties, constants):
        '''re-Formulation from explicit derivatives::
         https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
        return constants["kb"]*(sp.log(2*properties["degeneration"] + 1))*constants["JtoeV"]

    @staticmethod
    def svib3d(properties, constants):
        '''re-Formulation from explicit derivatives::
         https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
        temp = sp.symbols("temperature")
        qvib = 1
        if "freq3d" in properties:
            for freq in properties["freq3d"]:
                if freq > 0.0:
                    qvib *= ((sp.exp((-1/2*constants["hc"]*freq)/(constants["kb"]*temp))) /
                             (1-sp.exp((-constants["hc"]*freq)/(constants["kb"]*temp))))
        return (constants["kb"]*sp.log(qvib))*constants["JtoeV"] + properties["zpe3d"]/temp

    @staticmethod
    def svib2d(properties, constants):
        '''re-Formulation from explicit derivatives::
         https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html'''
        temp = sp.symbols("temperature")
        qvib = 1
        if "freq2d" in properties:
            for freq in properties["freq2d"]:
                if freq > 0.0:
                    qvib *= ((sp.exp((-1/2*constants["hc"]*freq)/(constants["kb"]*temp))) /
                             (1-sp.exp((-constants["hc"]*freq)/(constants["kb"]*temp))))
        return (constants["kb"]*sp.log(qvib))*constants["JtoeV"] + properties["zpe2d"]/temp


class Enthalpy:
    """ C.J. Cramer. Essentials of Computational Chemistry, Second Edition. Wiley, 2004.
    and Raymand Chang, "PHYSICAL CHEMISTRY for Chemical and Biological Science", ISBN: 1-891389-06-8,
        page 91 :: H=[dCp/dT](T1,T2)
        H =  E + ZPE + integral(Cp, 0 --> T) """
    def __init__(self, rconditions, systems, constants, restricted_arg):
        ''' Reaction conditions are set as symbols using SYMPY '''
        temp = sp.symbols("temperature")
        for name in systems.keys():     # species
            if systems[name]["kind"] == "molecule":
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["cp3d"] = self.cp3d(systems[name][nadsorbates])
                        systems[name][nadsorbates]["cp2d"] = self.cp2d(systems[name][nadsorbates])
                        enthalpy3d = systems[name][nadsorbates]["energy0"] + systems[name][nadsorbates]["zpe3d"]
                                # activate for production
                                # + sp.integrate(systems[name][nadsorbates]["cp3d"], temp))
                        systems[name][nadsorbates]["enthalpy3d"] = enthalpy3d
                        enthalpy2d = systems[name][nadsorbates]["energy0"] + systems[name][nadsorbates]["zpe2d"]
                                # activate for production
                                # + sp.integrate(systems[name][nadsorbates]["cp2d"], temp))
                        systems[name][nadsorbates]["enthalpy2d"] = enthalpy2d
                        datalabel = ["zpe3d", "cp3d", "enthalpy3d", "cp2d", "zpe2d", "enthalpy2d"]
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                              datalabel, "Enthalpy")
            else:
                for nadsorbates in systems[name].keys():    # number of species, i.e. "coverage"
                    if nadsorbates not in restricted_arg:       # only for nadsorbates
                        systems[name][nadsorbates]["cp3d"] = self.cp3d(systems[name][nadsorbates])
                        enthalpy3d = systems[name][nadsorbates]["energy0"] +  systems[name][nadsorbates]["zpe3d"]
                                    # activate for production
                                    # + sp.integrate(systems[name][nadsorbates]["cp3d"], temp))
                        systems[name][nadsorbates]["enthalpy3d"] = enthalpy3d
                        datalabel = ["zpe3d", "cp3d", "enthalpy3d"]
                        printdata(rconditions, name, nadsorbates, systems[name][nadsorbates], constants,
                              datalabel, "Enthalpy")
        self.systems = systems

    @staticmethod
    def cp3d(properties):
        '''Hans Kuhn, Horst-Dieter Försterling, David Hennessey Waldeck, "Principles of Physical Chemistry"
        ISBN: 9780470089644, page 551 :: Cp=T*[dS/dT](N,P)'''
        temp = sp.symbols("temperature")
        return temp*sp.diff(properties["sentropy3d"], temp)

    @staticmethod
    def cp2d(properties):
        ''' Hans Kuhn, Horst-Dieter Försterling, David Hennessey Waldeck, "Principles of Physical Chemistry"
        ISBN: 9780470089644, page 551 :: Cp=T*[dS/dT](N,P) '''
        temp = sp.symbols("temperature")
        return temp*sp.diff(properties["sentropy2d"], temp)





