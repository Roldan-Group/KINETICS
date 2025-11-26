
import sympy as sp


t, temp, cov = sp.symbols('time temperature coverage', positive=True, real=True)
h, kb, c, hc, R, Av, Fa, qelectron, JtoeV = sp.symbols(
	"h kb c hc R Av Fa qelectron JtoeV", real=True, positive=True)
qelectron = sp.symbols("qelectron", real=True)

constants = {
	h: sp.Float(6.62607588705515e-34),  # kg m^2 s^-1 == J s
	kb: sp.Float(1.380658045430573e-23),  # J K^-1
	c: sp.Float(299792458),  # m s^-1
	hc: sp.Float(1.98644586e-23),  # J cm (to balance the frequencies units of cm^-1)
	R: sp.Float(8.31446261815324),  # J⋅K−1⋅mol−1 ----------needed?
	Av: sp.Float(6.022139922973909e23),  # mols^-1
	Fa: sp.Float(96485.3321233100184),  # C⋅mol^−1
	qelectron: sp.Float(-1.60217648740e-19),  # C
	JtoeV: sp.Float(6.24150974e18)  # 1 J = 6.24..e18 eV
	}


def sym_equation(processes, name):  # process is processes[process]; i indicates the process number
    '''The rates expressions re-written using the symbolic rate constant instead of the long formula.'''
    equation  = 0
    for process in processes.keys():
        if name in processes[process]['products']:
            for r in range(len(processes[process]['products'])):
                if name == process['products'][r]:
                    equation += processes[process]['pstoichio'][r]  # rfactor
            equation *= sp.symbols(f'k_{process}')
            for r in range(len(processes[process]['reactants'])):
                equation *= sp.symbols(f"{processes[process]['reactants'][r]}") ** processes[process]['rstoichio'][r]
        elif name in processes[process]['reactants']:
            for r in range(len(processes[process]['reactants'])):
                if name == process['reactants'][r]:
                    equation -= processes[process]['rstoichio'][r]  # rfactor
            equation *= sp.symbols(f'k_{process}')
            for r in range(len(processes[process]['reactants'])):
                equation *= sp.symbols(f"{processes[process]['reactants'][r]}") ** porcesses[process]['rstoichio'][r]
    return equation
