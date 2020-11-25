'''
    Alberto Roldan

    STRUCTURE:
        - Library of interatomic distances (metals)
            and bulk lattices as distance between sites (oxides)

'''

import numpy as np


def opt_atom_distance(support, site, element):          # Defines the optimised ONE atom-support distance in Ansgtroms
    optdistances = {
                    ('MgO', 'O',    'Co',   ),    # 3rd row
                    ('MgO', 'Mg',   'Co',   ),
                    ('MgO', 'O',    'Ni',   ),
                    ('MgO', 'Mg',   'Ni',   ),            
                    ('MgO', 'O',    'Cu',   2.0184),
                    ('MgO', 'Mg',   'Cu',   2.7037),    # 3rd row
                    ('MgO', 'O',    'Ru',   2.0032),    # 4th row
                    ('MgO', 'Mg',   'Ru',   2.6564),
                    ('MgO', 'O',    'Rh',   2.0836),
                    ('MgO', 'Mg',   'Rh',   2.6412),
                    ('MgO', 'O',    'Pd',   2.1020),
                    ('MgO', 'Mg',   'Pd',   2.5317),
                    ('MgO', 'O',    'Ag',   2.4385),
                    ('MgO', 'Mg',   'Ag',   2.8251),    # 4th row
                    ('MgO', 'O',    'Ir',   2.0184),    # 5th row
                    ('MgO', 'Mg',   'Ir',   2.8757),
                    ('MgO', 'O',    'Pt',   1.9885),
                    ('MgO', 'Mg',   'Pt',   2.5948),
                    ('MgO', 'O',    'Au',   2.3179),
                    ('MgO', 'Mg',   'Au',   2.6857),
                    ('C', 'C',   'Au',   2.6857),    # 5th row
                    }
    for i, sys in enumerate(optdistances):
        if sys[0] == support:
            if sys[1] == site:
                if sys[2] == element:
                    optimum_distance = sys[3]
    return optimum_distance


def sites(support):         # defines the common adsorption sites on a support
    sites = {
        "MgO": ["Mg", "O"]
    }
    return sites[str(support)]


def isolated_atoms(element):                     # energies of isolated atoms in vaccuo (RPBE)
    elements = {
                 'Fe': -3.4949,                 # 3rd row
                 'Co': -2.1442,
                 'Ni': -0.8036,
                 'Cu': -0.0099,                 # 3rd row
                 'Ru': -2.5726,                 # 4th row
                 'Rh': -1.5157,
                 'Pd': -1.5043,
                 'Ag': -0.1921,                 # 4th row
                 'Ir': -1.2597,                 # 5th row
                 'Pt': -1.6193,
                 'Au': -0.1921,                 # 5th row
               }

    return elements[element]



def ecoh_bulk(element):                         # cohesion energies at bulk coordination (RPBE)
    ecoh = {
            'Co': [-7.1245, 12],    # hcp       # 3rd row
            'Ni': [-5.6273, 12],
            'Cu': [-4.0698, 12],                # 3rd row
            'Ru': [-9.4469, 12],    # hcp       # 4rd row
            'Rh': [-7.5247, 12],
            'Pd': [-5.5162, 12],
            'Ag': [-3.0364, 12],                # 4rd row
            'Ir': [-9.4589, 12],                # 5rd row
            'Pt': [-6.5738, 12],
            'Au': [-3.5650, 12],                # 5th row
           }
    return ecoh[element]


def ecoh_trend(element):                         # cohesion energies trend parameter (a in equations)
    coh_parameter = {
            'Au': 3.18740093,                # 5th row
           }
    return coh_parameter[element]


def areas(element, coordination):                # Atomic areas previously calculated from surfaces as a function of the atom's coordination [0-->12]
    area = {                                    #   interpolated using a  LORENTZIAN function(x,a,b,c,d):   a + b / (4 * ( x - c )**2 - d**2)
#           element || areas vs coordination [0-->12] in Angstroms || a,b,c,d and R^2 of interpolation
# 3rd row           0       1       2       3       4       5       6       7       8       9    10      11    12      a        b       c       d       R^2            
            'Co': [-7.1245, 12], # hcp
            'Ni': [-5.6273, 12],
            'Cu': [-4.0698, 12],
# 4rd row           0       1       2       3       4       5       6       7       8       9    10      11    12      a        b       c       d       R^2 
            'Ru': [-9.4469, 12], # hcp
            'Rh': [-7.5247, 12],
            'Pd': [-5.5162, 12],
            'Ag': [-3.0364, 12],
# 5rd row           0       1       2       3       4       5       6       7       8       9    10      11    12      a        b       c       d       R^2            
            'Ir': [-9.4589, 12],                
            'Pt': [-6.5738, 12],
            'Au': [15.033, -3126.937, 19.167, 0.491, 0.9665],
           }
    a, b, c, d, r2 = area[element]
    return a + b / (4 * (int(coordination) - c)**2 - d**2)


def surf_energies(element, coordination):        # Atomic Surface energies previously calculated from surfaces as a function of the atom's coordination [0-->12]
    surf_energy = {                                 #   interpolated using a straigh line a*x + b
#           element || energies vs coordination [0-->12] in Angstroms || a,b and R^2 of interpolation
# 3rd row           0       1       2       3       4       5       6       7       8       9    10      11    12      a        b       R^2 
            'Co': [-7.1245, 12], # hcp
            'Ni': [-5.6273, 12],
            'Cu': [-4.0698, 12],
# 4rd row           0       1       2       3       4       5       6       7       8       9    10      11    12      a        b        R^2 
            'Ru': [-9.4469, 12], # hcp
            'Rh': [-7.5247, 12],
            'Pd': [-5.5162, 12],
            'Ag': [-3.0364, 12],
# 5rd row           0       1       2       3       4    5     6      7      8      9     10      11    12      a      b      R^2
            'Ir': [-9.4589, 12],                
            'Pt': [-6.5738, 12],
            'Au': [2.001, 1.837, 1.673, 1.509, 1.345, 1.181, 1.067, 0.802, 0.656, 0.517, 0.361, 0.269, 0.000, -0.164, 2.001, 0.9108],
            }
    return surf_energy[element][coordination]


def morse_potential_depth(support, element, ics):                # considering as reference O sites and d_z
    mp_depth = {
            ('MgO', 'Au', -6.93621926, 2.96247442, 16.27817071),    # 5rd row
    }
    for i, sys in enumerate(mp_depth):
        if sys[0] is support:
            if sys[1] is element:
                a = sys[2]
                b = sys[3]
                c = sys[4]
    return a + b * np.log(ics + c)            # LOGARITHMIC


def morse_potential(support, site, element, distance, d_eq):
    mp = {
            ('MgO', 'O',    'Au',   1.853,  2.324),    # 5rd row
            ('MgO', 'Mg',   'Au',   1.139,  2.876),
        }
    for i, sys in enumerate(mp):
        if sys[0] is support:
            if sys[1] is site:
                if sys[2] is element:
                    a = sys[3]
                    r_eq = sys[4]
    return d_eq * (np.exp(-2 * a * (distance - r_eq)) - 2 * np.exp(-a * (distance - r_eq)))


