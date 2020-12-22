



import sys
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read 
from ase.visualize import view


#####################################################

def read_gulp(filename, atom_type_list=[]):
    """Import gulp output type files.
    atom_type_list is a list where each element is the desctiption of atom types
    example  [[O, O1, O2, O3], [H, H1, H2]]
    O is chemical symbols to be associated with O1, O2 and O3 types
    Reads unitcell, atom positions.
    """

    from ase import Atoms
    import numpy as np

    atoms = Atoms()
    fd = open(filename, 'r')
    lines = fd.readlines()
    fd.close()
    positions = []
    symbols = []
    a = (1, 0, 0)
    b = (0, 1, 0)
    c = (0, 0, 1)
    pbc = False

    for i, line in enumerate(lines):
# looking for the cell
        if line.find('Final Cartesian lattice vectors') != -1:
            s = i + 2
            a = (float(lines[s].split()[0]), float(lines[s].split()[1]), float(lines[s].split()[2]))
            s = s + 1
            b = (float(lines[s].split()[0]), float(lines[s].split()[1]), float(lines[s].split()[2]))
            s = s + 1
            c = (float(lines[s].split()[0]), float(lines[s].split()[1]), float(lines[s].split()[2]))
# looking for periodic conditions
        if line.find('Dimensionality') != -1:
            if int(line[19]) == 3:
                pbc = True
#                print(pbc)
# looking for atomic positions
        if line.find('Final fractional coordinates of atoms') != -1:
            s = i + 5
            while(True):
                s = s + 1
                if lines[s].find("------------") != -1:
                    break
                if lines[s].find(" s ") != -1:
                    continue
                x = float(lines[s].split()[3])
                y = float(lines[s].split()[4])
                z = float(lines[s].split()[5])
                positions.append([x, y, z])
                symbols.append(lines[s].split()[1])
            positions = np.array(positions)
# looking for chemical symbols
    try:
        atoms = Atoms(symbols)
    except KeyError:
        for i in range(len(symbols)):
            for types in atom_type_list:
                if symbols[i] in types[1:]:
                    symbols[i] = types[0]
        atoms = Atoms(symbols)

# setting Atoms
    atoms.set_cell([a, b, c])
    atoms.set_pbc((pbc, pbc, pbc))
    atoms.set_scaled_positions(positions)
#    atoms.set_positions(positions)

    return atoms
################################################################




filename = sys.argv[1]

atoms = read_gulp(filename)

print(atoms)
view(atoms)

