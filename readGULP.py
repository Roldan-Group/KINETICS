



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

    for i, line in enumerate(lines):


        if line.find("Final Cartesian lattice vectors") != -1:
            s = i + 1
            print(line)
            a = (float(lines[s].split()[0]), float(lines[s].split()[1]), float(lines[s].split()[2]))
            s = i + 1
            b = (float(lines[s].split()[0]), float(lines[s].split()[1]), float(lines[s].split()[2]))
            s = i + 1
            c = (float(lines[s].split()[0]), float(lines[s].split()[1]), float(lines[s].split()[2]))

            atoms.set_cell



        if line.find('Dimensionality') != -1:
            if int(line[19]) == 3:
                atoms.set_pbc((True, True, True))
            else:
                atoms.set_pbc((False, False, False))

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
            break

    atoms.set_scaled_positions(positions)

    try:
        atoms.set_chemical_symbols(symbols)
    except KeyError:
        for i in range(len(symbols)):
            for types in atom_type_list:
                if symbols[i] in types[1:]:
                    symbols[i] = types[0]
                    break
        atoms.set_chemical_symbols(symbols)

    return atoms
################################################################




filename = sys.argv[1]

atoms = read_gulp(filename)

print(atoms)
view(atoms)

