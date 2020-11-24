

from ase.build import bulk, surface, fcc111, fcc100, fcc110
from ase import Atoms
from ase.constraints import FixAtoms
import ase.io.vasp
from ase.visualize import view

miller = (0, 1, 4)
layers = 32       # when step, consider multi by miller
size = (1, 1, layers)
slab_expand = (3, 1, 1)
bulk_lattice = 4.171822846911442 # Au 3.56447997606178 # Cu 

name = str(layers)
vacuum = 15/bulk_lattice 	# 15A vacuum

bulk = bulk('Au', 'fcc', a=1)

surf = surface(bulk, miller, layers) * slab_expand
#surf = fcc111("Cu", a=1, size=size)

surf.center(vacuum=vacuum, axis=2)
surf.set_positions([[i.position[0], i.position[1], i.position[2]-vacuum] for i in surf])
z_max = max([i.position[2] for i in surf])
surf.set_constraint(FixAtoms(indices=[i.index for i in surf if i.position[2] < z_max/2]))
ase.io.vasp.write_vasp(name, surf, direct=False, vasp5=True, sort=True, ignore_constraints=False)

f = open(name)
new = f.readlines()
f.close()
xyz = [i.split() for i in new[9:]]
xyz = sorted(xyz, key=lambda x: x[2])

ifile = open(name, 'w+')
ifile.write("{} {} layers = {} vacuum = 15\n" .format(surf.symbols, miller, layers))
ifile.write(" {:.11f}\n" .format(bulk_lattice))
for i in range(2, 9):
    ifile.write(new[i])
for i in xyz:
    ifile.write(" {:.11f} {:.11f} {:.11f}  {:s} {:s} {:s}\n" .format(float(i[0]), float(i[1]), float(i[2]),
                                                                     i[3], i[4], i[5]))
ifile.close()

#view(ase.io.vasp.read_vasp(name))
