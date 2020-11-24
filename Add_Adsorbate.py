'''
		usage:	~.py adsorbent_file adsorbate_file
'''
import sys
import numpy as np
from ase.io import read
from ase.build import add_adsorbate
from ase import Atoms
from ase.constraints import FixAtoms
import ase.io.vasp
from ase.visualize import view


adsorbate_position = (0, 0)
adsorbate_offset = (.25, .85)
adsorbate_mol_index = 6
adsorbate_distance = 2.0

#-------------------------------------------------------------- Adsorbent
if len(sys.argv) > 1:
	adsorbent_file = sys.argv[1]
else:
	answer = input("Which file defines the adsorbent?\n")
	try:
		adsorbent_file = answer
	except IOError:
		raise Exception("   " + answer + " is not a valid!")
print(" The adsorbent is", read(adsorbent_file).symbols)
#-------------------------------------------------------------- Adsorbate
if len(sys.argv) > 2:
	adsorbate_file = sys.argv[2]
else:
	answer = input("Which file defines the adsorbate?\n")
	try:
		adsorbate_file = answer
	except IOError:
		pass
print(" The adsorbate is", read(adsorbate_file).symbols)


def Rewrite(file_name, adsorbate_name, adsorbate_position, adsorbate_offset, adsorbate_distance, adsorbate_mol_index):
	f = open(file_name)
	structure = f.readlines()
	f.close()
	ifile = open(file_name, 'w+')
	angs = "$\\gamma$"# "$\AA$"
	ifile.write("{} at position {}, offset {} and at {}{:s} from adsorbent through atom {}\n" .format(adsorbate_name,
																									 adsorbate_position,
																									 adsorbate_offset,
																									 adsorbate_distance,
																									 angs,
																								   adsorbate_mol_index))
	for i in range(1, 9):
		ifile.write(structure[i])
	elements = [int(i) for i in structure[6].split()]
	i = 9
	for n_atoms in elements:
		xyz = []
		for line in range(i, i+n_atoms):
			xyz.append([n for n in structure[line].split()])
		i += n_atoms
		xyz = sorted(xyz, key=lambda x: x[2], reverse=True)
		for a in xyz:
			ifile.write(" {:>15.11f} {:>15.11f} {:>15.11f}  {:s} {:s} {:s}\n" .format(float(a[0]), float(a[1]),
																				   float(a[2]), a[3], a[4], a[5]))
	ifile.close()


#------------------------------------------------------------- rotate adsorbate
for x in range(0, 91, 45):
	adsorbate = read(adsorbate_file)
	adsorbate.rotate(-x, 'x')
	for z in range(0, 91, 45):
		adsorbate.rotate(z, 'z')
		adsorbent = read(adsorbent_file)
		add_adsorbate(adsorbent, adsorbate, adsorbate_distance, position=adsorbate_position,
									  							offset=adsorbate_offset,
					  											mol_index=adsorbate_mol_index)  # use offset instead of position
		ase.io.vasp.write_vasp(str(adsorbent.symbols) + "_xRot" + str(x) + "_zRot" + str(z) + ".vasp",
									adsorbent, direct=False, vasp5=True, sort=True, ignore_constraints=False)
		Rewrite(str(adsorbent.symbols) + "_xRot" + str(x) + "_zRot" + str(z) + ".vasp", str(adsorbate.symbols),
								adsorbate_position, adsorbate_offset, adsorbate_distance, adsorbate_mol_index)

#		view(adsorbent)




