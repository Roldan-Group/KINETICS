


import sys
from ase.io import read, write
from ase.visualize import view

file_name = sys.argv[1]
system = read(file_name)

write(str(file_name) + ".cif", system, format="cif")
out = read(str(file_name) + ".cif")
#view(out)









