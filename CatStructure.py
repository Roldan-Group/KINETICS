




import os, sys, re, subprocess
from ase import Atoms
from ase.io import read, write



eleNames = ["Au"]                           # Elements in the Cluster
inputFiles = "OUTCAR"



# Check for a common folder to place the structures
#if not os.path.exists("../../Data2Structures"):
#    os.makedirs("../../Data2Structures")

# Finds the line number in DataLine2XYZ file
#OutputFile = str("../../Data2Structures/Data2XYZ.csv")
OutputFile = str("../Data2XYZ.csv")
try:
    ofile = open(OutputFile)
    nlines = len(ofile.readlines())
except:
    nlines = 0
    pass

# Get the information of the system
try:
    system = read(inputFiles)
except IOError:
    raise Exception("   The "+InputFile+" file for structure in line "+nlines+" is not a valid file!")


# identify the atoms in the cluster
XYZ = [system[i].position for i in range(len(system)) if system[i].symbol in eleNames]

# get the XYZ of each atom in the cluster
#XYZ = [ i.position for i in cluster]

# print XYZ in Data2XYZ
ofile = open(OutputFile,'a+')
ofile.write ("%d, " %(len(XYZ)))
for n in range(len(XYZ)-1):
    ofile.write ("%.5f, %.5f, %.5f,   " %(XYZ[n][0],XYZ[n][1],XYZ[n][2]))
ofile.write ("%.5f, %.5f, %.5f" %(XYZ[-1][0],XYZ[-1][1],XYZ[-1][2]))
ofile.write ("\n")
ofile.close



