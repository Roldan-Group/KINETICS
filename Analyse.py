

#####################################
#
#	Alberto 03/2020
#
#####################################

import os,sys,math, ase.io.vasp
from ase.io import read
import numpy as np
from ase import Atoms, neighborlist

inputFile = ["OUTCAR", "CONTCAR", "freq/OUTCAR", "ibrion/workfunction/Workfunction", "ibrion/charge/charge_SUM/ACF.dat", "ibrion/charge/charge_SUM/SpinDensity/ACF.dat"]
outputfile = "Analysed.dat"



def energy(inputFile,outputfile):
	output = read(inputFile,index=-1)
	E = output.get_total_energy()                                                    # total energy of the system
	a,b,c,alfa,beta,gamma = output.get_cell_lengths_and_angles()
	volume = output.get_volume()

	file = open(outputfile,'w+')
	file.write ("\n\nTotal Energy = %.2f\n" %(E))
	file.write ("Cell Lattices = %.3f %.3f %.3f\n" %(a,b,c))
	file.write ("Cell Angles = %.1f %.1f %.1f\n" %(alfa,beta,gamma))
	file.write ("Cell Volume = %.2f\n" %(volume))
	file.close()
	return

def structure(inputFile,outputfile,arguments):
	atoms = read(inputFile)
	allelements = atoms.get_chemical_symbols()
	symbols = []
	for i in atoms:
		if i.symbol not in list(symbols):
			symbols.append(i.symbol)

	file = open(outputfile,'a+')
	file.write ("Total number of atoms = {}\n". format(len(atoms)))
	file.close()

	if len(arguments) > 1:
		arguments.pop(0)
		elements = [ atoms[i] for i in range(len(atoms)) if atoms[i].symbol in arguments]
	else:
		arguments = []
		elements = []

	distances=[]
	for i in elements:                                                              # the distances
		tmp=allelements[0]
		dist=[];
		for j in atoms:
			if j.index != i.index:
				if j.symbol != tmp:
					distances.append(sorted(dist,key=lambda x : x[4])[0])
					dist=[]
					tmp=j.symbol
					dist.append((i.symbol,i.index+1,j.symbol,j.index+1,atoms.get_distance(i.index,j.index,mic=True,vector=False)))
				else:
					dist.append((i.symbol,i.index+1,j.symbol,j.index+1,atoms.get_distance(i.index,j.index,mic=True,vector=False)))
	file = open(outputfile,'a+')
	file.write("\n")
	for d in distances:
	    file.write ("The shortest distance between %s%s and %s%s in %.3f\n" %(d[0],d[1],d[2],d[3],d[4]))
	file.close()

	cutOff = neighborlist.natural_cutoffs(atoms,mult=1.2)							# average distances
	a,b,d = neighborlist.neighbor_list('ijd', atoms,cutOff)
	coord = np.bincount(a)

	average_distances=[]
	average_angles=[]
	for n in range(len(symbols)):
		imask=[(symbols[n]==atoms[i].symbol) for i in a]
		for m in range(len(symbols)):
			jmask=[(symbols[m]==atoms[j].symbol) for j in b]
			dist=[]
			angle=[]
			for i in range(len(d)):
				if imask[i] == True and jmask[i] == True:
					dist.append(d[i])
#					if coord[a[i]] > 1:
#						c=0
#						while a[i+c] and a[i] == a[i+c] and c < coord[a[i]]:
#							ang = atoms.get_angle(b[i],a[i],b[i+c],mic=True)
#							if ang > 1:
#								angle.append(ang)
#							c += 1
			if len(dist) > 0:
				average_distances.append((symbols[n],symbols[m],sum(dist)/len(dist)))
			if len(angle) > 0:
				average_angles.append((symbols[n],symbols[m],sum(angle)/len(angle)))

	file = open(outputfile,'a+')
	file.write("\n")
	for i,j,ave_d in average_distances:
		file.write ("The average distance of coordinated %s and %s is %.3f\n" %(i,j,ave_d))
	for i,j,ave_ang in average_angles:
		file.write ("The average angle of coordinated %s--%s--%s is %.1f\n" %(j,i,j,ave_ang))
	file.close()
	return

def Freq(inputFile,outputfile):
        lines = open(inputFile,"r").readlines()
        nlines = len(lines)
        nfreq=0; freq=[]
        for iline in range(len(lines)):
            line=lines[iline]
            linesplit = line.split()
            nwords = len(linesplit)
            if nwords > 7 and 'f' in linesplit[1] and 'THz' in linesplit[4] and 'cm-1' in linesplit[8]:
                if 'i' in linesplit[1]:
                    freq.append(-float(linesplit[7]))
                else:
                    freq.append(float(linesplit[7]))
                    nfreq=nfreq+1


        file = open(outputfile,'a+')
        file.write("\nReal frequencies are")
        for f in freq:
            file.write("   %.1f" %f)
        file.write("  cm-1\n")
        file.write("The classical ZPE contribution is %.2f eV\n" %(0.5*4.1357e-15*2.998e10*sum(freq)))
        file.close()

def WF(inputFile,outputfile):
	WorkF_line = open(inputFile,"r").readlines()[1]
	file = open(outputfile,'a+')
	file.write("\nThe energy at the vacuum is %.2f\n" %(float(WorkF_line.split()[0])))
	file.write("The Fermi energy is %.2f\n" %(float(WorkF_line.split()[1])))
	file.write("The WorkFunction is %.2f\n" %(float(WorkF_line.split()[2])))
	file.close()

def charges(inputFile, inputFile2, outputfile, arguments):
	atoms = read(inputFile)
	allelements = atoms.get_chemical_symbols()
	symbols=[]
	for i in atoms:
		if i.symbol not in list(symbols):
			symbols.append(i.symbol)

	if len(arguments) > 1:
		elements = [ atoms[i] for i in range(len(atoms)) if atoms[i].symbol in arguments]
	else:
		arguments = []
		elements = []

	valence = open(inputFile2,"r").readlines()[2:-4]
	charges = [ valence[i].split()[4] for i in range(len(valence))]

	e_transferred=0												# transferred charge
	for i in elements:
#		print (i.symbol,e_transferred)
		e_transferred += float(charges[i.index])
	file = open(outputfile,'a+')
	file.write("\nThe valence electrons trasferred to %s is %.1f\n" %(arguments,e_transferred))

	for n in symbols:											# average charge
		n_charges=[]
		for i in range(len(allelements)):
			if allelements[i] == n:
				n_charges.append(float(charges[i]))
		file.write("The average valence electrons of species %s is %.1f\n" %(n,sum(n_charges)/len(n_charges)))

	cutOff = neighborlist.natural_cutoffs(atoms,mult=1.2)                        
	a,b = neighborlist.neighbor_list('ij', atoms,cutOff)

	for e in elements:
		j_index = [ b[i] for i in range(len(a)) if e.index == a[i] ]
		for j in j_index:
			j_symbol = atoms[j].symbol
			file.write("The valence electrons of %s%s conected atom %s%s is %.1f\n" %(j_symbol,j+1,e.symbol,e.index + 1,float(charges[j])))
	file.close()



def spins(inputFile, inputFile2, outputfile, arguments):
	atoms = ase.io.read(inputFile)
	allelements = atoms.get_chemical_symbols()
	symbols=[]
	for i in atoms:
		if i.symbol not in list(symbols):
			symbols.append(i.symbol)

	if len(arguments) > 1:
		elements = [ atoms[i] for i in range(len(atoms)) if atoms[i].symbol in arguments]
	else:
		arguments = []
		elements = []

	valence = open(inputFile2,"r").readlines()[2:-4]
	spins = [ valence[i].split()[4] for i in range(len(valence))]

	file = open(outputfile,'a+')
	for n in symbols:                                                   # average spin
		n_spins=[]
		for i in range(len(allelements)):
			if allelements[i] == n:
				n_spins.append(float(spins[i]))
		file.write("The average spin density of species %s is %.1f\n" %(n,sum(n_spins)/len(n_spins)))

	cutOff = neighborlist.natural_cutoffs(atoms,mult=1.2)
	a,b = neighborlist.neighbor_list('ij', atoms,cutOff)

	for e in elements:
		j_index = [ b[i] for i in range(len(a)) if e.index == a[i] ]
		for j in j_index:
			j_symbol = atoms[j].symbol
			file.write("The spin density of %s%s conected atom %s%s is %.1f\n" %(j_symbol,j+1,e.symbol,e.index + 1,float(spins[j])))
	file.close()

def Bader_Volume(inputFile, inputFile2, outputfile):
        atoms = ase.io.read(inputFile)

        Bader = open(inputFile2,"r").readlines()[2:-4]
        Bvolume = sum([ float(Bader[i].split()[6]) for i in range(len(Bader))])
        volume = atoms.get_volume()
        mass = atoms.get_masses()

        file = open(outputfile,'a+')
        file.write("\nThe sum off all Bader atomic volumes is %.3G cm^3/g\n" %((volume-Bvolume)/sum(mass) *1E-24/6.02214086E-23))
        file.close()






###################################################################################################################################

try:
    fd = open(inputFile[0])
    if fd:
        energy(inputFile[0],outputfile)
except IOError: 
    print("   "+inputFile[0]+" is not a valid file!")
    pass
try:
    fd = open(inputFile[4])
    if fd:
        Bader_Volume(inputFile[1],inputFile[4],outputfile)
except IOError:
    print("   "+inputFile[4]+" is not a valid file!")
    pass
try:
    fd = open(inputFile[1])
    if fd:
        structure(inputFile[1], outputfile, sys.argv)
except IOError:
    print("   "+inputFile[1]+" is not a valid file!")
    pass
try:
    fd = open(inputFile[2])
    if fd:
        Freq(inputFile[2],outputfile)
except IOError:
    print("   "+inputFile[2]+" is not a valid file!")
    pass
try:
    fd = open(inputFile[3])
    if fd:
        WF(inputFile[3],outputfile)
except IOError: 
    print("   "+inputFile[3]+" is not a valid file!")
    pass
try:
    fd = open(inputFile[4])
    if fd:
        charges(inputFile[1], inputFile[4],outputfile, sys.argv)
except IOError:
    print("   "+inputFile[4]+" is not a valid file!")
    pass
try:
    fd = open(inputFile[5])
    if fd:
        spins(inputFile[1], inputFile[5],outputfile, sys.argv)
except IOError:
    print("   "+inputFile[5]+" is not a valid file!")
    pass

