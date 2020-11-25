'''
    Versions:
        Alberto: 08/2019

    STRUCTURE:
        - check structure and energy (OUTCAR) files
        - check boundary conditions
        - Import coordination data
        - Import energy from the OUTCAR

'''

from ase.io import read

class Write_labels:
    def __init__(self, outfile, data):
        line=[]; n=0
        for value in data:
            if type(value) is list:
                for label in value:
                    line.append(self.label_definition(label.strip(), n))
                    n += 1
            else:
                line.append(self.label_definition(value.strip(), n))
                n += 1

        heading=[]
        for value in data:
            if type(value) is list:
                [heading.append(i) for i in value]
            else:
                heading.append(value)

        output = open(outfile,"w+")
        output.write("#\n#\n#\n")
        for dat in line:
            output.write("#     %s\n" % dat)
        output.write("#\n#\n#\n#")
        for i, dat in enumerate(heading):
            if i <= 3:
                output.write(" {:>4}" .format(dat))
            else:
                output.write(" {:>9}" .format(dat))
        output.write("\n")

        output.close()
                   
    def label_definition(self, value, n):

        if value is "N":
            comment = str("Column {0:3d} = {1:^10s} = Total number of atoms forming the cluster" .format(n, value))
        elif value is 'i_c':
            comment = str("Column {0:3d} = {1:^10s} = Number of cluster atoms coordinating the surface" .format(n, value))
        elif value.startswith('cs_') is True:
            comment = str("Column {0:3d} = {1:^10s} = Number of surface sites coordinating with the cluster" .format(n, value))
        elif value is "cc":
            comment = str("Column {0:3d} = {1:^10s} = Average atomic coordination within the cluster" .format(n, value))
        elif value is "i_cc":
            comment = str("Column {0:3d} = {1:^10s} = Average coordination of cluster atoms at the interface within the cluster only" .format(n, value))
        elif value.startswith('dist_') is True:
            comment = str("Column {0:3d} = {1:^10s} = Average of minimum distances (in Å) between the surface sites and the clusters atoms" .format(n, value))
        elif value is 'cs_height':
            comment = str("Column {0:3d} = {1:^10s} = Distance (in Å) between the surface and the cluster" .format(n, value))
        elif value.startswith('Zdist') is True:
            comment = str("Column {0:3d} = {1:^10s} = Distance (in Å) between the average surface hight and the cluster's centre of mass" .format(n, value))

        elif value is "GCN":
            comment = str("Column {0:3d} = {1:^10s} = Average generalised coordination number for the atoms in the cluster excluding the coordination with the support" .format(n, value))
        elif value is "c_i_area":
            comment = str("Column {0:3d} = {1:^10s} = Cluster interface area (in Angstrom^2) -- check Library" .format(n, value))
        elif value is "c_s_area":
            comment = str("Column {0:3d} = {1:^10s} = Area exposed by the cluster excluding the interface (in Angstrom^2 per atom)" .format(n, value))

        elif value is 'Esurf':
            comment = str("Column {0:3d} = {1:^10s} = Exposed surface energy (in J/m^2) of the cluster (not interface) -- check Library" .format(n, value))
        elif value is 'Ecoh':
            comment = str("Column {0:3d} = {1:^10s} = Cohesion energy per cluster atom (in eV/atom) == (Ecluster -( N * Eatom))/N" .format(n, value))
        elif value is 'Eadh':
            comment = str("Column {0:3d} = {1:^10s} = Cluster adhesion energy (in eV) == Esystem - (Esurface + Ecluster)" .format(n, value))
        elif value is 'Eb':
            comment = str("Column {0:3d} = {1:^10s} = Binding energy per cluster atom (in eV/atom) == (Esystem -(Esurface + N * Eatom))/N" .format(n, value))
        elif value is 'Etotal':
            comment = str("Column {0:3d} = {1:^10s} = Total energy of the system (in eV)" .format(n, value))

        else:
            comment = str("Column {0:3d} = {1:^10s} = label not identified" .format(n, value))


        return comment



def write_results(outfile,data):
    line = []
    for value in data:
        if type(value) is list:
            [line.append(i) for i in value]
        elif type(value) is dict:
            [line.append(value[i]) for i in value]
        else:
            line.append(value)

    output = open(outfile, "a+")
    output.write("  ")
    for n, dat in enumerate(line):
        if type(dat) is int:
            output.write(" {:3d} " .format(dat))
        elif type(dat) is float:
            output.write(" {:> 9.3f}"  .format(dat))
        elif type(dat) is str:
            output.write(" # {:>s}" .format(dat))
        else:
            output.write(" {:> 9.3f}" .format(dat))

    output.write("\n")
    output.close()


def write_outcar(structure_file, energy):
    system = read(structure_file)
    output = open("OUTCAR", "w+")
    output.write(" energy  without entropy= {:> 12.6f} energy(sigma->0) = {:> 12.6f}\n" .format(energy, energy))
    output.write(" POSITION\n---------------------------\n")               # to adapt the reading from GA

    xyz = system.get_positions()
    for i in range(len(xyz)):
        x, y, z = xyz[i]
        output.write(" %.8f  %.8f  %.8f" %(x, y, z))
        output.write("\n")
    output.write("---------------------------\n")
    output.close()


