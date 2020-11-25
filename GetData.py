'''
	Versions:
		Alberto: 08/2019

	STRUCTURE:
		- check structure and energy (OUTCAR) files
		- check boundary conditions
		- Import coordination data
		- Import energy from the OUTCAR

'''

import os, sys
from CheckFiles import checkFiles
from Coordination import Coordination
from GCN import Generalised_coodination
from Areas import Areas
from Zdistance import Cluster_surface_distance
from Energies import Energies, Energy_prediction
from WriteData import Write_labels, write_results

if sys.argv:
	arguments=sys.argv
else:
	arguments[1]= "."

path = os.getcwd()
name = path.split("/")[-4]+"/"+path.split("/")[-3]+"/"+path.split("/")[-2]+"/"+path.split("/")[-1]

cluster_elements = ["Au"]                           # Elements in the Cluster
support = "MgO"                             # Surface name
inputfiles = ["OUTCAR", "CONTCAR"]
isolated_cluster = arguments[1]+"/gas/OUTCAR"           # ".../OTHER/Supported/MgO/Au/Basin_Hopping/1Au/gas/OUTCAR"
isolated_support = "/home/alberto/RESEARCH/OTHER/DATASET/rPBE/Supports/MgO/MgO/Surface/OUTCAR"


#from ase.io import read
#atoms = read(inputFiles[1])
#from ase.visualize import view
#view(atoms)
print(name)

#checkFiles([structureFile,energyFile,isolatedCluster])      	# it crashes with GenerationModel for Ecoh << eleNames
								# comment it the get isolated clusters Ecoh only
coordination = Coordination(inputfiles[1], cluster_elements, support)
gcn = Generalised_coodination(inputfiles[1], cluster_elements, support)
area = Areas(inputfiles[1], cluster_elements, support)
z_distance = Cluster_surface_distance(inputfiles[1], cluster_elements, support)
energy = Energies(inputfiles, isolated_support, isolated_cluster, cluster_elements, support)


labels = ["N", "i_c", coordination.site_cluster_coordination_label, "i_cc", coordination.cluster_coord_labels,
				coordination.support_cluster_min_distance_labels, "cs_height", z_distance.zlabels, "GCN", "c_i_area",
				"c_s_area", "Esurf", "Ecoh", "Eadh", "Eb", "Etotal", "  structure_path"]
values = [coordination.cluster_size, coordination.interface_cluster, coordination.site_cluster_coordination,
		  coordination.interface_cc_average, coordination.cluster_ave_coordination, coordination.support_cluster_min_distance,
		  z_distance.interface_height, z_distance.cluster_cm_surface_distance, float(gcn.gcn_average),
		  area.cluster_interface_area, area.cluster_surface_area, energy.e_cluster_surface, energy.cohesion,
		  energy.adhesion, energy.binding, energy.e_total, name]

Write_labels("labels.txt", labels)
write_results("data.dat", values)


energies = Energy_prediction(inputfiles[1], isolated_support, cluster_elements, support)
#print(round(energy.adhesion,2), round(energies.e_adh,2))
output = open("Eadh.dat", "w+")
output.write("%.3f %.3f   # %s\n" %(energy.adhesion, energies.e_adh, name))
output.close()
#output = open("cc_Ecoh.dat", "w+")
#output.write("%.3f %.3f  # %s\n" %(coordination.cluster_ave_coordination, energy.cohesion, name))
#output.close()
#output = open("Esurf.dat", "w+")
#output.write("%.3f %.3f  # %s\n" %(energy.e_cluster_surface, energies.e_cluster_surface, name))
#output.close()
output = open("Eb.dat", "w+")
output.write("%.3f %.3f  # %s\n" %(energy.binding, energies.e_binding, name))
output.close()
output = open("Etotal.dat", "w+")
output.write("%.3f %.3f  # %s\n" %(energy.e_total, energies.e_system, name))
output.close()

