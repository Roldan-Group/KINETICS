'''
    Versions:
        Alberto: 08/2019
        Alberto: 09/2020

'''

from Coordination import Coordination
from GCN import Generalised_coodination
from Areas import Areas
from Zdistance import Cluster_mass_centre_surface_distance
from Energies import Energy_prediction
from WriteData import Write_labels, write_results, write_outcar

#####################################################################################################
cluster_elements = ["Au"]                                                   # Elements in the Cluster
support = "MgO"                                                     # Surface name
structurefile = "POSCAR"
isolated_support = "/home/alberto/RESEARCH/OTHER/Supported/DataSet/MgO/MgO/Surface/OUTCAR"
####################################################################################################

coordination = Coordination(structurefile, cluster_elements, support)
gcn = Generalised_coodination(structurefile, cluster_elements, support)
area = Areas(structurefile, cluster_elements, support)
z_distance = Cluster_mass_centre_surface_distance(structurefile, cluster_elements, support)
energies = Energy_prediction(structurefile, isolated_support, cluster_elements, support)

labels = ["N", coordination.site_cluster_coordination_label, coordination.cluster_coord_labels,
				coordination.support_cluster_min_distance_labels, z_distance.zlabels, "GCN", "c_i_area",
				"c_s_area", " Esurf", "  Ecoh", "  Eadh", "   Eb  ", "Etotal", "   structure_path"]
values = [coordination.cluster_size, coordination.site_cluster_coordination, coordination.cluster_ave_coordination,
		  coordination.support_cluster_min_distance, z_distance.cluster_cm_surface_distance, gcn.gcn_average,
		  area.cluster_interface_area, area.cluster_surface_area, energy.e_cluster_surface, energy.e_cohesion,
		  energy.e_adhesion, energy.e_binding, energy.e_system,  name]

Write_labels("Predicted.dat", labels)
write_results("Predicted.dat", values)
write_outcar(structurefile, energies.e_system)

