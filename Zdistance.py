'''
        provides the distance from the cluster's centre of mass to the average top surface
'''


from ase.io import read
from Library import sites
from Coordination import Coordination



class Cluster_surface_distance:
    def __init__(self, inputfile, cluster_elements, support):
        system = read(inputfile)
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        support_index = [system[i].index for i in range(len(system)) if system[i].symbol in sites(support)]

        support_zmax = max([system[i].position[2] for i in support_index])
        sites_index = [i.index for i in system if i.symbol in sites(support) and i.position[2] >= support_zmax - 1]  # gets the site atoms index in the support
        pos = sum([(system[i].position[2] * system[i].mass) for i in cluster_index])
        self.cluster_mass_centre = float(pos/sum([system[i].mass for i in cluster_index]))

        self.cluster_cm_surface_distance = float(self.cluster_mass_centre - sum([system[i].position[2] for i in sites_index])\
                                           /len(sites_index))
        self.zlabels = str("Zdist")

        coordination = Coordination(inputfile, cluster_elements, support)
        c_interface = coordination.interface_cluster_index
        if len(c_interface) > 0:
            c_height_average = sum([system[i].position[2] for i in c_interface])/len(c_interface)
        else:
            cluster_zmin = min([system[i].position[2] for i in cluster_index])
            c_inter_index = [i for i in cluster_index if system[i].position[2] <= cluster_zmin + 1]
            c_height_average = sum([system[i].position[2] for i in c_inter_index])/len(c_inter_index)

        self.interface_height = float(c_height_average - sum([system[i].position[2] for i in sites_index])/len(sites_index))
#        print(c_height_average, sum([system[i].position[2] for i in sites_index])/len(sites_index), self.interface_height)




