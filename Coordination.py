'''
    by David Saadatmandi and Alberto Roldan

    Reads an input with geometries in XYZ format and calculates the
    the cluster coordination with the surface (cs)
    and coordination between metal atoms forming the cluster (cc).

'''


import numpy as np
from ase import neighborlist
from ase.io import read
from ase.build import bulk, molecule
from Library import opt_atom_distance, sites


class Coordination:
    def __init__(self, inputfile, cluster_elements, support):
        system = read(inputfile)

        self.cluster_coord_labels = ["cc"]
        self.site_cluster_coordination_label = [("cs_" + site) for site in sites(support)]
        self.support_cluster_min_distance_labels = [("dist_" + site) for site in sites(support)]
        self.support_cluster_min_distance, self.support_coordinating, self.sites_index_all = self.sites_coordination(
                                                                                    system, cluster_elements, support)
        self.interface_cluster = len(self.interface_cluster_index)
        self.cluster_size, self.cluster_ave_coordination, self.cluster_coordinating = self.cluster_coordination(system,
                                                                                                    cluster_elements)

    def cluster_coordination(self, system, cluster_elements):
        average_coordination = 0
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
#        del cluster_atoms[[atom.index for atom in cluster_atoms if atom.symbol not in cluster_elements]]
        coordinating = {}
        if len(cluster_index) > 1:
            cutoff = sum([self.bulk_distances(system[i].symbol)*1.2 for i in cluster_index])\
                     /len(cluster_index)
#            cutoff = neighborlist.natural_cutoffs(system, mult=1.2)
            a, b = neighborlist.neighbor_list('ij', system, cutoff)
            for i in cluster_index:
                coordinating[str(i)] = [b[n] for n in range(len(a)) if a[n] == i and b[n] in cluster_index]
            average_coordination = sum([len(coordinating[str(i)]) for i in cluster_index])/len(cluster_index)
        else:
            coordinating[str(cluster_index[0])] = []

        return int(len(cluster_index)), float(average_coordination), coordinating


    def sites_coordination(self, system, cluster_elements, support):

        self.interface_cluster_index = []
        self.interface_support_index = []
        self.site_cluster_coordination = {}
        self.interface_cc_average = 0
        s_sites = {}
        coordinating = {}
        cs_distance = {}
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]
        for site in sites(support):
            cs_distance[site] = 0
            distances = []
            self.site_cluster_coordination[site] = 0
            support_zmax = max([i.position[2] for i in system if i.symbol == site])
            sites_index = [i.index for i in system if i.symbol == site and i.position[2] >= support_zmax -1]  # gets the site atoms index in the support
            s_sites.update({site: sites_index})
            optimised_distance = [opt_atom_distance(support, site, i) for i in cluster_elements]
            distance_cutoff = sum(optimised_distance) / len(optimised_distance)
            cutoff = neighborlist.natural_cutoffs(system, distance_cutoff)
            a, b, d = neighborlist.neighbor_list('ijd', system, cutoff)

            for n in cluster_index:
                coord = [b[i] for i in range(len(a)) if a[i] == n and b[i] in sites_index and d[i] <= distance_cutoff*1.5]
                if n not in coordinating:
                    coordinating[n] = coord
                else:
                    for i in coord:
                        coordinating[n].append(i)
                if len(coord) > 0:
                    if n not in self.interface_cluster_index:
                        self.interface_cluster_index.append(n)
                    for j in coord:
                        if j not in self.interface_support_index:
                            self.interface_support_index.append(j)
                    distances.append(min([d[i] for i in range(len(a)) if a[i] == n and b[i] in coord]))

            if len(distances) > 0:
                cs_distance[site] = float(sum(distances)/len(distances))
                self.site_cluster_coordination[site] = int(len(distances))
#                print("coord>>", site, n, distances, cs_distance[site])
            else:
                dist_array = []
                for i in cluster_index:
                    for j in sites_index:
                        dist_array.append(system.get_distance(i, j, mic=True, vector=False))
                cs_distance[site] = float(min(dist_array))

        for n in cluster_index:
            self.interface_cc_average += float(len([b[i] for i in range(len(a)) if a[i] == n and b[i] in cluster_index
                                                   and d[i] <= distance_cutoff*1.5])/len(cluster_index))

        return cs_distance, coordinating, s_sites

    def bulk_distances(self, element):
        try:
            bulkdistance = sum(bulk(element).get_cell_lengths_and_angles()[
                               0:3]) / 3  # the average distance between atoms in the bulk
        except ValueError:
            molec = str(element.symbol + '2')
            bulkdistance = molecule(molec).get_distance(0, 1)  # interatomic distance in a diatomic molecule
        except ValueError:
            bulkdistance = 1  # minimum distance by default
        return bulkdistance


