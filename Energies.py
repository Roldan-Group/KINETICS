'''
    by Alberto Roldan

    Reads inputs containing the energies and calculate
        Ecohesion
        Eadhesion
        TotalEnergy

    Versions:
        Alberto: 08/2019


    STRUCTURE:

'''

import numpy as np
from ase.io import read
from Coordination import Coordination
from GCN import Generalised_coodination
from Areas import Areas
from Library import isolated_atoms, surf_energies, sites, morse_potential, ecoh_bulk, ecoh_trend, morse_potential_depth


class Energies:
    def __init__(self, inputfiles, isolated_support, isolated_cluster, cluster_elements, support):
        system = read(inputfiles[0], index=-1)
        cluster = read(isolated_cluster)
        slab = read(isolated_support)

        e_system = system.get_total_energy()
        e_cluster = cluster.get_total_energy()
        e_slab = slab.get_total_energy()
        e_atoms = sum([isolated_atoms(system[i].symbol) for i in range(len(system))
                       if system[i].symbol in cluster_elements])

        coordination = Coordination(inputfiles[1], cluster_elements, support)
        c_coord = coordination.cluster_coordinating
        gcn = Generalised_coodination(inputfiles[1], cluster_elements, support)
        c_surf = gcn.cluster_surface_index
        area = Areas(inputfiles[1], cluster_elements, support)
        c_surf_area = area.cluster_surface_area

        self.cohesion = float((e_cluster - e_atoms) / len(cluster))
        self.adhesion = float(e_system - (e_cluster + e_slab))
#        if sum(coordination.site_cluster_coordination) > 0:
#            self.normalised_adhesion = self.adhesion / sum(coordination.site_cluster_coordination)
#        else:
#            self.normalised_adhesion = self.adhesion
        self.binding = float((e_system - (e_slab + e_atoms))/len(cluster))
        self.e_total = float(e_system)

        self.e_cluster_surface = surface_energy(system, c_coord, c_surf, c_surf_area)


def surface_energy(system, c_coord, c_surf, c_surf_area):
    ev_to_joules = 1.60218E-19
    e_c_surf = sum([surf_energies(system[i].symbol, len(c_coord[str(i)])) for i in c_surf])
    return float(e_c_surf * ev_to_joules / (c_surf_area * 1e-20))                  # J/m^2


class Energy_prediction:
    def __init__(self, inputfile, isolated_support, cluster_elements, support):
        system = read(inputfile)
        slab = read(isolated_support)
        e_slab = slab.get_total_energy()

        coordination = Coordination(inputfile, cluster_elements, support)
        s_coord = coordination.support_coordinating
        s_index = coordination.sites_index_all
# D_eq for Morse Potential is defined by i_c and i_cc of the interface by a second order polynomial or logarithm function
        i_c = coordination.interface_cluster
        i_cc = coordination.interface_cc_average
        c_coord = coordination.cluster_coordinating

        gcn = Generalised_coodination(inputfile, cluster_elements, support)
        c_surf = gcn.cluster_surface_index

        area = Areas(inputfile, cluster_elements, support)
        c_surf_area = area.cluster_surface_area

        self.e_adhesion(system, support, s_coord, s_index, i_c, i_cc)
        self.e_cohesion(system, c_coord)
        self.e_system = self.e_adh + e_slab + (self.e_coh * int(coordination.cluster_size) + self.e_atom)
        self.e_binding = (self.e_system - (e_slab + self.e_atom))/int(coordination.cluster_size)
        self.e_cluster_surface = surface_energy(system, c_coord, c_surf, c_surf_area)

    def e_adhesion(self, system, support, s_coord, s_index, i_c, i_cc):
        self.e_adh = 0
        for site in sites(support):
            for i in s_coord:
                distances = []
# D_eq for Morse Potential is defined by i_c and i_cc of the interface by a second order polynomial or log function
                d_eq = morse_potential_depth(support, system[i].symbol, i_c*i_cc)         # R^2 = 0.94    <<< CHECK!
                if len(s_coord[i]) > 0:
                    for j in s_coord[i]:
                        if system[j].symbol is site:
                            distances.append([i, j, system.get_distance(i, j, mic=True, vector=False)])
                    else:
                        for j in s_index[site]:
                            distances.append([i, j, system.get_distance(i, j, mic=True, vector=False)])
                else:
                    for j in s_index[site]:
                        distances.append([i, j, system.get_distance(i, j, mic=True, vector=False)])

                distances.sort(key=lambda x: x[2])
                self.e_adh += morse_potential(support, site, system[distances[0][0]].symbol, distances[0][2], d_eq)\
                              / (len(s_coord) * len(sites(support)))    # WHAT IF s_coord = 0?!?!?

    def e_cohesion(self, system, c_coord):
        self.e_coh = 0
        self.e_atom = 0
        for i in c_coord:
            cc = len(c_coord[i])
            bulk_coh, bulk_coord = ecoh_bulk(system[int(i)].symbol)
            a = ecoh_trend(system[int(i)].symbol)
            self.e_coh += float(1/abs(bulk_coh) * (bulk_coh*np.log(a)/np.log(a/(a+bulk_coord)) -
                                                   (bulk_coh/np.log(a/(a+bulk_coord)))*np.log(a+cc)))
            self.e_atom += float(isolated_atoms(system[int(i)].symbol))

