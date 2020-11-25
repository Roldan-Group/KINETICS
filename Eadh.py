'''
    by Alberto Roldan

'''

from ase.io import read
from Coordination import Coordination
from Library import sites, morse_potential


class Adhesion_energy:
    def __init__(self, inputfiles, isolated_support, isolated_cluster, cluster_elements, support):
        system = read(inputfiles[1])
        coordination = Coordination(inputfiles[1], cluster_elements, support)
        s_coord = coordination.support_coordinating
        s_index = coordination.sites_index_all

        e_atom = 0
        for site in sites(support):
            for i in s_coord:
                distances = []
                if len(s_coord[i]) > 0:
                    for j in s_coord[i]:
                        if system[j].symbol is site:
                            distances.append([i, j, system.get_distance(i, j, mic=True, vector=False)])
                    distances.sort(key=lambda x: x[2])
                    if len(distances) > 0:
#                        print(len(distances),system[distances[0][1]].symbol, system[distances[0][0]].symbol)
                        e_atom += morse_potential(support, site, system[distances[0][0]].symbol, distances[0][2])
                    else:
                        for j in s_index[site]:
                            distances.append([i, j, system.get_distance(i, j, mic=True, vector=False)])
                        distances.sort(key=lambda x: x[2])
                        e_atom += morse_potential(support, site, system[distances[0][0]].symbol, distances[0][2])
                else:
                    for j in s_index[site]:
                        distances.append([i, j, system.get_distance(i, j, mic=True, vector=False)])
                    distances.sort(key=lambda x: x[2])
                    e_atom += morse_potential(support, site, system[distances[0][0]].symbol, distances[0][2])
#        print(d_range, e_0)
        self.e_adh_5 = e_atom
#        print(">>>>", round(energies.adhesion,2)," || method_5 = ",round(self.e_adh_5,2), "|| len s_coo=",s_coord)