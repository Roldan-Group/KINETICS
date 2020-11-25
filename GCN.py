'''
        provides information related to the Generalised Coordination Number
'''


from ase.io import read
from Coordination import Coordination
from Library import ecoh_bulk



class Generalised_coodination:
    def __init__(self, inputfile, cluster_elements, support):
        system = read(inputfile)
        coordination = Coordination(inputfile, cluster_elements, support)
        c_coord = coordination.cluster_coordinating
        s_coord = coordination.support_coordinating

        self.cluster_gcn(system, cluster_elements, c_coord, s_coord)

    def cluster_gcn(self, system, cluster_elements, c_coord, s_coord):
        cluster_index = [system[i].index for i in range(len(system)) if system[i].symbol in cluster_elements]

        self.gcn_average = 0
        self.gcn = {}
        self.cluster_bulk_index = []
        self.cluster_surface_index = []
        for n in cluster_index:
            i_gcn = 0
            e_coh, bulk_coordination = ecoh_bulk(system[n].symbol)
            if len(c_coord[str(n)]) > 0:
                for j in c_coord[str(n)]:
                    i_gcn += len(c_coord[str(j)])            # coordination of the coordinating atom of n
                self.gcn[n] = i_gcn/bulk_coordination
            if len(c_coord[str(n)]) < bulk_coordination and n not in self.cluster_surface_index:  # exposed atoms at the surface
                self.cluster_surface_index.append(n)
            if len(c_coord[str(n)]) == bulk_coordination and n not in self.cluster_bulk_index:
                self.cluster_bulk_index.append(n)

        if len([self.gcn[i] for i in self.gcn]) > 0:
            self.gcn_average = float(sum([self.gcn[i] for i in self.gcn])/len(self.gcn))

