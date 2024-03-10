"""
nvtmc
"""
import random
import numpy

import pathadd
import data_types as dt


class SilicaPolymerizationMC(dt.BccLattice3D):
    """ class containing polymerization reaction  """

    def __init__(self, cfg_obj: dt.CfgObj):
        """ Initialize the class
        Args:
            cfg_obj: the configuration obejct
        """
        self._cfg_obj = cfg_obj
        super().__init__(*self._cfg_obj.size)

        self.molecule_list = numpy.zeros(self._cfg_obj.nsio4, dtype=int)
        self.head_list = numpy.zeros(self.nsites, dtype=int)

    def are_2rings_formed(self, isite, ihead=None):
        """
        check if 2 rings are formed
        Args:
            isite: ID of the site
            ihead: head of isite. if ihead is None then the value of ihead is calculated from head_list

        Returns:
            true if it/will participate in 2 rings
        """
        if ihead is None:
            ihead = self.head_list[isite]
        ihead_sites = self.get_neighbors_isite(isite, ihead)
        for idx in range(6):
            jsite = self.n2list[isite, idx]
            if self.spin[jsite] == dt.Atom.SILICON.value:
                jhead = self.head_list[jsite]
                jhead_sites = self.get_neighbors_isite(jsite, jhead)
                if len(numpy.intersect1d(ihead_sites, jhead_sites)) > 1:
                    return True

        return False

    def can_place_molecule(self, isite, ihead):
        """ check if moldcule can be placed on isite with pointer as ihead """
        head_site_list = self.get_neighbors_isite(isite, ihead)
        can_place = False
        if self.spin[isite] == dt.Atom.VACANT.value:

            for jsite in head_site_list:
                jspin = self.spin[jsite]
                if jspin == dt.Atom.VACANT.value or jspin == dt.Atom.Oxygen.value:
                    two_rings_formed = self.are_2rings_formed(isite, ihead)
                    can_place = not two_rings_formed
                else:
                    can_place = False
        else:
            can_place = False

        return can_place

    def add_molecule_to_site(self, molecule_id, isite, ihead):
        """ add a molecule to the site """
        self.spin[isite] = dt.Atom.SILICON.value
        self.molecule_list[molecule_id] = isite
        head_site_list = self.get_neighbors_isite(isite, ihead)
        for neigh_site in head_site_list:
            self.spin[neigh_site] = self.spin[neigh_site] + dt.Atom.Oxygen.value

    def initialize_system(self):
        """ initialize the system """

        for molecule in range(len(self.molecule_list)):
            isite = random.randint(0, self.nsites - 1)
            ihead = random.randint(0, 1)

            attempts = 0
            can_place = self.can_place_molecule(isite, ihead)
            while not can_place and attempts < 100:
                can_place = self.can_place_molecule(isite, ihead)
                attempts += 1
                if can_place:
                    self.add_molecule_to_site(molecule, isite, ihead)
                    break

    def run(self):
        """ the MC running program """
        self.initialize_system()


def main():
    """ the main executable function """

    raw_data = {
        'temperature': 0.15,
        'lx': 8,
        'ly': 8,
        'lz': 8,
        'nsio4': 16,
        'nsteps': 10000,
        'neqsteps': 5000,
        'dsteps': 100
    }
    cfg_obj = dt.CfgObj(raw_data)
    SilicaPolymerizationMC(cfg_obj).run()


if __name__ == '__main__':
    main()
