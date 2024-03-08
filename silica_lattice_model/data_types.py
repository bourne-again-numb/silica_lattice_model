"""
The module containing custom datatypes
"""

import numpy
from enum import Enum


class Atom(Enum):
    """ enum types for atoms """
    VACANT = 0
    SILICON = 1
    OH = 2
    Oxygen = 4




class BccLattice3D:
    """ The BCC lattice class """

    def __init__(self, lx, ly, lz):
        """
        initialize the class
        Args:
            lx: size of the lattice in X
            ly: size of lattice in Y
            lz: size of the lattice in z
        """
        self._lx = lx
        self._ly = ly
        self._lz = lz
        self.nsites = 2 * lx * ly * lz

        self.rx = numpy.zeros(self.nsites, dtype=int)
        self.ry = numpy.zeros(self.nsites, dtype=int)
        self.rz = numpy.zeros(self.nsites, dtype=int)
        self.tag_to_site = numpy.zeros(8 * self._lx * self._ly * self._lz + 1, dtype=int)
        self.n1list = numpy.zeros((self.nsites, 8), dtype=int)
        self.n2list = numpy.zeros((self.nsites, 6), dtype=int)
        self.spin = numpy.zeros(self.nsites, dtype=int)
        for idx in range(len(self.spin)):
            self.spin[idx] = Atom.VACANT.value

        self._setup_lattice()

    def _tag(self, ix, iy, iz):
        """ get side ID """
        # itag = int(ix + 2 * self._lx * (iy + 2 * self._ly * iz))
        itag = (ix - 1) * 2 * self._ly * 2 * self._lz + (iy - 1) * 2 * self._lz + iz
        return itag

    def _setup_site_tags(self):
        """ setup site tags """
        isite = 0
        for ix in range(1, 2 * self._lx + 1):
            for iy in range(1, 2 * self._ly + 1):
                for iz in range(1, 2 * self._lz + 1):
                    icheck = (ix % 2) + (iy % 2) + (iz % 2)
                    if icheck == 0 or icheck == 3:
                        itag = self._tag(ix, iy, iz)
                        self.tag_to_site[itag] = isite
                        # print(isite, itag, '\t', ix, iy, iz)
                        self.rx[isite] = ix
                        self.ry[isite] = iy
                        self.rz[isite] = iz

                        isite = isite + 1
        assert isite == self.nsites  # sanity check the final size number needs to be nsites

    def _setup_n1list(self):
        """ setup lattice """
        isite = 0
        for ix in range(1, 2 * self._lx + 1):
            for iy in range(1, 2 * self._ly + 1):
                for iz in range(1, 2 * self._lz + 1):
                    icheck = (ix % 2) + (iy % 2) + (iz % 2)

                    if icheck == 0 or icheck == 3:

                        # periodic boundary conditions
                        ix1 = ix - 1
                        if ix1 < 1:
                            ix1 = 2 * self._lx
                        ix2 = ix + 1
                        if ix2 > 2 * self._lx:
                            ix2 = 1

                        iy1 = iy - 1
                        if iy1 < 1:
                            iy1 = 2 * self._ly
                        iy2 = iy + 1
                        if iy2 > 2 * self._ly:
                            iy2 = 1

                        iz1 = iz - 1
                        if iz1 < 1:
                            iz1 = 2 * self._lz
                        iz2 = iz + 1
                        if iz2 > 2 * self._lz:
                            iz2 = 1

                        self.n1list[isite, 0] = self.tag_to_site[self._tag(ix1, iy1, iz1)]
                        self.n1list[isite, 1] = self.tag_to_site[self._tag(ix2, iy1, iz1)]
                        self.n1list[isite, 2] = self.tag_to_site[self._tag(ix2, iy2, iz1)]
                        self.n1list[isite, 3] = self.tag_to_site[self._tag(ix1, iy2, iz1)]
                        self.n1list[isite, 4] = self.tag_to_site[self._tag(ix1, iy1, iz2)]
                        self.n1list[isite, 5] = self.tag_to_site[self._tag(ix2, iy1, iz2)]
                        self.n1list[isite, 6] = self.tag_to_site[self._tag(ix2, iy2, iz2)]
                        self.n1list[isite, 7] = self.tag_to_site[self._tag(ix1, iy2, iz2)]

                        isite = isite + 1

        assert isite == self.nsites  # sanity check the final size number needs to be nsites

    def _setup_n2list(self):
        """ setup lattice """
        isite = 0
        for ix in range(1, 2 * self._lx + 1):
            for iy in range(1, 2 * self._ly + 1):
                for iz in range(1, 2 * self._lz + 1):
                    icheck = (ix % 2) + (iy % 2) + (iz % 2)

                    if icheck == 0 or icheck == 3:

                        # periodic boundary conditions
                        ix1 = ix - 2
                        if ix1 < 1:
                            ix1 = ix1 + 2 * self._lx
                        ix2 = ix + 2
                        if ix2 > 2 * self._lx:
                            ix2 = ix2 - 2 * self._lx

                        iy1 = iy - 2
                        if iy1 < 1:
                            iy1 = iy1 + 2 * self._ly
                        iy2 = iy + 2
                        if iy2 > 2 * self._ly:
                            iy2 = iy2 - 2 * self._ly

                        iz1 = iz - 2
                        if iz1 < 1:
                            iz1 = iz1 + 2 * self._lz
                        iz2 = iz + 2
                        if iz2 > 2 * self._lz:
                            iz2 = iz2 - 2 * self._lz

                        self.n2list[isite, 0] = self.tag_to_site[self._tag(ix1, iy, iz)]
                        self.n2list[isite, 1] = self.tag_to_site[self._tag(ix2, iy, iz)]
                        self.n2list[isite, 2] = self.tag_to_site[self._tag(ix, iy1, iz)]
                        self.n2list[isite, 3] = self.tag_to_site[self._tag(ix, iy2, iz)]
                        self.n2list[isite, 4] = self.tag_to_site[self._tag(ix, iy, iz1)]
                        self.n2list[isite, 5] = self.tag_to_site[self._tag(ix, iy, iz2)]

                        isite = isite + 1

        assert isite == self.nsites  # sanity check the final size number needs to be nsites

    def _setup_lattice(self):
        """ setup lattice """
        # debine tag_to_site
        self._setup_site_tags()
        self._setup_n1list()
        self._setup_n2list()

    @staticmethod
    def neighbors_for_head(head):
        """
        get the neighbors index for a given head.
        Args:
            head: the value of head can be 0 or 1

        Returns:
            numpy.array of neighbors ids
        """
        if head == 0:
            result = numpy.array([0, 2, 5, 7], dtype=int)
        elif head == 1:
            result = numpy.array([1, 3, 4, 6], dtype=int)
        else:
            raise ValueError(f'Head can only take value 0 or 1: {head}')

        return result

    def get_neighbors_isite(self, isite, head):
        """ get the heighbors for isite given a head """
        neighs_list = numpy.zeros(4, dtype=int)
        neighs_id = self.neighbors_for_head(head)
        for idx, head_id in enumerate(neighs_id):
            neighs_list[idx] = self.n1list[isite, idx]

        return neighs_list


class CfgObj:
    """ the config object class """

    def __init__(self, raw_data):
        """ initialize the class """
        self._raw_data = raw_data

    @property
    def size(self):
        return self._raw_data['lx'], self._raw_data['ly'], self._raw_data['lz']

    @property
    def nsio4(self):
        return self._raw_data['nsio4']

    @property
    def nsda(self):
        return 0

    @property
    def temperature(self):
        return self._raw_data['temperature']

    @property
    def mc_steps(self):
        return self._raw_data['nsteps'], self._raw_data['neqsteps'], self._raw_data['dsteps']
