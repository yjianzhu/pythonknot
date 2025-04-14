import numpy as np

# design for generating knot configurations for diffrent simulation software
# lammps read_data, gromacs gmx editconf, openmm pdb
# get the coordinates of the knot, then write to the corresponding file
class KnotGenerator:
    def __init__(self, knot):
        """ knot should be a numpy array with shape (N_atoms, 3)"""
        self.knot = knot

    def write_lammps_data(self, filename):
        pass

    def write_gromacs_gro(self, filename):
        pass

    def write_openmm_pdb(self, filename):
        pass

    def write_xyz(self, filename):
        pass

