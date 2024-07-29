import numpy as np

def read_xyz(file):
    """Read xyz file and get all the conf in same file, so return N_frames*n_atoms*3 np array"""
    with open(file) as f:
        lines = f.readlines()
        N_frames = int(len(lines)/(int(lines[0])+2))
        N_atoms = int(lines[0])
        print(N_frames, N_atoms)
        coords = np.zeros((N_frames, N_atoms, 3))
        for i in range(N_frames):
            lines = lines[2:]
            for j in range(N_atoms):
                line = lines[j].split()
                coords[i,j,0] = float(line[1])
                coords[i,j,1] = float(line[2])
                coords[i,j,2] = float(line[3])
            lines = lines[N_atoms:]

    return coords