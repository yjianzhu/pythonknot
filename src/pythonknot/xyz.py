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

def write_xyz(file, coords):
    """Write xyz file with coords, coords is a N_frames*N_atoms*3 np array, all the frames have same number of atoms"""
    N_atoms = coords.shape[1]
    N_frames = coords.shape[0]
    with open(file, 'w') as f:
        for i in range(N_frames):
            f.write(str(N_atoms)+'\n')
            f.write('frame '+str(i)+'\n')
            for j in range(N_atoms):
                f.write('1 '+str(coords[i,j,0])+' '+str(coords[i,j,1])+' '+str(coords[i,j,2])+'\n')